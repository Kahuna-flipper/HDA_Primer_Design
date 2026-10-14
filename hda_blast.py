"""BLAST-based specificity checks for tHDA primer sets.

Two analyses, both run against either a local BLAST database or NCBI (remote):

1. Primer specificity (Primer-BLAST-like)
   All primers are BLASTed (blastn-short). A hit counts as a *binding site* when
   the alignment reaches the primer's 3' end and the total number of mismatches
   (mismatches + gaps + unaligned 5' bases) is <= ``max_mismatches``. Binding
   sites on the same subject that face each other within ``max_product`` bp are
   reported as predicted products (any primer combination: F/R, F/F, R/R).

2. Amplicon uniqueness across species
   Each amplicon is BLASTed (blastn) and hits are grouped by species. For each
   amplicon this reports how many species carry a similar sequence, the best
   identity/coverage per species, and which non-target species are
   "cross-reactive" (identity >= ``min_identity`` over >= ``min_coverage`` of the
   amplicon).

Species are taken from BLAST taxonomy (``sscinames``) when available (NCBI
remote, or local DBs built with taxonomy). Otherwise they are parsed from the
subject title: ``[species=Genus species]`` tags written by
``build_local_db`` are used first, then the first two words of the title.

Remote searches go to NCBI's public BLAST servers and can take many minutes;
restrict them with ``entrez_query`` (e.g. ``'Bacteria[organism]'``) where possible.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

OUTFMT_FIELDS = ['qseqid', 'sacc', 'stitle', 'staxids', 'sscinames', 'pident', 'length',
                 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'slen',
                 'qlen', 'evalue', 'bitscore']
_INT_FIELDS = {'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send',
               'slen', 'qlen'}
_FLOAT_FIELDS = {'pident', 'evalue', 'bitscore'}
_MIN_REMOTE_VERSION = (2, 10)


# ----------------------------------------------------------------------------
# BLAST execution
# ----------------------------------------------------------------------------

def _blast_version(exe: str) -> tuple[int, int] | None:
    try:
        out = subprocess.run([exe, '-version'], capture_output=True, text=True).stdout
    except OSError:
        return None
    m = re.search(r'blastn: (\d+)\.(\d+)', out)
    return (int(m.group(1)), int(m.group(2))) if m else None


def find_blastn(require_remote: bool = False) -> str:
    """Locate blastn: $HDA_BLASTN, then the active Python env, then PATH."""
    candidates = [os.environ.get('HDA_BLASTN'),
                  str(Path(sys.executable).parent / 'blastn'),
                  shutil.which('blastn')]
    for exe in filter(None, candidates):
        ver = _blast_version(exe)
        if ver and (not require_remote or ver >= _MIN_REMOTE_VERSION):
            return exe
    raise RuntimeError(
        'No suitable blastn found (remote BLAST needs BLAST+ >= 2.10). Install with\n'
        '  conda install -n primer_design -c bioconda blast\n'
        'or set HDA_BLASTN=/path/to/blastn')


@dataclass
class BlastConfig:
    """Where and how to BLAST.

    Give ``db`` as a local database path, or set ``remote=True`` to use NCBI
    (``db`` then names an NCBI database, default ``core_nt``).
    """
    db: str | None = None
    remote: bool = False
    entrez_query: str | None = None   # remote only, e.g. 'Bacteria[organism]'
    threads: int = 4                  # local only
    blastn: str | None = None
    max_target_seqs: int = 5000      # per query; species counts are a lower bound if hit

    def __post_init__(self):
        if self.remote and self.db is None:
            self.db = 'core_nt'
        if self.db is None:
            raise ValueError('BlastConfig needs a local db path or remote=True')
        if self.blastn is None:
            self.blastn = find_blastn(require_remote=self.remote)


def run_blast(queries: dict[str, str], cfg: BlastConfig, task: str,
              evalue: float, extra: Iterable[str] = ()) -> list[dict]:
    """BLAST {id: seq} queries, return tabular hits as dicts (OUTFMT_FIELDS)."""
    if not queries:
        return []
    with tempfile.TemporaryDirectory(prefix='hda_blast_') as tmp:
        qpath = Path(tmp) / 'query.fa'
        qpath.write_text(''.join(f'>{k}\n{v}\n' for k, v in queries.items()))
        cmd = [cfg.blastn, '-task', task, '-query', str(qpath), '-db', cfg.db,
               '-evalue', str(evalue), '-max_target_seqs', str(cfg.max_target_seqs),
               '-outfmt', '6 ' + ' '.join(OUTFMT_FIELDS), *extra]
        if cfg.remote:
            cmd.append('-remote')
            if cfg.entrez_query:
                cmd += ['-entrez_query', cfg.entrez_query]
        else:
            cmd += ['-num_threads', str(cfg.threads)]
        env = None
        if not cfg.remote:
            # Lets blastn find taxdb.bt* next to the database for sscinames.
            db_dir = str(Path(cfg.db).resolve().parent)
            env = dict(os.environ)
            env['BLASTDB'] = os.pathsep.join(filter(None, [db_dir, env.get('BLASTDB')]))
        proc = subprocess.run(cmd, capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        raise RuntimeError(f'blastn failed ({proc.returncode}):\n{proc.stderr.strip()}')

    hits = []
    for line in proc.stdout.splitlines():
        if not line or line.startswith('#'):
            continue
        row = dict(zip(OUTFMT_FIELDS, line.split('\t')))
        for k in _INT_FIELDS:
            row[k] = int(row[k])
        for k in _FLOAT_FIELDS:
            row[k] = float(row[k])
        row['species'] = species_of(row)
        hits.append(row)
    return hits


_SPECIES_TAG = re.compile(r'\[species=([^\]]+)\]')
_TITLE_SPECIES = re.compile(r'([A-Z][a-z]+ (?:sp\.|[a-z][a-z-]+))')


def species_level(name: str) -> str:
    """Drop strain designations: 'Citrobacter koseri ATCC BAA-895' -> 'Citrobacter koseri'.

    Unnamed species ('Bacillus sp. XYZ', 'Firmicutes bacterium ABC') are kept whole
    since the strain is the only thing distinguishing them.
    """
    words = name.split()
    if words and words[0] == 'Candidatus':
        return ' '.join(words[:3])
    if len(words) >= 2 and words[1] in ('sp.', 'bacterium', 'archaeon'):
        return name
    return ' '.join(words[:2])


def species_of(hit: dict) -> str:
    sci = hit.get('sscinames', '')
    if sci and sci != 'N/A':
        return species_level(sci.split(';')[0])
    title = hit.get('stitle', '')
    m = _SPECIES_TAG.search(title) or _TITLE_SPECIES.search(title)
    return species_level(m.group(1)) if m else (title.split()[0] if title else hit['sacc'])


def matches_species(species: str, targets: Iterable[str]) -> bool:
    """Case-insensitive prefix match, so 'Escherichia coli' matches strains/serovars."""
    s = species.lower()
    return any(s.startswith(t.lower()) for t in targets)


def build_local_db(fasta_to_species: dict[str, str], out_prefix: str,
                   makeblastdb: str | None = None) -> str:
    """Concatenate genome FASTAs into a local nucleotide BLAST DB.

    Args:
        fasta_to_species: {path_to_fasta: 'Genus species'}; headers are rewritten
            as ``>{file_stem}_{n} [genome={file_stem}] [species=Genus species] header`` so the
            genome and species can be recovered from hits without a taxonomy database.
        out_prefix: DB path prefix to create.
    Returns:
        out_prefix (pass as BlastConfig(db=...)).
    """
    exe = makeblastdb or str(Path(find_blastn()).with_name('makeblastdb'))
    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    combined = f'{out_prefix}.fa'
    with open(combined, 'w') as out:
        for path, species in fasta_to_species.items():
            stem, n = Path(path).stem, 0
            with open(path) as fh:
                for line in fh:
                    if line.startswith('>'):
                        n += 1
                        out.write(f'>{stem}_{n} [genome={stem}] [species={species}] '
                                  f'{line[1:].strip()}\n')
                    else:
                        out.write(line)
    subprocess.run([exe, '-in', combined, '-dbtype', 'nucl', '-out', out_prefix],
                   check=True, capture_output=True)
    os.remove(combined)
    return out_prefix


# ----------------------------------------------------------------------------
# 1. Primer specificity
# ----------------------------------------------------------------------------

@dataclass
class BindingSite:
    primer_id: str
    subject: str
    species: str
    strand: str          # '+' primer anneals to minus strand and extends rightward
    five_prime: int      # subject coordinate of the primer's (projected) 5' end
    mismatches: int


@dataclass
class PredictedProduct:
    pair_label: str
    subject: str
    species: str
    start: int
    end: int
    size: int
    fwd_primer: str
    rev_primer: str
    total_mismatches: int
    on_target: bool


def _binding_sites(hits: list[dict], max_mismatches: int) -> list[BindingSite]:
    sites = []
    for h in hits:
        if h['qend'] != h['qlen']:          # 3' end must be in the alignment
            continue
        mm = h['mismatch'] + h['gapopen'] + (h['qstart'] - 1)
        if mm > max_mismatches:
            continue
        if h['sstart'] <= h['send']:
            strand, five = '+', h['sstart'] - (h['qstart'] - 1)
        else:
            strand, five = '-', h['sstart'] + (h['qstart'] - 1)
        sites.append(BindingSite(h['qseqid'], h['sacc'], h['species'], strand, five, mm))
    return sites


def predict_products(sites: list[BindingSite], primer_ids: set[str], pair_label: str,
                     max_product: int,
                     target_species: list[str]) -> list[PredictedProduct]:
    by_subject: dict[str, list[BindingSite]] = defaultdict(list)
    for s in sites:
        if s.primer_id in primer_ids:
            by_subject[s.subject].append(s)
    products = []
    for subject, ss in by_subject.items():
        plus = [s for s in ss if s.strand == '+']
        minus = [s for s in ss if s.strand == '-']
        for f in plus:
            for r in minus:
                size = r.five_prime - f.five_prime + 1
                if 0 < size <= max_product:
                    mm = f.mismatches + r.mismatches
                    # Without a target species nothing can be called on-target.
                    on = bool(target_species) and matches_species(f.species, target_species)
                    products.append(PredictedProduct(
                        pair_label, subject, f.species, f.five_prime, r.five_prime, size,
                        f.primer_id, r.primer_id, mm, on))
    return products


# ----------------------------------------------------------------------------
# 2. Amplicon uniqueness across species
# ----------------------------------------------------------------------------

@dataclass
class SpeciesHit:
    amplicon_id: str
    species: str
    n_subjects: int
    best_identity: float
    best_coverage: float
    is_target: bool


@dataclass
class AmpliconUniqueness:
    amplicon_id: str
    n_species: int                        # species with a hit covering >= min_coverage
    crossreactive_species: list[str]      # non-target, >= identity & coverage thresholds
    max_offtarget_identity: float         # identity * coverage of closest non-target
    closest_offtarget: str
    target_found: bool | None
    per_species: list[SpeciesHit] = field(default_factory=list)


def amplicon_uniqueness(hits: list[dict], amplicon_id: str, target_species: list[str],
                        min_identity: float, min_coverage: float) -> AmpliconUniqueness:
    agg: dict[str, dict] = {}
    for h in hits:
        if h['qseqid'] != amplicon_id:
            continue
        cov = 100.0 * (h['qend'] - h['qstart'] + 1) / h['qlen']
        a = agg.setdefault(h['species'], {'subjects': set(), 'ident': 0.0, 'cov': 0.0,
                                          'score': 0.0})
        a['subjects'].add(h['sacc'])
        # Pick the single best alignment by identity x coverage.
        score = h['pident'] * cov / 100.0
        if score > a['score']:
            a.update(ident=h['pident'], cov=cov, score=score)

    per_species = sorted(
        (SpeciesHit(amplicon_id, sp, len(a['subjects']), a['ident'], a['cov'],
                    matches_species(sp, target_species) if target_species else False)
         for sp, a in agg.items()),
        key=lambda s: -s.best_identity * s.best_coverage)

    off = [s for s in per_species if not s.is_target]
    cross = [s.species for s in off
             if s.best_identity >= min_identity and s.best_coverage >= min_coverage]
    closest = off[0] if off else None
    return AmpliconUniqueness(
        amplicon_id=amplicon_id,
        n_species=sum(s.best_coverage >= min_coverage for s in per_species),
        crossreactive_species=cross,
        max_offtarget_identity=(closest.best_identity * closest.best_coverage / 100.0
                                if closest else 0.0),
        closest_offtarget=closest.species if closest else '',
        target_found=any(s.is_target for s in per_species) if target_species else None,
        per_species=per_species,
    )


# ----------------------------------------------------------------------------
# Driver
# ----------------------------------------------------------------------------

@dataclass
class SpecificityResult:
    products: list[PredictedProduct]
    uniqueness: dict[str, AmpliconUniqueness]   # keyed by pair label
    primer_hits: list[dict]
    amplicon_hits: list[dict]
    target_species: list[str]
    primer_ids: dict[str, str]       # BLAST query id -> primer sequence
    max_mismatches: int = 3


@dataclass
class PrimerSpecificity:
    primer: str
    pairs: str                       # pair labels using this primer
    n_sites: int                     # binding sites (<= max_mismatches, 3' end aligned)
    n_species: int                   # species with >= 1 binding site
    n_perfect_species: int           # species with a 0-mismatch site
    n_offtarget_species: int
    n_offtarget_perfect_species: int
    offtarget_perfect_species: str
    target_bound: bool | None


def primer_specificity(res: 'SpecificityResult', pairs: list) -> list[PrimerSpecificity]:
    """Per-primer view: how many species/sites each primer can bind on its own."""
    sites = _binding_sites(res.primer_hits, res.max_mismatches)
    by_primer: dict[str, list[BindingSite]] = defaultdict(list)
    for s in sites:
        by_primer[s.primer_id].append(s)
    used_in: dict[str, list[str]] = defaultdict(list)
    for p in pairs:
        for seq in (p.left_seq, p.right_seq):
            used_in[seq].append(pair_label(p))
    t = res.target_species
    out = []
    for pid, seq in res.primer_ids.items():
        ss = by_primer.get(pid, [])
        species = {s.species for s in ss}
        perfect = {s.species for s in ss if s.mismatches == 0}
        off = {sp for sp in species if not matches_species(sp, t)} if t else species
        off_perfect = {sp for sp in perfect if not matches_species(sp, t)} if t else perfect
        out.append(PrimerSpecificity(
            primer=seq, pairs=';'.join(dict.fromkeys(used_in[seq])),
            n_sites=len(ss), n_species=len(species), n_perfect_species=len(perfect),
            n_offtarget_species=len(off), n_offtarget_perfect_species=len(off_perfect),
            offtarget_perfect_species=';'.join(sorted(off_perfect)),
            target_bound=any(matches_species(sp, t) for sp in species) if t else None,
        ))
    return out


def infer_target_species(hits: list[dict]) -> list[str]:
    """The one species with full-length perfect hits for every query, if unambiguous.

    Returns [] (with a warning) when several species match perfectly, e.g. for
    mobile genes such as resistance cassettes, since any choice would be arbitrary.
    """
    perfect: dict[str, set[str]] = defaultdict(set)
    for h in hits:
        if h['pident'] == 100.0 and h['qstart'] == 1 and h['qend'] == h['qlen']:
            perfect[h['qseqid']].add(h['species'])
    if not perfect:
        return []
    common = set.intersection(*perfect.values())
    if len(common) == 1:
        return list(common)
    if common:
        print(f'WARNING: {len(common)} species match all queries perfectly '
              f'({", ".join(sorted(common)[:5])}{" ..." if len(common) > 5 else ""}); '
              'cannot infer target, pass target_species. All species treated as non-target.',
              file=sys.stderr)
    return []


def pair_label(p) -> str:
    return f'{p.seq_id}_R{p.region}.{p.rank_in_region}'


def check_specificity(
    pairs: list,
    cfg: BlastConfig,
    target_species: str | list[str] | None = None,
    max_mismatches: int = 3,
    max_product: int = 1000,
    min_identity: float = 90.0,
    min_coverage: float = 80.0,
    check_primers: bool = True,
    check_amplicons: bool = True,
) -> SpecificityResult:
    """BLAST primers and amplicons of ``pairs`` and annotate them in place.

    Adds to each PrimerPair:
        n_products, n_offtarget_products, offtarget_product_species,
        amp_n_species, amp_crossreactive_species, amp_closest_offtarget,
        amp_closest_offtarget_idcov, amp_target_found, specific
    """
    if isinstance(target_species, str):
        target_species = [target_species]
    target_species = list(target_species or [])

    primer_q, amp_q, seen_primer, seen_amp = {}, {}, {}, {}
    for p in pairs:
        for seq in (p.left_seq, p.right_seq):
            seen_primer.setdefault(seq, f'P{len(seen_primer) + 1}')
        seen_amp.setdefault(p.amplicon, f'A{len(seen_amp) + 1}')
    primer_q = {v: k for k, v in seen_primer.items()}
    amp_q = {v: k for k, v in seen_amp.items()}

    primer_hits = amplicon_hits = []
    if check_primers:
        print(f'BLAST: {len(primer_q)} primers vs {cfg.db}'
              f'{" (remote)" if cfg.remote else ""} ...', file=sys.stderr)
        primer_hits = run_blast(primer_q, cfg, 'blastn-short', evalue=1000,
                                extra=['-word_size', '7'])
    if check_amplicons:
        print(f'BLAST: {len(amp_q)} amplicons vs {cfg.db}'
              f'{" (remote)" if cfg.remote else ""} ...', file=sys.stderr)
        amplicon_hits = run_blast(amp_q, cfg, 'blastn', evalue=1e-5)

    if not target_species:
        target_species = infer_target_species(primer_hits + amplicon_hits)
        if target_species:
            print(f'Target species not given; inferred from top hits: {target_species[0]}',
                  file=sys.stderr)

    sites = _binding_sites(primer_hits, max_mismatches)
    products, uniq = [], {}
    for p in pairs:
        label = pair_label(p)
        if check_primers:
            ids = {seen_primer[p.left_seq], seen_primer[p.right_seq]}
            prods = predict_products(sites, ids, label, max_product,
                                     target_species)
            products.extend(prods)
            off = [x for x in prods if not x.on_target]
            p.n_products = len(prods)
            p.n_offtarget_products = len(off)
            p.offtarget_product_species = ';'.join(sorted({x.species for x in off}))
        if check_amplicons:
            u = amplicon_uniqueness(amplicon_hits, seen_amp[p.amplicon], target_species,
                                    min_identity, min_coverage)
            u.amplicon_id = label
            for s in u.per_species:
                s.amplicon_id = label
            uniq[label] = u
            p.amp_n_species = u.n_species
            p.amp_crossreactive_species = ';'.join(u.crossreactive_species)
            p.amp_closest_offtarget = u.closest_offtarget
            p.amp_closest_offtarget_idcov = u.max_offtarget_identity
            p.amp_target_found = u.target_found
        # None = undetermined: off-target cannot be told apart from target.
        p.specific = (None if not target_species else
                      getattr(p, 'n_offtarget_products', 0) == 0
                      and not getattr(p, 'amp_crossreactive_species', ''))
    return SpecificityResult(products, uniq, primer_hits, amplicon_hits, target_species,
                             primer_q, max_mismatches)


SPECIFICITY_COLUMNS = [
    'specific', 'n_products', 'n_offtarget_products', 'offtarget_product_species',
    'amp_n_species', 'amp_crossreactive_species', 'amp_closest_offtarget',
    'amp_closest_offtarget_idcov', 'amp_target_found',
]


def write_products(products: list[PredictedProduct], path: str) -> None:
    import csv
    cols = list(PredictedProduct.__dataclass_fields__)
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(cols)
        for x in products:
            w.writerow([getattr(x, c) for c in cols])


def write_species_table(uniq: dict[str, AmpliconUniqueness], path: str) -> None:
    import csv
    cols = list(SpeciesHit.__dataclass_fields__)
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(cols)
        for u in uniq.values():
            for s in u.per_species:
                w.writerow([f'{v:.2f}' if isinstance(v, float) else v
                            for v in (getattr(s, c) for c in cols)])


def write_primer_table(rows: list[PrimerSpecificity], path: str) -> None:
    import csv
    cols = list(PrimerSpecificity.__dataclass_fields__)
    with open(path, 'w', newline='') as fh:
        w = csv.writer(fh, delimiter='\t')
        w.writerow(cols)
        for r in rows:
            w.writerow([getattr(r, c) for c in cols])


def summarize_specificity(pairs: list, res: SpecificityResult, top_species: int = 5) -> str:
    lines = []
    for p in pairs:
        label = pair_label(p)
        status = {True: 'SPECIFIC', False: 'NOT SPECIFIC',
                  None: 'UNKNOWN (no target species)'}[p.specific]
        if not getattr(p, 'n_products', 0) and not getattr(p, 'amp_n_species', 0):
            status += ' (no hits at all: target not in database?)'
        parts = [f'  {label:24s} {status}']
        if hasattr(p, 'n_products'):
            if res.target_species:
                parts.append(f'products {p.n_products} ({p.n_offtarget_products} off-target'
                             + (f': {p.offtarget_product_species}'
                                if p.offtarget_product_species else '') + ')')
            else:
                parts.append(f'products {p.n_products} in '
                             f'{len({x.species for x in res.products if x.pair_label == label})}'
                             ' species')
        u = res.uniqueness.get(label)
        if u:
            parts.append(f'amplicon in {u.n_species} species')
            if u.target_found is not None:
                parts.append('target found' if u.target_found else 'TARGET NOT FOUND')
            if u.closest_offtarget:
                parts.append(f'closest non-target {u.closest_offtarget} '
                             f'({u.max_offtarget_identity:.0f}% id x cov)')
            if u.crossreactive_species:
                parts.append(f'cross-reactive: {", ".join(u.crossreactive_species[:top_species])}'
                             + (' ...' if len(u.crossreactive_species) > top_species else ''))
        lines.append(' | '.join(parts))
    return '\n'.join(lines)


# ----------------------------------------------------------------------------
# Single-genome scan: where else can a primer pair bind?
# ----------------------------------------------------------------------------

def blast_db_from_fasta(fasta: str, out_prefix: str | None = None,
                        makeblastdb: str | None = None, force: bool = False) -> str:
    """Build (or reuse) a BLAST database from one FASTA, e.g. a single genome.

    Args:
        fasta: genome/sequence FASTA; its own headers are kept (``-parse_seqids``),
            so hits report the original contig names.
        out_prefix: DB path prefix; default ``<fasta without extension>_blastdb``.
        force: rebuild even if the DB is newer than the FASTA.
    Returns:
        the DB prefix, to pass as ``BlastConfig(db=...)``.
    """
    src = Path(fasta)
    if not src.exists():
        raise FileNotFoundError(fasta)
    prefix = Path(out_prefix) if out_prefix else src.with_suffix('').with_name(
        src.with_suffix('').name + '_blastdb')
    built = next((p for p in (prefix.with_suffix('.nsq'), prefix.with_suffix('.nal'))
                  if p.exists()), None)
    if built and not force and built.stat().st_mtime >= src.stat().st_mtime:
        return str(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    exe = makeblastdb or str(Path(find_blastn()).with_name('makeblastdb'))
    proc = subprocess.run([exe, '-in', str(src), '-dbtype', 'nucl', '-parse_seqids',
                           '-out', str(prefix)], capture_output=True, text=True)
    if proc.returncode != 0:
        # -parse_seqids rejects duplicate/malformed ids; fall back to generic ids.
        proc = subprocess.run([exe, '-in', str(src), '-dbtype', 'nucl', '-out', str(prefix)],
                              capture_output=True, text=True)
        if proc.returncode != 0:
            raise RuntimeError(f'makeblastdb failed:\n{proc.stderr.strip()}')
    print(f'Built BLAST database {prefix}', file=sys.stderr)
    return str(prefix)


@dataclass
class GenomeSite:
    """One place in the genome where a primer can prime (3' end aligned)."""
    primer: str
    role: str            # 'forward' or 'reverse'
    contig: str
    position: int        # 1-based genome coordinate of the primer's 5' end
    strand: str          # '+' extends rightward, '-' leftward
    mismatches: int      # mismatches + gaps + unaligned 5' bases
    identity: float      # over the full primer length
    aligned: int         # matching primer bases


def _fetch_windows(cfg: BlastConfig, windows: list[tuple[str, int, int, str]]) -> list[str] | None:
    """blastdbcmd the given (contig, start, stop, strand) 1-based windows, in order."""
    exe = str(Path(cfg.blastn).with_name('blastdbcmd'))
    spec = '\n'.join(f'{c} {a}-{b} {st}' for c, a, b, st in windows) + '\n'
    with tempfile.TemporaryDirectory(prefix='hda_fetch_') as tmp:
        batch = Path(tmp) / 'batch.txt'
        batch.write_text(spec)
        proc = subprocess.run([exe, '-db', cfg.db, '-entry_batch', str(batch),
                               '-outfmt', '%s'], capture_output=True, text=True,
                              env={**os.environ,
                                   'BLASTDB': str(Path(cfg.db).resolve().parent)})
    if proc.returncode != 0:
        return None
    seqs = [l.strip().upper() for l in proc.stdout.splitlines() if l.strip()]
    return seqs if len(seqs) == len(windows) else None


def _exact_sites(hits: list[dict], queries: dict[str, str], cfg: BlastConfig,
                 max_mismatches: int, require_3p_match: int = 1
                 ) -> list[tuple[dict, int, int, str, str]]:
    """Re-count mismatches over the primer's full footprint.

    BLAST may report only part of a primer as aligned; counting the unaligned
    bases as mismatches would hide real binding sites (a 5' mismatch still primes).
    Every hit is therefore extended along its diagonal to the full primer length
    and compared base by base with the genome.

    A site is kept when the last ``require_3p_match`` bases of the primer match
    exactly (a 3'-terminal mismatch blocks extension, so such sites are not
    amplifiable) and the whole primer has <= ``max_mismatches`` mismatches.

    Returns [(hit, five_prime_pos, mismatches, strand, genomic_seq)], deduplicated
    per (primer, contig, 3' position).
    """
    cands: dict[tuple, dict] = {}
    for h in hits:
        # Project the hit along its diagonal to the primer's 3' end.
        overhang = h['qlen'] - h['qend']
        three = h['send'] + overhang if h['sstart'] <= h['send'] else h['send'] - overhang
        key = (h['qseqid'], h['sacc'], three)
        if key not in cands or h['bitscore'] > cands[key]['bitscore']:
            cands[key] = dict(h, _three=three)
    if not cands:
        return []

    ordered = list(cands.values())
    windows, meta = [], []
    for h in ordered:
        qlen = h['qlen']
        plus = h['sstart'] <= h['send']
        three = h['_three']
        five = three - (qlen - 1) if plus else three + (qlen - 1)
        lo, hi = (five, three) if plus else (three, five)
        clip_lo, clip_hi = max(1, lo), min(h['slen'], hi)
        windows.append((h['sacc'], clip_lo, clip_hi, 'plus' if plus else 'minus'))
        meta.append((h, five, plus, lo, hi, clip_lo, clip_hi))

    seqs = _fetch_windows(cfg, windows)
    out = []
    if seqs is None:                            # e.g. remote DB: conservative fallback
        for h in ordered:
            if h['qend'] != h['qlen']:          # cannot verify the 3' end without the sequence
                continue
            mm = h['mismatch'] + h['gapopen'] + (h['qstart'] - 1)
            if mm <= max_mismatches:
                plus = h['sstart'] <= h['send']
                five = h['send'] - (h['qlen'] - 1) if plus else h['send'] + (h['qlen'] - 1)
                out.append((h, five, mm, '+' if plus else '-', ''))
        return out

    for (h, five, plus, lo, hi, clip_lo, clip_hi), gseq in zip(meta, seqs):
        primer = queries[h['qseqid']]
        # Pad where the footprint ran off the contig; missing bases cannot pair.
        pad_5 = clip_lo - lo if plus else hi - clip_hi
        pad_3 = hi - clip_hi if plus else clip_lo - lo
        aligned = '-' * pad_5 + gseq + '-' * pad_3
        mm = sum(1 for a, b in zip(primer, aligned) if a != b) + abs(len(primer) - len(aligned))
        tail = require_3p_match
        if tail and (len(aligned) < len(primer) or primer[-tail:] != aligned[-tail:]):
            continue                            # 3'-terminal mismatch: will not extend
        if mm <= max_mismatches:
            out.append((h, five, mm, '+' if plus else '-', gseq))
    return out


def scan_genome(pair, cfg: BlastConfig, max_mismatches: int = 3, max_product: int = 3000,
                require_3p_match: int = 1
                ) -> tuple[list[GenomeSite], list[PredictedProduct]]:
    """Find every binding site of a pair's two primers in one genome, and the
    products they could make.

    A binding site is a place where the last ``require_3p_match`` bases of the
    primer match exactly (a 3'-terminal mismatch blocks extension) and the
    full-length primer footprint has <= ``max_mismatches`` mismatches. Two facing
    sites within ``max_product`` bp are a predicted product; the intended one is
    flagged via ``PredictedProduct.on_target``.

    Returns:
        (sites sorted by mismatches then position, products, intended first).
    """
    queries = {'F': pair.left_seq, 'R': pair.right_seq}
    hits = run_blast(queries, cfg, 'blastn-short', evalue=1000, extra=['-word_size', '7'])
    exact = _exact_sites(hits, queries, cfg, max_mismatches, require_3p_match)

    role = {'F': 'forward', 'R': 'reverse'}
    sites, bsites = [], []
    for h, five, mm, strand, gseq in exact:
        primer = queries[h['qseqid']]
        sites.append(GenomeSite(
            primer=primer, role=role[h['qseqid']], contig=h['sacc'], position=five,
            strand=strand, mismatches=mm,
            identity=100.0 * (len(primer) - mm) / len(primer), aligned=len(primer) - mm))
        bsites.append(BindingSite(h['qseqid'], h['sacc'], h['sacc'], strand, five, mm))
    sites.sort(key=lambda s: (s.mismatches, s.contig, s.position))

    products = predict_products(bsites, {'F', 'R'}, pair_label(pair), max_product, [])
    for x in products:
        x.on_target = x.total_mismatches == 0 and x.size == pair.product_size
    products.sort(key=lambda x: (not x.on_target, x.total_mismatches, x.size))
    return sites, products


def summarize_genome_scan(sites: list[GenomeSite], products: list[PredictedProduct],
                          max_mismatches: int = 3) -> str:
    lines = []
    for role in ('forward', 'reverse'):
        rs = [s for s in sites if s.role == role]
        perfect = [s for s in rs if s.mismatches == 0]
        lines.append(f'{role.capitalize():8s} {rs[0].primer if rs else "":30s} '
                     f'{len(rs)} binding site(s) with <= {max_mismatches} mismatches, '
                     f'{len(perfect)} perfect')
        for s in rs[:10]:
            lines.append(f'    {s.contig}:{s.position} ({s.strand}) '
                         f'{s.mismatches} mm, {s.aligned}/{len(s.primer)} nt aligned')
        if len(rs) > 10:
            lines.append(f'    ... {len(rs) - 10} more')
    intended = [x for x in products if x.on_target]
    extra = [x for x in products if not x.on_target]
    lines.append(f'Products: {len(intended)} intended, {len(extra)} additional')
    for x in products[:15]:
        tag = 'INTENDED' if x.on_target else 'extra   '
        lines.append(f'    {tag} {x.species}:{min(x.start, x.end)}-{max(x.start, x.end)} '
                     f'{x.size} bp, {x.fwd_primer}/{x.rev_primer}, {x.total_mismatches} mm')
    if len(products) > 15:
        lines.append(f'    ... {len(products) - 15} more')
    return '\n'.join(lines)
