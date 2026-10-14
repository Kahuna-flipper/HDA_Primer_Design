"""Strain-specific tHDA amplicons: design primers only where the input is unique.

Workflow (``design_strain_specific``):

1. Prevalence scan. The input sequence(s) are split into chunks and BLASTed
   (megablast) against one or more background databases, e.g. a panel of other
   strains of the same species plus ``ref_prok_rep_genomes`` for other species.
   For every stretch of the input, the number of background *genomes* carrying it
   (an HSP with >= ``min_identity`` %) is computed.
2. Unique regions. Stretches carried by <= ``max_background_genomes`` genomes and
   at least the minimum product size long are kept.
3. Primer design. NEB tHDA primers (hda_primer_design) are designed inside each
   unique region only.
4. Validation. The best candidates' primers and amplicons are BLASTed against
   the same backgrounds; pairs are ranked by the number of background genomes
   with a predicted product, then with a similar amplicon, then by design score.

Background genomes are identified by ``[genome=...]`` tags (written by
``hda_blast.build_local_db``), else by the ``{stem}_{n}`` sequence ids it writes,
else by species (adequate for one-genome-per-species DBs such as
ref_prok_rep_genomes).

If the input genome itself is in a background DB, it would make everything look
non-unique: genomes covering >= 99 % of the input at 100 % identity are excluded
automatically (with a warning), and others can be excluded via ``exclude``.
"""

from __future__ import annotations

import re
import sys
from collections import defaultdict
from dataclasses import dataclass, field

from hda_blast import (BlastConfig, _binding_sites, predict_products, run_blast,
                       species_of)
from hda_primer_design import (HDAParams, PrimerPair, clean_sequence, design_hda_primers,
                               gc_content)

_GENOME_TAG = re.compile(r'\[genome=([^\]]+)\]')
_STEM_ID = re.compile(r'^(\S+)_\d+$')


def genome_of(hit: dict) -> str:
    title = hit.get('stitle', '')
    m = _GENOME_TAG.search(title)
    if m:
        return m.group(1)
    first = title.split()[0] if title else ''
    m = _STEM_ID.match(first)
    if m and '[species=' in title:
        return m.group(1)
    return species_of(hit)


# ----------------------------------------------------------------------------
# 1-2. Prevalence scan and unique regions
# ----------------------------------------------------------------------------

def _merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge half-open intervals."""
    out: list[list[int]] = []
    for s, e in sorted(intervals):
        if out and s <= out[-1][1]:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(s, e) for s, e in out]


@dataclass
class Prevalence:
    """Background carriage along one input record."""
    record: str
    length: int
    segments: list[tuple[int, int, int]]          # (start, end, n_genomes), half-open
    genome_intervals: dict[str, list[tuple[int, int]]]

    def count_at(self, start: int, end: int) -> int:
        """Max number of background genomes carrying any base of [start, end)."""
        return max((n for s, e, n in self.segments if s < end and e > start), default=0)

    def genomes_overlapping(self, start: int, end: int) -> list[str]:
        return sorted(g for g, iv in self.genome_intervals.items()
                      if any(s < end and e > start for s, e in iv))


@dataclass
class UniqueRegion:
    record: str
    start: int          # 0-based, half-open
    end: int
    max_background: int

    @property
    def length(self) -> int:
        return self.end - self.start


def scan_prevalence(
    records: list[tuple[str, str]],
    backgrounds: list[BlastConfig],
    min_identity: float = 90.0,
    min_hsp_len: int = 50,
    exclude: list[str] | None = None,
    auto_exclude_self: bool = True,
    chunk: int = 20000,
    overlap: int = 200,
    task: str = 'megablast',
) -> tuple[dict[str, Prevalence], set[str]]:
    """BLAST chunks of every record against each background and sweep coverage.

    Returns ({record: Prevalence}, excluded genome names incl. auto-detected self).
    """
    exclude = set(exclude or [])
    queries, offsets = {}, {}
    for rid, seq in records:
        for s in range(0, len(seq), chunk - overlap):
            qid = f'q{len(queries)}'
            queries[qid] = seq[s:s + chunk]
            offsets[qid] = (rid, s)
            if s + chunk >= len(seq):
                break

    # genome -> record -> intervals; genome keys are namespaced per DB when >1 DB.
    raw: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    exact: dict[str, dict[str, list]] = defaultdict(lambda: defaultdict(list))
    for cfg in backgrounds:
        print(f'Prevalence scan: {len(queries)} chunks vs {cfg.db} ({task}) ...',
              file=sys.stderr)
        hits = run_blast(queries, cfg, task, evalue=1e-10,
                         extra=['-perc_identity', str(min_identity)])
        prefix = f'{cfg.db.rsplit("/", 1)[-1]}:' if len(backgrounds) > 1 else ''
        for h in hits:
            if h['pident'] < min_identity or h['length'] < min_hsp_len:
                continue
            rid, off = offsets[h['qseqid']]
            g = prefix + genome_of(h)
            iv = (off + h['qstart'] - 1, off + h['qend'])
            raw[g][rid].append(iv)
            if h['pident'] == 100.0:
                exact[g][rid].append(iv)

    total = sum(len(s) for _, s in records)
    if auto_exclude_self:
        for g, per_rec in exact.items():
            covered = sum(e - s for ivs in per_rec.values() for s, e in _merge(ivs))
            if covered >= 0.99 * total and g not in exclude:
                print(f'NOTE: excluding background genome {g!r}: identical to the input '
                      f'over {100 * covered / total:.1f}% (the input itself?)', file=sys.stderr)
                exclude.add(g)

    out = {}
    for rid, seq in records:
        genome_iv = {g: _merge(per_rec[rid]) for g, per_rec in raw.items()
                     if g not in exclude and per_rec.get(rid)}
        events = defaultdict(int)
        for ivs in genome_iv.values():
            for s, e in ivs:
                events[s] += 1
                events[e] -= 1
        segments, depth, prev = [], 0, 0
        for pos in sorted(events) + [len(seq)]:
            if pos > prev:
                segments.append((prev, min(pos, len(seq)), depth))
            depth += events.get(pos, 0)
            prev = max(prev, pos)
        out[rid] = Prevalence(rid, len(seq), segments, genome_iv)
    return out, exclude


def unique_regions(prev: Prevalence, max_background_genomes: int = 0,
                   min_len: int = 80) -> list[UniqueRegion]:
    runs, cur = [], None
    for s, e, n in prev.segments:
        if n <= max_background_genomes:
            if cur and cur[1] == s:
                cur = [cur[0], e, max(cur[2], n)]
            else:
                if cur:
                    runs.append(cur)
                cur = [s, e, n]
        elif cur:
            runs.append(cur)
            cur = None
    if cur:
        runs.append(cur)
    return [UniqueRegion(prev.record, s, e, n) for s, e, n in runs if e - s >= min_len]


# ----------------------------------------------------------------------------
# 3-4. Design inside unique regions and validate
# ----------------------------------------------------------------------------

@dataclass
class StrainSpecificPair:
    pair: PrimerPair
    record: str
    unique_region: UniqueRegion
    n_bg_genomes_product: int | None = None      # background genomes with a predicted product
    n_bg_genomes_amplicon: int | None = None     # ... carrying the amplicon (>= id & cov)
    bg_product_genomes: list[str] = field(default_factory=list)
    bg_amplicon_genomes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        d = self.pair.as_dict()
        d.update(record=self.record,
                 unique_region=f'{self.unique_region.start + 1}-{self.unique_region.end}',
                 unique_region_len=self.unique_region.length,
                 n_bg_genomes_product=self.n_bg_genomes_product,
                 n_bg_genomes_amplicon=self.n_bg_genomes_amplicon,
                 bg_product_genomes=';'.join(self.bg_product_genomes),
                 bg_amplicon_genomes=';'.join(self.bg_amplicon_genomes))
        return d


def _design_in_region(seq: str, region: UniqueRegion, params: HDAParams,
                      n_regions: int, pairs_per_region: int) -> list[PrimerPair]:
    sub = seq[region.start:region.end]
    try:
        pairs, _ = design_hda_primers(sub, seq_id=region.record, params=params,
                                      n_regions=n_regions, pairs_per_region=pairs_per_region,
                                      fast=True)
    except ValueError:
        return []
    for p in pairs:                      # shift back to record coordinates
        p.left_start += region.start
        p.right_end += region.start
    return pairs


def validate_pairs(cands: list[StrainSpecificPair], backgrounds: list[BlastConfig],
                   exclude: set[str], max_mismatches: int = 3, max_product: int = 1000,
                   min_identity: float = 90.0, min_coverage: float = 80.0) -> None:
    """BLAST primers/amplicons; fill the n_bg_* fields in place."""
    primers: dict[str, str] = {}
    amps: dict[str, str] = {}
    for c in cands:
        for s in (c.pair.left_seq, c.pair.right_seq):
            primers.setdefault(s, f'P{len(primers) + 1}')
        amps.setdefault(c.pair.amplicon, f'A{len(amps) + 1}')
    pq = {v: k for k, v in primers.items()}
    aq = {v: k for k, v in amps.items()}

    prod_genomes: dict[int, set] = defaultdict(set)
    amp_genomes: dict[int, set] = defaultdict(set)
    for cfg in backgrounds:
        prefix = f'{cfg.db.rsplit("/", 1)[-1]}:' if len(backgrounds) > 1 else ''
        print(f'Validation: {len(pq)} primers + {len(aq)} amplicons vs {cfg.db} ...',
              file=sys.stderr)
        p_hits = run_blast(pq, cfg, 'blastn-short', evalue=1000, extra=['-word_size', '7'])
        a_hits = run_blast(aq, cfg, 'blastn', evalue=1e-5)
        # Products are predicted per genome: reuse hda_blast by keying on genome.
        for h in p_hits + a_hits:
            h['species'] = prefix + genome_of(h)
        sites = _binding_sites(p_hits, max_mismatches)
        for i, c in enumerate(cands):
            ids = {primers[c.pair.left_seq], primers[c.pair.right_seq]}
            for x in predict_products(sites, ids, '', max_product, []):
                if x.species not in exclude:
                    prod_genomes[i].add(x.species)
            aid = amps[c.pair.amplicon]
            for h in a_hits:
                if h['qseqid'] != aid or h['species'] in exclude:
                    continue
                cov = 100.0 * (h['qend'] - h['qstart'] + 1) / h['qlen']
                if h['pident'] >= min_identity and cov >= min_coverage:
                    amp_genomes[i].add(h['species'])
    for i, c in enumerate(cands):
        c.bg_product_genomes = sorted(prod_genomes[i])
        c.bg_amplicon_genomes = sorted(amp_genomes[i])
        c.n_bg_genomes_product = len(c.bg_product_genomes)
        c.n_bg_genomes_amplicon = len(c.bg_amplicon_genomes)


def design_strain_specific(
    sequence: str | list[tuple[str, str]],
    backgrounds: BlastConfig | list[BlastConfig],
    seq_id: str = 'target',
    params: HDAParams | None = None,
    exclude: list[str] | None = None,
    max_background_genomes: int = 0,
    min_identity: float = 90.0,
    n_candidates: int = 30,
    pairs_per_region: int = 2,
    n_report: int = 10,
    validate: bool = True,
    max_mismatches: int = 3,
    auto_exclude_self: bool = True,
) -> tuple[list[StrainSpecificPair], dict]:
    """Design NEB-compliant tHDA primers in regions unique to the input.

    Args:
        sequence: a sequence, or [(record_id, sequence), ...] (e.g. genome contigs).
        backgrounds: BLAST DB(s) of genomes to be distinguished from, e.g. a panel of
            other strains (build with hda_blast.build_local_db) and/or
            ref_prok_rep_genomes.
        exclude: background genome names to ignore (e.g. the input strain itself).
        max_background_genomes: a region is "unique" if at most this many background
            genomes carry it (0 = strictly unique; >0 = "mostly unique").
        min_identity: % identity for a background genome to count as carrying a region.
        n_candidates: candidate pairs (best design score) to BLAST-validate.
        pairs_per_region: primer sets kept per unique region before validation.
        n_report: pairs returned after validation ranking.

    Returns:
        (ranked StrainSpecificPairs, info dict with unique regions and prevalence).
    """
    records = ([(seq_id, clean_sequence(sequence))] if isinstance(sequence, str)
               else [(r, clean_sequence(s)) for r, s in sequence])
    backgrounds = [backgrounds] if isinstance(backgrounds, BlastConfig) else list(backgrounds)
    exclude_set = set(exclude or [])
    if params is None:
        params = HDAParams.for_target(gc_content(''.join(s for _, s in records)))

    prev, exclude_set = scan_prevalence(records, backgrounds, min_identity=min_identity,
                                        exclude=list(exclude_set),
                                        auto_exclude_self=auto_exclude_self)
    seen = {g for p in prev.values() for g in p.genome_intervals}
    regions = []
    for rid, _ in records:
        regions += unique_regions(prev[rid], max_background_genomes, params.product_size[0])
    print(f'{len(regions)} unique regions ({sum(r.length for r in regions):,} bp) '
          f'with <= {max_background_genomes} background genomes; '
          f'{len(seen)} background genomes share some sequence', file=sys.stderr)

    seqs = dict(records)
    cands: list[StrainSpecificPair] = []
    for reg in regions:
        n_sub = max(1, min(5, reg.length // 150))
        for p in _design_in_region(seqs[reg.record], reg, params, n_sub, pairs_per_region):
            cands.append(StrainSpecificPair(p, reg.record, reg))
    cands.sort(key=lambda c: c.pair.score)
    print(f'{len(cands)} candidate primer pairs in unique regions', file=sys.stderr)

    # Spread validation over regions: best pair of each region first, then the rest.
    firsts = {}
    for c in cands:
        firsts.setdefault((c.record, c.unique_region.start, c.pair.region), c)
    ordered = list(firsts.values()) + [c for c in cands if c not in firsts.values()]
    chosen = ordered[:n_candidates]

    if validate and chosen:
        validate_pairs(chosen, backgrounds, exclude_set, max_mismatches=max_mismatches,
                       min_identity=min_identity)
        chosen.sort(key=lambda c: (c.n_bg_genomes_product, c.n_bg_genomes_amplicon,
                                   c.pair.score))
    else:
        chosen.sort(key=lambda c: c.pair.score)

    # Report distinct unique regions first (NEB: test several amplicon regions),
    # keeping the validation ranking within each pass.
    by_region: dict[tuple, list[StrainSpecificPair]] = defaultdict(list)
    for c in chosen:
        by_region[(c.record, c.unique_region.start)].append(c)
    chosen = [c for k in range(max((len(v) for v in by_region.values()), default=0))
              for c in chosen if by_region[(c.record, c.unique_region.start)].index(c) == k]
    for rank, c in enumerate(chosen, 1):
        c.pair.region, c.pair.rank_in_region = rank, 1
    info = dict(prevalence=prev, unique_regions=regions, params=params, excluded=exclude_set,
                n_candidates=len(cands), background_genomes_sharing=len(seen))
    return chosen[:n_report], info


def summarize_strain_specific(results: list[StrainSpecificPair], info: dict) -> str:
    regions = info['unique_regions']
    lines = [f"# {len(regions)} unique regions, {sum(r.length for r in regions):,} bp; "
             f"{info['n_candidates']} candidate pairs"]
    for i, c in enumerate(results, 1):
        p = c.pair
        val = ('' if c.n_bg_genomes_product is None else
               f'  background: product in {c.n_bg_genomes_product}, '
               f'amplicon in {c.n_bg_genomes_amplicon} genomes')
        lines.append(
            f'  #{i} {c.record}:{p.left_start + 1}-{p.right_end + 1} ({p.product_size} bp, '
            f'Tm {p.product_tm:.1f}, GC {p.amplicon_gc:.1f}%) in unique region '
            f'{c.unique_region.start + 1}-{c.unique_region.end} '
            f'({c.unique_region.length} bp){val}\n'
            f'      F {p.left_seq}  ({p.left_len} nt, Tm {p.left_tm:.1f}, GC {p.left_gc:.0f}%)\n'
            f'      R {p.right_seq}  ({p.right_len} nt, Tm {p.right_tm:.1f}, GC {p.right_gc:.0f}%)')
        if c.bg_product_genomes:
            lines.append(f'      products in: {", ".join(c.bg_product_genomes[:8])}'
                         + (' ...' if len(c.bg_product_genomes) > 8 else ''))
    return '\n'.join(lines)


STRAIN_COLUMNS = [
    'record', 'region', 'n_bg_genomes_product', 'n_bg_genomes_amplicon',
    'unique_region', 'unique_region_len', 'amplicon_start_1based', 'amplicon_end_1based',
    'product_size', 'product_tm', 'amplicon_gc', 'primer_mean_gc', 'score',
    'left_seq', 'left_len', 'left_tm', 'left_gc', 'right_seq', 'right_len', 'right_tm',
    'right_gc', 'left_hairpin_th', 'right_hairpin_th', 'pair_compl_any_th',
    'bg_product_genomes', 'bg_amplicon_genomes', 'amplicon',
]


def write_strain_table(results: list[StrainSpecificPair], path: str) -> None:
    import csv
    with open(path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=STRAIN_COLUMNS, delimiter='\t',
                           extrasaction='ignore', restval='')
        w.writeheader()
        for c in results:
            w.writerow({k: f'{v:.2f}' if isinstance(v, float) else v
                        for k, v in c.as_dict().items()})


def write_unique_regions(info: dict, path: str) -> None:
    with open(path, 'w') as fh:
        fh.write('record\tstart_1based\tend\tlength\tmax_background_genomes\n')
        for r in info['unique_regions']:
            fh.write(f'{r.record}\t{r.start + 1}\t{r.end}\t{r.length}\t{r.max_background}\n')


def main(argv: list[str] | None = None) -> int:
    import argparse
    from hda_primer_design import read_fasta
    ap = argparse.ArgumentParser(
        description='Design tHDA primers in regions unique to the input (strain-specific).')
    ap.add_argument('fasta', help='target sequence / genome (FASTA, multi-contig ok)')
    ap.add_argument('-b', '--background', action='append', required=True,
                    help='BLAST DB of genomes to distinguish from (repeatable)')
    ap.add_argument('-o', '--output', help='output TSV (also writes <stem>.unique_regions.tsv)')
    ap.add_argument('--exclude', action='append', default=[],
                    help='background genome to ignore, e.g. the input strain (repeatable)')
    ap.add_argument('--max-background-genomes', type=int, default=0,
                    help='region counts as unique if carried by <= N background genomes')
    ap.add_argument('--min-identity', type=float, default=90.0)
    ap.add_argument('--n-candidates', type=int, default=30,
                    help='candidate pairs to BLAST-validate (default 30)')
    ap.add_argument('--n-report', type=int, default=10)
    ap.add_argument('--threads', type=int, default=16)
    ap.add_argument('--tm-model', default='neb_legacy')
    args = ap.parse_args(argv)

    records = read_fasta(args.fasta)
    gc = gc_content(''.join(s for _, s in records))
    params = HDAParams.for_target(gc, tm_model=args.tm_model)
    bgs = [BlastConfig(db=db, threads=args.threads) for db in args.background]
    results, info = design_strain_specific(
        records, bgs, params=params, exclude=args.exclude,
        max_background_genomes=args.max_background_genomes, min_identity=args.min_identity,
        n_candidates=args.n_candidates, n_report=args.n_report)
    print(summarize_strain_specific(results, info), file=sys.stderr)
    print('\nNote: order desalted (salt-free) primers; tHDA is sensitive to salt.',
          file=sys.stderr)
    if args.output:
        write_strain_table(results, args.output)
        write_unique_regions(info, args.output.rsplit('.', 1)[0] + '.unique_regions.tsv')
    else:
        write_strain_table(results, '/dev/stdout')
    return 0


if __name__ == '__main__':
    sys.exit(main())
