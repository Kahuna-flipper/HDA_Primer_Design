"""Primer design for thermophilic helicase-dependent amplification (tHDA).

Wraps primer3 with the NEB IsoAmp / tHDA primer design recommendations:

    Product size : 80 - 120 bp (opt ~100)
    Product Tm   : min 68, opt 71, max 75 C
    Primer size  : min 24, opt 27, max 33 nt
    Primer Tm    : min 60, opt 68, max 74 C
    Primer GC%   : min 35, opt 44, max 60

plus the NEB footnotes / "other considerations":

    * opt primer size -> 26 nt when target GC > 45 %
    * opt primer Tm   -> 65 C when target GC < 37.5 %, 71 C when target GC > 45 %
    * amplicons with GC ~40 % are preferred                         (consideration 1)
    * primer sets whose mean GC is closest to amplicon GC preferred (consideration 2)
    * test several amplicon regions, several primer sets per region (consideration 3)
    * walk serial longer/shorter primers around a working set       (consideration 4)
    * order desalted / salt-free primers                            (consideration 5)

Tm model
--------
NEB's numbers were calibrated with the original Primer3 (e.g. "a ~100 bp amplicon
with ~40 % GC has a product Tm of ~71 C"). Modern primer3 defaults (SantaLucia
parameters, 1.5 mM Mg2+, 0.6 mM dNTP) give Tm values several degrees higher, which
would reject most good tHDA amplicons. The default ``tm_model='neb_legacy'``
therefore reproduces the original Primer3 settings (Breslauer table, Schildkraut
salt correction, 50 mM Na+, no Mg2+/dNTP). Use ``tm_model='primer3_default'``
to keep modern primer3 thermodynamics.

Usage
-----
    python hda_primer_design.py input.fasta -o primers.tsv
    python hda_primer_design.py --sequence ACGT... --n-regions 5 --pairs-per-region 3

    # + primer specificity and amplicon uniqueness across species (see hda_blast.py)
    python hda_primer_design.py input.fasta -o primers.tsv --blast-db /path/to/db \\
        --target-species "Escherichia coli"
    python hda_primer_design.py input.fasta -o primers.tsv \\
        --blast-remote --entrez-query "Bacteria[organism]"
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from dataclasses import dataclass, field, replace
from functools import lru_cache
from typing import Iterable

import primer3

COMPLEMENT = str.maketrans('ACGTNacgtn', 'TGCANtgcan')

TM_MODELS = {
    # Settings of the original Primer3 used to derive the NEB recommendations.
    'neb_legacy': dict(
        PRIMER_TM_FORMULA=0,          # Breslauer et al. 1986
        PRIMER_SALT_CORRECTIONS=0,    # Schildkraut & Lifson 1965
        PRIMER_SALT_MONOVALENT=50.0,
        PRIMER_SALT_DIVALENT=0.0,
        PRIMER_DNTP_CONC=0.0,
        PRIMER_DNA_CONC=50.0,
    ),
    # Current primer3 defaults.
    'primer3_default': dict(
        PRIMER_TM_FORMULA=1,          # SantaLucia 1998
        PRIMER_SALT_CORRECTIONS=1,    # SantaLucia 1998
        PRIMER_SALT_MONOVALENT=50.0,
        PRIMER_SALT_DIVALENT=1.5,
        PRIMER_DNTP_CONC=0.6,
        PRIMER_DNA_CONC=50.0,
    ),
}

_CALC_TM_METHOD = {0: 'breslauer', 1: 'santalucia'}
_CALC_SALT_METHOD = {0: 'schildkraut', 1: 'santalucia', 2: 'owczarzy'}


# ----------------------------------------------------------------------------
# Sequence helpers
# ----------------------------------------------------------------------------

def clean_sequence(seq: str) -> str:
    """Upper-case and strip whitespace/digits (e.g. GenBank-formatted input)."""
    return ''.join(c for c in seq.upper() if c.isalpha())


def reverse_complement(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def gc_content(seq: str) -> float:
    """GC percentage over unambiguous bases."""
    seq = seq.upper()
    acgt = sum(seq.count(b) for b in 'ACGT')
    if acgt == 0:
        return 0.0
    return 100.0 * (seq.count('G') + seq.count('C')) / acgt


def read_fasta(path: str) -> list[tuple[str, str]]:
    """Return [(id, sequence), ...] from a FASTA file (or a bare-sequence file)."""
    records: list[tuple[str, str]] = []
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if name is not None:
                    records.append((name, clean_sequence(''.join(chunks))))
                name, chunks = line[1:].split()[0] if len(line) > 1 else 'seq', []
            else:
                chunks.append(line)
    if name is None and chunks:
        name = 'seq'
    if name is not None:
        records.append((name, clean_sequence(''.join(chunks))))
    return records


# ----------------------------------------------------------------------------
# Parameters
# ----------------------------------------------------------------------------

@dataclass
class HDAParams:
    """NEB tHDA design parameters. Ranges are (min, opt, max)."""

    product_size: tuple[int, int, int] = (80, 100, 120)
    product_tm: tuple[float, float, float] = (68.0, 71.0, 75.0)
    primer_size: tuple[int, int, int] = (24, 27, 33)
    primer_tm: tuple[float, float, float] = (60.0, 68.0, 74.0)
    primer_gc: tuple[float, float, float] = (35.0, 44.0, 60.0)
    max_primer_tm_diff: float = 5.0
    max_secondary_tm: float = 47.0      # hairpin / self-dimer / primer-dimer Tm (primer3 default)
    preferred_amplicon_gc: float = 40.0
    tm_model: str = 'neb_legacy'

    # Ranking weights applied on top of the primer3 pair penalty.
    weight_amplicon_gc: float = 0.2     # per % away from preferred_amplicon_gc
    weight_primer_amplicon_gc: float = 0.2  # per % |mean primer GC - amplicon GC|

    # Extra primer3 tags (override anything above).
    extra_primer3_args: dict = field(default_factory=dict)

    @classmethod
    def for_target(cls, target_gc: float, **kwargs) -> 'HDAParams':
        """Apply the NEB GC-dependent adjustments to primer size / Tm optima."""
        p = cls(**kwargs)
        smin, sopt, smax = p.primer_size
        tmin, topt, tmax = p.primer_tm
        if target_gc > 45.0:
            sopt = 26
            topt = 71.0
        elif target_gc < 37.5:
            topt = 65.0
        return replace(p, primer_size=(smin, sopt, smax), primer_tm=(tmin, topt, tmax))

    def thermo_settings(self) -> dict:
        try:
            return dict(TM_MODELS[self.tm_model])
        except KeyError:
            raise ValueError(f'tm_model must be one of {list(TM_MODELS)}') from None

    def to_primer3(self, num_return: int) -> dict:
        args = {
            'PRIMER_TASK': 'generic',
            'PRIMER_PICK_LEFT_PRIMER': 1,
            'PRIMER_PICK_RIGHT_PRIMER': 1,
            'PRIMER_PICK_INTERNAL_OLIGO': 0,
            'PRIMER_NUM_RETURN': num_return,
            'PRIMER_MIN_SIZE': self.primer_size[0],
            'PRIMER_OPT_SIZE': self.primer_size[1],
            'PRIMER_MAX_SIZE': self.primer_size[2],
            'PRIMER_MIN_TM': self.primer_tm[0],
            'PRIMER_OPT_TM': self.primer_tm[1],
            'PRIMER_MAX_TM': self.primer_tm[2],
            'PRIMER_MIN_GC': self.primer_gc[0],
            'PRIMER_OPT_GC_PERCENT': self.primer_gc[1],
            'PRIMER_MAX_GC': self.primer_gc[2],
            'PRIMER_WT_GC_PERCENT_LT': 0.1,
            'PRIMER_WT_GC_PERCENT_GT': 0.1,
            'PRIMER_PAIR_MAX_DIFF_TM': self.max_primer_tm_diff,
            'PRIMER_PRODUCT_SIZE_RANGE': [[self.product_size[0], self.product_size[2]]],
            'PRIMER_PRODUCT_OPT_SIZE': self.product_size[1],
            'PRIMER_PAIR_WT_PRODUCT_SIZE_LT': 0.05,
            'PRIMER_PAIR_WT_PRODUCT_SIZE_GT': 0.05,
            'PRIMER_PRODUCT_MIN_TM': self.product_tm[0],
            'PRIMER_PRODUCT_OPT_TM': self.product_tm[1],
            'PRIMER_PRODUCT_MAX_TM': self.product_tm[2],
            'PRIMER_PAIR_WT_PRODUCT_TM_LT': 0.2,
            'PRIMER_PAIR_WT_PRODUCT_TM_GT': 0.2,
            'PRIMER_MAX_HAIRPIN_TH': self.max_secondary_tm,
            'PRIMER_MAX_SELF_ANY_TH': self.max_secondary_tm,
            'PRIMER_MAX_SELF_END_TH': self.max_secondary_tm,
            'PRIMER_PAIR_MAX_COMPL_ANY_TH': self.max_secondary_tm,
            'PRIMER_PAIR_MAX_COMPL_END_TH': self.max_secondary_tm,
            'PRIMER_MAX_POLY_X': 4,
            'PRIMER_MAX_NS_ACCEPTED': 0,
            'PRIMER_EXPLAIN_FLAG': 1,
        }
        args.update(self.thermo_settings())
        args.update(self.extra_primer3_args)
        return args


# ----------------------------------------------------------------------------
# Thermodynamic helpers (use the same model primer3 is configured with)
# ----------------------------------------------------------------------------

def primer_tm(seq: str, params: HDAParams) -> float:
    t = params.thermo_settings()
    return primer3.calc_tm(
        seq,
        mv_conc=t['PRIMER_SALT_MONOVALENT'],
        dv_conc=t['PRIMER_SALT_DIVALENT'],
        dntp_conc=t['PRIMER_DNTP_CONC'],
        dna_conc=t['PRIMER_DNA_CONC'],
        tm_method=_CALC_TM_METHOD[t['PRIMER_TM_FORMULA']],
        salt_corrections_method=_CALC_SALT_METHOD[t['PRIMER_SALT_CORRECTIONS']],
    )


def _thal_kwargs(params: HDAParams) -> dict:
    t = params.thermo_settings()
    return dict(mv_conc=t['PRIMER_SALT_MONOVALENT'], dv_conc=t['PRIMER_SALT_DIVALENT'],
                dntp_conc=t['PRIMER_DNTP_CONC'], dna_conc=t['PRIMER_DNA_CONC'])


def product_tm(amplicon: str, params: HDAParams) -> float:
    """Amplicon Tm, identical to primer3's PRIMER_PAIR_*_PRODUCT_TM (long_seq_tm)."""
    t = params.thermo_settings()
    na = t['PRIMER_SALT_MONOVALENT']
    dv, dntp = t['PRIMER_SALT_DIVALENT'], t['PRIMER_DNTP_CONC']
    if dv > dntp:
        na += 120.0 * math.sqrt(dv - dntp)
    gc = gc_content(amplicon)
    return 81.5 + 16.6 * math.log10(na / 1000.0) + 0.41 * gc - 600.0 / len(amplicon)


# ----------------------------------------------------------------------------
# Design
# ----------------------------------------------------------------------------

@dataclass
class PrimerPair:
    seq_id: str
    region: int
    rank_in_region: int
    left_seq: str
    right_seq: str
    left_start: int          # 0-based, 5' end of left primer on template
    left_len: int
    right_end: int           # 0-based, 5' end of right primer on template (amplicon end)
    right_len: int
    left_tm: float
    right_tm: float
    left_gc: float
    right_gc: float
    product_size: int
    product_tm: float
    amplicon_gc: float
    primer_mean_gc: float
    gc_match: float          # |mean primer GC - amplicon GC|  (consideration 2)
    primer3_penalty: float
    score: float             # lower is better
    left_hairpin_th: float
    right_hairpin_th: float
    pair_compl_any_th: float
    pair_compl_end_th: float
    amplicon: str

    @property
    def amplicon_start(self) -> int:
        return self.left_start

    @property
    def amplicon_end(self) -> int:
        return self.right_end

    def as_dict(self) -> dict:
        # vars() also picks up annotations added later, e.g. by hda_blast.check_specificity.
        d = dict(vars(self))
        # 1-based inclusive coordinates are friendlier for reporting.
        d['amplicon_start_1based'] = self.left_start + 1
        d['amplicon_end_1based'] = self.right_end + 1
        return d


@lru_cache(maxsize=200000)
def _secondary_tms(seq: str, other: str | None, tm_model: str) -> tuple[float, float]:
    """(hairpin Tm, self-dimer Tm) of seq, or (heterodimer Tm, 0) with other."""
    kw = _thal_kwargs(HDAParams(tm_model=tm_model))
    if other is not None:
        return primer3.calc_heterodimer(seq, other, **kw).tm, 0.0
    return primer3.calc_hairpin(seq, **kw).tm, primer3.calc_homodimer(seq, **kw).tm


def _parse_pairs(res: dict, template: str, params: HDAParams, seq_id: str,
                 thal_post: bool = False) -> list[PrimerPair]:
    """Parse primer3 pairs. With ``thal_post`` (primer3 run without thermodynamic
    alignment), hairpin/dimer Tms are computed here and failing pairs dropped."""
    pairs = []
    for i in range(res.get('PRIMER_PAIR_NUM_RETURNED', 0)):
        l_start, l_len = res[f'PRIMER_LEFT_{i}']
        r_end, r_len = res[f'PRIMER_RIGHT_{i}']
        amplicon = template[l_start:r_end + 1]
        l_seq = res[f'PRIMER_LEFT_{i}_SEQUENCE']
        r_seq = res[f'PRIMER_RIGHT_{i}_SEQUENCE']
        l_gc = res[f'PRIMER_LEFT_{i}_GC_PERCENT']
        r_gc = res[f'PRIMER_RIGHT_{i}_GC_PERCENT']
        amp_gc = gc_content(amplicon)
        mean_gc = (l_gc + r_gc) / 2.0
        gc_match = abs(mean_gc - amp_gc)
        penalty = res[f'PRIMER_PAIR_{i}_PENALTY']
        if thal_post:
            l_hp, l_sd = _secondary_tms(l_seq, None, params.tm_model)
            r_hp, r_sd = _secondary_tms(r_seq, None, params.tm_model)
            het, _ = _secondary_tms(*sorted((l_seq, r_seq)), params.tm_model)
            if max(l_hp, l_sd, r_hp, r_sd, het) > params.max_secondary_tm:
                continue
            # Only the any-position heterodimer Tm is computed; it bounds the 3'-end one.
            thal = (l_hp, r_hp, het, het)
        else:
            thal = (res[f'PRIMER_LEFT_{i}_HAIRPIN_TH'], res[f'PRIMER_RIGHT_{i}_HAIRPIN_TH'],
                    res[f'PRIMER_PAIR_{i}_COMPL_ANY_TH'], res[f'PRIMER_PAIR_{i}_COMPL_END_TH'])
        score = (penalty
                 + params.weight_amplicon_gc * abs(amp_gc - params.preferred_amplicon_gc)
                 + params.weight_primer_amplicon_gc * gc_match)
        pairs.append(PrimerPair(
            seq_id=seq_id, region=-1, rank_in_region=-1,
            left_seq=l_seq, right_seq=r_seq,
            left_start=l_start, left_len=l_len, right_end=r_end, right_len=r_len,
            left_tm=res[f'PRIMER_LEFT_{i}_TM'], right_tm=res[f'PRIMER_RIGHT_{i}_TM'],
            left_gc=l_gc, right_gc=r_gc,
            product_size=res[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
            product_tm=res[f'PRIMER_PAIR_{i}_PRODUCT_TM'],
            amplicon_gc=amp_gc, primer_mean_gc=mean_gc, gc_match=gc_match,
            primer3_penalty=penalty, score=score,
            left_hairpin_th=thal[0], right_hairpin_th=thal[1],
            pair_compl_any_th=thal[2], pair_compl_end_th=thal[3],
            amplicon=amplicon,
        ))
    return pairs


def _overlap_fraction(a: PrimerPair, b: PrimerPair) -> float:
    ov = min(a.right_end, b.right_end) - max(a.left_start, b.left_start) + 1
    return max(0, ov) / min(a.product_size, b.product_size)


def design_hda_primers(
    sequence: str,
    seq_id: str = 'target',
    params: HDAParams | None = None,
    n_regions: int = 5,
    pairs_per_region: int = 3,
    region_overlap: float = 0.5,
    included_region: tuple[int, int] | None = None,
    window_step: int = 20,
    candidates_per_window: int = 20,
    auto_adjust_gc: bool = True,
    fast: bool | None = None,
    target: tuple[int, int] | None = None,
    keep: 'callable | None' = None,
) -> tuple[list[PrimerPair], dict]:
    """Design tHDA primer sets for several distinct amplicon regions.

    Primer3 is run on overlapping windows tiled along the template so that
    candidates cover the whole sequence. Candidates are scored (primer3 penalty +
    amplicon-GC deviation from 40 % + primer/amplicon GC mismatch) and greedily
    clustered into regions: a pair joins a region if its amplicon overlaps the
    region's best amplicon by more than ``region_overlap``.

    Args:
        sequence: template DNA.
        params: HDAParams; if None, NEB defaults adjusted for target GC.
        n_regions: number of distinct amplicon regions to report.
        pairs_per_region: primer sets reported per region.
        region_overlap: fractional overlap above which two amplicons are the same region.
        included_region: optional (start, length), 0-based, to restrict amplicons.
        window_step: step between primer3 search windows (bp).
        candidates_per_window: PRIMER_NUM_RETURN for each window.
        auto_adjust_gc: apply the NEB GC-dependent size/Tm optima (only if params is None).
        target: (start, length) relative to ``sequence``, 0-based; the product must
            span it (primer3 SEQUENCE_TARGET), and only windows containing it are searched.
        keep: optional predicate applied to every candidate pair before ranking.
        fast: run primer3 without thermodynamic oligo alignment and apply the
            hairpin/self-dimer/primer-dimer Tm checks only to returned pairs (~50-100x
            faster; same final criteria). Default: on for targets > 2 kb.

    Returns:
        (pairs, info) where info contains the target GC and parameters used.
    """
    template = clean_sequence(sequence)
    lo, length = included_region if included_region else (0, len(template))
    hi = lo + length
    region_seq = template[lo:hi]
    target_gc = gc_content(region_seq)

    if params is None:
        params = HDAParams.for_target(target_gc) if auto_adjust_gc else HDAParams()

    if len(region_seq) < params.product_size[0]:
        raise ValueError(f'{seq_id}: target ({len(region_seq)} bp) shorter than minimum '
                         f'product size ({params.product_size[0]} bp)')

    if fast is None:
        fast = len(region_seq) > 2000
    if fast:
        global_args = params.to_primer3(candidates_per_window)
        global_args['PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT'] = 0
        # The score-based self-complementarity limits used without thermodynamic
        # alignment are much stricter for 24-33-mers; disable them and rely on the
        # thermodynamic post-check in _parse_pairs instead.
        for k in ('PRIMER_MAX_SELF_ANY', 'PRIMER_MAX_SELF_END',
                  'PRIMER_PAIR_MAX_COMPL_ANY', 'PRIMER_PAIR_MAX_COMPL_END'):
            global_args[k] = 1000
        win = min(len(region_seq), 250)
        window_step = max(window_step, 50)
    else:
        global_args = params.to_primer3(candidates_per_window)
        # Windows a bit larger than the max product so each window yields many options.
        win = min(len(region_seq), params.product_size[2] + 40)
    if target is not None:
        t0, tlen = target
        t1 = t0 + tlen
        if tlen > params.product_size[2]:
            raise ValueError(f'target ({tlen} bp) longer than the maximum product size '
                             f'({params.product_size[2]} bp)')
        win = max(win, tlen)
        # Only windows that fully contain the target can yield a spanning product.
        w_lo = max(lo, t1 - win)
        w_hi = min(t0, hi - win)
        starts = list(range(w_lo, w_hi + 1, max(1, window_step)))
        if not starts:
            raise ValueError('no search window contains the target; increase flank')
        if starts[-1] != w_hi:
            starts.append(w_hi)
    else:
        starts = list(range(lo, max(lo, hi - win) + 1, window_step))
        if starts[-1] != hi - win:
            starts.append(hi - win)

    candidates: dict[tuple, PrimerPair] = {}
    explain: list[str] = []
    for ws in starts:
        seq_args = {'SEQUENCE_ID': seq_id, 'SEQUENCE_TEMPLATE': template,
                    'SEQUENCE_INCLUDED_REGION': [ws, win]}
        if target is not None:
            seq_args['SEQUENCE_TARGET'] = [target[0], target[1]]
        res = primer3.bindings.design_primers(seq_args, global_args)
        if res.get('PRIMER_PAIR_NUM_RETURNED', 0) == 0:
            explain.append('; '.join(
                f'{k.lower()}: {res.get(f"PRIMER_{k}_EXPLAIN", "")}'
                for k in ('LEFT', 'RIGHT', 'PAIR')))
        for p in _parse_pairs(res, template, params, seq_id, thal_post=fast):
            if keep is not None and not keep(p):
                continue
            key = (p.left_start, p.left_len, p.right_end, p.right_len)
            candidates.setdefault(key, p)

    ranked = sorted(candidates.values(), key=lambda p: p.score)
    regions: list[list[PrimerPair]] = []
    for p in ranked:
        for members in regions:
            if _overlap_fraction(p, members[0]) > region_overlap:
                if len(members) < pairs_per_region:
                    members.append(p)
                break
        else:
            if len(regions) < n_regions:
                regions.append([p])

    out = []
    for r_idx, members in enumerate(regions, start=1):
        for k, p in enumerate(members, start=1):
            p.region, p.rank_in_region = r_idx, k
            out.append(p)

    info = {
        'seq_id': seq_id,
        'length': len(template),
        'target_gc': target_gc,
        'params': params,
        'n_candidates': len(candidates),
        'primer3_explain_when_empty': explain[:3],
    }
    return out, info


# ----------------------------------------------------------------------------
# Consideration 4: serial longer / shorter primer variants
# ----------------------------------------------------------------------------

@dataclass
class PrimerVariant:
    pair_label: str
    side: str                 # 'left' or 'right'
    delta: int                # nt added (+) / removed (-) at the 5' end
    sequence: str
    length: int
    tm: float
    gc: float
    hairpin_tm: float
    homodimer_tm: float
    product_size: int
    product_tm: float
    within_neb_ranges: bool


def serial_primer_variants(
    template: str,
    pair: PrimerPair,
    params: HDAParams,
    max_delta: int = 4,
) -> list[PrimerVariant]:
    """Generate primers 5'-extended / 5'-trimmed relative to a working set.

    The 3' end (which determines priming) is kept fixed; each primer is varied
    by -max_delta..+max_delta nt at its 5' end while its partner is unchanged.
    """
    template = clean_sequence(template)
    thal = _thal_kwargs(params)
    variants: list[PrimerVariant] = []
    label = f'{pair.seq_id}_R{pair.region}.{pair.rank_in_region}'

    for side in ('left', 'right'):
        for d in range(-max_delta, max_delta + 1):
            if side == 'left':
                start = pair.left_start - d
                length = pair.left_len + d
                if start < 0 or length < 15:
                    continue
                seq = template[start:start + length]
                amp_start, amp_end = start, pair.right_end
            else:
                end = pair.right_end + d
                length = pair.right_len + d
                if end >= len(template) or length < 15:
                    continue
                seq = reverse_complement(template[end - length + 1:end + 1])
                amp_start, amp_end = pair.left_start, end
            amplicon = template[amp_start:amp_end + 1]
            tm = primer_tm(seq, params)
            gc = gc_content(seq)
            ptm = product_tm(amplicon, params)
            ok = (params.primer_size[0] <= length <= params.primer_size[2]
                  and params.primer_tm[0] <= tm <= params.primer_tm[2]
                  and params.primer_gc[0] <= gc <= params.primer_gc[2]
                  and params.product_size[0] <= len(amplicon) <= params.product_size[2]
                  and params.product_tm[0] <= ptm <= params.product_tm[2])
            variants.append(PrimerVariant(
                pair_label=label, side=side, delta=d, sequence=seq, length=length,
                tm=tm, gc=gc,
                hairpin_tm=primer3.calc_hairpin(seq, **thal).tm,
                homodimer_tm=primer3.calc_homodimer(seq, **thal).tm,
                product_size=len(amplicon), product_tm=ptm, within_neb_ranges=ok,
            ))
    return variants


# ----------------------------------------------------------------------------
# Reporting
# ----------------------------------------------------------------------------

REPORT_COLUMNS = [
    'seq_id', 'region', 'rank_in_region', 'score',
    'amplicon_start_1based', 'amplicon_end_1based', 'product_size', 'product_tm',
    'amplicon_gc', 'primer_mean_gc', 'gc_match',
    'left_seq', 'left_len', 'left_tm', 'left_gc',
    'right_seq', 'right_len', 'right_tm', 'right_gc',
    'left_hairpin_th', 'right_hairpin_th', 'pair_compl_any_th', 'pair_compl_end_th',
    'primer3_penalty', 'amplicon',
]


def _fmt(v):
    return f'{v:.2f}' if isinstance(v, float) else v


def write_table(pairs: Iterable[PrimerPair], path_or_handle, delimiter: str = '\t') -> None:
    pairs = list(pairs)
    columns = list(REPORT_COLUMNS)
    if any(hasattr(p, 'specific') for p in pairs):
        from hda_blast import SPECIFICITY_COLUMNS
        columns[4:4] = SPECIFICITY_COLUMNS
    close = False
    fh = path_or_handle
    if isinstance(path_or_handle, str):
        fh = open(path_or_handle, 'w', newline='')
        close = True
    try:
        w = csv.DictWriter(fh, fieldnames=columns, delimiter=delimiter,
                           extrasaction='ignore', restval='')
        w.writeheader()
        for p in pairs:
            w.writerow({k: _fmt(v) for k, v in p.as_dict().items()})
    finally:
        if close:
            fh.close()


def summarize(pairs: list[PrimerPair], info: dict) -> str:
    p = info['params']
    lines = [
        f"# {info['seq_id']}: {info['length']} bp, target GC {info['target_gc']:.1f}%, "
        f"tm_model={p.tm_model}, primer opt size {p.primer_size[1]} nt, "
        f"primer opt Tm {p.primer_tm[1]:.0f} C, {info['n_candidates']} candidate pairs",
    ]
    if not pairs:
        lines.append('  No primer pairs satisfied the tHDA constraints.')
        for e in info['primer3_explain_when_empty']:
            lines.append(f'  primer3: {e}')
        return '\n'.join(lines)
    for q in pairs:
        lines.append(
            f"  R{q.region}.{q.rank_in_region}  {q.left_start + 1}-{q.right_end + 1} "
            f"({q.product_size} bp, Tm {q.product_tm:.1f}, GC {q.amplicon_gc:.1f}%)  "
            f"score {q.score:.2f}\n"
            f"      F {q.left_seq}  ({q.left_len} nt, Tm {q.left_tm:.1f}, GC {q.left_gc:.0f}%)\n"
            f"      R {q.right_seq}  ({q.right_len} nt, Tm {q.right_tm:.1f}, GC {q.right_gc:.0f}%)"
        )
    return '\n'.join(lines)


# ----------------------------------------------------------------------------
# CLI
# ----------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(
        description='Design tHDA primers following NEB recommendations.')
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument('fasta', nargs='?', help='FASTA file with one or more target sequences')
    src.add_argument('--sequence', help='target sequence given directly')
    ap.add_argument('--id', default='target', help='sequence id when using --sequence')
    ap.add_argument('-o', '--output', help='write TSV (or CSV if name ends in .csv)')
    ap.add_argument('--n-regions', type=int, default=5)
    ap.add_argument('--pairs-per-region', type=int, default=3)
    ap.add_argument('--region', help='restrict amplicons to START-END (1-based, inclusive)')
    ap.add_argument('--tm-model', choices=list(TM_MODELS), default='neb_legacy')
    ap.add_argument('--no-gc-adjust', action='store_true',
                    help='do not adjust primer size/Tm optima for target GC')
    ap.add_argument('--variants', type=int, default=0, metavar='N',
                    help='also report 5\' serial variants (+/-N nt) for each top pair')

    bl = ap.add_argument_group('BLAST specificity / cross-species uniqueness')
    bdb = bl.add_mutually_exclusive_group()
    bdb.add_argument('--blast-db', help='local BLAST nucleotide database')
    bdb.add_argument('--blast-remote', nargs='?', const='core_nt', metavar='DB',
                     help='BLAST at NCBI (default db core_nt); can take many minutes')
    bl.add_argument('--entrez-query', help="remote only, e.g. 'Bacteria[organism]'")
    bl.add_argument('--target-species', action='append',
                    help='intended target species (repeatable); inferred if omitted')
    bl.add_argument('--max-mismatches', type=int, default=3,
                    help='max primer mismatches for a binding site (default 3)')
    bl.add_argument('--max-offtarget-product', type=int, default=1000,
                    help='max size of predicted off-target products (default 1000)')
    bl.add_argument('--min-identity', type=float, default=90.0,
                    help='%% identity for a cross-reactive species (default 90)')
    bl.add_argument('--min-coverage', type=float, default=80.0,
                    help='%% amplicon coverage for a cross-reactive species (default 80)')
    bl.add_argument('--blast-threads', type=int, default=4)
    args = ap.parse_args(argv)
    if args.entrez_query and not args.blast_remote:
        ap.error('--entrez-query requires --blast-remote')

    records = ([(args.id, clean_sequence(args.sequence))] if args.sequence
               else read_fasta(args.fasta))
    included = None
    if args.region:
        s, e = (int(x) for x in args.region.split('-'))
        included = (s - 1, e - s + 1)

    all_pairs: list[PrimerPair] = []
    for seq_id, seq in records:
        gc = gc_content(seq[included[0]:sum(included)] if included else seq)
        params = (HDAParams(tm_model=args.tm_model) if args.no_gc_adjust
                  else HDAParams.for_target(gc, tm_model=args.tm_model))
        pairs, info = design_hda_primers(
            seq, seq_id=seq_id, params=params, n_regions=args.n_regions,
            pairs_per_region=args.pairs_per_region, included_region=included)
        print(summarize(pairs, info), file=sys.stderr)
        all_pairs.extend(pairs)

        if args.variants:
            for p in (q for q in pairs if q.rank_in_region == 1):
                print(f'  serial variants for R{p.region}.1:', file=sys.stderr)
                for v in serial_primer_variants(seq, p, params, args.variants):
                    flag = 'ok ' if v.within_neb_ranges else '-- '
                    print(f'    {flag}{v.side:5s} {v.delta:+d} {v.sequence} '
                          f'({v.length} nt, Tm {v.tm:.1f}, GC {v.gc:.0f}%, '
                          f'product {v.product_size} bp / {v.product_tm:.1f} C)',
                          file=sys.stderr)

    spec = None
    if (args.blast_db or args.blast_remote) and all_pairs:
        import hda_blast
        cfg = hda_blast.BlastConfig(
            db=args.blast_db or args.blast_remote, remote=bool(args.blast_remote),
            entrez_query=args.entrez_query, threads=args.blast_threads)
        spec = hda_blast.check_specificity(
            all_pairs, cfg, target_species=args.target_species,
            max_mismatches=args.max_mismatches, max_product=args.max_offtarget_product,
            min_identity=args.min_identity, min_coverage=args.min_coverage)
        print('\n# Specificity (target: '
              f'{", ".join(spec.target_species) or "unknown"})', file=sys.stderr)
        print(hda_blast.summarize_specificity(all_pairs, spec), file=sys.stderr)

    print('\nNote: order desalted (salt-free) primers; tHDA is sensitive to salt.',
          file=sys.stderr)
    if args.output:
        write_table(all_pairs, args.output, ',' if args.output.endswith('.csv') else '\t')
        if spec is not None:
            stem = args.output.rsplit('.', 1)[0]
            hda_blast.write_products(spec.products, f'{stem}.blast_products.tsv')
            hda_blast.write_species_table(spec.uniqueness, f'{stem}.amplicon_species.tsv')
            hda_blast.write_primer_table(hda_blast.primer_specificity(spec, all_pairs),
                                         f'{stem}.primer_specificity.tsv')
            print(f'Wrote {stem}.blast_products.tsv, {stem}.amplicon_species.tsv, '
                  f'{stem}.primer_specificity.tsv', file=sys.stderr)
    else:
        write_table(all_pairs, sys.stdout)
    return 0


if __name__ == '__main__':
    sys.exit(main())


def neb_check(pair: PrimerPair, params: HDAParams | None = None) -> list[tuple]:
    """Check one pair against the NEB tHDA ranges.

    Returns [(criterion, value, 'min-max', ok), ...]; ``ok`` is None for the
    preferences (NEB "other considerations") that are not hard limits.
    """
    p = params or HDAParams()
    rows = [
        ('Product size (bp)', pair.product_size, p.product_size, True),
        ('Product Tm (C)', pair.product_tm, p.product_tm, True),
        ('Forward length (nt)', pair.left_len, p.primer_size, True),
        ('Reverse length (nt)', pair.right_len, p.primer_size, True),
        ('Forward Tm (C)', pair.left_tm, p.primer_tm, True),
        ('Reverse Tm (C)', pair.right_tm, p.primer_tm, True),
        ('Forward GC (%)', pair.left_gc, p.primer_gc, True),
        ('Reverse GC (%)', pair.right_gc, p.primer_gc, True),
    ]
    out = []
    for name, value, (lo, _opt, hi), _hard in rows:
        out.append((name, value, f'{lo:g}-{hi:g}', lo <= value <= hi))
    out.append(('Primer Tm difference (C)', abs(pair.left_tm - pair.right_tm),
                f'<= {p.max_primer_tm_diff:g}',
                abs(pair.left_tm - pair.right_tm) <= p.max_primer_tm_diff))
    # NEB considerations 1 and 2: preferences, not limits.
    out.append(('Amplicon GC (%)', pair.amplicon_gc,
                f'~{p.preferred_amplicon_gc:g} preferred', None))
    out.append(('Mean primer GC vs amplicon GC (%)', pair.gc_match, 'small preferred', None))
    return out


def format_neb_check(rows: list[tuple]) -> str:
    lines = []
    for name, value, rng, ok in rows:
        mark = {True: 'ok', False: 'FAIL', None: '--'}[ok]
        v = f'{value:.1f}' if isinstance(value, float) else str(value)
        lines.append(f'  {mark:4s} {name:34s} {v:>8s}   ({rng})')
    return '\n'.join(lines)


def design_for_locus(
    genome: str,
    start: int,
    end: int,
    params: HDAParams | None = None,
    flank: int = 300,
    min_overlap: int = 1,
    n_regions: int = 8,
    pairs_per_region: int = 1,
    fast: bool | None = None,
) -> tuple[list[PrimerPair], dict]:
    """Design tHDA primers for an amplicon overlapping a locus in a larger sequence.

    Useful when a locus is too short, or too uneven in GC, to hold a compliant pair
    on its own: primers are searched in the locus plus ``flank`` bp either side, and
    only pairs whose amplicon overlaps the locus are kept.

    Args:
        genome: the full sequence (e.g. a chromosome).
        start, end: locus coordinates, 1-based inclusive.
        flank: bp of flanking sequence to search on each side.
        min_overlap: minimum bp of the locus the amplicon must cover.
        n_regions, pairs_per_region: as in design_hda_primers.

    Returns:
        (pairs sorted by locus coverage then design score, info). Coordinates are
        absolute (0-based) in ``genome``; each pair gains ``locus_overlap`` (bp) and
        ``locus_covered`` (% of the locus).
    """
    genome = clean_sequence(genome)
    lo0, hi0 = start - 1, end                      # half-open, 0-based
    if not 0 <= lo0 < hi0 <= len(genome):
        raise ValueError(f'locus {start}-{end} outside sequence (1-{len(genome)})')
    lo = max(0, lo0 - flank)
    hi = min(len(genome), hi0 + flank)
    sub = genome[lo:hi]
    if params is None:
        params = HDAParams.for_target(gc_content(sub))

    locus_len = hi0 - lo0
    if min_overlap >= locus_len:
        # Product must span the whole locus: let primer3 enforce it.
        kw = dict(target=(lo0 - lo, locus_len))
    else:
        want = min_overlap
        kw = dict(keep=lambda q: max(0, min(q.right_end + 1, hi0 - lo)
                                     - max(q.left_start, lo0 - lo)) >= want)
    pairs, info = design_hda_primers(sub, params=params, n_regions=n_regions,
                                     pairs_per_region=pairs_per_region, fast=fast, **kw)
    kept = []
    for q in pairs:
        q.left_start += lo
        q.right_end += lo
        q.locus_overlap = max(0, min(q.right_end + 1, hi0) - max(q.left_start, lo0))
        q.locus_covered = 100.0 * q.locus_overlap / (hi0 - lo0)
        if q.locus_overlap >= min_overlap:
            kept.append(q)
    kept.sort(key=lambda q: (-q.locus_overlap, q.score))
    for rank, q in enumerate(kept, 1):
        q.region, q.rank_in_region = rank, 1
    info.update(locus=(start, end), searched=(lo + 1, hi), n_overlapping=len(kept))
    return kept, info


def evaluate_primers(left_seq: str, right_seq: str, template: str | None = None,
                     params: HDAParams | None = None, seq_id: str = 'user') -> PrimerPair:
    """Score an existing primer pair against the NEB tHDA ranges.

    Both primers are given 5'->3'; the reverse primer is the reverse complement of
    the template's plus strand, as usual. If ``template`` is supplied (a region or a
    whole genome) and both primers match it exactly, the amplicon and its size, Tm
    and GC are filled in; otherwise those stay 0 and only per-primer values apply
    (use ``hda_blast.scan_genome`` to locate imperfect matches).

    Returns a PrimerPair, so ``neb_check`` and ``hda_blast.scan_genome`` accept it.
    """
    left_seq = clean_sequence(left_seq)
    right_seq = clean_sequence(right_seq)
    p = params or HDAParams()
    thal = _thal_kwargs(p)

    l_start = r_end = -1
    amplicon = ''
    if template:
        template = clean_sequence(template)
        l_start = template.find(left_seq)
        rc = reverse_complement(right_seq)
        r_start = template.find(rc, max(0, l_start))
        if l_start >= 0 and r_start >= 0:
            r_end = r_start + len(rc) - 1
            amplicon = template[l_start:r_end + 1]

    pair = PrimerPair(
        seq_id=seq_id, region=1, rank_in_region=1,
        left_seq=left_seq, right_seq=right_seq,
        left_start=l_start, left_len=len(left_seq),
        right_end=r_end, right_len=len(right_seq),
        left_tm=primer_tm(left_seq, p), right_tm=primer_tm(right_seq, p),
        left_gc=gc_content(left_seq), right_gc=gc_content(right_seq),
        product_size=len(amplicon),
        product_tm=product_tm(amplicon, p) if amplicon else 0.0,
        amplicon_gc=gc_content(amplicon) if amplicon else 0.0,
        primer_mean_gc=(gc_content(left_seq) + gc_content(right_seq)) / 2.0,
        gc_match=abs((gc_content(left_seq) + gc_content(right_seq)) / 2.0
                     - gc_content(amplicon)) if amplicon else 0.0,
        primer3_penalty=0.0, score=0.0,
        left_hairpin_th=primer3.calc_hairpin(left_seq, **thal).tm,
        right_hairpin_th=primer3.calc_hairpin(right_seq, **thal).tm,
        pair_compl_any_th=primer3.calc_heterodimer(left_seq, right_seq, **thal).tm,
        pair_compl_end_th=primer3.calc_end_stability(left_seq, right_seq, **thal).tm,
        amplicon=amplicon,
    )
    pair.left_homodimer_th = primer3.calc_homodimer(left_seq, **thal).tm
    pair.right_homodimer_th = primer3.calc_homodimer(right_seq, **thal).tm
    return pair
