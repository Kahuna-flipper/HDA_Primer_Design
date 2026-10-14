# tHDA Primer Design

Primer design for **thermophilic helicase-dependent amplification (tHDA)**, following the
NEB IsoAmp recommendations, with BLAST-based specificity checking.

Built on [primer3-py](https://github.com/libnano/primer3-py) and BLAST+.

## Why not plain Primer3

The NEB tHDA parameters were calibrated with the *original* Primer3 Tm model. Modern primer3
defaults (SantaLucia parameters, 1.5 mM Mg2+, 0.6 mM dNTP) give Tm values several degrees
higher, which rejects most good tHDA amplicons. The default `tm_model='neb_legacy'`
reproduces the original settings (Breslauer table, Schildkraut salt correction, 50 mM Na+,
no Mg2+/dNTP), so NEB's reference point holds: a ~100 bp amplicon at ~40 % GC gives a
product Tm of ~71 C. Pass `tm_model='primer3_default'` for modern thermodynamics.

## Install

```bash
conda create -n primer_design python=3.12
conda activate primer_design
conda install -c bioconda -c conda-forge primer3-py blast
```

Nothing else is required; pandas is optional (the notebook falls back to plain text).

## Quick start

`primer_design_analysis.ipynb` is the entry point: give it a template (or genome
coordinates) and it returns one NEB-compliant forward/reverse pair, then optionally
BLASTs the primers against a genome to show where else they could bind.

```python
from hda_primer_design import HDAParams, design_hda_primers, gc_content, neb_check, format_neb_check

params = HDAParams.for_target(gc_content(template))
pairs, info = design_hda_primers(template, params=params)
print(format_neb_check(neb_check(pairs[0], params)))
```

Command line:

```bash
# design
python hda_primer_design.py targets.fasta -o primers.tsv --n-regions 5 --pairs-per-region 3

# design + specificity across all bacteria
python hda_primer_design.py targets.fasta -o primers.tsv \
    --blast-db /path/to/ref_prok_rep_genomes --target-species "Escherichia coli"

# strain-specific amplicons: regions absent from a panel of other genomes
python hda_unique.py genome.fna -b /path/to/strain_panel -b /path/to/ref_prok_rep_genomes \
    -o strain_primers.tsv
```

## Modules

| file | contents |
|---|---|
| `hda_primer_design.py` | NEB parameters (`HDAParams`), design (`design_hda_primers`, `design_for_locus`), scoring of existing primers (`evaluate_primers`, `neb_check`), 5'-serial variants (`serial_primer_variants`) |
| `hda_blast.py` | BLAST databases (`blast_db_from_fasta`, `build_local_db`), single-genome binding sites (`scan_genome`), cross-species specificity and amplicon uniqueness (`check_specificity`) |
| `hda_unique.py` | strain-specific amplicons: prevalence scan, unique regions, design inside them, BLAST validation (`design_strain_specific`) |

### Notes on the specificity model

- A **binding site** is a place where the last `require_3p_match` bases of the primer match
  exactly (a 3'-terminal mismatch blocks extension) and the full-length primer footprint has
  at most `max_mismatches` mismatches. Mismatches are counted base by base against the
  genome, not read from the BLAST alignment, so partially aligned sites are not missed.
- A **predicted product** is two facing binding sites on one sequence within `max_product` bp,
  for any primer combination (F/R, F/F, R/R).
- For cross-species work, NCBI's `ref_prok_rep_genomes` (one reference genome per species,
  ~27 GB) plus `taxdb` is a good background:
  `update_blastdb.pl --decompress --source ncbi ref_prok_rep_genomes taxdb`.
  It tests *exclusivity* (other species), not *inclusivity* across strains of your target;
  for that, build a strain panel with `build_local_db`.

## Data

`e-coli-mg1655-genome.fasta` is the E. coli K-12 MG1655 reference, used by the notebook's
default example.

## NEB tHDA recommendations

The parameters implemented here, quoted from NEB:


Instructions

tHDA primers can be designed using either the PrimerQuest program or the Primer3 program. Recommended parameter settings:
Product size: 80 - 120 bp**
Product Tm: Min. 68; Opt. 71; Max. 75†
Primer size: Min. 24; Opt. 27; Max. 33‡
Primer Tm: Min. 60; Opt. 68; Max. 74§
Primer GC%: Min. 35; Opt. 44; Max. 60
(**) tHDA works most efficiently with a product size around 100 bp. Successful tHDA amplifications were achieved with a product as short as 85 bp and as long as 129 bp.
(†) The Tm of an amplicon with a product size around 100 bp and a G + C content around 40% is approximately 71°C from the calculation of Primer3. Successful tHDA amplifications were achieved with a product Tm as low as 68 °C and as high as 77°C.
(‡) The optimal primer size may be set at 26 bases when the G + C content of target sequence is larger than 45%. Successful tHDA amplifications were achieved with a primer size as short as 22 bases and as long as 32 bases.
(§) The optimal primer Tm may be set at 64 - 66°C when the G + C content of target sequence is smaller than 37.5% and at 70 - 72°C when the G + C content of target sequence is larger than 45%. Successful tHDA amplifications were achieved with a primer Tm as low as 60 °C and as high as 75 °C.

Other considerations:
1) Amplicons (regions to be amplified) containing a G + C content of approximately 40% are preferable (the G + C content of amplicons and target sequences can be calculated using the Oligonucleotide Properties Calculator program found at: (https://www.basic.northwestern.edu/biotools/oligocalc.html).
2) Primer sets obtained using Primer3 with an average G + C content closest to the G + C content of the corresponding amplicons are preferable. 
3) To obtain ideal amplification performance, test several amplicon regions for each target gene and design several primer sets for each selected region. The efficiency of the tHDA reactions may vary dramatically for different amplicons as well as for different primer sets within the specific region. Using primers optimized for PCR reactions may not guarantee similar performance in tHDA reactions.
4) Once a set of primers yields a positive amplification result, serial primer sets longer or shorter than this set may be analyzed to generate optimal primers for increased amplification efficiency.
5) Due to the sensitivity of tHDA to changes in salt concentrations in the reaction, the use of desalted primers is recommended. Successful tHDA reactions were achieved using salt-free primers synthesized by Operon Biotechnologies.
