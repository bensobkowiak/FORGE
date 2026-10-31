# FORGE

**Fast ONT Resistance and Genome Evaluation**

A Nextflow (DSL2) pipeline for Oxford Nanopore bacterial isolates: read QC, assembly,
plasmid reconstruction, taxonomic screening, typing and antimicrobial resistance
detection. One wide results table per collection, one directory per sample.

---

## What it does

| Stage | Tool |
|---|---|
| Read QC | NanoPlot |
| Read taxonomy | Kraken2 |
| Filtering & assembly | Filtlong → Flye (`--nano-hq`) |
| Plasmid reconstruction | MOB-suite (MOB-recon / MOB-typer) |
| Assembly taxonomy | Kraken2, per contig, length-weighted |
| Genome quality | CheckM2, BUSCO |
| Depth & circularity | minimap2 / samtools, Flye assembly graph |
| Typing | mlst, Kleborate (+ Kaptive) |
| AMR | AMRFinderPlus, ResFinder, RGI, CefiderocolFinder |

---

## Requirements

- Nextflow ≥ 23.04
- conda + mamba (environments are built automatically from `envs/*.yml`)
- Databases, supplied as paths:
  - **Kraken2** database directory (use the same build on every machine —
    results are not comparable across database versions)
  - **ResFinder** database directory
  - **AMRFinderPlus** database — *the conda package ships none*; create it once with
    `amrfinder_update -d <dir>`
  - **CheckM2** `.dmnd` (auto-downloaded if omitted)

---

## Quick start

### Single machine

```bash
./run_pipeline.sh --fastq /path/to/fastq_dir
```

Input is either barcode subdirectories or a flat directory of `fastq.gz`.
Defaults assume databases under `$PROGRAMS` (`~/Documents/Programs`):

```bash
./run_pipeline.sh --fastq reads/ \
  --maldi_csv Maldi_species.csv \
  --kraken_db    $PROGRAMS/KrakenDB \
  --resfinder_db $PROGRAMS/db_resfinder \
  --amrfinder_db $PROGRAMS/amrfinderplus_db/latest
```

| Option | Purpose |
|---|---|
| `--fastq <dir>` | input (required) |
| `--maldi_csv <file>` | MALDI-TOF species per isolate; enables contamination detection |
| `--results_file <file>` | wide results table (default `results.tsv`) |
| `--samples <file>` | process only these samples, one per line; **re-runs them even if already present** |
| `--outdir <dir>` | per-sample outputs (default `results/`) |
| `--fresh` | ignore the Nextflow cache |

### SGE cluster

```bash
qsub pipeline_sge.sh --fastq /path/to/fastq_dir [--maldi_csv <file>]
```

---

## Resuming

Samples already named in the results table are skipped, so re-running on a directory
only processes what is new. Most stages are additionally gated on their **published
output existing**, so work survives a deleted `work/` directory. `--samples` overrides
the skip, which is how you force a re-run.

---

## Output

```
results.tsv                     wide table: one row per isolate
results/<sample>/
  reads/nanoplot/               read QC
  assembly/flye/                assembly + assembly_info.txt
  assembly/mob_recon/           plasmid reconstruction
  assembly/contigs/             chromosome / plasmid / query FASTAs
  chromosome/                   kraken2, checkm2, busco, mlst, kleborate
  chromosome/{amrfinder,resfinder,rgi}/
  plasmids/{amrfinder,resfinder,rgi}/
```

`results.tsv` holds 64 fixed columns followed by one column per resistance gene, valued
by where the gene was found and which tool called it:

```
chromosome(AMRfinder:98.5%,RGI:91.0%);AA002(AMRfinder:97.2%)
```

Every column is defined in **`results_data_dictionary.xlsx`**.

---

## Behaviour worth knowing

**Genus filtering.** With a MALDI CSV, chromosomal contigs are filtered to the MALDI
genus. Filtering only applies where that genus accounts for ≥1 Mb **and** ≥25% of the
assembly; below either threshold the assembly is kept intact and `Maldi_contigs_used`
is set to `Mismatch`, rather than reducing the genome to a fragment.

**Mobile-element rescue.** MOB-recon calls a contig chromosomal whenever it finds no
replicon, so plasmids and ICEs carrying only a relaxase are misassigned. Chromosomal
contigs with a relaxase or *oriT* and length ≤ `mge_max_len` (1 Mb) are moved into the
plasmid set as `rescued_MGE_<n>`. Read depth is deliberately **not** a criterion -
plasmid and chromosome coverage are not reliably separated in ONT libraries.

**Query contigs.** Off-genus contigs below 1 Mb with no evidence of mobility are
neither chromosome nor plasmid. They are screened for AMR (reported with location
`query`) but excluded from chromosome and plasmid metrics. Off-genus contigs ≥1 Mb are
discarded and their lengths recorded in `Removed_Contig_Lengths`.

**Circularity** is judged on the **longest** chromosomal contig. An assembly can close
the chromosome while still carrying prophage or unplaced fragments, so extra contigs do
not make the answer `No`; `Chr_Contigs` reports fragmentation separately.

**Kleborate presets** are chosen from the chromosome's species — `kpsc`, `kosc` or
`escherichia`. Species outside all three (e.g. *K. aerogenes*) are not typed, since an
ST from an unrelated scheme is not interpretable. Kaptive K/O typing and clonal groups
are KpSC-only.

---

## Helper scripts

| Script | Purpose |
|---|---|
| `organise_by_genus.py` | reorganise `results/` into per-genus folders of assemblies and reports |
| `check_organised_fastas.py` | verify those FASTAs against the results table |

---

## Configuration

Tunable parameters live in `nextflow.config`; profiles are `standard`, `macstudio`,
`sge`, `lsf`, `docker`, `singularity` and `test`.

| Parameter | Default | Meaning |
|---|---|---|
| `filtlong_target_bases` | 750000000 | assembly input cap (~150× for 5 Mb); depth is measured on the full read set |
| `flye_genome_size` | `null` | let Flye estimate per sample |
| `mge_rescue` | `true` | enable the mobile-element rescue |
| `mge_max_len` | 1000000 | largest plausible non-chromosomal element |
| `maldi_min_genus_len` | 1000000 | genus filtering threshold (bp) |
| `maldi_min_genus_frac` | 0.25 | genus filtering threshold (fraction) |

---

## Notes

- `errorStrategy = 'ignore'`: a failing sample is dropped, not fatal to the run.
- CefiderocolFinder and ParEx are **not** installed by `envs/amr.yml`; without them
  those steps write a stub and their columns stay empty.
- `PDC_type` is currently unpopulated (see above); PDC variants are still reported in
  the AMR gene columns by AMRFinder/RGI.
