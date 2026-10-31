#!/usr/bin/env python3
"""
Build (or update) the wide results TSV.

One column per unique resistance gene across all samples and tools.
Column set grows automatically as new genes are found in new runs.

Usage:
  write_results.py --results_file <path>
                   --metrics   <s1.metrics.json> [<s2.metrics.json> ...]
                   --amr       <s1.amr.tsv>      [<s2.amr.tsv> ...]

Strategy:
  1. Load the existing results file (retains old sample rows and old AMR columns).
  2. Load per-sample fixed metrics JSON files from this run.
  3. Load per-sample AMR TSV files from this run.
  4. Union all AMR gene column names (old + new).
  5. Rewrite the entire file: fixed cols, then AMR cols (sorted).

AMR gene column values use the format produced by combine_amr.py, e.g.:
  chromosome(AMRfinder:98.5%,RGI:91.0%);pUTI89(AMRfinder:97.2%)
"""
import os
import sys
import json
import csv
import argparse
from collections import defaultdict

FIXED_COLS = [
    'Name', 'Data_GB', 'Num_Reads', 'Num_Reads_Q10', 'Median_Q', 'Median_Read_Length',
    'Major_Genus_Reads', 'Pct_Major_Genus_Reads', 'Genera_gt5pct_Reads', 'Major_Species_Reads', 'Pct_Major_Species_Reads', 'Species_gt5pct_Reads', 'Pct_Unclassified_Reads',
    'MALDI_Species_Match',
    'Maldi_genus_assembled',
    'Maldi_contigs_used',
    'Assembly_Major_Genus', 'Assembly_Pct_Major_Genus', 'Assembly_Genera_gt5pct', 'Assembly_Major_Species', 'Assembly_Pct_Major_Species', 'Assembly_Species_gt5pct', 'Assembly_Pct_Unclassified',
    'Removed_Contig_Lengths',
    'Chr_Contig_Species',
    'Circularised', 'Chr_Contigs', 'Completeness', 'Contamination',
    'Genome_Size', 'Contig_N50', 'Max_Contig_Length', 'Avg_Gene_Length', 'Total_CDS', 'GC_Content',
    'Chr_Mean_Depth', 'Coverage_Cap_Applied',
    'Num_Plasmids', 'Num_Rescued_MGE', 'Num_Query_Contigs', 'Plasmid_Names', 'Plasmid_Species',
    'Primary_Mob_Clusters', 'Secondary_Mob_Clusters', 'Inc_Replicons',
    'Plasmid_Mean_Depths',
    'BUSCO_Complete', 'BUSCO_Single_Copy', 'BUSCO_Duplicated', 'BUSCO_Fragmented', 'BUSCO_Missing',
    'MLST_ST', 'K_type', 'O_type', 'Kleborate_ST', 'Clonal_Group',
    'Virulence_Score', 'Yersiniabactin', 'Colibactin', 'Aerobactin', 'Salmochelin', 'RmpADC', 'RmpA2',
    'PDC_type',
]

def load_existing(results_file):
    """Return (rows_dict, amr_cols) from an existing results TSV.
    rows_dict: sample_name -> full row dict
    amr_cols:  ordered list of AMR column names (everything after FIXED_COLS)
    """
    rows = {}
    amr_cols = []
    if not os.path.exists(results_file) or os.path.getsize(results_file) == 0:
        return rows, amr_cols
    with open(results_file, newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        file_cols = reader.fieldnames or []
        amr_cols  = [c for c in file_cols if c not in FIXED_COLS]
        for row in reader:
            name = row.get('Name', '').strip()
            if name:
                rows[name] = dict(row)
    return rows, amr_cols

def load_amr_tsv(path):
    """Return {gene: value} dict from a combine_amr.py output TSV."""
    genes = {}
    if not path or not os.path.exists(path):
        return genes
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            gene  = row.get('gene', '').strip()
            value = row.get('value', '').strip()
            if gene and gene != 'gene':
                genes[gene] = value
    return genes

def load_sample_tsvs(paths):
    """Load per-sample TSV files (header + 1 data row each).
    Returns (rows_dict, extra_cols) where extra_cols are non-FIXED_COLS columns found."""
    rows = {}
    extra_cols = []
    seen_extra = set()
    fixed_set = set(FIXED_COLS)
    for path in paths:
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            continue
        with open(path, newline='', encoding='utf-8') as fh:
            reader = csv.DictReader(fh, delimiter='\t')
            for col in (reader.fieldnames or []):
                if col not in fixed_set and col not in seen_extra:
                    seen_extra.add(col)
                    extra_cols.append(col)
            for row in reader:
                name = row.get('Name', '').strip()
                if name:
                    rows[name] = dict(row)
    return rows, extra_cols


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--results_file',  required=True)
    p.add_argument('--sample_tsvs',   nargs='*', default=[],
                   help='Per-sample *_results.tsv files to merge into results_file')
    p.add_argument('--metrics',        nargs='*', default=[])
    p.add_argument('--amr',            nargs='*', default=[])
    args = p.parse_args()

    # ------------------------------------------------------------------
    # 1. Load existing rows and their extra (AMR) columns
    # ------------------------------------------------------------------
    old_rows, old_amr_cols = load_existing(args.results_file)

    # ------------------------------------------------------------------
    # 2. Load new data — either from per-sample TSVs or metrics + AMR
    # ------------------------------------------------------------------
    seen = set(old_amr_cols)
    new_rows = {}

    if args.sample_tsvs:
        # Simple path: each sample TSV already has all columns in one row
        tsv_rows, tsv_extra = load_sample_tsvs(args.sample_tsvs)
        for col in tsv_extra:
            if col not in seen:
                seen.add(col)
        new_rows = tsv_rows
    else:
        # Legacy path: separate metrics JSON + AMR TSV files
        new_metrics = {}
        for path in args.metrics:
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                data = json.load(fh)
            name = data.get('Name', '').strip()
            if name:
                new_metrics[name] = data

        new_amr = {}
        for path in args.amr:
            if not os.path.exists(path):
                continue
            basename = os.path.basename(path)
            for suffix in ('.amr.tsv', '.amr_combined.txt'):
                if basename.endswith(suffix):
                    basename = basename[:-len(suffix)]
                    break
            new_amr[basename] = load_amr_tsv(path)

        for sample_amr in new_amr.values():
            for gene in sample_amr:
                if gene not in seen:
                    seen.add(gene)

        for name, metrics in new_metrics.items():
            row = dict(metrics)
            row.update(new_amr.get(name, {}))
            new_rows[name] = row

    all_amr_genes = sorted(seen, key=lambda g: g.lower())

    # ------------------------------------------------------------------
    # 3. Merge: old rows first, new rows overwrite/append
    # ------------------------------------------------------------------
    merged = dict(old_rows)
    for name, row in new_rows.items():
        merged[name] = row

    # ------------------------------------------------------------------
    # 4. Write wide TSV
    # ------------------------------------------------------------------
    all_cols = FIXED_COLS + all_amr_genes
    with open(args.results_file, 'w', newline='', encoding='utf-8') as fh:
        writer = csv.DictWriter(fh, fieldnames=all_cols, delimiter='\t',
                                extrasaction='ignore', restval='',
                                lineterminator='\n')
        writer.writeheader()
        order = list(old_rows.keys())
        for name in merged:
            if name not in old_rows:
                order.append(name)
        for name in order:
            if name in merged:
                writer.writerow(merged[name])

    n_new = len(set(new_rows) - set(old_rows))
    n_amr = len(all_amr_genes)
    print(f"Results written: {len(merged)} samples, {n_amr} AMR gene columns "
          f"({n_new} new samples added).")

if __name__ == '__main__':
    main()
