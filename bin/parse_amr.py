#!/usr/bin/env python3
"""
Parse AMRfinder, ResFinder, RGI, and CefiderocolFinder outputs into a
structured TSV with columns: gene  identity  coverage  contig  tool

This intermediate format is consumed by combine_amr.py which merges
chromosome and plasmid hits into the final combined AMR column.

Usage: parse_amr.py --amrfinder <file> --resfinder <dir>
                    --rgi <file> --cefiderocol <file>
                    > <sample_source.amr_parsed.tsv>
"""
import sys
import os
import argparse
import csv

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--amrfinder')
    p.add_argument('--resfinder')
    p.add_argument('--rgi')
    p.add_argument('--cefiderocol')
    return p.parse_args()

def emit(gene, identity, coverage, contig, tool):
    gene     = (gene     or 'NA').strip()
    identity = (identity or 'NA').strip().replace('%', '')
    coverage = (coverage or 'NA').strip().replace('%', '')
    contig   = (contig   or 'NA').strip()
    if gene and gene not in ('NA', '-', '.'):
        print(f"{gene}\t{identity}\t{coverage}\t{contig}\t{tool}")

def parse_amrfinder(path):
    if not path or not os.path.exists(path):
        return
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            raw_id = (row.get('% Identity to reference sequence') or
                      row.get('% Identity to reference') or
                      row.get('% Identity') or '').strip()
            method = row.get('Method', '').strip()
            if not raw_id:
                raw_id = method if method else 'NA'
            emit(
                row.get('Gene symbol') or row.get('Element symbol', ''),
                raw_id,
                row.get('% Coverage of reference sequence', 'NA'),
                row.get('Contig id', 'NA'),
                'AMRfinder'
            )

def parse_resfinder(res_dir):
    if not res_dir or not os.path.exists(res_dir):
        return
    result_tab = os.path.join(res_dir, 'ResFinder_results_tab.txt')
    if not os.path.exists(result_tab):
        return
    with open(result_tab) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            # Coverage expressed as "matched_len/ref_len"
            cov_raw = row.get('Query / Template length', 'NA')
            try:
                a, b = cov_raw.split('/')
                cov = f"{int(a)/int(b)*100:.1f}"
            except Exception:
                cov = cov_raw
            emit(
                row.get('Resistance gene', ''),
                row.get('Identity', 'NA'),
                cov,
                row.get('Contig', 'NA'),
                'ResFinder'
            )

def parse_rgi(path):
    if not path or not os.path.exists(path):
        return
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            if row.get('Cut_Off', '').strip().lower() not in ('strict', 'perfect'):
                continue
            emit(
                row.get('Best_Hit_ARO') or row.get('ARO', ''),
                row.get('Best_Identities', 'NA'),
                row.get('Percentage Length of Reference Sequence', 'NA'),
                row.get('Contig', 'NA'),
                'RGI'
            )

def parse_cefiderocol(path):
    if not path or not os.path.exists(path):
        return
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            cols = line.split('\t')
            # Expected format: gene  identity  coverage  [contig]
            gene     = cols[0] if len(cols) > 0 else ''
            identity = cols[1] if len(cols) > 1 else 'NA'
            coverage = cols[2] if len(cols) > 2 else 'NA'
            contig   = cols[3] if len(cols) > 3 else 'NA'
            emit(gene, identity, coverage, contig, 'CefiderocolFinder')

def main():
    args = parse_args()
    print("gene\tidentity\tcoverage\tcontig\ttool")
    parse_amrfinder(args.amrfinder)
    parse_resfinder(args.resfinder)
    parse_rgi(args.rgi)
    parse_cefiderocol(args.cefiderocol)

if __name__ == '__main__':
    main()
