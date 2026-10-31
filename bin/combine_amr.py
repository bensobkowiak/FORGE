#!/usr/bin/env python3
"""
Merge chromosome and plasmid AMR hits into a per-gene TSV.

Output (two columns, one row per unique gene):
  gene_name<TAB>location_string

where location_string shows every location in which the gene was found:
  chromosome(AMRfinder:98.5%,RGI:91.0%);AA002(AMRfinder:97.2%)

Plasmid location labels are the clean primary cluster IDs from contig_to_plasmid.tsv
(e.g. AA002, novel_1) with no sample-name prefix.

Usage:
  combine_amr.py --chr_amr   <sample.chr.amr_parsed.tsv>
                 --pla_amr   <sample.pla.amr_parsed.tsv>
                 --contig_map <contig_to_plasmid.tsv>
                 > <sample.amr.tsv>
"""
import sys
import os
import argparse
import csv
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--chr_amr',    required=True)
    p.add_argument('--pla_amr',    required=True)
    p.add_argument('--contig_map', required=True)
    p.add_argument('--qry_amr',    default=None,
                   help='AMR parsed from off-genus query contigs; all hits are\n'
                        'labelled "query" since these are neither chromosome nor plasmid')
    return p.parse_args()

def load_contig_map(path):
    mapping = {}
    if not os.path.exists(path):
        return mapping
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                mapping[parts[0].strip()] = parts[1].strip()
    return mapping

def lookup_location(contig_map, contig):
    """
    Return the plasmid name for a contig.

    Tries exact match first. Falls back to longest-prefix match to handle
    tools like RGI that append positional info to contig IDs
    (e.g. NODE_5_length_..._3000_3500_- where the key is NODE_5_length_...).
    """
    if contig in contig_map:
        return contig_map[contig]
    best_key = None
    for key in contig_map:
        if contig.startswith(key) and (best_key is None or len(key) > len(best_key)):
            best_key = key
    return contig_map[best_key] if best_key else 'plasmid_unknown'

def load_amr_tsv(path, source_override=None, contig_map=None):
    hits = []
    if not path or not os.path.exists(path):
        return hits
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            gene     = row.get('gene', '').strip()
            identity = row.get('identity', 'NA').strip()
            contig   = row.get('contig', 'NA').strip()
            tool     = row.get('tool', 'NA').strip()
            if not gene or gene in ('NA', 'gene'):
                continue
            if source_override:
                location = source_override
            elif contig_map:
                location = lookup_location(contig_map, contig)
            else:
                location = 'plasmid_unknown'
            hits.append({'gene': gene, 'identity': identity, 'location': location, 'tool': tool})
    return hits

def main():
    args  = parse_args()
    cmap  = load_contig_map(args.contig_map)
    chr_h = load_amr_tsv(args.chr_amr, source_override='chromosome')
    pla_h = load_amr_tsv(args.pla_amr, contig_map=cmap)
    qry_h = load_amr_tsv(args.qry_amr, source_override='query') if args.qry_amr else []
    all_h = chr_h + pla_h + qry_h

    print("gene\tvalue")
    if not all_h:
        return

    canonical = {}
    grouped = defaultdict(lambda: defaultdict(dict))
    for h in all_h:
        try:
            pct_val = float(h['identity'])
        except (ValueError, TypeError):
            pct_val = None
        gene, loc, tool = h['gene'], h['location'], h['tool']
        norm = gene.lower()
        if norm not in canonical:
            canonical[norm] = gene
        prev = grouped[norm][loc].get(tool)
        if prev is None or (pct_val is not None and (prev[0] is None or pct_val > prev[0])):
            grouped[norm][loc][tool] = (pct_val, h['identity'])

    for norm in sorted(grouped):
        locs = grouped[norm]
        order = sorted(locs, key=lambda l: (0 if l == 'chromosome' else
                                           2 if l == 'query' else 1, l))
        parts = []
        for loc in order:
            hits_str = ','.join(
                f"{t}:{f'{v:.1f}%' if v is not None else raw}"
                for t, (v, raw) in sorted(locs[loc].items())
            )
            parts.append(f"{loc}({hits_str})")
        print(f"{canonical[norm]}\t{';'.join(parts)}")

if __name__ == '__main__':
    main()
