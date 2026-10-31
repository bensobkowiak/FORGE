#!/usr/bin/env python3
"""
Parse CheckM2, BUSCO and MLST outputs for a sample.
Usage: parse_genome_eval.py --checkm2 <quality_report.tsv>
                             --busco <short_summary.txt>
                             --mlst <mlst.tsv>
                             --assembly <assembly.fasta>
                             --assembly_info <assembly_info.txt>
Prints one TSV row:
  circular  total_contigs  completeness  contamination  genome_size
  n50  max_contig  avg_gene_len  total_cds  gc_content
  busco_complete  busco_single  busco_dup  busco_frag  busco_missing
  mlst_st
"""
import sys
import os
import argparse
import re
import gzip

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--checkm2')
    p.add_argument('--checkm2_full', default=None,
                   help='CheckM2 report on full assembly — contamination value overrides --checkm2')
    p.add_argument('--busco')
    p.add_argument('--mlst')
    p.add_argument('--assembly')
    p.add_argument('--assembly_info', help="Flye assembly_info.txt (per-contig circularity)")
    return p.parse_args()

def parse_checkm2(path):
    defaults = ['NA'] * 8
    if not path:
        return defaults
    with open(path) as fh:
        hdr = fh.readline().strip().split('\t')
        for line in fh:
            cols = line.strip().split('\t')
            row = dict(zip(hdr, cols))
            return [
                row.get('Additional_Notes', 'NA'),  # placeholder for genome size from fasta
                row.get('Completeness', 'NA'),
                row.get('Contamination', 'NA'),
                row.get('Genome_Size', 'NA'),
                row.get('Contig_N50', 'NA'),
                row.get('Max_Contig_Length', 'NA'),
                row.get('Average_Gene_Length', 'NA'),
                row.get('Total_Coding_Sequences', 'NA'),
                row.get('GC_Content', 'NA'),
                row.get('Total_Contigs', 'NA'),
            ]
    return defaults

def parse_busco(path):
    defaults = ['NA'] * 5
    if not path:
        return defaults
    patterns = {k: re.compile(rf'\b{k}:([\d.]+)%') for k in ['C', 'S', 'D', 'F', 'M']}
    vals = {k: 'NA' for k in patterns}
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            for k, pat in patterns.items():
                if vals[k] == 'NA':
                    m = pat.search(line)
                    if m:
                        vals[k] = m.group(1)
    return [vals[k] for k in ['C', 'S', 'D', 'F', 'M']]

def parse_mlst(path):
    if not path:
        return 'NA'
    with open(path) as fh:
        for line in fh:
            cols = line.strip().split('\t')
            if len(cols) >= 3:
                return cols[2]  # ST field
    return 'NA'

def parse_circular(chromosome_path, info_path):
    """
    Circularity of the chromosome, from Flye's assembly_info.txt (column 4,
    'circ.'), judged on the LONGEST chromosomal contig.

    An assembly can close the chromosome into a circle and still carry further
    chromosomal contigs — unplaced fragments, prophage, or a second replicon —
    so contig count alone does not decide it. Chr_Contigs reports the
    fragmentation separately; this column reports whether the chromosome itself
    was closed.

    This is the authoritative value; combine_results.py falls back to the
    alignment-based call from parse_depth.py only when this returns NA.
    """
    if not chromosome_path or not info_path or not os.path.isfile(info_path):
        return 'NA'
    try:
        with open(chromosome_path) as fh:
            names = [l[1:].split()[0] for l in fh if l.startswith('>')]
    except OSError:
        return 'NA'
    if not names:
        return 'NA'

    circ, length = {}, {}
    try:
        with open(info_path) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) >= 4:
                    name = cols[0].strip()
                    circ[name] = cols[3].strip().upper()
                    try:
                        length[name] = int(cols[1])
                    except ValueError:
                        pass
    except OSError:
        return 'NA'

    main = max(names, key=lambda c: length.get(c, 0))
    flag = circ.get(main)
    if flag is None:
        return 'NA'
    return 'Yes' if flag in ('Y', 'YES', '+', 'TRUE', '1') else 'No'

def main():
    args = parse_args()

    # CheckM2 output cols: Name, Completeness, Contamination, Completeness_Model_Used,
    #   Translation_Table_Used, Coding_Density, Contig_N50, Average_Gene_Length,
    #   Genome_Size, GC_Content, Total_Coding_Sequences, Additional_Notes, Total_Contigs, Max_Contig_Length
    checkm2_vals = ['NA'] * 9
    if args.checkm2:
        try:
            with open(args.checkm2) as fh:
                hdr = fh.readline().strip().split('\t')
                for line in fh:
                    cols = line.strip().split('\t')
                    row = dict(zip(hdr, cols))
                    checkm2_vals = [
                        row.get('Total_Contigs', 'NA'),
                        row.get('Completeness', 'NA'),
                        row.get('Contamination', 'NA'),
                        row.get('Genome_Size', 'NA'),
                        row.get('Contig_N50', 'NA'),
                        row.get('Max_Contig_Length', 'NA'),
                        row.get('Average_Gene_Length', 'NA'),
                        row.get('Total_Coding_Sequences', 'NA'),
                        row.get('GC_Content', 'NA'),
                    ]
                    break
        except Exception:
            pass

    # Override contamination with full-assembly CheckM2 value when available.
    # completeness/genome size/etc. come from the filtered chromosome CheckM2.
    if args.checkm2_full and os.path.isfile(args.checkm2_full):
        try:
            with open(args.checkm2_full) as fh:
                hdr = fh.readline().strip().split('\t')
                for line in fh:
                    cols = line.strip().split('\t')
                    row = dict(zip(hdr, cols))
                    checkm2_vals[2] = row.get('Contamination', checkm2_vals[2])
                    break
        except Exception:
            pass

    busco_vals = parse_busco(args.busco)
    mlst_st    = parse_mlst(args.mlst)
    circular   = parse_circular(args.assembly, args.assembly_info)

    out = [circular] + checkm2_vals + busco_vals + [mlst_st]
    print('\t'.join(str(v) for v in out))

if __name__ == '__main__':
    main()
