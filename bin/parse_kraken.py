#!/usr/bin/env python3
"""
Parse Kraken2 output and emit:
  major_genus  pct_major_genus  genera_gt5pct \
  major_species  pct_major_species  species_gt5pct  pct_unclassified

genera_gt5pct / species_gt5pct are semicolon-separated "Name:pct" pairs.

With only a report (reads mode), percentages are Kraken2's own, i.e. the
fraction of SEQUENCES assigned to each taxon.

For an assembly, --kraken_out and --assembly switch to percentages weighted by
contig LENGTH. Sequence-weighted figures are badly misleading for assemblies:
a sample whose 4.4 Mb chromosome is one contig and which carries a single 5 kb
foreign fragment reports as 50/50 by contig count, but is 99.9/0.1 by bases.
"""
import argparse
import sys


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('report')
    p.add_argument('--kraken_out', default=None,
                   help='Kraken2 --output per-sequence file (assembly mode)')
    p.add_argument('--assembly',   default=None,
                   help='FASTA that was classified; used for contig lengths')
    return p.parse_args()


def emit(genera, species, unclassified_pct):
    if not genera:
        print(f"NA\tNA\tNA\tNA\tNA\tNA\t{unclassified_pct:.1f}")
        return
    major_g = max(genera, key=genera.get)
    gt5 = [f"{g}:{p:.1f}%" for g, p in sorted(genera.items(), key=lambda x: -x[1])
           if p >= 5.0 and g != major_g]
    g_str = ";".join(gt5) if gt5 else "NA"
    if species:
        major_s = max(species, key=species.get)
        s_gt5 = [f"{s}:{p:.1f}%" for s, p in sorted(species.items(), key=lambda x: -x[1])
                 if p >= 5.0 and s != major_s]
        s_str = f"{major_s}\t{species[major_s]:.1f}\t{';'.join(s_gt5) if s_gt5 else 'NA'}"
    else:
        s_str = "NA\tNA\tNA"
    print(f"{major_g}\t{genera[major_g]:.1f}\t{g_str}\t{s_str}\t{unclassified_pct:.1f}")


def parse_report(report_file):
    """Kraken2's own sequence-weighted percentages."""
    genera, species, unclassified_pct = {}, {}, 0.0
    with open(report_file) as fh:
        for line in fh:
            cols = line.strip().split('\t')
            if len(cols) < 6:
                continue
            pct, rank, name = float(cols[0]), cols[3], cols[5].strip()
            if rank == 'U':
                unclassified_pct = pct
            elif rank == 'G':
                genera[name] = pct
            elif rank == 'S':
                species[name] = pct
    emit(genera, species, unclassified_pct)


def taxid_lineage(report_file):
    """
    taxid -> (genus, species) from the report's depth-first structure.
    Indentation in the name column encodes depth, so the nearest enclosing
    G / S ancestor gives each taxid its genus and species.
    """
    lineage = {}
    stack = []          # [(indent, rank, name)]
    with open(report_file) as fh:
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 6:
                continue
            rank, taxid, name = cols[3].strip(), cols[4].strip(), cols[5]
            indent = len(name) - len(name.lstrip(' '))
            name = name.strip()
            while stack and stack[-1][0] >= indent:
                stack.pop()
            stack.append((indent, rank, name))
            genus = next((n for _i, r, n in reversed(stack) if r == 'G'), None)
            sp    = next((n for _i, r, n in reversed(stack) if r == 'S'), None)
            lineage[taxid] = (genus, sp)
    return lineage


def contig_lengths(fasta):
    lengths, name, n = {}, None, 0
    with open(fasta) as fh:
        for line in fh:
            if line.startswith('>'):
                if name:
                    lengths[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line.strip())
    if name:
        lengths[name] = n
    return lengths


def parse_length_weighted(report_file, kraken_out, fasta):
    lineage = taxid_lineage(report_file)
    lengths = contig_lengths(fasta)
    g_bases, s_bases, unclassified, total = {}, {}, 0, 0
    with open(kraken_out) as fh:
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 3:
                continue
            status, seq_id, taxid = cols[0], cols[1].split()[0], cols[2].strip()
            ln = lengths.get(seq_id, 0)
            if not ln:
                continue
            total += ln
            if status == 'U':
                unclassified += ln
                continue
            genus, sp = lineage.get(taxid, (None, None))
            if genus:
                g_bases[genus] = g_bases.get(genus, 0) + ln
            if sp:
                s_bases[sp] = s_bases.get(sp, 0) + ln
    if not total:
        print("NA\tNA\tNA\tNA\tNA\tNA\t100.0")
        return
    pct = lambda d: {k: 100.0 * v / total for k, v in d.items()}
    emit(pct(g_bases), pct(s_bases), 100.0 * unclassified / total)


def main():
    args = parse_args()
    if args.kraken_out and args.assembly:
        parse_length_weighted(args.report, args.kraken_out, args.assembly)
    else:
        parse_report(args.report)


if __name__ == '__main__':
    main()
