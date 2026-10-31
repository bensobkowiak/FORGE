#!/usr/bin/env python3
"""
Filter assembly FASTA to contigs classified as the MALDI genus by Kraken2.

Uses Kraken2 per-sequence output (--output) for per-contig taxid assignment
and the Kraken2 report (--report) to resolve which taxids fall within the
target genus clade (genus + all descendant species/strains).

Output flag:
  "Yes" - MALDI genus contigs are being used:
           either the assembly was already single-genus (MALDI),
           or non-MALDI contigs were stripped.
  "No"  - MALDI genus not found in assembly; unfiltered assembly kept.
  "Mismatch" - the MALDI genus is present but accounts for too little of the
           assembly to be the chromosome (below --min_genus_len or
           --min_genus_frac). Filtering to it would discard the genome and
           leave a fragment, so the unfiltered assembly is kept and the sample
           flagged for review: MALDI and sequence disagree.

Unclassified contigs are kept regardless (ambiguous, conservatively retained).

Off-genus contigs are split by length.  Below --max_len they are written to
the query-contig FASTA: too small to be a second chromosome, they may carry
resistance genes worth reporting, so they are screened for AMR but excluded
from both chromosome and plasmid metrics.  At or above --max_len they are
discarded as before, their lengths recorded in --out_removed.
"""
import sys
import os
import argparse
import shutil


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--assembly',      required=True)
    p.add_argument('--kraken_out',    required=True, help='Kraken2 --output per-contig file')
    p.add_argument('--kraken_report', required=True, help='Kraken2 --report file')
    p.add_argument('--maldi_genus',   required=True)
    p.add_argument('--out_fasta',     required=True)
    p.add_argument('--out_flag',      required=True)
    p.add_argument('--out_removed',   required=True,
                   help='File for lengths of removed non-MALDI contigs, ";"-separated')
    p.add_argument('--out_query',     required=True,
                   help='FASTA of off-genus contigs below --max_len (AMR screening only)')
    p.add_argument('--out_query_info', required=True,
                   help='Query contig count and lengths: "<n>\\t<len;len>"')
    p.add_argument('--max_len',       type=int, default=1000000,
                   help='Off-genus contigs at or above this length are discarded, '
                        'not treated as query contigs')
    p.add_argument('--min_genus_len',  type=int, default=1000000,
                   help='Minimum total length of MALDI-genus contigs for filtering '
                        'to proceed; below this the assembly is kept unfiltered '
                        'and flagged Mismatch')
    p.add_argument('--out_species',   required=True,
                   help='Species of the chromosome contigs carried forward, as '
                        '"Name" or "Name(n);Name(n)" when several are present')
    p.add_argument('--min_genus_frac', type=float, default=0.25,
                   help='Minimum fraction of assembly length that must be MALDI '
                        'genus for filtering to proceed')
    return p.parse_args()


def read_fasta(path):
    seqs = {}    # id → sequence
    headers = {} # id → full header line (without '>')
    name = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if name:
                    seqs[name] = ''.join(buf)
                full = line[1:]
                name = full.split()[0]
                headers[name] = full
                buf = []
            else:
                buf.append(line)
    if name:
        seqs[name] = ''.join(buf)
    return seqs, headers


def write_fasta(path, seqs, headers, order):
    with open(path, 'w') as fh:
        for name in order:
            if name not in seqs:
                continue
            fh.write(f">{headers.get(name, name)}\n")
            seq = seqs[name]
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + '\n')


def get_genus_taxids(report_path, maldi_genus):
    """
    Parse the Kraken2 report to collect all taxids in the MALDI genus clade.
    The report is depth-first; genus lines are rank 'G'/'G1'/'G2'.
    Descendants have greater indentation in the name column.
    """
    genus_taxids = set()
    in_target = False
    target_indent = -1
    maldi_lower = maldi_genus.strip().lower()

    with open(report_path) as fh:
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 6:
                continue
            rank  = cols[3].strip()
            taxid = cols[4].strip()
            name  = cols[5]
            indent = len(name) - len(name.lstrip(' '))

            if in_target:
                if indent <= target_indent:
                    in_target = False  # left the genus subtree
                else:
                    genus_taxids.add(taxid)
                    continue

            if rank.startswith('G') and maldi_lower in name.strip().lower():
                in_target = True
                target_indent = indent
                genus_taxids.add(taxid)

    return genus_taxids


def taxid_species(report_path):
    """
    taxid -> species name, using the report's indentation to find each taxid's
    enclosing species. A contig assigned to a strain (rank S1/S2) rolls up to
    its species rather than reporting the strain.
    """
    out, stack = {}, []
    with open(report_path) as fh:
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) < 6:
                continue
            rank, taxid, name = cols[3].strip(), cols[4].strip(), cols[5]
            indent = len(name) - len(name.lstrip(' '))
            while stack and stack[-1][0] >= indent:
                stack.pop()
            stack.append((indent, rank, name.strip()))
            sp = next((n for _i, r, n in reversed(stack) if r == 'S'), None)
            if sp:
                out[taxid] = sp
    return out


def summarise_species(contig_ids, contig_tax, tax_sp):
    """
    "Klebsiella pneumoniae" for a single species, or
    "Klebsiella pneumoniae(3);Escherichia coli(1)" when several are present.
    """
    if not contig_ids:
        return 'NA'
    counts = {}
    for c in contig_ids:
        taxid = contig_tax.get(c)
        name = tax_sp.get(taxid, 'Unclassified') if taxid else 'Unclassified'
        counts[name] = counts.get(name, 0) + 1
    if len(counts) == 1:
        return next(iter(counts))
    ordered = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return ';'.join(f"{n}({c})" for n, c in ordered)


def load_contig_taxids(kraken_out_path):
    """
    Parse Kraken2 per-sequence output.
    Returns dict: contig_id → taxid string, or None if unclassified.
    """
    contig_tax = {}
    with open(kraken_out_path) as fh:
        for line in fh:
            line = line.rstrip()
            if not line:
                continue
            cols = line.split('\t')
            if len(cols) < 3:
                continue
            status = cols[0]
            seq_id = cols[1].split()[0]
            taxid  = cols[2].strip()
            contig_tax[seq_id] = None if status == 'U' else taxid
    return contig_tax


def main():
    args = parse_args()

    def write_outputs(flag, removed='NA', query_ids=None, seqs=None, headers=None,
                      kept_ids=None, contig_tax=None, tax_sp=None):
        with open(args.out_flag, 'w') as fh:
            fh.write(flag + '\n')
        with open(args.out_removed, 'w') as fh:
            fh.write(removed + '\n')
        # Always write both query outputs, empty if there are none, so the
        # process has a stable output contract.
        query_ids = query_ids or []
        write_fasta(args.out_query, seqs or {}, headers or {}, order=query_ids)
        lens = ';'.join(str(len(seqs[c])) for c in query_ids) if query_ids else 'NA'
        with open(args.out_query_info, 'w') as fh:
            fh.write(f"{len(query_ids)}\t{lens}\n")
        with open(args.out_species, 'w') as fh:
            fh.write(summarise_species(kept_ids or [], contig_tax or {}, tax_sp or {}) + '\n')

    seqs, headers = read_fasta(args.assembly)
    contig_tax    = load_contig_taxids(args.kraken_out)
    tax_sp        = taxid_species(args.kraken_report)
    all_ids       = list(seqs)

    if not args.maldi_genus or args.maldi_genus in ('NA', '', 'None'):
        shutil.copy(args.assembly, args.out_fasta)
        write_outputs('No', kept_ids=all_ids, contig_tax=contig_tax, tax_sp=tax_sp)
        return

    genus_taxids = get_genus_taxids(args.kraken_report, args.maldi_genus)

    maldi_contigs = []
    other_contigs = []

    for cid in seqs:
        taxid = contig_tax.get(cid)
        if taxid is None:
            maldi_contigs.append(cid)  # unclassified — keep conservatively
        elif taxid in genus_taxids:
            maldi_contigs.append(cid)
        else:
            other_contigs.append(cid)

    if not maldi_contigs:
        # No MALDI genus contigs at all — can't filter, use full assembly
        shutil.copy(args.assembly, args.out_fasta)
        write_outputs('No', kept_ids=all_ids, contig_tax=contig_tax, tax_sp=tax_sp)
    elif not other_contigs:
        # Already single genus (MALDI) — no filtering needed
        shutil.copy(args.assembly, args.out_fasta)
        write_outputs('Yes', kept_ids=all_ids, contig_tax=contig_tax, tax_sp=tax_sp)
    elif (sum(len(seqs[c]) for c in maldi_contigs) < args.min_genus_len or
          sum(len(seqs[c]) for c in maldi_contigs) <
              args.min_genus_frac * sum(len(v) for v in seqs.values())):
        # The MALDI genus is present but trivial — e.g. one 5 kb contig beside a
        # 4.4 Mb chromosome of another genus. Filtering to it would throw the
        # genome away and leave a fragment on which completeness, depth and MLST
        # are meaningless. Keep the assembly intact and flag the disagreement.
        shutil.copy(args.assembly, args.out_fasta)
        write_outputs('Mismatch', kept_ids=all_ids, contig_tax=contig_tax, tax_sp=tax_sp)
    else:
        # Multi-species: keep only MALDI genus contigs on the chromosome.
        # Off-genus contigs below max_len become query contigs rather than
        # being discarded, so any resistance genes they carry are still found.
        write_fasta(args.out_fasta, seqs, headers, order=maldi_contigs)
        query_ids = [c for c in other_contigs if len(seqs[c]) < args.max_len]
        drop_ids  = [c for c in other_contigs if len(seqs[c]) >= args.max_len]
        removed_lens = sorted((len(seqs[c]) for c in drop_ids), reverse=True)
        write_outputs('Yes',
                      ';'.join(str(n) for n in removed_lens) if removed_lens else 'NA',
                      query_ids, seqs, headers,
                      kept_ids=maldi_contigs, contig_tax=contig_tax, tax_sp=tax_sp)


if __name__ == '__main__':
    main()
