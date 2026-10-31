#!/usr/bin/env python3
"""
Check the genus-organised FASTA tree against the results table.

For every sample, compares what is actually in the copied FASTA files with what
the table claims, so a truncated copy or a stale file is caught rather than
assumed correct:

  Chromosome_assemblies/<sample>.chromosome.fasta
      contig count   vs Chr_Contigs
      total bases    vs Genome_Size
      longest contig vs Max_Contig_Length

  Plasmid_assemblies/<sample>.plasmids.fasta
      element names  vs Plasmid_Names   (plasmids and rescued MGEs)
      rescued MGEs   vs Num_Rescued_MGE

  and that each file sits under the genus folder its row implies, using the
  same genus resolution organise_by_genus.py used to place it.

Contig headers are matched to element names from the table, which resolves the
one ambiguity in the naming: AA123_1/AA123_2 are one plasmid across two contigs,
whereas rescued_MGE_1/rescued_MGE_2 are two separate elements.

Usage:
  ./check_organised_fastas.py --by_genus by_genus --results_file results.tsv
"""
import argparse
import csv
import importlib.util
import os
import re
import sys
from collections import Counter


def load_safe_genus():
    """
    Reuse organise_by_genus.py's genus normalisation so 'correct folder' here
    means exactly what that script would have produced, not a second guess.
    """
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'organise_by_genus.py')
    if not os.path.isfile(path):
        return None
    spec = importlib.util.spec_from_file_location("obg", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.safe_genus


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--by_genus',     required=True, help='Genus-organised tree')
    p.add_argument('--results_file', required=True, help='Table to check against')
    p.add_argument('--report',       default='fasta_check_report.tsv')
    p.add_argument('--tolerance',    type=float, default=0.0,
                   help='Allowed fractional difference in total bases [0 = exact]')
    p.add_argument('--genus_from', default='assembly', choices=['assembly', 'reads'],
                   help='Genus column the tree was built from [assembly]')
    p.add_argument('--no_genus_check', action='store_true',
                   help='Skip checking that files sit in the right genus folder')
    p.add_argument('--quiet', action='store_true', help='Only print the summary')
    return p.parse_args()


def read_table(path):
    with open(path, newline='', encoding='utf-8') as fh:
        head = fh.readline()
        d = '\t' if head.count('\t') >= head.count(',') else ','
        fh.seek(0)
        return list(csv.DictReader(fh, delimiter=d))


def fasta_stats(path):
    """(n_records, total_bases, max_len, [record_ids]) without loading sequences."""
    n = total = longest = cur = 0
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith('>'):
                if n:
                    longest = max(longest, cur)
                    total += cur
                n += 1
                cur = 0
                ids.append(line[1:].split()[0])
            else:
                cur += len(line.strip())
    if n:
        longest = max(longest, cur)
        total += cur
    return n, total, longest, ids


def g(row, col):
    v = (row.get(col) or '').strip()
    return '' if v in ('NA', '') else v


def num(row, col):
    try:
        return float(g(row, col))
    except ValueError:
        return None


def match_elements(ids, expected, sample):
    """
    Map each contig header to one of the table's element names.
    A header is <sample>_<name>, <sample>_<name>_<n>, or the bare forms.
    Returns (found_names, unmatched_ids).
    """
    found, unmatched = set(), []
    for cid in ids:
        base = cid[len(sample) + 1:] if cid.startswith(sample + '_') else cid
        if base in expected:
            found.add(base)
            continue
        m = re.match(r'^(.*)_\d+$', base)          # AA123_2 -> AA123
        if m and m.group(1) in expected:
            found.add(m.group(1))
            continue
        unmatched.append(cid)
    return found, unmatched


def main():
    args = parse_args()
    if not os.path.isdir(args.by_genus):
        sys.exit(f"ERROR: not a directory: {args.by_genus}")
    rows = {r['Name'].strip(): r for r in read_table(args.results_file) if r.get('Name')}

    # locate every copied FASTA, whichever genus folder it landed in
    chrom, plasmid = {}, {}
    found_genus = {}
    for genus in sorted(os.listdir(args.by_genus)):
        gdir = os.path.join(args.by_genus, genus)
        if not os.path.isdir(gdir):
            continue
        for sub, store, suffix in (('Chromosome_assemblies', chrom, '.chromosome.fasta'),
                                   ('Plasmid_assemblies',   plasmid, '.plasmids.fasta')):
            d = os.path.join(gdir, sub)
            if not os.path.isdir(d):
                continue
            for fn in os.listdir(d):
                if fn.endswith(suffix):
                    sample_name = fn[:-len(suffix)]
                    store[sample_name] = os.path.join(d, fn)
                    found_genus.setdefault(sample_name, set()).add(genus)

    problems = Counter()
    lines = []
    checked = ok = 0

    safe_genus = None if args.no_genus_check else load_safe_genus()
    primary  = 'Assembly_Major_Genus' if args.genus_from == 'assembly' else 'Major_Genus_Reads'
    fallback = 'Major_Genus_Reads'    if args.genus_from == 'assembly' else 'Assembly_Major_Genus'

    for sample, row in sorted(rows.items()):
        issues = []

        if safe_genus and sample in found_genus:
            want = safe_genus(g(row, primary))
            if want == 'Unassigned':
                want = safe_genus(g(row, fallback))
            wrong = sorted(found_genus[sample] - {want})
            if wrong:
                issues.append(f"wrong_genus_folder:found={','.join(wrong)},expected={want}")

        cpath = chrom.get(sample)

        if cpath:
            n, total, longest, _ = fasta_stats(cpath)
            exp_n   = num(row, 'Chr_Contigs')
            exp_sz  = num(row, 'Genome_Size')
            exp_max = num(row, 'Max_Contig_Length')
            if exp_n is not None and n != int(exp_n):
                issues.append(f"chr_contigs:file={n},table={int(exp_n)}")
            if exp_sz is not None and exp_sz > 0:
                diff = abs(total - exp_sz) / exp_sz
                if diff > args.tolerance:
                    issues.append(f"genome_size:file={total},table={int(exp_sz)}")
            if exp_max is not None and longest != int(exp_max):
                issues.append(f"max_contig:file={longest},table={int(exp_max)}")
        elif g(row, 'Genome_Size'):
            issues.append("chromosome_fasta_missing")

        ppath = plasmid.get(sample)
        expected = {e.strip() for e in g(row, 'Plasmid_Names').split(';') if e.strip()}
        exp_mge = num(row, 'Num_Rescued_MGE') or 0
        if ppath:
            _, _, _, ids = fasta_stats(ppath)
            found, unmatched = match_elements(ids, expected, sample)
            missing = expected - found
            if missing:
                issues.append(f"plasmids_missing:{';'.join(sorted(missing))}")
            if unmatched:
                issues.append(f"unexpected_contigs:{';'.join(unmatched[:3])}")
            n_mge = sum(1 for e in found if e.startswith('rescued_MGE'))
            if n_mge != int(exp_mge):
                issues.append(f"rescued_mge:file={n_mge},table={int(exp_mge)}")
        elif expected:
            issues.append("plasmid_fasta_missing")

        checked += 1
        if issues:
            for i in issues:
                problems[i.split(':')[0]] += 1
            lines.append(f"{sample}\tFAIL\t{'; '.join(issues)}")
            if not args.quiet:
                print(f"  {sample}: {'; '.join(issues)}")
        else:
            ok += 1
            lines.append(f"{sample}\tOK\t")

    orphans = sorted((set(chrom) | set(plasmid)) - set(rows))
    with open(args.report, 'w') as fh:
        fh.write("Sample\tStatus\tIssues\n")
        fh.write("\n".join(lines) + "\n")
        for o in orphans:
            fh.write(f"{o}\tORPHAN\tfasta present but no row in the table\n")

    print(f"\nchecked {checked} samples: {ok} OK, {checked - ok} with issues")
    print(f"  chromosome FASTAs found: {len(chrom)}   plasmid FASTAs found: {len(plasmid)}")
    if problems:
        print("  issues by type:")
        for k, v in problems.most_common():
            print(f"     {k:28} {v}")
    if orphans:
        print(f"  FASTAs with no table row: {len(orphans)} {orphans[:5]}")
    print(f"  report: {args.report}")


if __name__ == '__main__':
    main()
