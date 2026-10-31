#!/usr/bin/env python3
"""
Reorganise the pipeline's per-sample results into genus-level folders.

The pipeline writes results/<sample>/..., with many files carrying generic
names (assembly.fasta, contigs_filtered.fasta, contig_report.txt). This copies
them into a per-genus tree, renaming every generic file with the sample name so
the files stay identifiable once pooled:

    <outdir>/<Genus>/
        Full_assembly/          <sample>.assembly.fasta
        Chromosome_assemblies/  <sample>.chromosome.fasta
        Plasmid_assemblies/     <sample>.plasmids.fasta
                                (all plasmids and rescued MGEs in one file,
                                 each contig renamed to its plasmid/MGE)
        Reports/                <sample>.mob_contig_report.txt
                                <sample>.mobtyper_results.txt
                                <sample>.amr.tsv
                                <sample>.kleborate.txt      (Klebsiella only)

Samples are whatever directories exist under --results. Genus is read from each
sample's own files, so no results table is needed: <sample>.metrics.json first,
then the Kraken2 parsed files. Pass --results_file to take genus from a wide
results table instead (tab- or comma-separated; detected from the header).
Files that a sample does not have are skipped and counted, so a partial run
reorganises cleanly rather than failing.

Usage:
    ./organise_by_genus.py --results results --results_file results.tsv \
                           --outdir by_genus [--dry-run] [--link] [--include-query]
"""
import argparse
import csv
import json
import os
import shutil
import sys
from collections import Counter, defaultdict

SUBDIRS = ['Full_assembly', 'Chromosome_assemblies', 'Plasmid_assemblies', 'Reports']


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--results',     required=True, help='Pipeline results directory')
    p.add_argument('--results_file', '--results_tsv', dest='results_file', default=None,
                   metavar='FILE',
                   help='Optional wide results table (.tsv or .csv) to take the genus '
                        'from. Without it, genus is read from each sample directory. '
                        '--results_tsv is accepted as an alias.')
    p.add_argument('--outdir',      required=True, help='Destination for the genus tree')
    p.add_argument('--genus_from',  default='assembly', choices=['assembly', 'reads'],
                   help='Genus source column: Assembly_Major_Genus (default) or '
                        'Major_Genus_Reads. Falls back to the other if blank.')
    p.add_argument('--include_query', action='store_true',
                   help='Also copy query contigs into a Query_contigs/ subfolder. '
                        'These are neither chromosome nor plasmid, so they are '
                        'excluded unless asked for.')
    p.add_argument('--bare_headers', action='store_true',
                   help='Name contigs AA123 / rescued_MGE_1 without the sample '
                        'prefix. Default keeps <sample>_<element> so headers stay '
                        'unique if the per-genus FASTAs are concatenated.')
    p.add_argument('--overwrite', action='store_true',
                   help='Replace files that already exist in the destination. '
                        'Without it existing files are left alone and counted, '
                        'so re-running only fills in what is missing.')
    p.add_argument('--link', action='store_true',
                   help='Hard-link instead of copying (same filesystem only; '
                        'saves space for large assemblies)')
    p.add_argument('--dry-run', dest='dry_run', action='store_true',
                   help='Report what would happen without writing anything')
    return p.parse_args()


def read_table(path):
    """
    Read the wide results table as dicts. The delimiter is taken from the header
    line rather than the file extension, since data rows legitimately contain
    commas (AMR columns look like "chromosome(AMRfinder:98.5%,RGI:91.0%)") and
    would mislead a naive sniff.
    """
    with open(path, newline='', encoding='utf-8') as fh:
        header = fh.readline()
        if not header:
            return []
        delim = '\t' if header.count('\t') >= header.count(',') else ','
        fh.seek(0)
        return list(csv.DictReader(fh, delimiter=delim)), delim


def genus_from_sample_dir(sdir, sample, prefer_reads=False):
    """
    Genus for one sample from its own published files, in order of reliability:
    the metrics JSON, then the Kraken2 parsed summaries (field 1 is the major
    genus). Returns '' if nothing is readable.
    """
    keys = (['Major_Genus_Reads', 'Assembly_Major_Genus'] if prefer_reads
            else ['Assembly_Major_Genus', 'Major_Genus_Reads'])
    mj = os.path.join(sdir, f'{sample}.metrics.json')
    if os.path.isfile(mj):
        try:
            with open(mj) as fh:
                data = json.load(fh)
            for k in keys:
                v = (data.get(k) or '').strip()
                if v and v.upper() != 'NA':
                    return v
        except (OSError, ValueError):
            pass
    parsed = ([os.path.join(sdir, 'reads', 'kraken2', f'{sample}_reads.kraken2.parsed'),
               os.path.join(sdir, 'chromosome', 'kraken2', f'{sample}_assembly.kraken2.parsed')]
              if prefer_reads else
              [os.path.join(sdir, 'chromosome', 'kraken2', f'{sample}_assembly.kraken2.parsed'),
               os.path.join(sdir, 'reads', 'kraken2', f'{sample}_reads.kraken2.parsed')])
    for path in parsed:
        if os.path.isfile(path):
            try:
                with open(path) as fh:
                    first = fh.readline().split('\t')
                if first and first[0].strip() and first[0].strip().upper() != 'NA':
                    return first[0].strip()
            except OSError:
                pass
    return ''


def safe_genus(value):
    """Folder-safe genus name; anything missing becomes Unassigned."""
    v = (value or '').strip()
    if not v or v.upper() in ('NA', 'NAN', 'NONE', '-'):
        return 'Unassigned'
    # Kraken2 can return multi-word or punctuated names
    v = v.split()[0]
    return ''.join(c if (c.isalnum() or c in '-_') else '_' for c in v)


def read_fasta(path):
    """Yield (contig_id, description, sequence_lines) in file order."""
    cid, desc, buf = None, '', []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if cid:
                    yield cid, desc, buf
                parts = line[1:].split(None, 1)
                cid  = parts[0] if parts else ''
                desc = parts[1] if len(parts) > 1 else ''
                buf = []
            else:
                buf.append(line)
    if cid:
        yield cid, desc, buf


def load_contig_plasmid_map(path):
    """contig_id -> plasmid/MGE name, from contig_to_plasmid.tsv ('chromosome' skipped)."""
    out = {}
    if not os.path.isfile(path):
        return out
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if len(cols) >= 2 and cols[1].strip() and cols[1].strip() != 'chromosome':
                out[cols[0].strip()] = cols[1].strip()
    return out


def write_renamed_plasmids(sdir, sample, dst, args, stats):
    """
    One FASTA per sample holding every plasmid and rescued MGE contig, with each
    contig renamed after the element it belongs to (AA123, rescued_MGE_1, ...)
    rather than keeping Flye's contig_N. Where an element spans several contigs
    they are numbered. The original contig id is kept in the description so a
    sequence can still be traced back to the assembly.
    """
    contigs_dir = os.path.join(sdir, 'assembly', 'contigs')
    src = os.path.join(contigs_dir, 'plasmids_all.fasta')
    if not os.path.isfile(src):
        stats['missing']['plasmids'] += 1
        return
    if os.path.getsize(src) == 0:
        stats['empty']['plasmids'] += 1
        return

    cmap = load_contig_plasmid_map(os.path.join(contigs_dir, 'contig_to_plasmid.tsv'))
    if not cmap:
        # No map (older run): copy unchanged rather than invent names
        stats['copied']['plasmids_unrenamed'] += 1
        transfer(src, dst, args, stats, 'plasmids')
        return

    records = list(read_fasta(src))
    counts = Counter(cmap.get(cid, 'unassigned_plasmid') for cid, _, _ in records)
    seen = Counter()
    lines, renamed = [], 0
    for cid, _desc, seq in records:
        element = cmap.get(cid, 'unassigned_plasmid')
        seen[element] += 1
        name = element if counts[element] == 1 else f"{element}_{seen[element]}"
        if not args.bare_headers:
            name = f"{sample}_{name}"
        lines.append(f">{name} {cid}")
        lines.extend(seq)
        renamed += 1
    if not renamed:
        stats['empty']['plasmids'] += 1
        return
    if os.path.exists(dst) and not args.overwrite:
        stats['skipped']['plasmids'] += 1
        return
    stats['copied']['plasmids'] += 1
    if args.dry_run:
        return
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, 'w') as fh:
        fh.write('\n'.join(lines) + '\n')


def transfer(src, dst, args, stats, label):
    """Copy/link one file, skipping missing and empty sources."""
    if not os.path.isfile(src):
        stats['missing'][label] += 1
        return False
    if os.path.getsize(src) == 0:
        stats['empty'][label] += 1
        return False
    if os.path.exists(dst) and not args.overwrite:
        stats['skipped'][label] += 1
        return False
    if args.dry_run:
        stats['copied'][label] += 1
        return True
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    if os.path.exists(dst):
        os.remove(dst)
    if args.link:
        try:
            os.link(src, dst)
        except OSError:
            shutil.copy2(src, dst)   # cross-device: fall back to a copy
    else:
        shutil.copy2(src, dst)
    stats['copied'][label] += 1
    return True


def main():
    args = parse_args()
    if not os.path.isdir(args.results):
        sys.exit(f"ERROR: not a directory: {args.results}")
    if args.results_file and not os.path.isfile(args.results_file):
        sys.exit(f"ERROR: results file not found: {args.results_file}")

    primary  = 'Assembly_Major_Genus' if args.genus_from == 'assembly' else 'Major_Genus_Reads'
    fallback = 'Major_Genus_Reads'    if args.genus_from == 'assembly' else 'Assembly_Major_Genus'

    stats = {'copied': Counter(), 'missing': Counter(), 'empty': Counter(),
             'skipped': Counter()}
    per_genus = Counter()
    no_dir = []

    if args.results_file:
        rows, delim = read_table(args.results_file)
        if not rows:
            sys.exit(f"ERROR: no data rows in {args.results_file}")
        if 'Name' not in (rows[0].keys()):
            sys.exit(f"ERROR: no 'Name' column in {args.results_file} "
                     f"(read as {'tab' if delim == chr(9) else 'comma'}-separated); "
                     f"found: {', '.join(list(rows[0].keys())[:6])}")
        print(f"Samples from {args.results_file} ({len(rows)})")
    else:
        # No table: every directory under --results is a sample, and genus is
        # read from that sample's own files.
        rows = [{'Name': d} for d in sorted(os.listdir(args.results))
                if os.path.isdir(os.path.join(args.results, d))]
        if not rows:
            sys.exit(f"ERROR: no sample directories under {args.results}")
        print(f"Samples from {args.results} ({len(rows)} directories)")

    for row in rows:
        sample = (row.get('Name') or '').strip()
        if not sample:
            continue
        sdir = os.path.join(args.results, sample)
        if not os.path.isdir(sdir):
            no_dir.append(sample)
            continue

        if args.results_file:
            genus = safe_genus(row.get(primary))
            if genus == 'Unassigned':
                genus = safe_genus(row.get(fallback))
        else:
            genus = safe_genus(genus_from_sample_dir(
                sdir, sample, prefer_reads=(args.genus_from == 'reads')))
        per_genus[genus] += 1
        gdir = os.path.join(args.outdir, genus)

        # --- Full assembly (generic name) ---
        transfer(os.path.join(sdir, 'assembly', 'flye', 'assembly.fasta'),
                 os.path.join(gdir, 'Full_assembly', f'{sample}.assembly.fasta'),
                 args, stats, 'full_assembly')

        # --- Chromosome: prefer the genus-filtered version actually analysed ---
        chrom = os.path.join(sdir, 'chromosome', 'contigs_filtered.fasta')
        if not os.path.isfile(chrom):
            chrom = os.path.join(sdir, 'assembly', 'contigs', 'chromosome.fasta')
        transfer(chrom,
                 os.path.join(gdir, 'Chromosome_assemblies', f'{sample}.chromosome.fasta'),
                 args, stats, 'chromosome')

        # --- Plasmids + rescued MGEs: one file, contigs named by element ---
        write_renamed_plasmids(
            sdir, sample,
            os.path.join(gdir, 'Plasmid_assemblies', f'{sample}.plasmids.fasta'),
            args, stats)

        # --- Optional: query contigs (neither chromosome nor plasmid) ---
        if args.include_query:
            transfer(os.path.join(sdir, 'chromosome', 'query_contigs.fasta'),
                     os.path.join(gdir, 'Query_contigs', f'{sample}.query_contigs.fasta'),
                     args, stats, 'query_contigs')

        # --- Reports ---
        mob = os.path.join(sdir, 'assembly', 'mob_recon', 'mob_recon_out')
        transfer(os.path.join(mob, 'contig_report.txt'),
                 os.path.join(gdir, 'Reports', f'{sample}.mob_contig_report.txt'),
                 args, stats, 'mob_contig_report')
        transfer(os.path.join(mob, 'mobtyper_results.txt'),
                 os.path.join(gdir, 'Reports', f'{sample}.mobtyper_results.txt'),
                 args, stats, 'mobtyper_results')
        transfer(os.path.join(sdir, f'{sample}.amr.tsv'),
                 os.path.join(gdir, 'Reports', f'{sample}.amr.tsv'),
                 args, stats, 'amr_report')
        # Kleborate exists only for Klebsiella samples
        transfer(os.path.join(sdir, 'chromosome', 'kleborate', f'{sample}.kleborate.txt'),
                 os.path.join(gdir, 'Reports', f'{sample}.kleborate.txt'),
                 args, stats, 'kleborate')

    # Create the standard subfolders even where a genus had nothing for one
    if not args.dry_run:
        for genus in per_genus:
            for sub in SUBDIRS:
                os.makedirs(os.path.join(args.outdir, genus, sub), exist_ok=True)

    print(f"{'DRY RUN - nothing written' if args.dry_run else 'Done'}: "
          f"{sum(per_genus.values())} samples into {len(per_genus)} genus folders\n")
    print("Samples per genus:")
    for g, n in sorted(per_genus.items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"  {g:28} {n}")
    print("\nFiles transferred:")
    for k, n in sorted(stats['copied'].items()):
        print(f"  {k:28} {n}")
    if stats['missing']:
        print("\nNot found (sample lacked the file):")
        for k, n in sorted(stats['missing'].items()):
            print(f"  {k:28} {n}")
    if stats['empty']:
        print("\nSkipped because empty:")
        for k, n in sorted(stats['empty'].items()):
            print(f"  {k:28} {n}")
    if stats['skipped']:
        print("\nAlready present, left unchanged (use --overwrite to replace):")
        for k, n in sorted(stats['skipped'].items()):
            print(f"  {k:28} {n}")
    if no_dir:
        print(f"\nIn the results file but no results/ directory ({len(no_dir)}): "
              f"{', '.join(no_dir[:10])}{' ...' if len(no_dir) > 10 else ''}")


if __name__ == '__main__':
    main()
