#!/usr/bin/env python3
"""
Parse samtools coverage output + contig_to_plasmid.tsv to report
length-weighted mean depth per location.

Optionally detects chromosome circularity from a BAM by scanning primary
alignments for SA tags that span the chromosome start/end junction.
No BAM index required (full scan; bacterial genomes are small enough).

Output: three tab-separated fields on one line:
  chr_depth   plasmid_depths   circularised

where plasmid_depths is "name1:depth1;name2:depth2" sorted by depth desc, or NA,
and circularised is Yes/No/NA.  Sample-name prefix is stripped from plasmid
names when --sample is supplied.
"""
import sys
import subprocess
import argparse


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--coverage',   required=True, help='samtools coverage output')
    p.add_argument('--contig_map', required=True, help='contig_to_plasmid.tsv')
    p.add_argument('--bam',        default=None,  help='BAM for circularity detection')
    p.add_argument('--sample',     default=None,  help='Sample name; prefix stripped from plasmid names')
    return p.parse_args()


def detect_circular(bam_path, contig_loc, junction_window=1000, min_reads=3,
                    min_chr_length=200_000):
    """
    Scan all primary alignments for SA tags where one part is near the
    chromosome end and the supplementary part is near the chromosome start
    (or vice versa).  No index required.

    Judged on the LONGEST chromosomal contig only: an assembly can close the
    chromosome and still carry further chromosomal contigs (unplaced fragments,
    prophage, a second replicon), so contig count alone does not decide it.
    """
    chrom_contigs = {c for c, loc in contig_loc.items() if loc == 'chromosome'}
    if not chrom_contigs:
        return 'NA'

    try:
        header_out = subprocess.run(
            ['samtools', 'view', '-H', bam_path],
            capture_output=True, text=True, check=True
        ).stdout
    except Exception:
        return 'NA'

    ctg_len = {}
    for line in header_out.split('\n'):
        if not line.startswith('@SQ'):
            continue
        parts = dict(p.split(':', 1) for p in line.split('\t')[1:] if ':' in p)
        name = parts.get('SN', '')
        if name:
            try:
                ctg_len[name] = int(parts.get('LN', 0))
            except (ValueError, TypeError):
                pass

    # Restrict to the chromosome proper
    main_ctg = max(chrom_contigs, key=lambda c: ctg_len.get(c, 0))
    chrom_contigs = {main_ctg}
    junction_counts = {main_ctg: 0}

    # -F 2048 excludes supplementary alignments; we want primaries which carry SA tags
    try:
        view_out = subprocess.run(
            ['samtools', 'view', '-F', '2048', bam_path],
            capture_output=True, text=True, check=True
        ).stdout
    except Exception:
        return 'NA'

    for aln_line in view_out.split('\n'):
        if not aln_line.strip():
            continue
        fields = aln_line.split('\t')
        if len(fields) < 12:
            continue

        primary_ctg = fields[2]
        if primary_ctg not in chrom_contigs:
            continue

        try:
            primary_pos = int(fields[3])
        except ValueError:
            continue

        ctg_length = ctg_len.get(primary_ctg, 0)
        if ctg_length < max(2 * junction_window, min_chr_length):
            continue

        near_end   = primary_pos >= ctg_length - junction_window
        near_start = primary_pos <= junction_window
        if not (near_end or near_start):
            continue

        for tag in fields[11:]:
            if not tag.startswith('SA:Z:'):
                continue
            for sa in tag[5:].rstrip(';').split(';'):
                sa_parts = sa.split(',')
                if len(sa_parts) < 2:
                    continue
                sa_ctg = sa_parts[0]
                if sa_ctg != primary_ctg:
                    continue
                try:
                    sa_pos = int(sa_parts[1])
                except ValueError:
                    continue
                if ((near_end   and sa_pos <= junction_window) or
                        (near_start and sa_pos >= ctg_length - junction_window)):
                    junction_counts[primary_ctg] += 1
                    break

        if junction_counts.get(primary_ctg, 0) >= min_reads:
            return 'Yes'

    return 'No'


def main():
    args = parse_args()
    prefix = f"{args.sample}_" if args.sample else None

    contig_loc = {}
    with open(args.contig_map) as fh:
        for line in fh:
            parts = line.rstrip('\n').split('\t')
            if len(parts) >= 2:
                contig_loc[parts[0].strip()] = parts[1].strip()

    loc_bases     = {}
    loc_depth_sum = {}

    with open(args.coverage) as fh:
        for line in fh:
            if line.startswith('#'):
                continue
            cols = line.strip().split('\t')
            if len(cols) < 7 or cols[0] == 'rname':
                continue
            contig = cols[0]
            try:
                length    = int(cols[2])
                meandepth = float(cols[6])
            except (ValueError, IndexError):
                continue
            loc = contig_loc.get(contig, 'unknown')
            loc_bases[loc]      = loc_bases.get(loc, 0) + length
            loc_depth_sum[loc]  = loc_depth_sum.get(loc, 0.0) + meandepth * length

    loc_depth = {
        loc: loc_depth_sum[loc] / loc_bases[loc]
        for loc in loc_bases if loc_bases[loc] > 0
    }

    chr_depth = f"{loc_depth['chromosome']:.1f}" if 'chromosome' in loc_depth else 'NA'

    def display_name(loc):
        if prefix and loc.startswith(prefix):
            return loc[len(prefix):]
        return loc

    pla_parts = [
        f"{display_name(loc)}:{depth:.1f}"
        for loc, depth in sorted(
            ((l, d) for l, d in loc_depth.items() if l != 'chromosome'),
            key=lambda x: -x[1]
        )
    ]
    pla_str = ';'.join(pla_parts) if pla_parts else 'NA'

    circular = detect_circular(args.bam, contig_loc) if args.bam else 'NA'

    print(f"{chr_depth}\t{pla_str}\t{circular}")


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print(f"parse_depth.py error: {exc}", file=sys.stderr)
        import traceback; traceback.print_exc(file=sys.stderr)
        print("NA\tNA\tNA")
