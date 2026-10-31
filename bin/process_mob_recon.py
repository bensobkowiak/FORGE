#!/usr/bin/env python3
"""
Post-process mob_recon output:
  - Read contig_report.txt to map contigs to chromosome/plasmid
  - Write chromosome.fasta
  - Write one plasmid fasta per unique plasmid cluster (named or novel)
  - Emit plasmid summary: count, names, species, secondary mob clusters, inc/replicon types
Usage: process_mob_recon.py <mob_recon_dir> <assembly_fasta> <out_dir> [sample_name]
"""
import sys
import os
import re
from collections import defaultdict

def read_fasta(path):
    seqs = {}     # contig_id -> sequence
    headers = {}  # contig_id -> full header (without '>')
    contig_id = None
    buf = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith('>'):
                if contig_id:
                    seqs[contig_id] = ''.join(buf)
                full = line[1:]
                contig_id = full.split()[0]
                headers[contig_id] = full
                buf = []
            else:
                buf.append(line)
    if contig_id:
        seqs[contig_id] = ''.join(buf)
    return seqs, headers

def write_fasta(path, seqs, headers):
    with open(path, 'w') as fh:
        for name, seq in seqs.items():
            fh.write(f">{headers.get(name, name)}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + '\n')

def extract_species(mash_id):
    """Return genus+species from a mob_suite mash_neighbor_identification string."""
    if not mash_id or mash_id.strip().upper() in ('NA', '-', '.', ''):
        return 'NA'
    parts = mash_id.strip().split()
    if len(parts) >= 2:
        return f"{parts[0]} {parts[1]}"
    return parts[0] if parts else 'NA'

def clean_mob_val(v):
    """Return v unchanged if meaningful, else '-'."""
    if not v or v.strip().upper() in ('NA', '-', '.', ''):
        return '-'
    return v.strip()

def load_mobtyper(mob_dir):
    """
    Read mobtyper_results.txt and return a dict keyed by primary_cluster_id:
        { primary_cluster_id: {'secondary': str, 'replicons': str} }

    secondary_cluster_id and rep_type(s) are in this file, not in contig_report.txt.
    """
    path = os.path.join(mob_dir, 'mobtyper_results.txt')
    lookup = {}
    if not os.path.exists(path):
        return lookup
    with open(path) as fh:
        header = fh.readline().strip().split('\t')
        cli = next((header.index(c) for c in ('primary_cluster_id', 'cluster_id')
                    if c in header), None)
        sci = next((header.index(c) for c in ('secondary_cluster_id', 'mob_cluster_id')
                    if c in header), None)
        rpi = next((header.index(c) for c in
                    ('rep_type(s)', 'rep_types', 'replicon_type(s)', 'replicon_types')
                    if c in header), None)
        if cli is None:
            return lookup
        for line in fh:
            cols = line.strip().split('\t')
            if not cols:
                continue
            cluster = cols[cli].strip() if cli < len(cols) else ''
            if not cluster or cluster.upper() in ('NA', '-', '.', ''):
                continue
            sec = clean_mob_val(cols[sci].strip() if sci is not None and sci < len(cols) else '-')
            rep = clean_mob_val(cols[rpi].strip() if rpi is not None and rpi < len(cols) else '-')
            if cluster not in lookup:
                lookup[cluster] = {'secondary': sec, 'replicons': rep}
    return lookup

def load_contig_coverage(info_path):
    """
    Per-contig coverage and length from Flye's assembly_info.txt.
    Columns: #seq_name  length  cov.  circ.  repeat  mult.  alt_group  graph_path
    """
    cov, length = {}, {}
    if not info_path or not os.path.isfile(info_path):
        return cov, length
    try:
        with open(info_path) as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                cols = line.rstrip('\n').split('\t')
                if len(cols) < 3:
                    continue
                name = cols[0].strip()
                try:
                    length[name] = int(cols[1])
                    cov[name] = float(cols[2])
                except ValueError:
                    continue
    except OSError:
        pass
    return cov, length

def load_mobility(report_path):
    """
    contig_id -> {'relaxase': str, 'orit': str} from mob_recon's contig_report.txt.
    MOB-recon calls a contig 'chromosome' when it finds no replicon, even if the
    contig carries a relaxase or oriT — those are the ones worth re-examining.
    """
    out = {}
    if not os.path.isfile(report_path):
        return out
    with open(report_path) as fh:
        header = fh.readline().strip().split('\t')
        try:
            ci = header.index('contig_id')
        except ValueError:
            return out
        rli = next((header.index(c) for c in
                    ('relaxase_type(s)', 'relaxase_types', 'relaxase_type') if c in header), None)
        oti = next((header.index(c) for c in
                    ('orit_type(s)', 'orit_types', 'orit_type') if c in header), None)
        for line in fh:
            cols = line.rstrip('\n').split('\t')
            if ci >= len(cols):
                continue
            out[cols[ci]] = {
                'relaxase': clean_mob_val(cols[rli]) if rli is not None and rli < len(cols) else '-',
                'orit':     clean_mob_val(cols[oti]) if oti is not None and oti < len(cols) else '-',
            }
    return out

def rescue_mge_contigs(chrom_contigs, mobility, cov, clen, max_len):
    """
    Find contigs MOB-recon assigned to the chromosome that are really mobile
    elements: they carry a relaxase or oriT.  Returns a list of
    (contig, length, coverage, relaxase, orit); callers move these out of the
    chromosome so the genus filter cannot discard them.

    Read depth is deliberately NOT a criterion.  Multi-copy elements should sit
    above chromosome depth in principle, but in these ONT libraries plasmid and
    chromosome coverage are not reliably separated, so a depth test discards
    genuine elements.  Coverage is still reported for inspection.

    Contigs longer than max_len are never rescued: beyond the size of a real
    plasmid or ICE, a relaxase-bearing contig is more likely a collapsed repeat
    or a second chromosome than a mobile element.
    """
    rescued = []
    if not chrom_contigs:
        return rescued
    main_ctg = max(chrom_contigs, key=lambda c: clen.get(c, 0), default=None)
    for c in sorted(chrom_contigs, key=lambda x: -clen.get(x, 0)):
        if c == main_ctg:          # never reclassify the chromosome itself
            continue
        if max_len and clen.get(c, 0) > max_len:
            continue
        mob = mobility.get(c, {})
        relax, orit = mob.get('relaxase', '-'), mob.get('orit', '-')
        if relax == '-' and orit == '-':
            continue
        rescued.append((c, clen.get(c, 0), cov.get(c), relax, orit))
    return rescued

def main():
    argv = sys.argv[1:]
    assembly_info = None
    max_len = 1000000
    mge_rescue = True
    positional = []
    i = 0
    while i < len(argv):
        if argv[i] == '--assembly_info' and i + 1 < len(argv):
            assembly_info = argv[i + 1]; i += 2
        elif argv[i] == '--mge_max_len' and i + 1 < len(argv):
            max_len = int(float(argv[i + 1])); i += 2
        elif argv[i] == '--no_mge_rescue':
            mge_rescue = False; i += 1
        else:
            positional.append(argv[i]); i += 1
    mob_dir, asm_fasta, out_dir = positional[0], positional[1], positional[2]
    sample_name = positional[3] if len(positional) > 3 else ''
    os.makedirs(out_dir, exist_ok=True)

    # Load per-plasmid typing info from mobtyper_results.txt (has secondary_cluster_id
    # and rep_type(s); these columns are absent or empty in contig_report.txt).
    mobtyper = load_mobtyper(mob_dir)

    # Parse contig_report.txt
    report = os.path.join(mob_dir, 'contig_report.txt')
    chrom_contigs = set()
    plasmid_map = defaultdict(list)  # cluster_id -> [contig_ids]
    plasmid_names = {}               # cluster_id -> display_name (with sample prefix)
    plasmid_primary_clean = {}       # cluster_id -> clean primary mob cluster (no prefix)
    plasmid_species = {}             # cluster_id -> species
    plasmid_secondary_mob = {}       # cluster_id -> secondary_cluster_id
    plasmid_replicons = {}           # cluster_id -> rep_type(s) string

    novel_counter = 0
    with open(report) as fh:
        header = fh.readline().strip().split('\t')
        try:
            ci  = header.index('contig_id')
            mti = header.index('molecule_type')
            cli = next((header.index(c) for c in ('primary_cluster_id', 'cluster_id') if c in header), None)
        except ValueError as e:
            sys.exit(f"ERROR: unexpected contig_report.txt columns: {e}\nHeaders: {header}")

        # Species column (contig_report.txt)
        spi = next((header.index(c) for c in
                    ('mash_neighbor_identification', 'predicted_host_range_overall_name')
                    if c in header), None)

        for line in fh:
            cols = line.strip().split('\t')
            if not cols:
                continue
            contig  = cols[ci]  if ci  < len(cols) else None
            mol     = cols[mti].lower() if mti < len(cols) else ''
            cluster = cols[cli].strip()  if cli is not None and cli < len(cols) else None

            if not contig:
                continue

            if 'chromosome' in mol:
                chrom_contigs.add(contig)
            elif 'plasmid' in mol and cluster:
                plasmid_map[cluster].append(contig)

                if cluster not in plasmid_names:
                    if cluster and cluster.upper() not in ('NA', '-', '.', ''):
                        display = f"{sample_name}_{cluster}" if sample_name else cluster
                        plasmid_names[cluster] = display
                        plasmid_primary_clean[cluster] = cluster
                    else:
                        novel_counter += 1
                        suffix = f"novel_{novel_counter}"
                        plasmid_names[cluster] = f"{sample_name}_{suffix}" if sample_name else suffix
                        plasmid_primary_clean[cluster] = suffix

                # Species — first non-NA contig in cluster
                if cluster not in plasmid_species and spi is not None and spi < len(cols):
                    sp = extract_species(cols[spi])
                    if sp != 'NA':
                        plasmid_species[cluster] = sp

                # Secondary mob + replicons — prefer mobtyper_results.txt lookup
                if cluster not in plasmid_secondary_mob:
                    plasmid_secondary_mob[cluster] = mobtyper.get(cluster, {}).get('secondary', '-')
                if cluster not in plasmid_replicons:
                    plasmid_replicons[cluster] = mobtyper.get(cluster, {}).get('replicons', '-')

    # ---------------------------------------------------------------------
    # Rescue misassigned mobile elements.
    # MOB-recon calls a contig 'chromosome' whenever it finds no replicon, so
    # plasmids and other MGEs carrying only a relaxase land in the chromosome.
    # Kraken2 then often assigns them to a different genus (they originate
    # elsewhere) and FILTER_GENUS_CONTIGS deletes them, losing any AMR genes
    # they carry.  Moving them to the plasmid set here — upstream of that
    # filter — keeps them, and sends them through the plasmid AMR tools.
    # ---------------------------------------------------------------------
    rescued = []
    if mge_rescue:
        cov, clen = load_contig_coverage(assembly_info)
        mobility = load_mobility(report)
        rescued = rescue_mge_contigs(chrom_contigs, mobility, cov, clen, max_len)
        for n, (contig, ln, ccov, relax, orit) in enumerate(rescued, 1):
            chrom_contigs.discard(contig)
            cluster = f"rescued_MGE_{n}"
            while cluster in plasmid_map:   # never collide with a real cluster
                n += 1
                cluster = f"rescued_MGE_{n}"
            plasmid_map[cluster].append(contig)
            display = f"{sample_name}_{cluster}" if sample_name else cluster
            plasmid_names[cluster] = display
            plasmid_primary_clean[cluster] = cluster
            plasmid_species[cluster] = 'NA'
            plasmid_secondary_mob[cluster] = '-'
            # Record why it was rescued, so AMR hits stay traceable to the evidence
            plasmid_replicons[cluster] = (f"relaxase:{relax}" if relax != '-'
                                          else f"oriT:{orit}")

    # Audit trail — written even when empty so the output always exists
    with open(os.path.join(out_dir, 'rescued_mge_contigs.tsv'), 'w') as fh:
        fh.write("contig\tlength\tcoverage\trelaxase\torit\n")
        for contig, ln, ccov, relax, orit in rescued:
            cov_s = f"{ccov:.1f}" if ccov is not None else 'NA'
            fh.write(f"{contig}\t{ln}\t{cov_s}\t{relax}\t{orit}\n")

    # Read full assembly
    all_seqs, all_headers = read_fasta(asm_fasta)

    # Write chromosome fasta
    chrom_seqs = {c: s for c, s in all_seqs.items() if c in chrom_contigs}
    write_fasta(os.path.join(out_dir, 'chromosome.fasta'), chrom_seqs, all_headers)

    # Write per-plasmid fastas; build summary lists
    plasmid_display_names = []
    plasmid_species_list  = []
    secondary_mob_list    = []
    inc_replicon_list     = []
    for cluster, contigs in plasmid_map.items():
        display = plasmid_names.get(cluster, f"novel_{cluster}")
        plasmid_display_names.append(display)
        plasmid_species_list.append(plasmid_species.get(cluster, 'NA'))
        secondary_mob_list.append(plasmid_secondary_mob.get(cluster, '-'))
        inc_replicon_list.append(plasmid_replicons.get(cluster, '-'))

        pla_seqs = {c: all_seqs[c] for c in contigs if c in all_seqs}
        safe_name = re.sub(r'[^A-Za-z0-9_.-]', '_', display)
        write_fasta(os.path.join(out_dir, f'{safe_name}.fasta'), pla_seqs, all_headers)

    num_plasmids = len(plasmid_map)
    names_str    = ";".join(plasmid_display_names) if plasmid_display_names else "NA"
    species_str  = ";".join(plasmid_species_list)  if plasmid_species_list  else "NA"
    sec_str      = ";".join(secondary_mob_list)    if secondary_mob_list    else "NA"
    inc_str      = ";".join(inc_replicon_list)     if inc_replicon_list     else "NA"
    print(f"{num_plasmids}\t{names_str}\t{species_str}\t{sec_str}\t{inc_str}")

    # Always write a combined plasmid fasta (empty if no plasmids)
    all_plasmid_seqs = {}
    for cluster, contigs in plasmid_map.items():
        for c in contigs:
            if c in all_seqs:
                all_plasmid_seqs[c] = all_seqs[c]
    write_fasta(os.path.join(out_dir, 'plasmids_all.fasta'), all_plasmid_seqs, all_headers)

    # Write contig → plasmid name map (used by parse_depth.py and combine_amr.py).
    # Uses the clean primary cluster ID (no sample prefix) so downstream tools
    # display names like AA002 directly without any prefix-stripping step.
    with open(os.path.join(out_dir, 'contig_to_plasmid.tsv'), 'w') as fh:
        for c in chrom_contigs:
            fh.write(f"{c}\tchromosome\n")
        for cluster, contigs in plasmid_map.items():
            clean = plasmid_primary_clean.get(cluster, cluster)
            for c in contigs:
                fh.write(f"{c}\t{clean}\n")

    # Write contig → mob annotation map (used by combine_amr.py for AMR location labels)
    # Format: chromosome  OR  IncType;PrimaryMobCluster;SecondaryMobCluster
    with open(os.path.join(out_dir, 'contig_to_mob_annotation.tsv'), 'w') as fh:
        for c in chrom_contigs:
            fh.write(f"{c}\tchromosome\n")
        for cluster, contigs in plasmid_map.items():
            inc = plasmid_replicons.get(cluster, '-')
            pri = plasmid_primary_clean.get(cluster, cluster)
            sec = plasmid_secondary_mob.get(cluster, '-')
            ann = f"{inc};{pri};{sec}"
            for c in contigs:
                fh.write(f"{c}\t{ann}\n")

if __name__ == '__main__':
    main()
