#!/usr/bin/env python3
"""
Collect per-sample fixed metrics and emit a JSON file.
AMR columns are handled separately by write_results.py.

Usage: combine_results.py --sample <name>
         --nanoplot_stats <NanoStats.txt>
         --kraken_reads <kraken_reads.parsed>
         --kraken_assembly <kraken_assembly.parsed>
         --genome_eval <genome_eval.parsed>
         --plasmid_summary <mob_plasmids.txt>
         [--kleborate <kleborate.txt>]
         [--parex <parex.txt>]
         > <sample.metrics.json>
"""
import os
import sys
import json
import argparse

# K. pneumoniae ST → Clonal Group lookup (SLV/DLV-based, following Pasteur MLST scheme).
# Sources: Wyres & Holt (2018) Curr Opin Microbiol; Navon-Venezia et al (2017) FEMS Microbiol Rev;
#          Gorrie et al (2022) Nat Commun; Pitout et al (2019) J Med Microbiol;
#          Lam et al (2021) Nat Commun; Long et al (2017) Lancet Infect Dis (hypervirulent);
#          EARS-Net/EUSCAPE surveillance reports; Rodrigues et al (2019) mBio (CG231).
_KLEBS_CG = {
    # ── Carbapenem-resistant / globally dominant MDR lineages ──────────────────
    # CG258 — dominant KPC/NDM/OXA-48 lineage worldwide
    '11': 'CG258', '258': 'CG258', '340': 'CG258', '437': 'CG258',
    '512': 'CG258', '1519': 'CG258', '2342': 'CG258',

    # CG15 — globally disseminated ESBL/carbapenemase (CTX-M-15, OXA-48, NDM)
    '15': 'CG15', '39': 'CG15', '1764': 'CG15', '2500': 'CG15', '3567': 'CG15',

    # CG147 — NDM/OXA-48; rapidly expanding globally
    '147': 'CG147', '1760': 'CG147', '3569': 'CG147',

    # CG307 — OXA-48/NDM; dominant in Europe and Middle East
    '307': 'CG307', '364': 'CG307', '2502': 'CG307',

    # CG101 — CTX-M-15/NDM; Europe and South Asia
    '101': 'CG101', '4': 'CG101',

    # CG37 — NDM/OXA-48; major lineage in Asia, Europe, Middle East
    '37': 'CG37', '2587': 'CG37',

    # CG231 — NDM/OXA-232/OXA-48; dominant in South Asia (India, Bangladesh), Turkey, spreading globally
    '231': 'CG231',

    # CG395 — OXA-48; documented across Europe
    '395': 'CG395',

    # CG323 — NDM-producing; documented in global surveillance
    '323': 'CG323', '491': 'CG323',

    # ── Prevalent MDR lineages ─────────────────────────────────────────────────
    # CG14 — MDR; ESBL/carbapenemase, documented globally
    '14': 'CG14',

    # CG16 — MDR lineage
    '16': 'CG16',

    # CG17 — MDR lineage; documented in European surveillance
    '17': 'CG17',

    # CG20 — MDR lineage
    '20': 'CG20',

    # CG29 — MDR lineage
    '29': 'CG29',

    # CG35 — MDR lineage
    '35': 'CG35',

    # CG43 — MDR; documented in European and Asian surveillance
    '43': 'CG43',

    # CG45 — MDR lineage
    '45': 'CG45',

    # CG48 — MDR lineage
    '48': 'CG48',

    # CG116 — MDR; documented in European and Asian surveillance
    '116': 'CG116',

    # CG133 — MDR lineage; documented in multiple regions
    '133': 'CG133',

    # CG152 — MDR lineage; documented in surveillance studies
    '152': 'CG152',

    # CG218 — carbapenem-resistant; documented in European surveillance
    '218': 'CG218',

    # CG380 — convergent AMR + virulence
    '380': 'CG380',

    # ── Hypervirulent lineages ─────────────────────────────────────────────────
    # CG23 — hypervirulent (KL1); community-acquired liver abscess; Asia and globally
    '23': 'CG23', '57': 'CG23', '25': 'CG23',

    # CG82 — hypervirulent (KL2); SLV of ST23; distinct clinical lineage
    '82': 'CG82',

    # CG86 — hypervirulent (KL2); liver abscess and metastatic infection
    '86': 'CG86', '693': 'CG86',

    # CG65 — hypervirulent (KL3)
    '65': 'CG65',

    # CG66 — hypervirulent
    '66': 'CG66',
}

def assign_clonal_group(st_str):
    """Map a Kleborate MLST ST string (e.g. 'ST39', '39') to its clonal group, or 'NA'."""
    if not st_str or st_str == 'NA':
        return 'NA'
    # strip leading 'ST' and any '-NLV' novel-variant suffix
    st = st_str.lstrip('S').lstrip('T').split('-')[0].strip()
    return _KLEBS_CG.get(st, 'NA')

FIXED_COLS = [
    'Name', 'Data_GB', 'Num_Reads', 'Num_Reads_Q10', 'Median_Q', 'Median_Read_Length',
    'Major_Genus_Reads', 'Pct_Major_Genus_Reads', 'Genera_gt5pct_Reads', 'Major_Species_Reads', 'Pct_Major_Species_Reads', 'Species_gt5pct_Reads', 'Pct_Unclassified_Reads',
    'Assembly_Major_Genus', 'Assembly_Pct_Major_Genus', 'Assembly_Genera_gt5pct', 'Assembly_Major_Species', 'Assembly_Pct_Major_Species', 'Assembly_Species_gt5pct', 'Assembly_Pct_Unclassified',
    'Circularised', 'Chr_Contigs', 'Completeness', 'Contamination',
    'Genome_Size', 'Contig_N50', 'Max_Contig_Length', 'Avg_Gene_Length', 'Total_CDS', 'GC_Content',
    'Chr_Mean_Depth',
    'Num_Plasmids', 'Num_Rescued_MGE', 'Num_Query_Contigs', 'Plasmid_Names', 'Plasmid_Species',
    'Primary_Mob_Clusters', 'Secondary_Mob_Clusters', 'Inc_Replicons',
    'Plasmid_Mean_Depths',
    'BUSCO_Complete', 'BUSCO_Single_Copy', 'BUSCO_Duplicated', 'BUSCO_Fragmented', 'BUSCO_Missing',
    'MLST_ST', 'K_type', 'O_type', 'Kleborate_ST', 'Clonal_Group',
    'Virulence_Score', 'Yersiniabactin', 'Colibactin', 'Aerobactin', 'Salmochelin', 'RmpADC', 'RmpA2',
    'PDC_type',
]

def read_first_line(path):
    if not path or not os.path.exists(path):
        return ''
    with open(path) as fh:
        return fh.readline().strip()

def pad(lst, n):
    return (lst + ['NA'] * n)[:n]

def read_kraken_line(path):
    """Read a kraken parsed file, handling both old (4-field) and new (7-field) formats."""
    line = read_first_line(path)
    if not line:
        return ['NA'] * 7
    fields = line.split('\t')
    if len(fields) <= 4:
        # Old format: major_genus, pct, genera_gt5, pct_unclassified
        # Insert NA placeholders for the three species columns before unclassified
        base = pad(fields[:3], 3)
        unclassified = fields[3].strip() if len(fields) >= 4 else 'NA'
        return base + ['NA', 'NA', 'NA', unclassified]
    return pad(fields, 7)

def parse_kleborate(path):
    if not path or not os.path.exists(path):
        return ['NA'] * 11
    with open(path, encoding='utf-8') as fh:
        hdr = fh.readline().strip().split('\t')
        for line in fh:
            row = dict(zip(hdr, line.strip().split('\t')))
            def get(col, *alts):
                v = row.get(col) or ''
                for a in alts:
                    if not v:
                        v = row.get(a) or ''
                if not v and '__' in col:
                    # Column names are prefixed with the preset's module, so a
                    # kosc or escherichia run writes klebsiella_oxytoca_complex__
                    # mlst__ST where kpsc writes klebsiella_pneumo_complex__
                    # mlst__ST. Match on the module/field suffix so every preset
                    # resolves. Preset-specific modules (kaptive) simply stay
                    # absent and fall through to NA.
                    tail = col.split('__', 1)[1]
                    for c2, v2 in row.items():
                        if c2.endswith(tail) and (v2 or '') not in ('', '-'):
                            v = v2
                            break
                return v if v not in ('', '-') else 'NA'
            k   = get('klebsiella_pneumo_complex__kaptive__K_type',   'K_type')
            o   = get('klebsiella_pneumo_complex__kaptive__O_type',   'O_type')
            st  = get('klebsiella_pneumo_complex__mlst__ST',          'ST')
            # Clonal groups are defined on the KpSC MLST scheme. A KoSC ST is a
            # different scheme whose numbers would map onto the wrong CG, so
            # only assign one when this really is a KpSC result.
            is_kpsc = any(c.startswith('klebsiella_pneumo_complex__') for c in hdr)
            cg  = assign_clonal_group(st) if is_kpsc else 'NA'
            vs  = get('klebsiella_pneumo_complex__virulence_score__virulence_score', 'Virulence score')
            ybt = get('klebsiella__ybst__Yersiniabactin',  'Yersiniabactin')
            clb = get('klebsiella__cbst__Colibactin',      'Colibactin')
            aer = get('klebsiella__abst__Aerobactin',      'Aerobactin')
            sal = get('klebsiella__smst__Salmochelin',     'Salmochelin')
            rmp = get('klebsiella__rmst__RmpADC',          'RmpADC')
            r2  = get('klebsiella__rmpa2__rmpA2',          'RmpA2')
            return [k, o, st, cg, vs, ybt, clb, aer, sal, rmp, r2]
    return ['NA'] * 11

def parse_parex(path):
    if not path or not os.path.exists(path):
        return 'NA'
    with open(path) as fh:
        for line in fh:
            if 'PDC' in line or 'pdc' in line.lower():
                cols = line.strip().split('\t')
                if len(cols) >= 2:
                    return cols[1]
    return 'NA'

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--sample',          required=True)
    p.add_argument('--nanoplot_stats',  required=True)
    p.add_argument('--kraken_reads',    required=True)
    p.add_argument('--kraken_assembly', default=None)
    p.add_argument('--genome_eval',     default=None)
    p.add_argument('--plasmid_summary', required=True)
    p.add_argument('--chr_contig_species', default=None,
                   help='Species of the chromosome contigs carried forward')
    p.add_argument('--num_query_contigs', default=None,
                   help='Count of off-genus query contigs (from FILTER_GENUS_CONTIGS)')
    p.add_argument('--kleborate',       default=None)
    p.add_argument('--parex',           default=None)
    p.add_argument('--depth',           default=None)
    p.add_argument('--maldi_species',       default=None,
                   help='Expected species from MALDI-TOF (e.g. "Klebsiella pneumoniae")')
    p.add_argument('--maldi_contigs_used',  default=None,
                   help='Yes/No/NA string indicating whether MALDI genus contigs were used')
    p.add_argument('--removed_contig_lengths', default=None,
                   help='";"-separated lengths of contigs removed by MALDI genus filtering')
    p.add_argument('--coverage_cap_bases', default=None,
                   help='filtlong --target_bases value; used to flag whether this '
                        'sample was actually subsampled before assembly')
    args = p.parse_args()

    # genome_eval has 16 fields: [circular, 9 checkm2, 5 busco, mlst]
    # Depth must go between GC_Content (index 9) and BUSCO_Complete (index 10)
    if args.genome_eval:
        geval_all  = pad(read_first_line(args.genome_eval).split('\t'), 16)
        geval_pre  = geval_all[:10]   # Unicycler_Circular through GC_Content
        geval_post = geval_all[10:]   # BUSCO_Complete through MLST_ST
    else:
        geval_pre  = ['NA'] * 10
        geval_post = ['NA'] * 6

    # depth file now has 3 fields: chr_depth, plasmid_depths, circular
    # Flye's assembly-graph call (geval_pre[0]) is authoritative: it is available
    # for every assembly, whereas the alignment-based call needs CALC_DEPTH to have
    # produced a BAM with enough junction-spanning reads.  Alignment is the fallback.
    if args.depth:
        depth_raw = pad(read_first_line(args.depth).split('\t'), 3)
    else:
        depth_raw = ['NA', 'NA', 'NA']
    if geval_pre[0] in ('NA', '') and depth_raw[2] not in ('NA', ''):
        geval_pre[0] = depth_raw[2]   # Flye said nothing; fall back to alignment
    depth_vals = depth_raw[:2]        # only chr_depth + plasmid_depths go into vals

    # Strip sample-name prefix from plasmid names (e.g. GRC852_AA338 -> AA338)
    prefix = f"{args.sample}_"
    def strip_prefix(val):
        if not val or val == 'NA':
            return val
        parts = []
        for part in val.split(';'):
            if ':' in part:
                name, rest = part.split(':', 1)
                name = name[len(prefix):] if name.startswith(prefix) else name
                parts.append(f"{name}:{rest}")
            else:
                parts.append(part[len(prefix):] if part.startswith(prefix) else part)
        return ';'.join(parts)

    depth_vals[1] = strip_prefix(depth_vals[1])  # Plasmid_Mean_Depths

    # plasmid_summary.txt: num  names  species  secondary_mob_clusters  inc_replicons
    # Older files may have only 3 columns; pad to 5 with 'NA' for backward compat.
    pla_summary = pad(read_first_line(args.plasmid_summary).split('\t'), 5)
    pla_summary[1] = strip_prefix(pla_summary[1])  # Plasmid_Names (= Primary_Mob_Clusters)

    # Contigs rescued from the chromosome by process_mob_recon.py are named
    # rescued_MGE_<n>, so they can be counted straight off the plasmid name
    # list — no extra input file needs threading through COMPILE_SAMPLE.
    _names = pla_summary[1]
    rescued_n = str(len([x for x in _names.split(';')
                         if x.strip().startswith('rescued_MGE_')])) \
        if _names and _names != 'NA' else '0'

    # Query contigs: off-genus, below the element-size limit, screened for AMR
    # but excluded from chromosome and plasmid metrics.
    query_n = (args.num_query_contigs or '0').strip()
    if not query_n or query_n == 'NA':
        query_n = '0'

    vals = (
        [args.sample]
        + pad(read_first_line(args.nanoplot_stats).split('\t'), 5)
        + read_kraken_line(args.kraken_reads)
        + read_kraken_line(args.kraken_assembly)
        + geval_pre
        + [depth_vals[0]]               # Chr_Mean_Depth
        + [pla_summary[0], rescued_n, query_n]  # Num_Plasmids, Num_Rescued_MGE, Num_Query_Contigs
        + pla_summary[1:3]              # Plasmid_Names, Plasmid_Species
        + [pla_summary[1]]              # Primary_Mob_Clusters  (same as Plasmid_Names)
        + pla_summary[3:]               # Secondary_Mob_Clusters, Inc_Replicons
        + [depth_vals[1]]               # Plasmid_Mean_Depths
        + geval_post
        + parse_kleborate(args.kleborate)
        + [parse_parex(args.parex)]
    )
    vals = pad(vals, len(FIXED_COLS))

    row = dict(zip(FIXED_COLS, vals))

    if args.maldi_species:
        maldi_sp    = args.maldi_species.strip()
        maldi_genus = maldi_sp.split()[0].lower()

        # MALDI_Species_Match: was the MALDI species the dominant species in reads?
        kraken_sp = row.get('Major_Species_Reads', '')
        if kraken_sp and kraken_sp != 'NA':
            row['MALDI_Species_Match'] = (
                'True' if kraken_sp.lower().strip() == maldi_sp.lower() else 'False'
            )
        else:
            row['MALDI_Species_Match'] = 'NA'

        # Maldi_genus_assembled: was the chromosome assembled from the MALDI genus?
        # Derived from Kraken2 on whichever chromosome was used (original or genus re-assembly).
        asm_genus = row.get('Assembly_Major_Genus', '')
        if asm_genus and asm_genus not in ('NA', ''):
            row['Maldi_genus_assembled'] = (
                'Yes' if asm_genus.lower().startswith(maldi_genus) else 'No'
            )
        else:
            row['Maldi_genus_assembled'] = 'NA'
    else:
        row['MALDI_Species_Match']    = 'NA'
        row['Maldi_genus_assembled'] = 'NA'

    row['Maldi_contigs_used'] = args.maldi_contigs_used or 'NA'
    row['Chr_Contig_Species'] = args.chr_contig_species or 'NA'
    row['Removed_Contig_Lengths'] = args.removed_contig_lengths or 'NA'

    # True only when this sample had more data than the cap, i.e. Filtlong
    # actually subsampled it. Assembly used a subset; Chr_Mean_Depth did not.
    if args.coverage_cap_bases:
        try:
            row['Coverage_Cap_Applied'] = (
                'True' if float(row.get('Data_GB')) * 1e9 > float(args.coverage_cap_bases)
                else 'False')
        except (TypeError, ValueError):
            row['Coverage_Cap_Applied'] = 'NA'
    else:
        row['Coverage_Cap_Applied'] = 'False'

    # Format GC_Content as a percentage to 2 decimal places
    gc = row.get('GC_Content', 'NA')
    if gc not in ('NA', ''):
        try:
            gc_f = float(gc)
            if gc_f <= 1.0:
                gc_f *= 100
            row['GC_Content'] = f"{gc_f:.2f}"
        except ValueError:
            pass

    print(json.dumps(row))

if __name__ == '__main__':
    main()
