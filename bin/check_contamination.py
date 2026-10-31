#!/usr/bin/env python3
"""
Detect potential contamination using MALDI-TOF species and Kraken2 reads output.

Contamination is flagged when ANY genus at >5% of reads does not match the MALDI genus,
i.e. the sample is not cleanly the expected species. This covers two cases:
  - The dominant genus is not MALDI genus (MALDI organism may be a minor component)
  - The dominant genus IS the MALDI genus but other genera are also at >5%

When flagged, the MALDI genus taxid is looked up so its reads can be extracted
for genus-specific re-assembly of the chromosomal contigs.

Output: one TSV line — is_contaminated (true/false) TAB maldi_genus TAB taxid
taxid is NA if contamination not detected or MALDI genus not found in Kraken2 report.
"""
import sys
import argparse


def find_genus_taxid(report_path, maldi_genus):
    """Return the taxid for the genus-rank entry matching maldi_genus, or None."""
    mg_lower = maldi_genus.lower()
    with open(report_path) as fh:
        for line in fh:
            cols = line.split('\t')
            if len(cols) < 6:
                continue
            rank  = cols[3].strip()
            taxid = cols[4].strip()
            name  = cols[5].strip()
            if rank == 'G' and name.lower().startswith(mg_lower):
                return taxid
    return None


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--parsed', required=True, help='Kraken2 reads .parsed file')
    p.add_argument('--report', required=True, help='Kraken2 reads .report file')
    p.add_argument('--genus',  required=True, help='MALDI genus name (e.g. Klebsiella)')
    args = p.parse_args()

    maldi_genus = args.genus.strip()
    mg_lower    = maldi_genus.lower()

    try:
        with open(args.parsed) as fh:
            line = fh.readline().strip()
        fields = line.split('\t') if line else []

        major_genus   = fields[0].strip() if len(fields) > 0 else 'NA'
        genera_gt5pct = fields[2].strip() if len(fields) > 2 else 'NA'

        maldi_is_major   = major_genus.lower().startswith(mg_lower)
        has_other_genera = genera_gt5pct != 'NA'

        # Trigger if any genus at >5% is NOT the MALDI genus:
        #   - major genus doesn't match MALDI, OR
        #   - there are other genera (non-major) at >5% regardless of what the major is
        is_contaminated = (not maldi_is_major) or has_other_genera

        if is_contaminated:
            taxid = find_genus_taxid(args.report, maldi_genus) or 'NA'
            print("true\t" + maldi_genus + "\t" + taxid)
        else:
            print("false\t" + maldi_genus + "\tNA")

    except Exception as e:
        sys.stderr.write("check_contamination.py ERROR for genus '{}': {}\n".format(maldi_genus, e))
        # Treat as no contamination so the pipeline falls back to the original assembly
        print("false\t" + maldi_genus + "\tNA")


if __name__ == '__main__':
    main()
