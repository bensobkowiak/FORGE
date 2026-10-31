#!/usr/bin/env python3
"""Parse NanoPlot NanoStats.txt and emit key read metrics as TSV."""
import sys
import re

def parse_nanostat(stat_file):
    metrics = {
        'total_bases': None,
        'num_reads': None,
        'median_q': None,
        'reads_q10': None,
        'median_read_length': None,
    }
    with open(stat_file) as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('Total bases'):
                m = re.search(r'([\d,]+)', line)
                if m:
                    metrics['total_bases'] = int(m.group(1).replace(',', ''))
            elif line.startswith('Number of reads'):
                m = re.search(r'([\d,]+)', line)
                if m:
                    metrics['num_reads'] = int(m.group(1).replace(',', ''))
            elif line.startswith('Median read quality'):
                m = re.search(r'([\d.]+)', line)
                if m:
                    metrics['median_q'] = float(m.group(1))
            elif line.startswith('Median read length'):
                m = re.search(r'([\d,.]+)', line)
                if m:
                    metrics['median_read_length'] = int(float(m.group(1).replace(',', '')))
            elif re.match(r'^>Q10', line):
                # ">Q10:  45,876 (91.0%) 224.3Mb"
                m = re.search(r'([\d,]+)\s*\(', line)
                if m:
                    metrics['reads_q10'] = int(m.group(1).replace(',', ''))

    data_gb           = round(metrics['total_bases'] / 1e9, 4) if metrics['total_bases'] else 'NA'
    num_reads         = metrics['num_reads']         if metrics['num_reads']         else 'NA'
    reads_q10         = metrics['reads_q10']         if metrics['reads_q10']         else 'NA'
    median_q          = metrics['median_q']          if metrics['median_q']          else 'NA'
    median_read_length = metrics['median_read_length'] if metrics['median_read_length'] else 'NA'

    print(f"{data_gb}\t{num_reads}\t{reads_q10}\t{median_q}\t{median_read_length}")

if __name__ == '__main__':
    parse_nanostat(sys.argv[1])
