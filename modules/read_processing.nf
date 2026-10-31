process CONCAT_READS {
    tag "$sample"
    label 'process_low'

    input:
    // Stage inputs into their own subdirectory so no input filename can ever
    // collide with an output name. Writing to a colliding name would follow
    // the staged symlink and destroy the original sequencing data.
    tuple val(sample), path(fastqs, stageAs: 'input/*')

    output:
    tuple val(sample), path("${sample}.clean.fastq.gz"), emit: reads
    tuple val(sample), path("${sample}.raw.fastq.gz"),   emit: raw_reads

    script:
    def fq_list = fastqs instanceof List ? fastqs.join(' ') : fastqs
    """
    cat ${fq_list} > ${sample}.raw.fastq.gz

    # Sanitise reads: drop records with seq/qual length mismatch or truncated
    # records at EOF. Output to ${sample}.clean.fastq.gz so the name never
    # collides with the staged input symlink (${sample}.fastq.gz), which would
    # cause Python to follow the symlink and overwrite the original input file.
    #
    # Scan first, rewrite only if something is actually wrong: recompressing a
    # multi-GB file costs ~20x more than reading it, and most runs are clean.
    # Binary mode throughout — TextIOWrapper decoding dominates otherwise.
    python3 - <<'PYEOF'
import gzip, os, shutil, sys

RAW   = '${sample}.raw.fastq.gz'
CLEAN = '${sample}.clean.fastq.gz'

def valid(header, seq, plus, qual):
    return (header.startswith(b'@') and plus.startswith(b'+')
            and seq.strip() and qual.strip()
            and len(seq.strip()) == len(qual.strip()))

# --- pass 1: scan only (no compression) -----------------------------------
total = bad = 0
broken_stream = False
try:
    with gzip.open(RAW, 'rb') as fin:
        while True:
            header = fin.readline()
            if not header:
                break
            seq, plus, qual = fin.readline(), fin.readline(), fin.readline()
            total += 1
            if not qual or not valid(header, seq, plus, qual):
                bad += 1
except (EOFError, OSError) as e:
    broken_stream = True
    sys.stderr.write('Truncated gzip stream: ' + str(e) + '\\n')

if bad == 0 and not broken_stream:
    # Nothing to fix — hard-link instead of recompressing.
    try:
        os.link(RAW, CLEAN)
    except OSError:
        shutil.copy(RAW, CLEAN)
    sys.stderr.write(str(total) + ' reads, no corruption; linked without rewrite.\\n')
    raise SystemExit(0)

# --- pass 2: rewrite, dropping the bad records ----------------------------
sys.stderr.write(str(bad) + ' corrupted record(s) found; rewriting.\\n')
kept = skipped = 0
with gzip.open(RAW, 'rb') as fin, \
     gzip.open(CLEAN, 'wb', compresslevel=1) as fout:
    while True:
        try:
            header = fin.readline()
            if not header:
                break
            seq, plus, qual = fin.readline(), fin.readline(), fin.readline()
        except (EOFError, OSError) as e:
            sys.stderr.write('Truncated gzip stream, stopping: ' + str(e) + '\\n')
            break
        if not qual:
            skipped += 1
            continue
        if valid(header, seq, plus, qual):
            fout.write(header + seq + plus + qual)
            kept += 1
        else:
            skipped += 1

sys.stderr.write(str(skipped) + ' corrupted read(s) removed, ' + str(kept) + ' kept.\\n')
PYEOF
    """
}

process NANOPLOT {
    tag "$sample"
    label 'process_medium'
    conda "${projectDir}/envs/qc.yml"
    publishDir { "${params.outdir}/${sample}/reads/nanoplot" }, mode: 'copy'

    input:
    tuple val(sample), path(reads)

    output:
    tuple val(sample), path("NanoStats.txt"), emit: stats
    path "*"

    script:
    """
    NanoPlot --fastq ${reads} \\
             --outdir . \\
             --threads ${task.cpus} \\
             --prefix "" \\
             --N50
    """
}

process PARSE_READ_STATS {
    tag "$sample"
    label 'process_low'
    conda "${projectDir}/envs/qc.yml"
    publishDir { "${params.outdir}/${sample}/reads" }, mode: 'copy'

    input:
    tuple val(sample), path(nanostats)

    output:
    tuple val(sample), path("${sample}.read_stats.tsv")

    script:
    """
    parse_nanoplot.py ${nanostats} > ${sample}.read_stats.tsv
    """
}
