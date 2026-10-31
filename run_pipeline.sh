#!/usr/bin/env bash


set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ── Defaults — override on the command line ───────────────────────────────────
INPUT_DIR=""
MALDI_CSV="${HERE}/Maldi_species.csv"
OUTDIR="${HERE}/results"
RESULTS_FILE="${HERE}/results.tsv"
WORK_DIR="${HERE}/work"
PROGRAMS="${PROGRAMS:-${HOME}/Documents/Programs}"
KRAKEN_DB="${PROGRAMS}/KrakenDB"
CHECKM2_DB=""          
RESFINDER_DB="${PROGRAMS}/db_resfinder"
AMRFINDER_DB="${PROGRAMS}/amrfinderplus_db/latest"
SAMPLE_LIST=""
FLYE_BIN=""
RESUME="-resume"
EXTRA=()

usage() {
    cat <<EOF
Usage: ./run_pipeline.sh --fastq <dir> [options]

  --fastq        <dir>   Input directory (barcode subdirs or flat fastq.gz)   [required]
  --maldi_csv    <file>  MALDI species CSV            [${MALDI_CSV}]
  --outdir       <dir>   Per-sample results           [${OUTDIR}]
  --results_file <file>  Wide results TSV             [${RESULTS_FILE}]
  --work         <dir>   Nextflow work directory      [${WORK_DIR}]
  --kraken_db    <dir>   Kraken2 database             [${KRAKEN_DB}]
  --checkm2_db   <file>  CheckM2 .dmnd file           [auto-found in ${PROGRAMS}/checkm2_db]
  --resfinder_db <dir>   ResFinder database           [${RESFINDER_DB}]
  --amrfinder_db <dir>   AMRFinder+ database          [${AMRFINDER_DB}]
  --samples      <file>  Only process samples named in this file (one per
                         line). Listed samples are reprocessed even if they
                         are already in the results file.
  --flye_bin     <path>  Specific flye executable     [whatever conda provides]
  --fresh                Ignore the Nextflow cache (omit -resume)
  --                     Everything after this is passed through to nextflow
EOF
    exit 1
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq)        INPUT_DIR="$2";    shift 2 ;;
        --maldi_csv)    MALDI_CSV="$2";    shift 2 ;;
        --outdir)       OUTDIR="$2";       shift 2 ;;
        --results_file) RESULTS_FILE="$2"; shift 2 ;;
        --work)         WORK_DIR="$2";     shift 2 ;;
        --kraken_db)    KRAKEN_DB="$2";    shift 2 ;;
        --checkm2_db)   CHECKM2_DB="$2";   shift 2 ;;
        --resfinder_db) RESFINDER_DB="$2"; shift 2 ;;
        --amrfinder_db) AMRFINDER_DB="$2"; shift 2 ;;
        --samples)      SAMPLE_LIST="$2";  shift 2 ;;
        --flye_bin)     FLYE_BIN="$2";     shift 2 ;;
        --fresh)        RESUME="";         shift   ;;
        --)             shift; EXTRA=("$@"); break ;;
        -h|--help)      usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done

[[ -n "${INPUT_DIR}" ]] || { echo "Error: --fastq is required" >&2; usage; }
[[ -d "${INPUT_DIR}" ]] || { echo "Error: input dir not found: ${INPUT_DIR}" >&2; exit 1; }
[[ -d "${KRAKEN_DB}" ]] || { echo "Error: kraken_db not found: ${KRAKEN_DB}" >&2; exit 1; }
[[ -z "${SAMPLE_LIST}" || -f "${SAMPLE_LIST}" ]] || \
    { echo "Error: --samples file not found: ${SAMPLE_LIST}" >&2; exit 1; }

command -v nextflow >/dev/null || { echo "Error: nextflow not on PATH" >&2; exit 1; }

# CheckM2 
if [[ -z "${CHECKM2_DB}" ]]; then
    CHECKM2_DB=$(find "${PROGRAMS}/checkm2_db" -name "*.dmnd" -type f 2>/dev/null | head -1)
fi
[[ -n "${CHECKM2_DB}"   && -f "${CHECKM2_DB}"   ]] || echo "Warning: no CheckM2 .dmnd found under ${PROGRAMS}/checkm2_db — CheckM2 will fail" >&2
[[ -d "${RESFINDER_DB}" ]] || echo "Warning: ResFinder db not found at ${RESFINDER_DB} — ResFinder will fail" >&2
if [[ ! -d "${AMRFINDER_DB}" ]]; then
    echo "Warning: AMRFinder db not found at ${AMRFINDER_DB} — AMRFinder will fail." >&2
    echo "         The conda package ships no database. Download it once with:" >&2
    echo "           amrfinder_update -d ${AMRFINDER_DB%/latest}" >&2
fi

# ── Size the run to this machine ──────────────────────────────────────────────

PERF_CORES=$(sysctl -n hw.perflevel0.physicalcpu 2>/dev/null || sysctl -n hw.ncpu)
TOTAL_GB=$(( $(sysctl -n hw.memsize) / 1073741824 ))
USABLE_GB=$(( TOTAL_GB * 90 / 100 ))

echo "Machine : $(sysctl -n machdep.cpu.brand_string 2>/dev/null || uname -m)"
echo "Cores   : ${PERF_CORES} performance"
echo "Memory  : ${TOTAL_GB} GB total, allowing ${USABLE_GB} GB"
echo "Input   : ${INPUT_DIR}"
echo "Kraken2 : ${KRAKEN_DB}"
echo "CheckM2 : ${CHECKM2_DB:-<missing>}"
echo "ResFinder: ${RESFINDER_DB}"
echo "AMRFinder: ${AMRFINDER_DB}"
[[ -n "${SAMPLE_LIST}" ]] && echo "Samples  : ${SAMPLE_LIST} ($(grep -cve '^[[:space:]]*$' -e '^#' "${SAMPLE_LIST}") listed)"
echo

mkdir -p "${OUTDIR}" "${WORK_DIR}"

OPT=()
[[ -n "${MALDI_CSV}"    && -f "${MALDI_CSV}"    ]] && OPT+=(--maldi_csv    "${MALDI_CSV}")
[[ -n "${CHECKM2_DB}"   ]] && OPT+=(--checkm2_db   "${CHECKM2_DB}")
[[ -n "${RESFINDER_DB}" ]] && OPT+=(--resfinder_db "${RESFINDER_DB}")
[[ -d "${AMRFINDER_DB}" ]] && OPT+=(--amrfinder_db "${AMRFINDER_DB}")
[[ -n "${SAMPLE_LIST}" ]] && OPT+=(--sample_list  "${SAMPLE_LIST}")
[[ -n "${FLYE_BIN}"     ]] && OPT+=(--flye_bin     "${FLYE_BIN}")

nextflow run "${HERE}/main.nf" \
    -profile macstudio \
    -w "${WORK_DIR}" \
    --input_dir    "${INPUT_DIR}" \
    --outdir       "${OUTDIR}" \
    --results_file "${RESULTS_FILE}" \
    --kraken_db    "${KRAKEN_DB}" \
    --max_cpus     4 \
    --max_memory   40.GB \
    -with-report   "${HERE}/pipeline_report.html" \
    -with-timeline "${HERE}/pipeline_timeline.html" \
    -with-trace    "${HERE}/pipeline_trace.txt" \
    "${OPT[@]}" \
    ${RESUME} \
    "${EXTRA[@]+"${EXTRA[@]}"}"

echo
echo "Done. Results: ${RESULTS_FILE}"
echo "Timing report: ${HERE}/pipeline_report.html"
