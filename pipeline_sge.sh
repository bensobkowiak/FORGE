#!/bin/bash
#$ -S /bin/bash
#$ -N forge
#$ -cwd
#$ -pe smp 1
#$ -R y
#$ -l h_rt=24:00:00
#
# FORGE — SGE array launcher. One array task per sample.
#
#   qsub -t 1-<N> [-tc 200] pipeline_sge.sh --fastq <dir> [options]
#

set -euo pipefail

usage() {
    cat <<EOF
Usage: qsub -t 1-<N> pipeline_sge.sh --fastq <dir> [options]

  --fastq        <dir>   input: barcode subdirs or flat fastq(.gz)   [required]
  --data_dir     <dir>   work/, conda_envs/, results/, results.tsv   [\$PWD]
  --pipeline_dir <dir>   location of main.nf and bin/                [submit dir]
  --maldi_csv    <file>  MALDI species CSV
  --outdir       <dir>   per-sample results        [<data_dir>/results]
  --results_file <file>  wide results table        [<data_dir>/results.tsv]
  --kraken_db    <dir>   Kraken2 database                            [required]
  --resfinder_db <dir>   ResFinder database
  --amrfinder_db <dir>   AMRFinderPlus database
  --checkm2_db   <file>  CheckM2 .dmnd
  --checkm2_env  <dir>   pre-built CheckM2 conda env
  --flye_bin     <path>  system flye (avoids SIGILL on older nodes)
  --nf_env       <name>  conda env providing nextflow + java   [none]
  --conda_sh     <path>  conda profile.d/conda.sh                [auto-detect]
  --nextflow     <path>  nextflow executable                     [auto-detect]
  --cpus         <n>     cpus per sample                              [NSLOTS]
  --memory       <n>GB   memory per sample                              [16.GB]
  --config       <file>  shell file of VAR=value defaults

Each option also reads an environment variable of the same name in upper case
(FASTQ, DATA_DIR, KRAKEN_DB, ...).
EOF
    exit 1
}

# ── Site config ───────────────────────────────────────────────────────────────
CONFIG="${CONFIG:-}"
args=("$@")
for ((i=0; i<${#args[@]}; i++)); do
    [[ "${args[$i]}" == "--config" ]] && CONFIG="${args[$((i+1))]}"
done
SCRIPT_DIR="${SGE_O_WORKDIR:-$PWD}"
[[ -z "${CONFIG}" && -f "${SCRIPT_DIR}/pipeline_sge.conf" ]] && CONFIG="${SCRIPT_DIR}/pipeline_sge.conf"
# shellcheck disable=SC1090
[[ -n "${CONFIG}" && -f "${CONFIG}" ]] && source "${CONFIG}"

# ── Defaults (environment wins over these, flags win over everything) ─────────
INPUT_DIR="${FASTQ:-}"
DATA_DIR="${DATA_DIR:-$PWD}"
PIPELINE_DIR="${PIPELINE_DIR:-$SCRIPT_DIR}"
MALDI_CSV="${MALDI_CSV:-}"
OUTDIR="${OUTDIR:-}"
RESULTS_FILE="${RESULTS_FILE:-}"
KRAKEN_DB="${KRAKEN_DB:-}"
RESFINDER_DB="${RESFINDER_DB:-}"
AMRFINDER_DB="${AMRFINDER_DB:-}"
CHECKM2_DB="${CHECKM2_DB:-}"
CHECKM2_ENV="${CHECKM2_ENV:-}"
FLYE_BIN="${FLYE_BIN:-}"
NF_ENV="${NF_ENV:-}"
CONDA_SH="${CONDA_SH:-}"
NEXTFLOW="${NEXTFLOW:-}"
CPUS="${CPUS:-}"
MEMORY="${MEMORY:-16.GB}"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --fastq)        INPUT_DIR="$2";    shift 2 ;;
        --data_dir)     DATA_DIR="$2";     shift 2 ;;
        --pipeline_dir) PIPELINE_DIR="$2"; shift 2 ;;
        --maldi_csv)    MALDI_CSV="$2";    shift 2 ;;
        --outdir)       OUTDIR="$2";       shift 2 ;;
        --results_file) RESULTS_FILE="$2"; shift 2 ;;
        --kraken_db)    KRAKEN_DB="$2";    shift 2 ;;
        --resfinder_db) RESFINDER_DB="$2"; shift 2 ;;
        --amrfinder_db) AMRFINDER_DB="$2"; shift 2 ;;
        --checkm2_db)   CHECKM2_DB="$2";   shift 2 ;;
        --checkm2_env)  CHECKM2_ENV="$2";  shift 2 ;;
        --flye_bin)     FLYE_BIN="$2";     shift 2 ;;
        --nf_env)       NF_ENV="$2";       shift 2 ;;
        --conda_sh)     CONDA_SH="$2";     shift 2 ;;
        --nextflow)     NEXTFLOW="$2";     shift 2 ;;
        --cpus)         CPUS="$2";         shift 2 ;;
        --memory)       MEMORY="$2";       shift 2 ;;
        --config)       shift 2 ;;
        -h|--help)      usage ;;
        *) echo "Unknown option: $1" >&2; usage ;;
    esac
done

[[ -n "${INPUT_DIR}" ]] || { echo "Error: --fastq is required" >&2; usage; }
[[ -d "${INPUT_DIR}" ]] || { echo "Error: input dir not found: ${INPUT_DIR}" >&2; exit 1; }

OUTDIR="${OUTDIR:-${DATA_DIR}/results}"
RESULTS_FILE="${RESULTS_FILE:-${DATA_DIR}/results.tsv}"
PIPELINE="${PIPELINE_DIR}/main.nf"
PIPELINE_BIN="${PIPELINE_DIR}/bin"
[[ -f "${PIPELINE}" ]] || { echo "Error: main.nf not found at ${PIPELINE} (use --pipeline_dir)" >&2; exit 1; }

# ── Conda / Nextflow ──────────────────────────────────────────────────────────
if [[ -z "${CONDA_SH}" ]]; then
    for c in "${CONDA_EXE:+$(dirname "$(dirname "${CONDA_EXE}")")/etc/profile.d/conda.sh}" \
             "${HOME}/miniforge3/etc/profile.d/conda.sh" \
             "${HOME}/miniconda3/etc/profile.d/conda.sh" \
             /opt/conda/etc/profile.d/conda.sh; do
        [[ -n "${c}" && -f "${c}" ]] && { CONDA_SH="${c}"; break; }
    done
fi
if [[ -n "${CONDA_SH}" ]]; then
    # shellcheck disable=SC1090
    source "${CONDA_SH}"
    [[ -f "${CONDA_SH%conda.sh}mamba.sh" ]] && source "${CONDA_SH%conda.sh}mamba.sh"
    [[ -n "${NF_ENV}" ]] && { conda activate "${NF_ENV}" 2>/dev/null \
        || echo "Note: conda env '${NF_ENV}' not activated" >&2; }
    export JAVA_HOME="${CONDA_PREFIX:-}"
    unset JAVA_CMD
fi

NEXTFLOW="${NEXTFLOW:-$(command -v nextflow || true)}"
[[ -n "${NEXTFLOW}" ]] || { echo "Error: nextflow not found (use --nextflow)" >&2; exit 1; }

export NXF_WORK="${NXF_WORK:-${DATA_DIR}/work}"
export NXF_CONDA_CACHEDIR="${NXF_CONDA_CACHEDIR:-${DATA_DIR}/conda_envs}"
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-${DATA_DIR}/conda_pkgs}"
export NXF_JVM_ARGS="${NXF_JVM_ARGS:--Xmx2g -Xss512k}"
ulimit -u unlimited 2>/dev/null || true

CPUS="${CPUS:-${NSLOTS:-1}}"
SAMPLE_LIST="${OUTDIR}/.samples.txt"
DONE_DIR="${OUTDIR}/.sge_done"
FAILED_TSV="${DATA_DIR}/failed.tsv"
mkdir -p "${OUTDIR}" "${DONE_DIR}"

# ── Sample list: first task builds it, the rest wait on the lock ─────────────
if [[ ! -s "${SAMPLE_LIST}" ]]; then
    (
        flock -w 60 9
        if [[ ! -s "${SAMPLE_LIST}" ]]; then
            SAMPLES_A=$(find "${INPUT_DIR}" -mindepth 1 -maxdepth 1 -type d -print0 \
                | xargs -0 -I{} bash -c \
                    'ls "{}/"*.fastq.gz "{}/"*.fastq 2>/dev/null | head -1 | grep -q . && basename "{}"' \
                2>/dev/null | sort || true)
            SAMPLES_B=$(find "${INPUT_DIR}" -maxdepth 1 \
                \( -name "*.fastq.gz" -o -name "*.fastq" \) \
                | xargs -I{} basename {} \
                | sed -E 's/\.(fastq\.gz|fastq)$//' | sed -E 's/_[0-9]+$//' | sort -u || true)
            ALL_SAMPLES=$(printf '%s\n%s\n' "${SAMPLES_A}" "${SAMPLES_B}" | sort -u | grep -v '^$' || true)
            [[ -n "${ALL_SAMPLES}" ]] || { echo "ERROR: no samples found in ${INPUT_DIR}" >&2; exit 1; }

            DONE_SAMPLES=""
            [[ -f "${RESULTS_FILE}" ]] && \
                DONE_SAMPLES+=$(awk -F'\t' 'NR>1 && $1!="" {gsub(/\r/,"",$1); print $1}' "${RESULTS_FILE}")$'\n'
            [[ -f "${FAILED_TSV}" ]] && \
                DONE_SAMPLES+=$(awk -F'\t' 'NR>1 && $1!="" {gsub(/\r/,"",$1); print $1}' "${FAILED_TSV}")$'\n'
            DONE_SET=$(printf '%s' "${DONE_SAMPLES}" | sort -u)

            MISSING=""
            while IFS= read -r sample; do
                [[ -n "${sample}" ]] || continue
                printf '%s\n' "${DONE_SET}" | grep -qxF "${sample}" || MISSING+="${sample}"$'\n'
            done <<< "${ALL_SAMPLES}"

            printf '%s' "${MISSING}" > "${SAMPLE_LIST}.tmp"
            mv "${SAMPLE_LIST}.tmp" "${SAMPLE_LIST}"
            echo "Sample list created: $(wc -l < "${SAMPLE_LIST}") outstanding samples"
        fi
    ) 9>"${SAMPLE_LIST}.lock"
fi

N_TOTAL=$(wc -l < "${SAMPLE_LIST}" | tr -d ' ')
TASK_ID="${SGE_TASK_ID:-1}"
if [[ "${TASK_ID}" -gt "${N_TOTAL}" ]]; then
    echo "Task ${TASK_ID}: no sample (total ${N_TOTAL}) — exiting."
    exit 0
fi
SAMPLE=$(sed -n "${TASK_ID}p" "${SAMPLE_LIST}")

# Nodes differ in instruction-set support, and a binary built for a newer CPU
# dies with SIGILL on an older one — an opaque tool failure without this.
echo "Host: $(hostname)  task=${TASK_ID}  sample=${SAMPLE}"
[[ -r /proc/cpuinfo ]] && echo "CPU flags: $(grep -m1 ^flags /proc/cpuinfo | tr ' ' '\n' \
    | grep -xE 'sse4_2|popcnt|avx|avx2|avx512f|bmi2' | paste -sd, -)"

# ── Which step failed, inferred from published outputs ───────────────────────
detect_failure() {
    if   [[ -s "${OUTDIR}/${SAMPLE}/${SAMPLE}.amr.tsv" ]] && \
         [[ -s "${OUTDIR}/${SAMPLE}/${SAMPLE}.metrics.json" ]];                                  then echo "Complete"
    elif [[ -s "${OUTDIR}/${SAMPLE}/${SAMPLE}.amr.tsv" ]];                                       then echo "Metrics compile failed"
    elif [[ -s "${OUTDIR}/${SAMPLE}/${SAMPLE}.metrics.json" ]];                                  then echo "AMR combine failed"
    elif [[ -f "${OUTDIR}/${SAMPLE}/chromosome/amrfinder/${SAMPLE}_chromosome.amrfinder.tsv" ]]; then echo "AMR compile failed"
    elif [[ -f "${OUTDIR}/${SAMPLE}/chromosome/mlst/${SAMPLE}.mlst.tsv" ]];                      then echo "AMR analysis failed"
    elif [[ -d "${OUTDIR}/${SAMPLE}/chromosome/busco/busco_out" ]];                              then echo "MLST failed"
    elif [[ -f "${OUTDIR}/${SAMPLE}/chromosome/checkm2/checkm2_output/quality_report.tsv" ]];    then echo "BUSCO failed"
    elif [[ -f "${OUTDIR}/${SAMPLE}/assembly/contigs/chromosome.fasta" ]];                       then echo "CheckM2 failed"
    elif [[ -d "${OUTDIR}/${SAMPLE}/assembly/mob_recon/mob_recon_out" ]];                        then echo "Contig separation failed"
    elif [[ -f "${OUTDIR}/${SAMPLE}/assembly/flye/assembly.fasta" ]];                            then echo "Mob recon failed"
    elif [[ -f "${OUTDIR}/${SAMPLE}/reads/nanoplot/NanoStats.txt" ]];                            then echo "Flye assembly failed"
    elif [[ -f "${OUTDIR}/${SAMPLE}/reads/${SAMPLE}.read_stats.tsv" ]];                          then echo "Kraken2 / read QC failed"
    else                                                                                              echo "Read QC / NanoPlot failed"
    fi
}

on_exit() {
    local exit_code=$?
    local reason
    reason=$(detect_failure)
    if [[ "${reason}" != "Complete" ]]; then
        echo "FAILED: ${SAMPLE} — ${reason}"
        local rs="${OUTDIR}/${SAMPLE}/reads/${SAMPLE}.read_stats.tsv"
        local kp="${OUTDIR}/${SAMPLE}/reads/kraken2/${SAMPLE}_reads.kraken2.parsed"
        local data_gb num_reads num_reads_q10 median_q median_len
        if [[ -s "${rs}" ]]; then
            IFS=$'\t' read -r data_gb num_reads num_reads_q10 median_q median_len < "${rs}"
        else
            data_gb=''; num_reads=''; num_reads_q10=''; median_q=''; median_len=''
        fi
        local mg pmg g5 ms pms s5 pu
        if [[ -s "${kp}" ]]; then
            IFS=$'\t' read -r mg pmg g5 ms pms s5 pu < "${kp}"
        else
            mg=''; pmg=''; g5=''; ms=''; pms=''; s5=''; pu=''
        fi
        (
            flock -w 30 8
            if [[ ! -f "${FAILED_TSV}" ]]; then
                printf 'Sample\tFailure_Reason\tData_GB\tNum_Reads\tNum_Reads_Q10\tMedian_Q\tMedian_Read_Length\tMajor_Genus_Reads\tPct_Major_Genus_Reads\tGenera_gt5pct_Reads\tMajor_Species_Reads\tPct_Major_Species_Reads\tSpecies_gt5pct_Reads\tPct_Unclassified_Reads\n' \
                    > "${FAILED_TSV}"
            fi
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "${SAMPLE}" "${reason}" "${data_gb}" "${num_reads}" "${num_reads_q10}" \
                "${median_q}" "${median_len}" "${mg}" "${pmg}" "${g5}" "${ms}" "${pms}" "${s5}" "${pu}" \
                >> "${FAILED_TSV}"
        ) 8>"${FAILED_TSV}.lock"
    fi
    # Failed samples must still count, so the combine step fires once every
    # submitted task has finished.
    touch "${DONE_DIR}/${SAMPLE}.done"
    exit ${exit_code}
}
trap on_exit EXIT

for f in "${RESULTS_FILE}" "${FAILED_TSV}"; do
    if awk -F'\t' 'NR>1 {gsub(/\r/,"",$1); print $1}' "${f}" 2>/dev/null | grep -qxF "${SAMPLE}"; then
        echo "Task ${TASK_ID}: ${SAMPLE} already in $(basename "${f}") — skipping."
        exit 0
    fi
done

# Nextflow needs POSIX file locks for .nextflow, which many shared filesystems
# do not provide. Launch from node-local storage, keep work and results shared.
LAUNCH_DIR="${TMPDIR:-/tmp}/nf_${SAMPLE}_${JOB_ID:-0}.${TASK_ID}"
mkdir -p "${LAUNCH_DIR}"
cd "${LAUNCH_DIR}"

OPT=()
[[ -n "${MALDI_CSV}"    ]] && OPT+=(--maldi_csv    "${MALDI_CSV}")
[[ -n "${KRAKEN_DB}"    ]] && OPT+=(--kraken_db    "${KRAKEN_DB}")
[[ -n "${RESFINDER_DB}" ]] && OPT+=(--resfinder_db "${RESFINDER_DB}")
[[ -n "${AMRFINDER_DB}" ]] && OPT+=(--amrfinder_db "${AMRFINDER_DB}")
[[ -n "${CHECKM2_DB}"   ]] && OPT+=(--checkm2_db   "${CHECKM2_DB}")
[[ -n "${CHECKM2_ENV}"  ]] && OPT+=(--checkm2_env  "${CHECKM2_ENV}")
[[ -n "${FLYE_BIN}"     ]] && OPT+=(--flye_bin     "${FLYE_BIN}")

"${NEXTFLOW}" run "${PIPELINE}" \
    -profile standard \
    -w "${NXF_WORK}/${SAMPLE}" \
    --sample_name           "${SAMPLE}" \
    --input_dir             "${INPUT_DIR}" \
    --outdir                "${OUTDIR}" \
    --results_file          /dev/null \
    --skip_write_results    true \
    --max_cpus              "${CPUS}" \
    --max_memory            "${MEMORY}" \
    --kraken_memory_mapping true \
    "${OPT[@]+"${OPT[@]}"}" \
    -resume

echo "Done: ${SAMPLE} at $(date)"

# ── Append this sample to the results table ──────────────────────────────────
METRICS_FILE="${OUTDIR}/${SAMPLE}/${SAMPLE}.metrics.json"
AMR_FILE="${OUTDIR}/${SAMPLE}/${SAMPLE}.amr.tsv"
if [[ -s "${METRICS_FILE}" ]]; then
    AMR_ARG=()
    [[ -s "${AMR_FILE}" ]] && AMR_ARG=(--amr "${AMR_FILE}")
    (
        flock -w 120 7
        python3 "${PIPELINE_BIN}/write_results.py" \
            --results_file "${RESULTS_FILE}" \
            --metrics      "${METRICS_FILE}" \
            "${AMR_ARG[@]+"${AMR_ARG[@]}"}"
    ) 7>"${RESULTS_FILE}.lock"
    echo "Results updated: ${SAMPLE}"
    rm -rf "${NXF_WORK:?}/${SAMPLE}"
else
    echo "WARNING: ${METRICS_FILE} missing — leaving work dir in place for ${SAMPLE}" >&2
fi
