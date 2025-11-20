#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
METADATA_FILE="${PROJECT_ROOT}/metadata/samples.csv"
DATA_DIR="${PROJECT_ROOT}/data"
THREADS="${THREADS:-4}"
PARALLEL="${PARALLEL:-2}"
DRY_RUN="${DRY_RUN:-0}"

if ! command -v fasterq-dump >/dev/null 2>&1; then
    echo "fasterq-dump not found. Install the SRA Toolkit and ensure fasterq-dump is on PATH." >&2
    exit 1
fi

if ! command -v pigz >/dev/null 2>&1; then
    echo "pigz not found. Install pigz (parallel gzip) or adjust the script." >&2
    exit 1
fi

if [[ ! -f "${METADATA_FILE}" ]]; then
    echo "Metadata file not found: ${METADATA_FILE}" >&2
    exit 1
fi

if ! [[ "${THREADS}" =~ ^[0-9]+$ ]] || (( THREADS < 1 )); then
    echo "THREADS must be a positive integer (received: ${THREADS})" >&2
    exit 1
fi

if ! [[ "${PARALLEL}" =~ ^[0-9]+$ ]] || (( PARALLEL < 1 )); then
    echo "PARALLEL must be a positive integer (received: ${PARALLEL})" >&2
    exit 1
fi

mkdir -p "${DATA_DIR}"

download_sample() {
    local sample_id="$1"
    local fastq_path="$2"
    local dest="${PROJECT_ROOT}/${fastq_path}"
    local dest_dir
    dest_dir="$(dirname "${dest}")"
    mkdir -p "${dest_dir}"

    if [[ -f "${dest}" ]]; then
        echo "[$(date +%T)] Skipping ${sample_id}; FASTQ already exists at ${dest}"
        return 0
    fi

    echo "[$(date +%T)] Downloading ${sample_id} -> ${dest}"
    if [[ "${DRY_RUN}" == "1" ]]; then
        return 0
    fi

    fasterq-dump "${sample_id}" --outdir "${dest_dir}" --threads "${THREADS}" --progress

    local raw_fastq="${dest_dir}/${sample_id}.fastq"
    local gz_fastq="${dest_dir}/${sample_id}.fastq.gz"

    if [[ -f "${raw_fastq}" ]]; then
        pigz -f "${raw_fastq}"
    fi

    if [[ ! -f "${gz_fastq}" ]]; then
        echo "Unable to locate FASTQ emitted by fasterq-dump for ${sample_id}" >&2
        return 1
    fi

    if [[ "${gz_fastq}" != "${dest}" ]]; then
        mv "${gz_fastq}" "${dest}"
    else
        echo "[$(date +%T)] Output already at ${dest}"
    fi
}

declare -a JOBS=()

while IFS=',' read -r sample_id condition fastq_path; do
    fastq_path="${fastq_path//$'\r'/}"
    [[ -z "${sample_id}" || -z "${fastq_path}" ]] && continue

    download_sample "${sample_id}" "${fastq_path}" &
    JOBS+=($!)

    if (( ${#JOBS[@]} >= PARALLEL )); then
        wait "${JOBS[0]}" || exit 1
        JOBS=("${JOBS[@]:1}")
    fi
done < <(tail -n +2 "${METADATA_FILE}")

for pid in "${JOBS[@]}"; do
    wait "${pid}" || exit 1
done

if command -v sha256sum >/dev/null 2>&1; then
    shopt -s nullglob
    fastqs=("${PROJECT_ROOT}"/data/*.fastq.gz)
    if (( ${#fastqs[@]} )); then
        echo "Computing SHA256 checksums..."
        (cd "${PROJECT_ROOT}" && sha256sum data/*.fastq.gz > metadata/fastq_checksums.txt)
    fi
    shopt -u nullglob
fi

echo "All FASTQs downloaded to ${DATA_DIR}"
