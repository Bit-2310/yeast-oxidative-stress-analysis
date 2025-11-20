#!/usr/bin/env bash

set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF_DIR="${PROJECT_ROOT}/data/genome"
INDEX_NAME="${INDEX_NAME:-sacCer3}"
FASTA_URL="${FASTA_URL:-https://ftp.ensembl.org/pub/release-111/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz}"
GTF_URL="${GTF_URL:-https://ftp.ensembl.org/pub/release-111/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.111.gtf.gz}"
BUILD_THREADS="${BUILD_THREADS:-8}"

if ! command -v hisat2-build >/dev/null 2>&1; then
    echo "hisat2-build not found. Install HISAT2 and ensure hisat2-build is on PATH." >&2
    exit 1
fi

if ! [[ "${BUILD_THREADS}" =~ ^[0-9]+$ ]] || (( BUILD_THREADS < 1 )); then
    echo "BUILD_THREADS must be a positive integer (received: ${BUILD_THREADS})" >&2
    exit 1
fi

download_with_resume() {
    local url="$1"
    local output="$2"
    if command -v aria2c >/dev/null 2>&1; then
        aria2c -x 8 -s 8 -o "${output}" "${url}"
    else
        curl -L -C - "${url}" -o "${output}"
    fi
}

verify_not_html() {
    local file="$1"
    if file "${file}" | grep -qi 'HTML'; then
        echo "Downloaded file ${file} appears to be HTML (download likely failed). Check the URL or network access." >&2
        exit 1
    fi
}
mkdir -p "${REF_DIR}"
cd "${REF_DIR}"

GENOME_FASTA="${REF_DIR}/${INDEX_NAME}.fa"
ANNOTATION_GTF="${PROJECT_ROOT}/data/Saccharomyces_cerevisiae.R64-1-1.gtf"

if [[ ! -f "${GENOME_FASTA}" ]]; then
    echo "Downloading genome FASTA..."
    download_with_resume "${FASTA_URL}" "${INDEX_NAME}.fa.gz"
    verify_not_html "${INDEX_NAME}.fa.gz"
    if command -v pigz >/dev/null 2>&1; then
        pigz -df "${INDEX_NAME}.fa.gz"
    else
        gunzip -f "${INDEX_NAME}.fa.gz"
    fi
fi

if [[ ! -f "${ANNOTATION_GTF}" ]]; then
    echo "Downloading annotation GTF..."
    tmp_gtf="${PROJECT_ROOT}/data/$(basename "${GTF_URL}")"
    download_with_resume "${GTF_URL}" "${tmp_gtf}"
    verify_not_html "${tmp_gtf}"
    if command -v pigz >/dev/null 2>&1; then
        pigz -df "${tmp_gtf}"
    else
        gunzip -f "${tmp_gtf}"
    fi
    mv "${tmp_gtf%.gz}" "${ANNOTATION_GTF}"
fi

if ls "${REF_DIR}/${INDEX_NAME}".*.ht2 >/dev/null 2>&1; then
    echo "HISAT2 index already present in ${REF_DIR}"
else
    echo "Building HISAT2 index..."
    hisat2-build -p "${BUILD_THREADS}" "${GENOME_FASTA}" "${INDEX_NAME}"
fi

echo "Reference genome prepared under ${REF_DIR}"
echo "Annotation saved to ${ANNOTATION_GTF}"
