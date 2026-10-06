#!/bin/bash

# ============================================================
# ASTRA-Drug: Adapter trimming + post-trimming QC
# ============================================================
#
#
#
#
#
# ============================================================


# ------------------------------------------------------------
# Path / Parameter
# ------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

FASTQ_DIR="${PROJECT_DIR}/data/raw/fastq"
METADATA="${PROJECT_DIR}/metadata/sample_metadata.tsv"

TRIM_DIR="${PROJECT_DIR}/data/processed/trimmed"
PAIRED_DIR="${TRIM_DIR}/paired"
UNPAIRED_DIR="${TRIM_DIR}/unpaired"

QC_DIR="${PROJECT_DIR}/results/qc_trimmed"
FASTQC_DIR="${QC_DIR}/fastqc"
MULTIQC_DIR="${QC_DIR}/multiqc"

LOG_DIR="${PROJECT_DIR}/logs/trimmomatic"

THREADS="${N_THREADS:-8}"

MINLEN=36

# Trimmomatic ILLUMINACLIP parameters
SEED_MISMATCHES=2
PALINDROME_CLIP_THRESHOLD=30
SIMPLE_CLIP_THRESHOLD=10


# ------------------------------------------------------------
# ------------------------------------------------------------

mkdir -p "${PAIRED_DIR}"
mkdir -p "${UNPAIRED_DIR}"
mkdir -p "${FASTQC_DIR}"
mkdir -p "${MULTIQC_DIR}"
mkdir -p "${LOG_DIR}"


# ------------------------------------------------------------
# ------------------------------------------------------------

ADAPTER_FA=$(find "${CONDA_PREFIX}/share" \
    -type f \
    -name "TruSeq3-PE-2.fa" \
    2>/dev/null | head -n 1)


if [[ -z "${ADAPTER_FA}" ]]; then

    echo "ERROR: TruSeq3-PE-2.fa was not found."
    echo "CONDA_PREFIX=${CONDA_PREFIX}"
    exit 1

fi


echo "========================================"
echo "ASTRA-Drug adapter trimming"
echo "========================================"
echo "FASTQ_DIR    : ${FASTQ_DIR}"
echo "METADATA     : ${METADATA}"
echo "TRIM_DIR     : ${TRIM_DIR}"
echo "ADAPTER_FA   : ${ADAPTER_FA}"
echo "THREADS      : ${THREADS}"
echo "MINLEN       : ${MINLEN}"
echo "========================================"


# ------------------------------------------------------------
# ------------------------------------------------------------

if [[ ! -f "${METADATA}" ]]; then
    echo "ERROR: metadata file not found:"
    echo "${METADATA}"
    exit 1
fi


SAMPLE_COUNT=$(tail -n +2 "${METADATA}" | cut -f1 | wc -l)

echo ""
echo "Samples in metadata: ${SAMPLE_COUNT}"


if [[ "${SAMPLE_COUNT}" -ne 24 ]]; then

    echo "ERROR: Expected 24 samples."
    echo "Found: ${SAMPLE_COUNT}"
    exit 1

fi


# ------------------------------------------------------------
# Step 1: Trimmomatic
# ------------------------------------------------------------

echo ""
echo "========================================"
echo "Step 1: Trimmomatic"
echo "========================================"


while IFS=$'\t' read -r SAMPLE_ID REST
do

    if [[ "${SAMPLE_ID}" == "sample_id" ]]; then
        continue
    fi


    R1="${FASTQ_DIR}/${SAMPLE_ID}_combined_R1.fastq.gz"
    R2="${FASTQ_DIR}/${SAMPLE_ID}_combined_R2.fastq.gz"

    P1="${PAIRED_DIR}/${SAMPLE_ID}_R1.trimmed.fastq.gz"
    U1="${UNPAIRED_DIR}/${SAMPLE_ID}_R1.unpaired.fastq.gz"

    P2="${PAIRED_DIR}/${SAMPLE_ID}_R2.trimmed.fastq.gz"
    U2="${UNPAIRED_DIR}/${SAMPLE_ID}_R2.unpaired.fastq.gz"

    SAMPLE_LOG="${LOG_DIR}/${SAMPLE_ID}.trimmomatic.log"


    echo ""
    echo "----------------------------------------"
    echo "Sample: ${SAMPLE_ID}"
    echo "----------------------------------------"


    if [[ ! -f "${R1}" ]]; then
        echo "ERROR: R1 not found: ${R1}"
        exit 1
    fi

    if [[ ! -f "${R2}" ]]; then
        echo "ERROR: R2 not found: ${R2}"
        exit 1
    fi


    trimmomatic PE \
        -threads "${THREADS}" \
        -phred33 \
        "${R1}" \
        "${R2}" \
        "${P1}" \
        "${U1}" \
        "${P2}" \
        "${U2}" \
        "ILLUMINACLIP:${ADAPTER_FA}:${SEED_MISMATCHES}:${PALINDROME_CLIP_THRESHOLD}:${SIMPLE_CLIP_THRESHOLD}" \
        "MINLEN:${MINLEN}" \
        2>&1 | tee "${SAMPLE_LOG}"


    if [[ ${PIPESTATUS[0]} -ne 0 ]]; then

        echo ""
        echo "ERROR: Trimmomatic failed for ${SAMPLE_ID}"
        exit 1

    fi


done < "${METADATA}"


echo ""
echo "Trimmomatic completed for all samples."


# ------------------------------------------------------------
# ------------------------------------------------------------

PAIRED_COUNT=$(find "${PAIRED_DIR}" \
    -maxdepth 1 \
    -type f \
    -name "*.trimmed.fastq.gz" | wc -l)


echo "Trimmed paired FASTQ files: ${PAIRED_COUNT}"


if [[ "${PAIRED_COUNT}" -ne 48 ]]; then

    echo "ERROR: Expected 48 paired trimmed FASTQ files."
    exit 1

fi


# ------------------------------------------------------------
# Step 2: FastQC
# ------------------------------------------------------------

echo ""
echo "========================================"
echo "Step 2: FastQC of trimmed paired reads"
echo "========================================"


fastqc \
    --threads "${THREADS}" \
    --outdir "${FASTQC_DIR}" \
    "${PAIRED_DIR}"/*.trimmed.fastq.gz \
    2>&1 | tee "${LOG_DIR}/posttrim_fastqc.log"


if [[ ${PIPESTATUS[0]} -ne 0 ]]; then

    echo "ERROR: FastQC failed."
    exit 1

fi


# ------------------------------------------------------------
# Step 3: MultiQC
# ------------------------------------------------------------

echo ""
echo "========================================"
echo "Step 3: MultiQC"
echo "========================================"


multiqc \
    "${FASTQC_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "ASTRA_Drug_trimmed_MultiQC.html" \
    --force \
    2>&1 | tee "${LOG_DIR}/posttrim_multiqc.log"


if [[ ${PIPESTATUS[0]} -ne 0 ]]; then

    echo "ERROR: MultiQC failed."
    exit 1

fi


echo ""
echo "========================================"
echo "02 trimming + QC completed successfully"
echo ""
echo "Trimmed paired reads:"
echo "${PAIRED_DIR}"
echo ""
echo "MultiQC report:"
echo "${MULTIQC_DIR}/ASTRA_Drug_trimmed_MultiQC.html"
echo "========================================"
