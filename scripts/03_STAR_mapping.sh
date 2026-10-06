#!/bin/bash

# ============================================================
# ASTRA-Drug
# 03: STAR mapping + mapping QC
# ============================================================
#
#
#
#   - coordinate-sorted BAM
#   - BAM index (.bai)
#   - STAR Log.final.out
#   - samtools flagstat
#
#
# ============================================================


# ------------------------------------------------------------
# Path / Parameter
# ------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

METADATA="${PROJECT_DIR}/metadata/sample_metadata.tsv"

TRIMMED_DIR="${PROJECT_DIR}/data/processed/trimmed/paired"

STAR_INDEX="${ASTRA_STAR_INDEX:-${PROJECT_DIR}/reference/star_index}"

STAR_OUT_DIR="${PROJECT_DIR}/data/processed/star"

QC_DIR="${PROJECT_DIR}/results/qc_mapping"
MULTIQC_DIR="${QC_DIR}/multiqc"

LOG_DIR="${PROJECT_DIR}/logs/star"

THREADS="${N_THREADS:-8}"


# ------------------------------------------------------------
# ------------------------------------------------------------

mkdir -p "${STAR_OUT_DIR}"
mkdir -p "${MULTIQC_DIR}"
mkdir -p "${LOG_DIR}"


echo "============================================================"
echo "ASTRA-Drug 03: STAR mapping"
echo "============================================================"
echo "PROJECT_DIR  : ${PROJECT_DIR}"
echo "METADATA     : ${METADATA}"
echo "TRIMMED_DIR  : ${TRIMMED_DIR}"
echo "STAR_INDEX   : ${STAR_INDEX}"
echo "STAR_OUT_DIR : ${STAR_OUT_DIR}"
echo "THREADS      : ${THREADS}"
echo "============================================================"


# ------------------------------------------------------------
# ------------------------------------------------------------

if [[ ! -f "${METADATA}" ]]; then
    echo "ERROR: metadata not found:"
    echo "${METADATA}"
    exit 1
fi


if [[ ! -f "${STAR_INDEX}/Genome" ]]; then
    echo "ERROR: STAR index not found:"
    echo "${STAR_INDEX}"
    exit 1
fi


SAMPLE_COUNT=$(tail -n +2 "${METADATA}" | cut -f1 | wc -l)

echo ""
echo "Samples in metadata: ${SAMPLE_COUNT}"


if [[ "${SAMPLE_COUNT}" -ne 24 ]]; then
    echo "ERROR: Expected 24 samples."
    exit 1
fi


# ------------------------------------------------------------
# Software version
# ------------------------------------------------------------

echo ""
echo "---- Software versions ----"

echo -n "STAR: "
STAR --version

echo ""
samtools --version | head -n 2

echo ""


# ------------------------------------------------------------
# Step 1: STAR mapping
# ------------------------------------------------------------

while IFS=$'\t' read -r SAMPLE_ID REST
do

    if [[ "${SAMPLE_ID}" == "sample_id" ]]; then
        continue
    fi


    R1="${TRIMMED_DIR}/${SAMPLE_ID}_R1.trimmed.fastq.gz"
    R2="${TRIMMED_DIR}/${SAMPLE_ID}_R2.trimmed.fastq.gz"

    SAMPLE_DIR="${STAR_OUT_DIR}/${SAMPLE_ID}"

    PREFIX="${SAMPLE_DIR}/${SAMPLE_ID}."

    BAM="${PREFIX}Aligned.sortedByCoord.out.bam"
    BAI="${BAM}.bai"

    mkdir -p "${SAMPLE_DIR}"


    echo ""
    echo "============================================================"
    echo "Sample: ${SAMPLE_ID}"
    echo "============================================================"
    echo "R1: ${R1}"
    echo "R2: ${R2}"


    # --------------------------------------------------------
    # FASTQ existence check
    # --------------------------------------------------------

    if [[ ! -f "${R1}" ]]; then
        echo "ERROR: R1 not found:"
        echo "${R1}"
        exit 1
    fi

    if [[ ! -f "${R2}" ]]; then
        echo "ERROR: R2 not found:"
        echo "${R2}"
        exit 1
    fi


    # --------------------------------------------------------
    #
    # --------------------------------------------------------

    if [[ -s "${BAM}" && -s "${PREFIX}Log.final.out" ]]; then

        echo ""
        echo "Existing STAR output detected."
        echo "STAR mapping will be skipped for ${SAMPLE_ID}."

    else

        echo ""
        echo "Running STAR..."
        echo ""

        STAR \
            --runThreadN "${THREADS}" \
            --genomeDir "${STAR_INDEX}" \
            --readFilesIn "${R1}" "${R2}" \
            --readFilesCommand zcat \
            --outFileNamePrefix "${PREFIX}" \
            --outSAMtype BAM SortedByCoordinate \
            --quantMode GeneCounts \
            2>&1 | tee "${LOG_DIR}/${SAMPLE_ID}.STAR.log"


        if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
            echo ""
            echo "ERROR: STAR failed for ${SAMPLE_ID}"
            exit 1
        fi

    fi


    # --------------------------------------------------------
    # --------------------------------------------------------

    if [[ ! -s "${BAM}" ]]; then
        echo "ERROR: BAM was not created:"
        echo "${BAM}"
        exit 1
    fi


    echo ""
    echo "Checking BAM integrity..."

    samtools quickcheck -v "${BAM}"


    if [[ $? -ne 0 ]]; then
        echo "ERROR: samtools quickcheck failed for ${SAMPLE_ID}"
        exit 1
    fi


    # --------------------------------------------------------
    # BAM index
    # --------------------------------------------------------

    if [[ ! -s "${BAI}" ]]; then

        echo "Creating BAM index..."

        samtools index \
            -@ "${THREADS}" \
            "${BAM}"

    fi


    # --------------------------------------------------------
    # flagstat
    # --------------------------------------------------------

    echo "Running samtools flagstat..."

    samtools flagstat \
        -@ "${THREADS}" \
        "${BAM}" \
        > "${SAMPLE_DIR}/${SAMPLE_ID}.flagstat.txt"


    echo "Completed: ${SAMPLE_ID}"


done < "${METADATA}"


# ------------------------------------------------------------
# ------------------------------------------------------------

BAM_COUNT=$(find "${STAR_OUT_DIR}" \
    -type f \
    -name "*.Aligned.sortedByCoord.out.bam" | wc -l)


echo ""
echo "============================================================"
echo "STAR mapping completed"
echo "BAM count: ${BAM_COUNT}"
echo "============================================================"


if [[ "${BAM_COUNT}" -ne 24 ]]; then
    echo "ERROR: Expected 24 BAM files."
    exit 1
fi


# ------------------------------------------------------------
# Step 3: MultiQC
# ------------------------------------------------------------

echo ""
echo "============================================================"
echo "Running MultiQC for STAR + samtools"
echo "============================================================"


multiqc \
    "${STAR_OUT_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "ASTRA_Drug_STAR_MultiQC.html" \
    --force \
    2>&1 | tee "${LOG_DIR}/STAR_multiqc.log"


if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    echo "ERROR: MultiQC failed."
    exit 1
fi


echo ""
echo "============================================================"
echo "03 STAR mapping successfully completed"
echo ""
echo "BAM directory:"
echo "${STAR_OUT_DIR}"
echo ""
echo "MultiQC report:"
echo "${MULTIQC_DIR}/ASTRA_Drug_STAR_MultiQC.html"
echo "============================================================"
