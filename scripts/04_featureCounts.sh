#!/bin/bash

# ============================================================
# ASTRA-Drug
# 04: featureCounts
# ============================================================
#
#
#
# 1. published_like
#    -M -O -p, unstranded (-s 0)
#
# 2. strand_aware
#    -M -O -p, reverse stranded (-s 2)
#
#
#   - clean gene count matrix
#   - ERCC total counts
#   - featureCounts summary
#   - MultiQC
#
# ============================================================


# ------------------------------------------------------------
# Path / Parameter
# ------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

METADATA="${PROJECT_DIR}/metadata/sample_metadata.tsv"

STAR_DIR="${PROJECT_DIR}/data/processed/star"

GTF="${ASTRA_GTF:-${PROJECT_DIR}/reference/mm10_ERCC_EGFP_mcherry.gtf}"

COUNT_BASE="${PROJECT_DIR}/data/processed/counts"

PUB_DIR="${COUNT_BASE}/published_like"
STRAND_DIR="${COUNT_BASE}/strand_aware"

QC_BASE="${PROJECT_DIR}/results/qc_featureCounts"

PUB_QC_DIR="${QC_BASE}/published_like"
STRAND_QC_DIR="${QC_BASE}/strand_aware"

LOG_DIR="${PROJECT_DIR}/logs/featureCounts"

THREADS="${N_THREADS:-8}"


# ------------------------------------------------------------
# ------------------------------------------------------------

mkdir -p "${PUB_DIR}"
mkdir -p "${STRAND_DIR}"

mkdir -p "${PUB_QC_DIR}"
mkdir -p "${STRAND_QC_DIR}"

mkdir -p "${LOG_DIR}"


echo "============================================================"
echo "ASTRA-Drug 04: featureCounts"
echo "============================================================"
echo "PROJECT_DIR : ${PROJECT_DIR}"
echo "STAR_DIR    : ${STAR_DIR}"
echo "GTF         : ${GTF}"
echo "THREADS     : ${THREADS}"
echo "============================================================"


# ------------------------------------------------------------
# ------------------------------------------------------------

if [[ ! -f "${METADATA}" ]]; then
    echo "ERROR: metadata not found:"
    echo "${METADATA}"
    exit 1
fi


if [[ ! -f "${GTF}" ]]; then
    echo "ERROR: GTF not found:"
    echo "${GTF}"
    exit 1
fi


# ------------------------------------------------------------
# ------------------------------------------------------------

BAMS=()
SAMPLE_IDS=()


while IFS=$'\t' read -r SAMPLE_ID REST
do

    if [[ "${SAMPLE_ID}" == "sample_id" ]]; then
        continue
    fi

    BAM="${STAR_DIR}/${SAMPLE_ID}/${SAMPLE_ID}.Aligned.sortedByCoord.out.bam"

    if [[ ! -s "${BAM}" ]]; then
        echo "ERROR: BAM not found or empty:"
        echo "${BAM}"
        exit 1
    fi

    BAMS+=("${BAM}")
    SAMPLE_IDS+=("${SAMPLE_ID}")

done < "${METADATA}"


BAM_COUNT=${#BAMS[@]}

echo ""
echo "BAM count: ${BAM_COUNT}"


if [[ "${BAM_COUNT}" -ne 24 ]]; then
    echo "ERROR: Expected 24 BAM files."
    exit 1
fi


echo ""
echo "BAM order:"
printf '%s\n' "${SAMPLE_IDS[@]}"


# ------------------------------------------------------------
# ------------------------------------------------------------

ERCC_GTF_LINES=$(grep -c "ERCC" "${GTF}" || true)

echo ""
echo "ERCC-containing GTF lines: ${ERCC_GTF_LINES}"

if [[ "${ERCC_GTF_LINES}" -eq 0 ]]; then
    echo "ERROR: No ERCC annotation found in GTF."
    exit 1
fi


# ============================================================
# Step 1
# published-like counts
# ============================================================

PUB_RAW="${PUB_DIR}/featureCounts_published_like.txt"

echo ""
echo "============================================================"
echo "Step 1: published-like featureCounts"
echo "============================================================"
echo ""

featureCounts \
    -T "${THREADS}" \
    -a "${GTF}" \
    -o "${PUB_RAW}" \
    -t exon \
    -g gene_id \
    -M \
    -O \
    -p \
    --countReadPairs \
    -s 0 \
    "${BAMS[@]}" \
    2>&1 | tee "${LOG_DIR}/published_like.featureCounts.log"


if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    echo "ERROR: published-like featureCounts failed."
    exit 1
fi


# ============================================================
# Step 2
# strand-aware counts
# ============================================================

STRAND_RAW="${STRAND_DIR}/featureCounts_strand_aware.txt"

echo ""
echo "============================================================"
echo "Step 2: strand-aware featureCounts"
echo "============================================================"
echo ""

featureCounts \
    -T "${THREADS}" \
    -a "${GTF}" \
    -o "${STRAND_RAW}" \
    -t exon \
    -g gene_id \
    -M \
    -O \
    -p \
    --countReadPairs \
    -s 2 \
    "${BAMS[@]}" \
    2>&1 | tee "${LOG_DIR}/strand_aware.featureCounts.log"


if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    echo "ERROR: strand-aware featureCounts failed."
    exit 1
fi


# ============================================================
# Function:
# ============================================================

make_clean_matrix () {

    INPUT_FILE="$1"
    OUTPUT_FILE="$2"

    awk '
    BEGIN {
        FS = OFS = "\t"
    }

    /^#/ {
        next
    }

    $1 == "Geneid" {

        printf "gene_id"

        for (i = 7; i <= NF; i++) {

            x = $i

            gsub(/^.*\//, "", x)

            sub(/\.Aligned\.sortedByCoord\.out\.bam$/, "", x)

            printf OFS x
        }

        printf "\n"

        next
    }

    {
        printf $1

        for (i = 7; i <= NF; i++) {
            printf OFS $i
        }

        printf "\n"
    }

    ' "${INPUT_FILE}" > "${OUTPUT_FILE}"
}


# ============================================================
# Function:
# ============================================================

make_ercc_totals () {

    INPUT_FILE="$1"
    OUTPUT_FILE="$2"

    awk '
    BEGIN {
        FS = OFS = "\t"
    }

    NR == 1 {

        NCOLS = NF

        for (i = 2; i <= NF; i++) {
            SAMPLE[i] = $i
        }

        next
    }

    $1 ~ /^ERCC/ {

        for (i = 2; i <= NF; i++) {
            ERCC[i] += $i
        }
    }

    END {

        print "sample_id", "ERCC_total"

        for (i = 2; i <= NCOLS; i++) {
            print SAMPLE[i], ERCC[i] + 0
        }
    }

    ' "${INPUT_FILE}" > "${OUTPUT_FILE}"
}


# ============================================================
# Step 3
# clean matrices
# ============================================================

echo ""
echo "============================================================"
echo "Step 3: Creating clean count matrices"
echo "============================================================"


PUB_MATRIX="${PUB_DIR}/gene_counts_published_like.tsv"

STRAND_MATRIX="${STRAND_DIR}/gene_counts_strand_aware.tsv"


make_clean_matrix \
    "${PUB_RAW}" \
    "${PUB_MATRIX}"


make_clean_matrix \
    "${STRAND_RAW}" \
    "${STRAND_MATRIX}"


# ============================================================
# Step 4
# ERCC totals
# ============================================================

echo ""
echo "============================================================"
echo "Step 4: Calculating ERCC totals"
echo "============================================================"


make_ercc_totals \
    "${PUB_MATRIX}" \
    "${PUB_DIR}/ercc_totals.tsv"


make_ercc_totals \
    "${STRAND_MATRIX}" \
    "${STRAND_DIR}/ercc_totals.tsv"


# ------------------------------------------------------------
#
# gene_id + 24 samples = 25 columns
# ------------------------------------------------------------

PUB_COLS=$(head -n 1 "${PUB_MATRIX}" | awk -F'\t' '{print NF}')

STRAND_COLS=$(head -n 1 "${STRAND_MATRIX}" | awk -F'\t' '{print NF}')


echo ""
echo "published-like matrix columns : ${PUB_COLS}"
echo "strand-aware matrix columns   : ${STRAND_COLS}"


if [[ "${PUB_COLS}" -ne 25 ]]; then
    echo "ERROR: published-like matrix should have 25 columns."
    exit 1
fi


if [[ "${STRAND_COLS}" -ne 25 ]]; then
    echo "ERROR: strand-aware matrix should have 25 columns."
    exit 1
fi


# ============================================================
# Step 5
# MultiQC
# ============================================================

echo ""
echo "============================================================"
echo "Step 5: MultiQC"
echo "============================================================"


multiqc \
    "${PUB_DIR}" \
    --outdir "${PUB_QC_DIR}" \
    --filename "ASTRA_Drug_featureCounts_published_like_MultiQC.html" \
    --force \
    2>&1 | tee "${LOG_DIR}/published_like.MultiQC.log"


if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    echo "ERROR: published-like MultiQC failed."
    exit 1
fi


multiqc \
    "${STRAND_DIR}" \
    --outdir "${STRAND_QC_DIR}" \
    --filename "ASTRA_Drug_featureCounts_strand_aware_MultiQC.html" \
    --force \
    2>&1 | tee "${LOG_DIR}/strand_aware.MultiQC.log"


if [[ ${PIPESTATUS[0]} -ne 0 ]]; then
    echo "ERROR: strand-aware MultiQC failed."
    exit 1
fi


# ============================================================
# Finished
# ============================================================

echo ""
echo "============================================================"
echo "04 featureCounts successfully completed"
echo "============================================================"

echo ""
echo "Published-like matrix:"
echo "${PUB_MATRIX}"

echo ""
echo "Published-like ERCC totals:"
echo "${PUB_DIR}/ercc_totals.tsv"

echo ""
echo "Strand-aware matrix:"
echo "${STRAND_MATRIX}"

echo ""
echo "Strand-aware ERCC totals:"
echo "${STRAND_DIR}/ercc_totals.tsv"

echo ""
echo "Published-like MultiQC:"
echo "${PUB_QC_DIR}/ASTRA_Drug_featureCounts_published_like_MultiQC.html"

echo ""
echo "Strand-aware MultiQC:"
echo "${STRAND_QC_DIR}/ASTRA_Drug_featureCounts_strand_aware_MultiQC.html"

echo "============================================================"
