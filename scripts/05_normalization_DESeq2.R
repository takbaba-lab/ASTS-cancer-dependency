#!/usr/bin/env Rscript

# ============================================================
# ASTRA-Drug
# 05: Counting-method comparison and dual normalization
# ============================================================
#
#
#
#
#
#    A. ERCC-sum normalization
#
#    B. standard DESeq2 normalization
#
#
#    M_ORX_DHT vs M_ORX     : male DHT effect
#    F_OVX_DHT vs F_OVX     : female DHT effect
#    M_ORX vs M_SHAM        : androgen withdrawal
#    M_ORX_E2 vs M_ORX      : male E2 effect
#    F_OVX_E2 vs F_OVX      : female E2 effect
#    F_SHAM vs M_SHAM       : physiological sex difference
#    F_OVX vs F_SHAM        : ovarian hormone withdrawal
#
#
# ============================================================


# ============================================================
# Path / Parameter
# ============================================================

PROJECT_DIR <- normalizePath(
    Sys.getenv("ASTRA_PROJECT_DIR", unset = getwd()),
    mustWork = FALSE
)

METADATA_FILE <- file.path(
    PROJECT_DIR,
    "metadata/sample_metadata.tsv"
)

PUB_COUNT_FILE <- file.path(
    PROJECT_DIR,
    "data/processed/counts/published_like/gene_counts_published_like.tsv"
)

STRAND_COUNT_FILE <- file.path(
    PROJECT_DIR,
    "data/processed/counts/strand_aware/gene_counts_strand_aware.tsv"
)

PUB_ERCC_FILE <- file.path(
    PROJECT_DIR,
    "data/processed/counts/published_like/ercc_totals.tsv"
)

STRAND_ERCC_FILE <- file.path(
    PROJECT_DIR,
    "data/processed/counts/strand_aware/ercc_totals.tsv"
)


OUT_DIR <- file.path(
    PROJECT_DIR,
    "results/05_normalization"
)

NORM_DIR <- file.path(
    PROJECT_DIR,
    "data/processed/normalized"
)

STANDARD_DE_DIR <- file.path(
    OUT_DIR,
    "DESeq2_standard"
)

ERCC_DE_DIR <- file.path(
    OUT_DIR,
    "DESeq2_ERCC"
)


MIN_TOTAL_COUNT <- 10


GROUP_LEVELS <- c(
    "M_SHAM",
    "M_ORX",
    "M_ORX_DHT",
    "M_ORX_E2",
    "F_SHAM",
    "F_OVX",
    "F_OVX_DHT",
    "F_OVX_E2"
)


# ============================================================
# Package
# ============================================================

suppressPackageStartupMessages({
    library(DESeq2)
})


# ============================================================
# Output directory
# ============================================================

dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    NORM_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    STANDARD_DE_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    ERCC_DE_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# Helper functions
# ============================================================

read_count_matrix <- function(file) {

    x <- read.delim(
        file,
        check.names = FALSE,
        stringsAsFactors = FALSE
    )

    if (!"gene_id" %in% colnames(x)) {
        stop("gene_id column not found: ", file)
    }

    if (anyDuplicated(x$gene_id) > 0) {
        stop("Duplicated gene_id found: ", file)
    }

    rownames(x) <- x$gene_id
    x$gene_id <- NULL

    x <- as.matrix(x)

    storage.mode(x) <- "integer"

    return(x)
}


read_ercc_totals <- function(file) {

    x <- read.delim(
        file,
        stringsAsFactors = FALSE
    )

    if (!all(c("sample_id", "ERCC_total") %in% colnames(x))) {
        stop(
            "ERCC file must contain sample_id and ERCC_total: ",
            file
        )
    }

    return(x)
}


is_external_gene <- function(gene_id) {

    grepl(
        "^ERCC",
        gene_id,
        ignore.case = TRUE
    ) |
    grepl(
        "^(EGFP|mCherry)$",
        gene_id,
        ignore.case = TRUE
    )
}


write_matrix <- function(mat, file) {

    out <- data.frame(
        gene_id = rownames(mat),
        mat,
        check.names = FALSE
    )

    write.table(
        out,
        file = file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}


# ============================================================
# Step 1
# Read metadata
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 1: Reading metadata\n")
cat("============================================================\n")


metadata <- read.delim(
    METADATA_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


if (!all(c("sample_id", "group") %in% colnames(metadata))) {
    stop("metadata must contain sample_id and group columns.")
}


if (nrow(metadata) != 24) {
    stop(
        "Expected 24 samples in metadata, found ",
        nrow(metadata)
    )
}


metadata$group <- factor(
    metadata$group,
    levels = GROUP_LEVELS
)


if (any(is.na(metadata$group))) {
    stop("Unknown group found in metadata.")
}


rownames(metadata) <- metadata$sample_id


cat("Samples:", nrow(metadata), "\n")


# ============================================================
# Step 2
# Read count matrices
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 2: Reading count matrices\n")
cat("============================================================\n")


counts_pub <- read_count_matrix(
    PUB_COUNT_FILE
)

counts_strand <- read_count_matrix(
    STRAND_COUNT_FILE
)


cat(
    "published-like genes:",
    nrow(counts_pub),
    "\n"
)

cat(
    "strand-aware genes:",
    nrow(counts_strand),
    "\n"
)


if (!all(metadata$sample_id %in% colnames(counts_pub))) {
    stop("Some metadata samples are missing in published-like matrix.")
}

if (!all(metadata$sample_id %in% colnames(counts_strand))) {
    stop("Some metadata samples are missing in strand-aware matrix.")
}


counts_pub <- counts_pub[
    ,
    metadata$sample_id,
    drop = FALSE
]

counts_strand <- counts_strand[
    ,
    metadata$sample_id,
    drop = FALSE
]


# ============================================================
# Step 3
# Compare published-like vs strand-aware counts
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 3: Comparing counting methods\n")
cat("============================================================\n")


common_genes <- intersect(
    rownames(counts_pub),
    rownames(counts_strand)
)


pub_common <- counts_pub[
    common_genes,
    ,
    drop = FALSE
]

strand_common <- counts_strand[
    common_genes,
    ,
    drop = FALSE
]


comparison_list <- lapply(
    metadata$sample_id,
    function(sample_id) {

        x <- pub_common[, sample_id]
        y <- strand_common[, sample_id]

        data.frame(
            sample_id = sample_id,

            published_total = sum(x),

            strand_total = sum(y),

            strand_over_published =
                sum(y) / sum(x),

            pearson_log1p = cor(
                log1p(x),
                log1p(y),
                method = "pearson"
            ),

            spearman = cor(
                x,
                y,
                method = "spearman"
            )
        )
    }
)


count_method_comparison <- do.call(
    rbind,
    comparison_list
)


write.table(
    count_method_comparison,
    file = file.path(
        OUT_DIR,
        "counting_method_sample_comparison.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


cat(
    "Mean Pearson(log1p):",
    mean(count_method_comparison$pearson_log1p),
    "\n"
)

cat(
    "Mean Spearman:",
    mean(count_method_comparison$spearman),
    "\n"
)


# ============================================================
# Step 4
# Compare ERCC totals
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 4: Comparing ERCC totals\n")
cat("============================================================\n")


ercc_pub <- read_ercc_totals(
    PUB_ERCC_FILE
)

ercc_strand <- read_ercc_totals(
    STRAND_ERCC_FILE
)


ercc_compare <- merge(
    ercc_pub,
    ercc_strand,
    by = "sample_id",
    suffixes = c(
        "_published",
        "_strand"
    )
)


ercc_compare$strand_over_published <-
    ercc_compare$ERCC_total_strand /
    ercc_compare$ERCC_total_published


ercc_compare <- ercc_compare[
    match(
        metadata$sample_id,
        ercc_compare$sample_id
    ),
]


write.table(
    ercc_compare,
    file = file.path(
        OUT_DIR,
        "ERCC_total_comparison.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# Step 5
# Define primary endogenous count matrix
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 5: Preparing primary endogenous matrix\n")
cat("============================================================\n")


external_gene <- is_external_gene(
    rownames(counts_strand)
)


cat(
    "External genes removed:",
    sum(external_gene),
    "\n"
)


counts_endogenous <- counts_strand[
    !external_gene,
    ,
    drop = FALSE
]


write_matrix(
    counts_endogenous,
    file.path(
        NORM_DIR,
        "counts_strand_aware_endogenous_raw.tsv"
    )
)


# ============================================================
# Step 6
# Standard DESeq2 normalization
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 6: Standard DESeq2 normalization\n")
cat("============================================================\n")


dds_standard_all <- DESeqDataSetFromMatrix(
    countData = counts_endogenous,
    colData = metadata,
    design = ~ group
)


dds_standard_all <- estimateSizeFactors(
    dds_standard_all
)


standard_sf <- sizeFactors(
    dds_standard_all
)


standard_norm <- counts(
    dds_standard_all,
    normalized = TRUE
)


write_matrix(
    standard_norm,
    file.path(
        NORM_DIR,
        "counts_standard_DESeq2_normalized.tsv"
    )
)


# ============================================================
# Step 7
# ERCC-sum normalization
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 7: ERCC-sum normalization\n")
cat("============================================================\n")


ercc_primary <- ercc_strand[
    match(
        metadata$sample_id,
        ercc_strand$sample_id
    ),
]


if (
    !identical(
        ercc_primary$sample_id,
        metadata$sample_id
    )
) {
    stop("ERCC sample order does not match metadata.")
}


if (any(ercc_primary$ERCC_total <= 0)) {
    stop("ERCC_total must be > 0.")
}


# ------------------------------------------------------------
# ERCC size factor
#
# size factor =
# ERCC_total / geometric mean(ERCC_total)
#
# ------------------------------------------------------------

ercc_geometric_mean <- exp(
    mean(
        log(
            ercc_primary$ERCC_total
        )
    )
)


ercc_sf <-
    ercc_primary$ERCC_total /
    ercc_geometric_mean


names(ercc_sf) <-
    ercc_primary$sample_id


# normalized count
ercc_norm <- sweep(
    counts_endogenous,
    2,
    ercc_sf,
    "/"
)


write_matrix(
    ercc_norm,
    file.path(
        NORM_DIR,
        "counts_ERCC_sum_normalized.tsv"
    )
)


# ============================================================
# Step 8
# Sample-level scaling summary
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 8: Sample scaling summary\n")
cat("============================================================\n")


raw_endogenous_total <-
    colSums(counts_endogenous)

standard_norm_total <-
    colSums(standard_norm)

ercc_norm_total <-
    colSums(ercc_norm)


sample_scaling <- data.frame(
    sample_id = metadata$sample_id,
    group = as.character(metadata$group),

    raw_endogenous_total =
        raw_endogenous_total[metadata$sample_id],

    ERCC_total =
        ercc_primary$ERCC_total,

    ERCC_size_factor =
        ercc_sf[metadata$sample_id],

    DESeq2_size_factor =
        standard_sf[metadata$sample_id],

    ERCC_normalized_endogenous_total =
        ercc_norm_total[metadata$sample_id],

    DESeq2_normalized_endogenous_total =
        standard_norm_total[metadata$sample_id]
)


write.table(
    sample_scaling,
    file = file.path(
        OUT_DIR,
        "sample_scaling_factors.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# Step 9
# Group-level global output summary
# ============================================================

group_list <- split(
    sample_scaling,
    sample_scaling$group
)


group_summary_list <- lapply(
    group_list,
    function(x) {

        data.frame(
            group = unique(x$group),

            n = nrow(x),

            mean_ERCC_normalized_total =
                mean(
                    x$ERCC_normalized_endogenous_total
                ),

            sd_ERCC_normalized_total =
                sd(
                    x$ERCC_normalized_endogenous_total
                ),

            CV_ERCC_normalized_total =
                sd(
                    x$ERCC_normalized_endogenous_total
                ) /
                mean(
                    x$ERCC_normalized_endogenous_total
                ),

            mean_standard_normalized_total =
                mean(
                    x$DESeq2_normalized_endogenous_total
                )
        )
    }
)


group_summary <- do.call(
    rbind,
    group_summary_list
)


group_summary <- group_summary[
    match(
        GROUP_LEVELS,
        group_summary$group
    ),
]


write.table(
    group_summary,
    file = file.path(
        OUT_DIR,
        "group_global_output_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# Step 10
# Prepare DESeq2 datasets
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 10: Running DESeq2 models\n")
cat("============================================================\n")


keep <- rowSums(
    counts_endogenous
) >= MIN_TOTAL_COUNT


cat(
    "Genes retained for DESeq2:",
    sum(keep),
    "/",
    length(keep),
    "\n"
)


# standard
dds_standard <- dds_standard_all[
    keep,
]


dds_standard <- DESeq(
    dds_standard,
    quiet = TRUE
)


# ERCC-normalized
dds_ercc <- DESeqDataSetFromMatrix(
    countData = counts_endogenous[
        keep,
        ,
        drop = FALSE
    ],
    colData = metadata,
    design = ~ group
)


sizeFactors(dds_ercc) <-
    ercc_sf[
        colnames(dds_ercc)
    ]


dds_ercc <- DESeq(
    dds_ercc,
    quiet = TRUE
)


# ============================================================
# Step 11
# Contrasts
# ============================================================

contrasts <- list(

    M_DHT = c(
        "M_ORX_DHT",
        "M_ORX"
    ),

    F_DHT = c(
        "F_OVX_DHT",
        "F_OVX"
    ),

    M_androgen_withdrawal = c(
        "M_ORX",
        "M_SHAM"
    ),

    M_E2 = c(
        "M_ORX_E2",
        "M_ORX"
    ),

    F_E2 = c(
        "F_OVX_E2",
        "F_OVX"
    ),

    sex_control = c(
        "F_SHAM",
        "M_SHAM"
    ),

    F_ovarian_withdrawal = c(
        "F_OVX",
        "F_SHAM"
    )
)


extract_result <- function(
    dds,
    contrast_name,
    numerator,
    denominator,
    out_dir
) {

    res <- results(
        dds,
        contrast = c(
            "group",
            numerator,
            denominator
        )
    )

    df <- as.data.frame(res)

    df$gene_id <- rownames(df)

    df <- df[
        ,
        c(
            "gene_id",
            "baseMean",
            "log2FoldChange",
            "lfcSE",
            "stat",
            "pvalue",
            "padj"
        )
    ]


    ord <- order(
        is.na(df$padj),
        df$padj,
        df$pvalue
    )

    df <- df[ord, ]


    outfile <- file.path(
        out_dir,
        paste0(
            contrast_name,
            ".tsv"
        )
    )


    write.table(
        df,
        file = outfile,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )


    finite_lfc <- df$log2FoldChange[
        is.finite(df$log2FoldChange)
    ]


    sig <- !is.na(df$padj) &
        df$padj < 0.05


    data.frame(
        contrast = contrast_name,

        numerator = numerator,

        denominator = denominator,

        genes_tested =
            length(finite_lfc),

        median_log2FC =
            median(
                finite_lfc
            ),

        fraction_log2FC_negative =
            mean(
                finite_lfc < 0
            ),

        fraction_log2FC_positive =
            mean(
                finite_lfc > 0
            ),

        significant_up =
            sum(
                sig &
                df$log2FoldChange > 0,
                na.rm = TRUE
            ),

        significant_down =
            sum(
                sig &
                df$log2FoldChange < 0,
                na.rm = TRUE
            )
    )
}


standard_summary <- list()
ercc_summary <- list()


for (contrast_name in names(contrasts)) {

    numerator <-
        contrasts[[contrast_name]][1]

    denominator <-
        contrasts[[contrast_name]][2]


    cat(
        "Contrast:",
        numerator,
        "vs",
        denominator,
        "\n"
    )


    standard_summary[[contrast_name]] <-
        extract_result(
            dds_standard,
            contrast_name,
            numerator,
            denominator,
            STANDARD_DE_DIR
        )


    ercc_summary[[contrast_name]] <-
        extract_result(
            dds_ercc,
            contrast_name,
            numerator,
            denominator,
            ERCC_DE_DIR
        )
}


standard_summary <- do.call(
    rbind,
    standard_summary
)

standard_summary$normalization <-
    "standard_DESeq2"


ercc_summary <- do.call(
    rbind,
    ercc_summary
)

ercc_summary$normalization <-
    "ERCC_sum"


contrast_summary <- rbind(
    standard_summary,
    ercc_summary
)


contrast_summary <- contrast_summary[
    ,
    c(
        "normalization",
        "contrast",
        "numerator",
        "denominator",
        "genes_tested",
        "median_log2FC",
        "fraction_log2FC_negative",
        "fraction_log2FC_positive",
        "significant_up",
        "significant_down"
    )
]


write.table(
    contrast_summary,
    file = file.path(
        OUT_DIR,
        "contrast_summary.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# Step 12
# Simple text report
# ============================================================

report_file <- file.path(
    OUT_DIR,
    "05_summary.txt"
)


sink(report_file)


cat("ASTRA-Drug Step 05 summary\n")
cat("==========================\n\n")


cat("Counting method comparison\n")
cat("--------------------------\n")

cat(
    "Mean Pearson correlation (log1p counts): ",
    mean(
        count_method_comparison$pearson_log1p
    ),
    "\n",
    sep = ""
)

cat(
    "Mean Spearman correlation: ",
    mean(
        count_method_comparison$spearman
    ),
    "\n\n",
    sep = ""
)


cat("ERCC size factor\n")
cat("----------------\n")

cat(
    "ERCC geometric mean: ",
    ercc_geometric_mean,
    "\n\n",
    sep = ""
)


cat("Group global-output summary\n")
cat("---------------------------\n")

print(
    group_summary,
    row.names = FALSE
)


cat("\n\nContrast summary\n")
cat("----------------\n")

print(
    contrast_summary,
    row.names = FALSE
)


sink()


# ============================================================
# Finished
# ============================================================

cat("\n")
cat("============================================================\n")
cat("05 successfully completed\n")
cat("============================================================\n")

cat("\nOutput directory:\n")
cat(OUT_DIR, "\n")

cat("\nNormalized matrices:\n")
cat(NORM_DIR, "\n")

cat("\nKey files:\n")

cat(
    file.path(
        OUT_DIR,
        "counting_method_sample_comparison.tsv"
    ),
    "\n"
)

cat(
    file.path(
        OUT_DIR,
        "ERCC_total_comparison.tsv"
    ),
    "\n"
)

cat(
    file.path(
        OUT_DIR,
        "sample_scaling_factors.tsv"
    ),
    "\n"
)

cat(
    file.path(
        OUT_DIR,
        "group_global_output_summary.tsv"
    ),
    "\n"
)

cat(
    file.path(
        OUT_DIR,
        "contrast_summary.tsv"
    ),
    "\n"
)

cat(
    file.path(
        OUT_DIR,
        "05_summary.txt"
    ),
    "\n"
)
