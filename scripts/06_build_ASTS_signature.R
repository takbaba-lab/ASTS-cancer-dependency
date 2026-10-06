#!/usr/bin/env Rscript

# ============================================================
# ASTRA-Drug
# 06: Build ASTS candidate signature
# ============================================================
#
#
# ASTS (Androgen-Suppressed Transcriptional State)
#
#
#   1. M_ORX / M_ORX_DHT
#      F_OVX / F_OVX_DHT
#
#   2. additive model
#         ~ sex + dht
#
#   3. interaction model
#         ~ sex + dht + sex:dht
#
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


RAW_COUNT_FILE <- file.path(
    PROJECT_DIR,
    "data/processed/normalized/counts_strand_aware_endogenous_raw.tsv"
)


SCALING_FILE <- file.path(
    PROJECT_DIR,
    "results/05_normalization/sample_scaling_factors.tsv"
)


DE05_DIR <- file.path(
    PROJECT_DIR,
    "results/05_normalization/DESeq2_standard"
)


OUT_DIR <- file.path(
    PROJECT_DIR,
    "results/06_ASTS_signature"
)


FIG_DIR <- file.path(
    OUT_DIR,
    "figures"
)


# ------------------------------------------------------------
# ASTS selection parameters
# ------------------------------------------------------------

POOLED_PADJ_MAX <- 0.05

MIN_ABS_POOLED_LFC <- 0.25

MIN_ABS_INDIVIDUAL_LFC <- 0.10


INTERACTION_PADJ_MAX <- 0.05

INTERACTION_MIN_ABS_LFC <- 0.50


MIN_TOTAL_COUNT <- 10


TOP_N_VALUES <- c(
    25,
    50,
    100
)


# ============================================================
# Packages
# ============================================================

suppressPackageStartupMessages({
    library(DESeq2)
})


# ============================================================
# Output directories
# ============================================================

dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)

dir.create(
    FIG_DIR,
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
        stop("Duplicated gene_id found.")
    }

    rownames(x) <- x$gene_id
    x$gene_id <- NULL

    x <- as.matrix(x)

    storage.mode(x) <- "integer"

    return(x)
}


read_DE_result <- function(
    file,
    prefix
) {

    x <- read.delim(
        file,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    keep <- c(
        "gene_id",
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj"
    )

    x <- x[, keep]

    colnames(x)[
        colnames(x) != "gene_id"
    ] <- paste0(
        prefix,
        "_",
        colnames(x)[
            colnames(x) != "gene_id"
        ]
    )

    return(x)
}


write_tsv <- function(
    x,
    file
) {

    write.table(
        x,
        file = file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}


safe_cor <- function(
    x,
    y
) {

    ok <- is.finite(x) &
        is.finite(y)

    cor(
        x[ok],
        y[ok],
        method = "pearson"
    )
}


# ============================================================
# Step 1
# Read metadata and counts
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 1: Reading data\n")
cat("============================================================\n")


metadata <- read.delim(
    METADATA_FILE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


counts_all <- read_count_matrix(
    RAW_COUNT_FILE
)


scaling <- read.delim(
    SCALING_FILE,
    stringsAsFactors = FALSE
)


# ============================================================
# Step 2
# Select the four core DHT groups
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 2: Selecting DHT discovery samples\n")
cat("============================================================\n")


CORE_GROUPS <- c(
    "M_ORX",
    "M_ORX_DHT",
    "F_OVX",
    "F_OVX_DHT"
)


meta <- metadata[
    metadata$group %in% CORE_GROUPS,
    ,
    drop = FALSE
]


if (nrow(meta) != 12) {

    stop(
        "Expected 12 DHT discovery samples, found ",
        nrow(meta)
    )
}


if (!all(meta$sample_id %in% colnames(counts_all))) {
    stop("Some samples are missing from count matrix.")
}


counts <- counts_all[
    ,
    meta$sample_id,
    drop = FALSE
]


# ------------------------------------------------------------
# ------------------------------------------------------------

meta$sex_model <- ifelse(
    grepl("^M_", meta$group),
    "male",
    "female"
)


meta$dht <- ifelse(
    grepl("_DHT$", meta$group),
    "present",
    "absent"
)


meta$sex_model <- factor(
    meta$sex_model,
    levels = c(
        "male",
        "female"
    )
)


meta$dht <- factor(
    meta$dht,
    levels = c(
        "absent",
        "present"
    )
)


rownames(meta) <- meta$sample_id


cat("\nSample table:\n")

print(
    meta[
        ,
        c(
            "sample_id",
            "group",
            "sex_model",
            "dht"
        )
    ],
    row.names = FALSE
)


# ============================================================
# Step 3
# Use the DESeq2 size factors already determined in Step 05
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 3: Loading Step-05 DESeq2 size factors\n")
cat("============================================================\n")


if (!all(meta$sample_id %in% scaling$sample_id)) {
    stop("Missing sample in sample_scaling_factors.tsv")
}


standard_sf <- scaling$DESeq2_size_factor[
    match(
        meta$sample_id,
        scaling$sample_id
    )
]


names(standard_sf) <- meta$sample_id


if (any(!is.finite(standard_sf))) {
    stop("Invalid DESeq2 size factor.")
}


print(
    data.frame(
        sample_id = names(standard_sf),
        size_factor = standard_sf
    ),
    row.names = FALSE
)


# ============================================================
# Step 4
# Gene filter
# ============================================================

keep_gene <- rowSums(
    counts
) >= MIN_TOTAL_COUNT


counts_use <- counts[
    keep_gene,
    ,
    drop = FALSE
]


cat(
    "\nGenes retained:",
    nrow(counts_use),
    "/",
    nrow(counts),
    "\n"
)


# ============================================================
# Step 5
# Additive model: common DHT effect
#
# design = ~ sex + DHT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 5: Pooled DHT additive model\n")
cat("============================================================\n")


dds_add <- DESeqDataSetFromMatrix(
    countData = counts_use,
    colData = meta,
    design = ~ sex_model + dht
)


sizeFactors(dds_add) <-
    standard_sf[
        colnames(dds_add)
    ]


dds_add <- DESeq(
    dds_add,
    quiet = TRUE
)


cat("\nResult names:\n")
print(
    resultsNames(dds_add)
)


pooled_name <- "dht_present_vs_absent"


if (!pooled_name %in% resultsNames(dds_add)) {

    stop(
        "Expected coefficient not found: ",
        pooled_name
    )
}


res_pooled <- results(
    dds_add,
    name = pooled_name
)


pooled_df <- as.data.frame(
    res_pooled
)


pooled_df$gene_id <-
    rownames(pooled_df)


pooled_df <- pooled_df[
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


colnames(pooled_df)[
    -1
] <- paste0(
    "pooled_",
    colnames(pooled_df)[
        -1
    ]
)


write_tsv(
    pooled_df,
    file.path(
        OUT_DIR,
        "pooled_DHT_additive.tsv"
    )
)


# ============================================================
# Step 6
# Interaction model
#
# design = ~ sex + DHT + sex:DHT
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 6: sex x DHT interaction model\n")
cat("============================================================\n")


dds_int <- DESeqDataSetFromMatrix(
    countData = counts_use,
    colData = meta,
    design = ~ sex_model + dht + sex_model:dht
)


sizeFactors(dds_int) <-
    standard_sf[
        colnames(dds_int)
    ]


dds_int <- DESeq(
    dds_int,
    quiet = TRUE
)


int_names <- resultsNames(
    dds_int
)


cat("\nInteraction model result names:\n")
print(int_names)


interaction_name <- grep(
    "sex_model.*dht|dht.*sex_model",
    int_names,
    value = TRUE
)


if (length(interaction_name) != 1) {

    stop(
        "Could not uniquely identify interaction coefficient."
    )
}


cat(
    "\nInteraction coefficient: ",
    interaction_name,
    "\n",
    sep = ""
)


res_interaction <- results(
    dds_int,
    name = interaction_name
)


interaction_df <- as.data.frame(
    res_interaction
)


interaction_df$gene_id <-
    rownames(interaction_df)


interaction_df <- interaction_df[
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


colnames(interaction_df)[
    -1
] <- paste0(
    "interaction_",
    colnames(interaction_df)[
        -1
    ]
)


write_tsv(
    interaction_df,
    file.path(
        OUT_DIR,
        "sex_DHT_interaction.tsv"
    )
)


# ============================================================
# Step 7
# Load Step-05 individual DHT and E2 results
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 7: Loading individual-sex and E2 effects\n")
cat("============================================================\n")


M_DHT <- read_DE_result(
    file.path(
        DE05_DIR,
        "M_DHT.tsv"
    ),
    "M_DHT"
)


F_DHT <- read_DE_result(
    file.path(
        DE05_DIR,
        "F_DHT.tsv"
    ),
    "F_DHT"
)


M_E2 <- read_DE_result(
    file.path(
        DE05_DIR,
        "M_E2.tsv"
    ),
    "M_E2"
)


F_E2 <- read_DE_result(
    file.path(
        DE05_DIR,
        "F_E2.tsv"
    ),
    "F_E2"
)


# ============================================================
# Step 8
# Merge all information
# ============================================================

annotation <- Reduce(
    function(x, y) {

        merge(
            x,
            y,
            by = "gene_id",
            all = FALSE
        )

    },
    list(
        pooled_df,
        interaction_df,
        M_DHT,
        F_DHT,
        M_E2,
        F_E2
    )
)


# ============================================================
# Step 9
# Define ASTS criteria
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 9: Defining ASTS broad and core\n")
cat("============================================================\n")


# ------------------------------------------------------------
# ------------------------------------------------------------

annotation$same_direction_DHT <-
    (
        annotation$pooled_log2FoldChange *
        annotation$M_DHT_log2FoldChange
    ) > 0 &
    (
        annotation$pooled_log2FoldChange *
        annotation$F_DHT_log2FoldChange
    ) > 0


# ------------------------------------------------------------
# strong interaction
# ------------------------------------------------------------

annotation$strong_sex_DHT_interaction <-
    !is.na(
        annotation$interaction_padj
    ) &
    annotation$interaction_padj <
        INTERACTION_PADJ_MAX &
    abs(
        annotation$interaction_log2FoldChange
    ) >=
        INTERACTION_MIN_ABS_LFC


# ------------------------------------------------------------
# ASTS broad
# ------------------------------------------------------------

annotation$ASTS_broad <-
    !is.na(
        annotation$pooled_padj
    ) &
    annotation$pooled_padj <
        POOLED_PADJ_MAX &
    annotation$same_direction_DHT


# ------------------------------------------------------------
# ASTS core
# ------------------------------------------------------------

annotation$ASTS_core <-
    annotation$ASTS_broad &
    abs(
        annotation$pooled_log2FoldChange
    ) >=
        MIN_ABS_POOLED_LFC &
    abs(
        annotation$M_DHT_log2FoldChange
    ) >=
        MIN_ABS_INDIVIDUAL_LFC &
    abs(
        annotation$F_DHT_log2FoldChange
    ) >=
        MIN_ABS_INDIVIDUAL_LFC &
    !annotation$strong_sex_DHT_interaction


# ------------------------------------------------------------
# DHT direction
# ------------------------------------------------------------

annotation$ASTS_direction <- ifelse(
    annotation$pooled_log2FoldChange > 0,
    "DHT_up",
    "DHT_down"
)


# ============================================================
# Step 10
# E2 annotation
# ============================================================

annotation$mean_E2_log2FC <-
    rowMeans(
        cbind(
            annotation$M_E2_log2FoldChange,
            annotation$F_E2_log2FoldChange
        ),
        na.rm = TRUE
    )


annotation$E2_same_direction_male <-
    (
        annotation$pooled_log2FoldChange *
        annotation$M_E2_log2FoldChange
    ) > 0


annotation$E2_same_direction_female <-
    (
        annotation$pooled_log2FoldChange *
        annotation$F_E2_log2FoldChange
    ) > 0


annotation$E2_sig_male <-
    !is.na(
        annotation$M_E2_padj
    ) &
    annotation$M_E2_padj < 0.05


annotation$E2_sig_female <-
    !is.na(
        annotation$F_E2_padj
    ) &
    annotation$F_E2_padj < 0.05


annotation$E2_class <- "mixed_or_weak"


annotation$E2_class[
    annotation$E2_same_direction_male &
    annotation$E2_same_direction_female &
    (
        annotation$E2_sig_male |
        annotation$E2_sig_female
    )
] <- "E2_concordant"


annotation$E2_class[
    !annotation$E2_same_direction_male &
    !annotation$E2_same_direction_female &
    (
        annotation$E2_sig_male |
        annotation$E2_sig_female
    )
] <- "E2_opposite"


annotation$E2_class[
    !annotation$E2_sig_male &
    !annotation$E2_sig_female
] <- "no_significant_E2"


# ============================================================
# Step 11
# Save full annotation
# ============================================================

annotation <- annotation[
    order(
        is.na(annotation$pooled_padj),
        annotation$pooled_padj,
        -abs(
            annotation$pooled_log2FoldChange
        )
    ),
]


write_tsv(
    annotation,
    file.path(
        OUT_DIR,
        "ASTS_all_gene_annotation.tsv"
    )
)


# ============================================================
# Step 12
# Broad / Core tables
# ============================================================

broad <- annotation[
    annotation$ASTS_broad,
]


core <- annotation[
    annotation$ASTS_core,
]


core_up <- core[
    core$pooled_log2FoldChange > 0,
]


core_down <- core[
    core$pooled_log2FoldChange < 0,
]


write_tsv(
    broad,
    file.path(
        OUT_DIR,
        "ASTS_broad.tsv"
    )
)


write_tsv(
    core,
    file.path(
        OUT_DIR,
        "ASTS_core.tsv"
    )
)


write_tsv(
    core_up,
    file.path(
        OUT_DIR,
        "ASTS_core_DHT_up.tsv"
    )
)


write_tsv(
    core_down,
    file.path(
        OUT_DIR,
        "ASTS_core_DHT_down.tsv"
    )
)


# ============================================================
# Step 13
# Ranked Top-N signatures
# ============================================================

rank_signature <- function(
    x,
    n,
    direction
) {

    x <- x[
        order(
            -abs(
                x$pooled_log2FoldChange
            ),
            x$pooled_padj
        ),
        ,
        drop = FALSE
    ]


    n_use <- min(
        n,
        nrow(x)
    )


    if (n_use == 0) {
        return(x[FALSE, ])
    }


    out <- x[
        seq_len(n_use),
        ,
        drop = FALSE
    ]


    out$signature_direction <-
        direction


    out$score_weight <- ifelse(
        direction == "DHT_up",
        1,
        -1
    )


    out$signature_rank <-
        seq_len(
            nrow(out)
        )


    return(out)
}


for (n in TOP_N_VALUES) {

    up_n <- rank_signature(
        core_up,
        n,
        "DHT_up"
    )


    down_n <- rank_signature(
        core_down,
        n,
        "DHT_down"
    )


    sig <- rbind(
        up_n,
        down_n
    )


    write_tsv(
        sig,
        file.path(
            OUT_DIR,
            paste0(
                "ASTS_candidate_top",
                n,
                "_each_direction.tsv"
            )
        )
    )
}


# ============================================================
# Step 14
# Diagnostic figures
# ============================================================

# ------------------------------------------------------------
# Figure 1
# Male DHT vs Female DHT
# ------------------------------------------------------------

png(
    file.path(
        FIG_DIR,
        "01_Male_vs_Female_DHT_log2FC.png"
    ),
    width = 1800,
    height = 1800,
    res = 200
)


plot(
    annotation$M_DHT_log2FoldChange,
    annotation$F_DHT_log2FoldChange,
    pch = 16,
    cex = 0.35,
    xlab = "Male DHT log2FC",
    ylab = "Female DHT log2FC",
    main = "DHT response: male vs female"
)


abline(
    h = 0,
    v = 0,
    lty = 2
)


abline(
    a = 0,
    b = 1,
    lty = 3
)


points(
    core$M_DHT_log2FoldChange,
    core$F_DHT_log2FoldChange,
    pch = 16,
    cex = 0.55
)


dev.off()


# ------------------------------------------------------------
# Figure 2
# pooled DHT vs mean E2
# ------------------------------------------------------------

png(
    file.path(
        FIG_DIR,
        "02_Pooled_DHT_vs_mean_E2_log2FC.png"
    ),
    width = 1800,
    height = 1800,
    res = 200
)


plot(
    core$pooled_log2FoldChange,
    core$mean_E2_log2FC,
    pch = 16,
    cex = 0.6,
    xlab = "Pooled DHT log2FC",
    ylab = "Mean E2 log2FC",
    main = "ASTS core: DHT vs E2"
)


abline(
    h = 0,
    v = 0,
    lty = 2
)


dev.off()


# ============================================================
# Step 15
# Summary
# ============================================================

summary_file <- file.path(
    OUT_DIR,
    "06_summary.txt"
)


sink(
    summary_file
)


cat("ASTRA-Drug Step 06 summary\n")
cat("==========================\n\n")


cat("Parameters\n")
cat("----------\n")

cat(
    "Pooled padj max:",
    POOLED_PADJ_MAX,
    "\n"
)

cat(
    "Min |pooled log2FC|:",
    MIN_ABS_POOLED_LFC,
    "\n"
)

cat(
    "Min |individual-sex log2FC|:",
    MIN_ABS_INDIVIDUAL_LFC,
    "\n"
)

cat(
    "Strong interaction: padj <",
    INTERACTION_PADJ_MAX,
    "and |LFC| >=",
    INTERACTION_MIN_ABS_LFC,
    "\n\n"
)


cat("Gene counts\n")
cat("-----------\n")

cat(
    "Genes tested:",
    nrow(annotation),
    "\n"
)

cat(
    "ASTS broad:",
    nrow(broad),
    "\n"
)

cat(
    "ASTS core:",
    nrow(core),
    "\n"
)

cat(
    "ASTS core DHT-up:",
    nrow(core_up),
    "\n"
)

cat(
    "ASTS core DHT-down:",
    nrow(core_down),
    "\n\n"
)


cat("Cross-sex DHT effect correlation\n")
cat("--------------------------------\n")

cat(
    "All genes:",
    safe_cor(
        annotation$M_DHT_log2FoldChange,
        annotation$F_DHT_log2FoldChange
    ),
    "\n"
)

cat(
    "ASTS broad:",
    safe_cor(
        broad$M_DHT_log2FoldChange,
        broad$F_DHT_log2FoldChange
    ),
    "\n"
)

cat(
    "ASTS core:",
    safe_cor(
        core$M_DHT_log2FoldChange,
        core$F_DHT_log2FoldChange
    ),
    "\n\n"
)


cat("E2 classification among ASTS core\n")
cat("---------------------------------\n")

print(
    table(
        core$E2_class,
        useNA = "ifany"
    )
)


cat("\n\nTop-N files\n")
cat("-----------\n")

for (n in TOP_N_VALUES) {

    cat(
        "Top",
        n,
        "per direction:",
        file.path(
            OUT_DIR,
            paste0(
                "ASTS_candidate_top",
                n,
                "_each_direction.tsv"
            )
        ),
        "\n"
    )
}


sink()


cat("\n")
cat("============================================================\n")
cat("06 completed successfully\n")
cat("============================================================\n")

cat(
    "\nASTS broad genes:",
    nrow(broad),
    "\n"
)

cat(
    "ASTS core genes:",
    nrow(core),
    "\n"
)

cat(
    "  DHT-up:",
    nrow(core_up),
    "\n"
)

cat(
    "  DHT-down:",
    nrow(core_down),
    "\n"
)

cat(
    "\nSummary:\n",
    summary_file,
    "\n"
)
