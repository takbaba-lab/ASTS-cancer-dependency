#!/usr/bin/env Rscript

# ============================================================
# ASTRA-Drug
# 06b: Freeze ASTS v1.0 signature
# ============================================================
#
#
#
#
# Ranking:
#
# Primary ASTS v1.0:
#   DHT-relative-up   Top 50
#   DHT-relative-down Top 50
#
# Sensitivity signatures:
#   Top 25 + Top 25   = 50 genes
#   Top 100 + Top 100 = 200 genes
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

INPUT_CORE <- file.path(
    PROJECT_DIR,
    "results/06_ASTS_signature/ASTS_core.tsv"
)

OUT_DIR <- file.path(
    PROJECT_DIR,
    "results/06b_ASTS_v1"
)


# Primary signature
PRIMARY_N_PER_DIRECTION <- 50


# Sensitivity analysis
SENSITIVITY_N_VALUES <- c(
    25,
    100
)


# ASTS version
ASTS_VERSION <- "ASTS_v1.0"


# ============================================================
# Output directory
# ============================================================

dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# Helper functions
# ============================================================

write_tsv <- function(x, file) {

    write.table(
        x,
        file = file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}


rank_genes <- function(x) {

    x$rank_score_abs_wald <-
        abs(x$pooled_stat)


    x <- x[
        order(
            -x$rank_score_abs_wald,
            -abs(x$pooled_log2FoldChange),
            x$pooled_padj,
            x$gene_id
        ),
        ,
        drop = FALSE
    ]


    x$direction_rank <-
        seq_len(
            nrow(x)
        )


    return(x)
}


make_signature <- function(
    ranked_up,
    ranked_down,
    n_per_direction,
    signature_name
) {

    if (nrow(ranked_up) < n_per_direction) {

        stop(
            "Not enough DHT-up genes. Required: ",
            n_per_direction,
            ", available: ",
            nrow(ranked_up)
        )
    }


    if (nrow(ranked_down) < n_per_direction) {

        stop(
            "Not enough DHT-down genes. Required: ",
            n_per_direction,
            ", available: ",
            nrow(ranked_down)
        )
    }


    up <- ranked_up[
        seq_len(n_per_direction),
        ,
        drop = FALSE
    ]


    down <- ranked_down[
        seq_len(n_per_direction),
        ,
        drop = FALSE
    ]


    up$score_weight <- 1
    down$score_weight <- -1


    up$score_component <- "DHT_relative_up"
    down$score_component <- "DHT_relative_down"


    sig <- rbind(
        up,
        down
    )


    sig$signature_name <-
        signature_name


    sig$n_per_direction <-
        n_per_direction


    first_columns <- c(
        "gene_id",
        "signature_name",
        "score_component",
        "score_weight",
        "direction_rank",
        "rank_score_abs_wald",
        "pooled_log2FoldChange",
        "pooled_stat",
        "pooled_padj",
        "M_DHT_log2FoldChange",
        "F_DHT_log2FoldChange",
        "interaction_log2FoldChange",
        "interaction_padj",
        "E2_class"
    )


    remaining_columns <- setdiff(
        colnames(sig),
        first_columns
    )


    sig <- sig[
        ,
        c(
            first_columns,
            remaining_columns
        ),
        drop = FALSE
    ]


    return(sig)
}


# ============================================================
# Step 1
# Read ASTS core
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 1: Reading ASTS core\n")
cat("============================================================\n")


core <- read.delim(
    INPUT_CORE,
    stringsAsFactors = FALSE,
    check.names = FALSE
)


required_columns <- c(
    "gene_id",
    "pooled_log2FoldChange",
    "pooled_stat",
    "pooled_padj",
    "M_DHT_log2FoldChange",
    "F_DHT_log2FoldChange",
    "interaction_log2FoldChange",
    "interaction_padj",
    "ASTS_core",
    "ASTS_direction",
    "E2_class"
)


missing_columns <- setdiff(
    required_columns,
    colnames(core)
)


if (length(missing_columns) > 0) {

    stop(
        "Missing required columns: ",
        paste(
            missing_columns,
            collapse = ", "
        )
    )
}


if (anyDuplicated(core$gene_id) > 0) {
    stop("Duplicated gene_id found in ASTS_core.tsv.")
}


if (!all(core$ASTS_core)) {
    stop("ASTS_core.tsv contains non-core genes.")
}


cat(
    "ASTS core genes:",
    nrow(core),
    "\n"
)


# ============================================================
# Step 2
# Split by DHT-relative direction
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 2: Splitting directions\n")
cat("============================================================\n")


core_up <- core[
    core$ASTS_direction == "DHT_up",
    ,
    drop = FALSE
]


core_down <- core[
    core$ASTS_direction == "DHT_down",
    ,
    drop = FALSE
]


cat(
    "DHT-relative-up:",
    nrow(core_up),
    "\n"
)

cat(
    "DHT-relative-down:",
    nrow(core_down),
    "\n"
)


# ============================================================
# Step 3
# Wald-statistic ranking
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 3: Wald-statistic ranking\n")
cat("============================================================\n")


ranked_up <- rank_genes(
    core_up
)


ranked_down <- rank_genes(
    core_down
)


write_tsv(
    ranked_up,
    file.path(
        OUT_DIR,
        "ASTS_core_DHT_relative_up_Wald_ranked.tsv"
    )
)


write_tsv(
    ranked_down,
    file.path(
        OUT_DIR,
        "ASTS_core_DHT_relative_down_Wald_ranked.tsv"
    )
)


# ============================================================
# Step 4
# Primary ASTS v1.0
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 4: Creating primary ASTS v1.0\n")
cat("============================================================\n")


primary_name <- paste0(
    ASTS_VERSION,
    "_primary_",
    PRIMARY_N_PER_DIRECTION,
    "up_",
    PRIMARY_N_PER_DIRECTION,
    "down"
)


primary <- make_signature(
    ranked_up = ranked_up,
    ranked_down = ranked_down,
    n_per_direction = PRIMARY_N_PER_DIRECTION,
    signature_name = primary_name
)


primary_file <- file.path(
    OUT_DIR,
    "ASTS_v1.0_primary_100genes.tsv"
)


write_tsv(
    primary,
    primary_file
)


primary_simple <- primary[
    ,
    c(
        "gene_id",
        "score_component",
        "score_weight",
        "direction_rank"
    )
]


write_tsv(
    primary_simple,
    file.path(
        OUT_DIR,
        "ASTS_v1.0_primary_100genes_simple.tsv"
    )
)


# ============================================================
# Step 5
# Sensitivity signatures
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 5: Creating sensitivity signatures\n")
cat("============================================================\n")


for (n in SENSITIVITY_N_VALUES) {

    signature_name <- paste0(
        ASTS_VERSION,
        "_sensitivity_",
        n,
        "up_",
        n,
        "down"
    )


    sig <- make_signature(
        ranked_up = ranked_up,
        ranked_down = ranked_down,
        n_per_direction = n,
        signature_name = signature_name
    )


    total_genes <- 2 * n


    outfile <- file.path(
        OUT_DIR,
        paste0(
            "ASTS_v1.0_sensitivity_",
            total_genes,
            "genes.tsv"
        )
    )


    write_tsv(
        sig,
        outfile
    )


    simple <- sig[
        ,
        c(
            "gene_id",
            "score_component",
            "score_weight",
            "direction_rank"
        )
    ]


    write_tsv(
        simple,
        file.path(
            OUT_DIR,
            paste0(
                "ASTS_v1.0_sensitivity_",
                total_genes,
                "genes_simple.tsv"
            )
        )
    )
}


# ============================================================
# Step 6
# E2 composition of primary signature
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 6: E2 annotation summary\n")
cat("============================================================\n")


e2_summary <- as.data.frame(
    table(
        primary$score_component,
        primary$E2_class
    ),
    stringsAsFactors = FALSE
)


colnames(e2_summary) <- c(
    "score_component",
    "E2_class",
    "n_genes"
)


write_tsv(
    e2_summary,
    file.path(
        OUT_DIR,
        "ASTS_v1.0_primary_E2_composition.tsv"
    )
)


# ============================================================
# Step 7
# Freeze manifest
# ============================================================

cat("\n")
cat("============================================================\n")
cat("Step 7: Writing freeze manifest\n")
cat("============================================================\n")


manifest_file <- file.path(
    OUT_DIR,
    "ASTS_v1.0_FREEZE_MANIFEST.txt"
)


sink(
    manifest_file
)


cat("ASTRA-Drug ASTS v1.0 freeze manifest\n")
cat("====================================\n\n")


cat("Version\n")
cat("-------\n")
cat(ASTS_VERSION, "\n\n")


cat("Source\n")
cat("------\n")
cat("Input: Step 06 ASTS_core.tsv\n")
cat("Normalization: standard DESeq2\n")
cat("ERCC-normalized DE results were NOT used for gene selection.\n\n")


cat("ASTS core definition inherited from Step 06\n")
cat("-------------------------------------------\n")
cat("Pooled DHT padj < 0.05\n")
cat("|pooled DHT log2FC| >= 0.25\n")
cat("|male DHT log2FC| >= 0.10\n")
cat("|female DHT log2FC| >= 0.10\n")
cat("Male and female DHT effects have the same direction.\n")
cat("Strong sex x DHT interaction excluded:\n")
cat("  interaction padj < 0.05 AND |interaction log2FC| >= 0.50\n\n")


cat("Ranking rule\n")
cat("------------\n")
cat("Genes are ranked separately within DHT-relative-up and DHT-relative-down.\n")
cat("Primary ranking: descending |pooled Wald statistic|.\n")
cat("Tie breakers:\n")
cat("  1. descending |pooled log2FC|\n")
cat("  2. ascending pooled padj\n")
cat("  3. gene_id\n\n")


cat("Primary signature\n")
cat("-----------------\n")
cat(
    PRIMARY_N_PER_DIRECTION,
    " DHT-relative-up genes\n",
    sep = ""
)

cat(
    PRIMARY_N_PER_DIRECTION,
    " DHT-relative-down genes\n",
    sep = ""
)

cat(
    "Total: ",
    2 * PRIMARY_N_PER_DIRECTION,
    " genes\n\n",
    sep = ""
)


cat("Sensitivity signatures\n")
cat("----------------------\n")

for (n in SENSITIVITY_N_VALUES) {

    cat(
        n,
        " up + ",
        n,
        " down = ",
        2 * n,
        " genes\n",
        sep = ""
    )
}


cat("\nInterpretation\n")
cat("--------------\n")
cat(
    "DHT-relative-up/down describe relative expression changes after ",
    "conventional compositional normalization.\n",
    sep = ""
)

cat(
    "They must not be interpreted as absolute RNA-per-cell increases ",
    "or decreases.\n\n",
    sep = ""
)


cat("Independence rule\n")
cat("-----------------\n")
cat(
    "ASTS v1.0 is frozen before examining DepMap, PRISM, GDSC, ",
    "or other drug-response associations.\n"
)

cat(
    "Drug-response results must not be used to change the primary ",
    "signature definition.\n"
)


sink()


# ============================================================
# Step 8
# Summary
# ============================================================

summary_file <- file.path(
    OUT_DIR,
    "06b_summary.txt"
)


sink(
    summary_file
)


cat("ASTRA-Drug Step 06b summary\n")
cat("===========================\n\n")


cat(
    "Input ASTS core genes:",
    nrow(core),
    "\n"
)

cat(
    "DHT-relative-up:",
    nrow(core_up),
    "\n"
)

cat(
    "DHT-relative-down:",
    nrow(core_down),
    "\n\n"
)


cat("Primary ASTS v1.0\n")
cat("-----------------\n")

cat(
    "Up genes:",
    sum(
        primary$score_component ==
            "DHT_relative_up"
    ),
    "\n"
)

cat(
    "Down genes:",
    sum(
        primary$score_component ==
            "DHT_relative_down"
    ),
    "\n"
)

cat(
    "Total:",
    nrow(primary),
    "\n\n"
)


cat("Top 10 DHT-relative-up genes\n")
cat("----------------------------\n")

print(
    ranked_up[
        seq_len(
            min(
                10,
                nrow(ranked_up)
            )
        ),
        c(
            "gene_id",
            "pooled_log2FoldChange",
            "pooled_stat",
            "pooled_padj",
            "M_DHT_log2FoldChange",
            "F_DHT_log2FoldChange",
            "E2_class"
        )
    ],
    row.names = FALSE
)


cat("\n\nTop 10 DHT-relative-down genes\n")
cat("------------------------------\n")

print(
    ranked_down[
        seq_len(
            min(
                10,
                nrow(ranked_down)
            )
        ),
        c(
            "gene_id",
            "pooled_log2FoldChange",
            "pooled_stat",
            "pooled_padj",
            "M_DHT_log2FoldChange",
            "F_DHT_log2FoldChange",
            "E2_class"
        )
    ],
    row.names = FALSE
)


cat("\n\nE2 composition of primary ASTS v1.0\n")
cat("-----------------------------------\n")

print(
    e2_summary,
    row.names = FALSE
)


cat("\n\nPrimary file\n")
cat("------------\n")
cat(primary_file, "\n")


sink()


# ============================================================
# Finished
# ============================================================

cat("\n")
cat("============================================================\n")
cat("06b completed successfully\n")
cat("============================================================\n")

cat(
    "\nPrimary ASTS v1.0:",
    primary_file,
    "\n"
)

cat(
    "\nFreeze manifest:",
    manifest_file,
    "\n"
)

cat(
    "\nSummary:",
    summary_file,
    "\n"
)
