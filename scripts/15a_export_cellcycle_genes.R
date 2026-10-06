#!/usr/bin/env Rscript

# ============================================================
# ASTRA-Drug
# Step 15A: Export predefined Hallmark cell-cycle gene sets
# ============================================================
#
#
#   HALLMARK_E2F_TARGETS
#   HALLMARK_G2M_CHECKPOINT
#   HALLMARK_MYC_TARGETS_V1
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


HALLMARK_RDS <- file.path(
    PROJECT_DIR,
    "data/interim/msigdb/Hallmark_Homo_sapiens.rds"
)


OUT_DIR <- file.path(
    PROJECT_DIR,
    "results/15_cellcycle_confounder"
)


OUT_FILE <- file.path(
    OUT_DIR,
    "15_hallmark_cellcycle_genes.tsv"
)


TARGET_SETS <- c(
    "HALLMARK_E2F_TARGETS",
    "HALLMARK_G2M_CHECKPOINT",
    "HALLMARK_MYC_TARGETS_V1"
)


# ============================================================
# Main
# ============================================================

if (!file.exists(HALLMARK_RDS)) {
    stop(
        "Hallmark cache was not found:\n",
        HALLMARK_RDS
    )
}


dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


hallmark <- readRDS(
    HALLMARK_RDS
)


required <- c(
    "gs_name",
    "gene_symbol"
)


if (!all(required %in% colnames(hallmark))) {
    stop(
        "Hallmark RDS does not contain gs_name/gene_symbol."
    )
}


result <- hallmark[
    hallmark$gs_name %in% TARGET_SETS,
    c(
        "gs_name",
        "gene_symbol"
    )
]


result <- unique(
    result
)


missing_sets <- setdiff(
    TARGET_SETS,
    unique(result$gs_name)
)


if (length(missing_sets) > 0) {

    stop(
        "Missing Hallmark sets: ",
        paste(
            missing_sets,
            collapse = ", "
        )
    )
}


write.table(
    result,
    OUT_FILE,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


cat(
    "Step 15A completed\n\n"
)


print(
    table(
        result$gs_name
    )
)


cat(
    "\nOutput:\n",
    OUT_FILE,
    "\n"
)
