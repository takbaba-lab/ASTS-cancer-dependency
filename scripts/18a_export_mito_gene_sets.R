#!/usr/bin/env Rscript

# ============================================================
# ASTRA-Drug
# Step 18A: Export mitochondrial gene sets
# ============================================================
#
#
#
# Primary:
#   HALLMARK_OXIDATIVE_PHOSPHORYLATION
#   REACTOME_MITOCHONDRIAL_TRANSLATION
#   REACTOME_RESPIRATORY_ELECTRON_TRANSPORT
#
#   MITO_CORE_UNION
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


REACTOME_RDS <- file.path(
    PROJECT_DIR,
    "data/interim/msigdb/Reactome_Homo_sapiens.rds"
)


OUT_DIR <- file.path(
    PROJECT_DIR,
    "results/18_mito_dependency_specificity"
)


OUT_FILE <- file.path(
    OUT_DIR,
    "18_mitochondrial_gene_sets.tsv"
)


PRIMARY_HALLMARK <- c(
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION"
)


PRIMARY_REACTOME <- c(
    "REACTOME_MITOCHONDRIAL_TRANSLATION",
    "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT"
)


SECONDARY_REACTOME <- c(
    "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",
    "REACTOME_COMPLEX_I_BIOGENESIS",
    "REACTOME_COMPLEX_IV_ASSEMBLY",
    "REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT",
    "REACTOME_MITOCHONDRIAL_PROTEIN_DEGRADATION",
    "REACTOME_MITOCHONDRIAL_TRNA_AMINOACYLATION",
    "REACTOME_RRNA_PROCESSING_IN_THE_MITOCHONDRION",
    "REACTOME_CRISTAE_FORMATION",
    "REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING"
)


# ============================================================
# Main
# ============================================================

dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


for (path in c(HALLMARK_RDS, REACTOME_RDS)) {

    if (!file.exists(path)) {
        stop(
            "Required file not found:\n",
            path
        )
    }
}


hallmark <- readRDS(
    HALLMARK_RDS
)


reactome <- readRDS(
    REACTOME_RDS
)


required <- c(
    "gs_name",
    "gene_symbol"
)


if (
    !all(required %in% colnames(hallmark))
    ||
    !all(required %in% colnames(reactome))
) {

    stop(
        "MSigDB RDS lacks gs_name or gene_symbol."
    )
}


# ------------------------------------------------------------
# Select predefined sets
# ------------------------------------------------------------

hallmark_out <- hallmark[
    hallmark$gs_name %in% PRIMARY_HALLMARK,
    c("gs_name", "gene_symbol")
]


hallmark_out$source <- "Hallmark"

hallmark_out$analysis_level <- "primary"


reactome_primary <- reactome[
    reactome$gs_name %in% PRIMARY_REACTOME,
    c("gs_name", "gene_symbol")
]


reactome_primary$source <- "Reactome"

reactome_primary$analysis_level <- "primary"


reactome_secondary <- reactome[
    reactome$gs_name %in% SECONDARY_REACTOME,
    c("gs_name", "gene_symbol")
]


reactome_secondary$source <- "Reactome"

reactome_secondary$analysis_level <- "secondary"


result <- rbind(
    hallmark_out,
    reactome_primary,
    reactome_secondary
)


result <- unique(
    result
)


# ------------------------------------------------------------
# Check expected primary sets
# ------------------------------------------------------------

expected_primary <- c(
    PRIMARY_HALLMARK,
    PRIMARY_REACTOME
)


missing_primary <- setdiff(
    expected_primary,
    unique(result$gs_name)
)


if (length(missing_primary) > 0) {

    stop(
        "Missing primary gene sets: ",
        paste(
            missing_primary,
            collapse = ", "
        )
    )
}


# ------------------------------------------------------------
# Mito core union
# ------------------------------------------------------------

core_genes <- unique(
    result$gene_symbol[
        result$gs_name %in% expected_primary
    ]
)


union_rows <- data.frame(
    gs_name = "MITO_CORE_UNION",
    gene_symbol = core_genes,
    source = "Derived_union",
    analysis_level = "primary",
    stringsAsFactors = FALSE
)


result <- rbind(
    result,
    union_rows
)


write.table(
    result,
    OUT_FILE,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


cat(
    "Step 18A completed successfully\n\n"
)


cat(
    "Gene-set sizes:\n"
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
