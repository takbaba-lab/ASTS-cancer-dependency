#!/usr/bin/env Rscript

# ============================================================
# ASTRA-Drug
# Step 18C: Global-burden-adjusted CRISPR dependency GSEA
# ============================================================
#
#
#
#   HALLMARK_OXIDATIVE_PHOSPHORYLATION
#   REACTOME_MITOCHONDRIAL_TRANSLATION
#   REACTOME_RESPIRATORY_ELECTRON_TRANSPORT
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


INPUT_FILE <- file.path(
    PROJECT_DIR,
    "results/18_mito_dependency_specificity",
    "18_CRISPR_genomewide_burden_adjusted.tsv"
)


HALLMARK_RDS <- file.path(
    PROJECT_DIR,
    "data/interim/msigdb",
    "Hallmark_Homo_sapiens.rds"
)


REACTOME_RDS <- file.path(
    PROJECT_DIR,
    "data/interim/msigdb",
    "Reactome_Homo_sapiens.rds"
)


OUT_DIR <- file.path(
    PROJECT_DIR,
    "results/18_mito_dependency_specificity/GSEA"
)


MIN_SIZE <- 15

MAX_SIZE <- 500

N_PERM_SIMPLE <- 10000

FDR_THRESHOLD <- 0.05

TOP_N_PLOT <- 25


FOCUS_PATHWAYS <- c(
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "REACTOME_MITOCHONDRIAL_TRANSLATION",
    "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
    "REACTOME_AEROBIC_RESPIRATION_AND_RESPIRATORY_ELECTRON_TRANSPORT",
    "REACTOME_COMPLEX_I_BIOGENESIS",
    "REACTOME_MITOCHONDRIAL_PROTEIN_IMPORT",
    "REACTOME_FORMATION_OF_ATP_BY_CHEMIOSMOTIC_COUPLING"
)


# ============================================================
# Packages
# ============================================================

suppressPackageStartupMessages({

    library(data.table)
    library(fgsea)
    library(ggplot2)

})


dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


# ============================================================
# Input
# ============================================================

for (
    path
    in c(
        INPUT_FILE,
        HALLMARK_RDS,
        REACTOME_RDS
    )
) {

    if (!file.exists(path)) {

        stop(
            "Required file not found:\n",
            path
        )
    }
}


x <- fread(
    INPUT_FILE
)


required <- c(
    "gene_symbol",
    "t_stat_burden_adjusted",
    "is_burden_reference_gene"
)


if (
    !all(
        required
        %in%
        colnames(x)
    )
) {

    stop(
        "Required columns missing."
    )
}


# ------------------------------------------------------------
# Remove burden reference genes
# ------------------------------------------------------------

x <- x[
    is_burden_reference_gene == FALSE
]


x <- x[
    is.finite(
        t_stat_burden_adjusted
    )
]


x <- x[
    !is.na(gene_symbol)
    &
    gene_symbol != ""
]


x[
    ,
    abs_t :=
        abs(
            t_stat_burden_adjusted
        )
]


setorder(
    x,
    -abs_t
)


x <- x[
    !duplicated(
        gene_symbol
    )
]


ranks <- x$t_stat_burden_adjusted

names(
    ranks
) <- x$gene_symbol


ranks <- sort(
    ranks,
    decreasing = TRUE
)


cat(
    "Genes in adjusted ranking: ",
    length(ranks),
    "\n",
    sep = ""
)


# ============================================================
# Gene sets
# ============================================================

hallmark_df <- readRDS(
    HALLMARK_RDS
)


reactome_df <- readRDS(
    REACTOME_RDS
)


hallmark_pathways <- split(
    hallmark_df$gene_symbol,
    hallmark_df$gs_name
)


reactome_pathways <- split(
    reactome_df$gene_symbol,
    reactome_df$gs_name
)


hallmark_pathways <- lapply(
    hallmark_pathways,
    unique
)


reactome_pathways <- lapply(
    reactome_pathways,
    unique
)


# ============================================================
# GSEA helper
# ============================================================

run_gsea <- function(
    pathways,
    label
) {

    cat(
        "\nRunning ",
        label,
        "...\n",
        sep = ""
    )


    result <- fgseaMultilevel(

        pathways =
            pathways,

        stats =
            ranks,

        minSize =
            MIN_SIZE,

        maxSize =
            MAX_SIZE,

        eps =
            0,

        nPermSimple =
            N_PERM_SIMPLE
    )


    result <- as.data.table(
        result
    )


    result <- result[
        is.finite(NES)
        &
        is.finite(pval)
        &
        is.finite(padj)
    ]


    result[
        ,
        abs_NES :=
            abs(NES)
    ]


    setorder(
        result,
        padj,
        -abs_NES
    )


    result[
        ,
        abs_NES :=
            NULL
    ]


    result[
        ,
        direction :=
            fifelse(
                NES > 0,
                "ASTS_high_less_dependent",
                "ASTS_high_more_dependent"
            )
    ]


    return(
        result
    )
}


hallmark <- run_gsea(
    hallmark_pathways,
    "Hallmark"
)


reactome <- run_gsea(
    reactome_pathways,
    "Reactome"
)


# ============================================================
# Flatten and save
# ============================================================

flatten <- function(result) {

    out <- copy(
        result
    )


    if (
        "leadingEdge"
        %in%
        colnames(out)
    ) {

        out[
            ,
            leadingEdge :=
                vapply(
                    leadingEdge,
                    function(z) {

                        paste(
                            z,
                            collapse = ";"
                        )
                    },
                    character(1)
                )
        ]
    }


    return(out)
}


fwrite(
    flatten(hallmark),
    file.path(
        OUT_DIR,
        "18c_Hallmark_burden_adjusted.tsv"
    ),
    sep = "\t"
)


fwrite(
    flatten(reactome),
    file.path(
        OUT_DIR,
        "18c_Reactome_burden_adjusted.tsv"
    ),
    sep = "\t"
)


# ============================================================
# Focus table
# ============================================================

combined_focus <- rbind(

    hallmark[
        pathway %in% FOCUS_PATHWAYS,
        .(
            collection = "Hallmark",
            pathway,
            NES,
            pval,
            padj,
            direction
        )
    ],

    reactome[
        pathway %in% FOCUS_PATHWAYS,
        .(
            collection = "Reactome",
            pathway,
            NES,
            pval,
            padj,
            direction
        )
    ]
)


fwrite(
    combined_focus,
    file.path(
        OUT_DIR,
        "18c_mito_focus_pathways.tsv"
    ),
    sep = "\t"
)


# ============================================================
# Plot helper
# ============================================================

make_plot <- function(
    result,
    title_text,
    output_file
) {

    plot_data <- copy(
        result
    )


    plot_data[
        ,
        abs_NES :=
            abs(NES)
    ]


    setorder(
        plot_data,
        padj,
        -abs_NES
    )


    plot_data <- plot_data[
        seq_len(
            min(
                TOP_N_PLOT,
                .N
            )
        )
    ]


    plot_data[
        ,
        pathway :=
            factor(
                pathway,
                levels =
                    rev(pathway)
            )
    ]


    p <- ggplot(
        plot_data,
        aes(
            x = NES,
            y = pathway
        )
    ) +

        geom_col() +

        geom_vline(
            xintercept = 0,
            linetype = "dashed"
        ) +

        labs(

            title =
                title_text,

            subtitle =
                "Global CRISPR dependency burden adjusted",

            x =
                paste0(
                    "Normalized enrichment score\n",
                    "positive = ASTS-high less dependent; ",
                    "negative = ASTS-high more dependent"
                ),

            y = NULL
        ) +

        theme_bw(
            base_size = 10
        )


    ggsave(
        output_file,
        p,
        width = 11,
        height = 8,
        dpi = 200
    )
}


make_plot(

    hallmark,

    "ASTS-associated CRISPR dependency: Hallmark",

    file.path(
        OUT_DIR,
        "18c_Hallmark_burden_adjusted.png"
    )
)


make_plot(

    reactome,

    "ASTS-associated CRISPR dependency: Reactome",

    file.path(
        OUT_DIR,
        "18c_Reactome_burden_adjusted.png"
    )
)


# ============================================================
# Summary
# ============================================================

summary_file <- file.path(
    OUT_DIR,
    "18c_summary.txt"
)


sink(
    summary_file
)


cat(
    "ASTRA-Drug Step 18C summary\n"
)

cat(
    "============================\n\n"
)


cat(
    "Ranking:\n"
)

cat(
    "Genome-wide CRISPR ASTS t-statistic after adjustment for\n"
)

cat(
    "global common-essential burden and global centered shift.\n\n"
)


cat(
    "Hallmark pathways tested: ",
    nrow(hallmark),
    "\n",
    sep = ""
)


cat(
    "Hallmark FDR < 0.05: ",
    sum(
        hallmark$padj < FDR_THRESHOLD,
        na.rm = TRUE
    ),
    "\n\n",
    sep = ""
)


cat(
    "Reactome pathways tested: ",
    nrow(reactome),
    "\n",
    sep = ""
)


cat(
    "Reactome FDR < 0.05: ",
    sum(
        reactome$padj < FDR_THRESHOLD,
        na.rm = TRUE
    ),
    "\n\n",
    sep = ""
)


cat(
    "Mitochondrial focus pathways\n"
)

cat(
    "----------------------------\n"
)


print(
    combined_focus
)


sink()


cat(
    "\nStep 18C completed successfully\n"
)

cat(
    "Summary:\n",
    summary_file,
    "\n"
)
