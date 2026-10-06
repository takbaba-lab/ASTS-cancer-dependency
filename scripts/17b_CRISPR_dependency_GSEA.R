#!/usr/bin/env Rscript

# ============================================================
# ASTRA-Drug
# Step 17B: Genome-wide CRISPR dependency GSEA
# ============================================================
#
#
#
#   Hallmark
#   Reactome
#
#
#
#   NES > 0
#
#   NES < 0
#
#
# ============================================================


# ============================================================
# 1. Path / Parameter
# ============================================================

PROJECT_DIR <- normalizePath(
    Sys.getenv("ASTRA_PROJECT_DIR", unset = getwd()),
    mustWork = FALSE
)


INPUT_FILE <- file.path(
    PROJECT_DIR,
    "results/17_CRISPR_genomewide",
    "17_CRISPR_genomewide_ASTS_association.tsv"
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
    "results/17_CRISPR_genomewide/GSEA"
)


MIN_SIZE <- 15

MAX_SIZE <- 500

N_PERM_SIMPLE <- 10000

FDR_THRESHOLD <- 0.05

TOP_N_PLOT <- 25



# ============================================================
# 2. Packages
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
# 3. Input checks
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



# ============================================================
# 4. Load genome-wide ranking
# ============================================================

cat(
    "Loading genome-wide CRISPR ranking...\n"
)


x <- fread(
    INPUT_FILE
)


required_columns <- c(
    "gene_symbol",
    "t_stat"
)


if (
    !all(
        required_columns
        %in%
        colnames(x)
    )
) {

    stop(
        "Required columns are missing."
    )
}


x[
    ,
    t_stat :=
        as.numeric(
            t_stat
        )
]


x <- x[
    is.finite(
        t_stat
    )
]


x <- x[
    !is.na(
        gene_symbol
    )
    &
    gene_symbol != ""
]


# duplicate protection
x[
    ,
    abs_t :=
        abs(
            t_stat
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


x[
    ,
    abs_t :=
        NULL
]


ranks <- x$t_stat

names(
    ranks
) <- x$gene_symbol


ranks <- sort(
    ranks,
    decreasing = TRUE
)


cat(
    "Genes in ranking: ",
    length(
        ranks
    ),
    "\n",
    sep = ""
)



# ============================================================
# 5. Load cached MSigDB
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
# 6. GSEA helper
# ============================================================

run_gsea <- function(
    pathways,
    label
) {

    cat(
        "\nRunning ",
        label,
        " GSEA...\n",
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


    cat(
        "Pathways returned: ",
        nrow(
            result
        ),
        "\n",
        sep = ""
    )


    # remove uncalculable pathways
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
            abs(
                NES
            )
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



# ============================================================
# 7. Run GSEA
# ============================================================

hallmark <- run_gsea(
    hallmark_pathways,
    "Hallmark"
)


reactome <- run_gsea(
    reactome_pathways,
    "Reactome"
)



# ============================================================
# 8. Prepare leadingEdge for TSV
# ============================================================

flatten_leading_edge <- function(result) {

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


    return(
        out
    )
}



# ============================================================
# 9. Save tables
# ============================================================

fwrite(
    flatten_leading_edge(
        hallmark
    ),
    file.path(
        OUT_DIR,
        "17b_Hallmark_CRISPR_dependency.tsv"
    ),
    sep = "\t"
)


fwrite(
    flatten_leading_edge(
        reactome
    ),
    file.path(
        OUT_DIR,
        "17b_Reactome_CRISPR_dependency.tsv"
    ),
    sep = "\t"
)



# ============================================================
# 10. Plot helper
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
            abs(
                NES
            )
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
                    rev(
                        pathway
                    )
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

            x =
                paste0(
                    "Normalized enrichment score\n",
                    "positive = ASTS-high less dependent; ",
                    "negative = ASTS-high more dependent"
                ),

            y =
                NULL
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



# ============================================================
# 11. Plots
# ============================================================

make_plot(

    hallmark,

    "ASTS-associated CRISPR dependency: Hallmark",

    file.path(
        OUT_DIR,
        "17b_Hallmark_CRISPR_dependency.png"
    )
)


make_plot(

    reactome,

    "ASTS-associated CRISPR dependency: Reactome",

    file.path(
        OUT_DIR,
        "17b_Reactome_CRISPR_dependency.png"
    )
)



# ============================================================
# 12. Summary
# ============================================================

summary_file <- file.path(
    OUT_DIR,
    "17b_summary.txt"
)


sink(
    summary_file
)


cat(
    "ASTRA-Drug Step 17B summary\n"
)

cat(
    "============================\n\n"
)


cat(
    "Ranking:\n"
)

cat(
    "Genome-wide CRISPR Gene Effect ~ ASTS stringent-model t statistic\n\n"
)


cat(
    "Interpretation:\n"
)

cat(
    "Positive NES = ASTS-high less dependent\n"
)

cat(
    "Negative NES = ASTS-high more dependent\n\n"
)


cat(
    "Hallmark\n"
)

cat(
    "--------\n"
)


cat(
    "Pathways tested: ",
    nrow(
        hallmark
    ),
    "\n",
    sep = ""
)


cat(
    "FDR < 0.05: ",
    sum(
        hallmark$padj
        <
        FDR_THRESHOLD,
        na.rm = TRUE
    ),
    "\n\n",
    sep = ""
)


print(
    hallmark[
        seq_len(
            min(
                30,
                .N
            )
        ),
        .(
            pathway,
            NES,
            pval,
            padj,
            direction
        )
    ]
)


cat(
    "\n\nReactome\n"
)

cat(
    "--------\n"
)


cat(
    "Pathways tested: ",
    nrow(
        reactome
    ),
    "\n",
    sep = ""
)


cat(
    "FDR < 0.05: ",
    sum(
        reactome$padj
        <
        FDR_THRESHOLD,
        na.rm = TRUE
    ),
    "\n\n",
    sep = ""
)


print(
    reactome[
        seq_len(
            min(
                50,
                .N
            )
        ),
        .(
            pathway,
            NES,
            pval,
            padj,
            direction
        )
    ]
)


sink()


cat(
    "\nStep 17B completed successfully\n"
)

cat(
    "Summary:\n",
    summary_file,
    "\n"
)
