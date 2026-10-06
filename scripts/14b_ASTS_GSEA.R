#!/usr/bin/env Rscript

# ============================================================
# ASTRA-Drug
# Step 14B: ASTS pathway characterization by GSEA
# ============================================================
#
#
#
# Primary ranking:
#
# Primary gene sets:
#     MSigDB Hallmark
#
# Secondary gene sets:
#     MSigDB Reactome
#
#
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


INPUT_DIR <- file.path(
    PROJECT_DIR,
    "results/14_ASTS_transcriptome"
)


PRIMARY_FILE <- file.path(
    INPUT_DIR,
    "14_gene_association_primary.tsv"
)


DRIVER_FILE <- file.path(
    INPUT_DIR,
    "14_gene_association_driver_adjusted.tsv"
)


OUT_DIR <- file.path(
    PROJECT_DIR,
    "results/14_ASTS_transcriptome/GSEA"
)

# ------------------------------------------------------------
# fgseaMultilevel initial permutation number
#
# ------------------------------------------------------------

N_PERM_SIMPLE <- 10000

# ------------------------------------------------------------
# Gene-set cache
#
# ------------------------------------------------------------

GENESET_CACHE_DIR <- file.path(
    PROJECT_DIR,
    "data/interim/msigdb"
)


HALLMARK_CACHE_FILE <- file.path(
    GENESET_CACHE_DIR,
    "Hallmark_Homo_sapiens.rds"
)


REACTOME_CACHE_FILE <- file.path(
    GENESET_CACHE_DIR,
    "Reactome_Homo_sapiens.rds"
)


# ------------------------------------------------------------
# GSEA parameters
# ------------------------------------------------------------

MIN_SIZE <- 15

MAX_SIZE <- 500

FDR_THRESHOLD <- 0.05

TOP_N_PLOT <- 20



# ============================================================
# 2. Packages
# ============================================================

suppressPackageStartupMessages({

    library(data.table)
    library(msigdbr)
    library(fgsea)
    library(ggplot2)

})


cat(
    "\n============================================================\n"
)

cat(
    "ASTRA-Drug Step 14B\n"
)

cat(
    "ASTS pathway characterization by GSEA\n"
)

cat(
    "============================================================\n\n"
)


cat(
    "R version: ",
    R.version.string,
    "\n",
    sep = ""
)


cat(
    "msigdbr version: ",
    as.character(
        packageVersion("msigdbr")
    ),
    "\n",
    sep = ""
)


cat(
    "fgsea version: ",
    as.character(
        packageVersion("fgsea")
    ),
    "\n\n",
    sep = ""
)



# ============================================================
# 3. Create directories
# ============================================================

dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


dir.create(
    GENESET_CACHE_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)



# ============================================================
# 4. Input check
# ============================================================

for (
    path
    in c(
        PRIMARY_FILE,
        DRIVER_FILE
    )
) {

    if (!file.exists(path)) {

        stop(
            "Required input file was not found:\n",
            path
        )
    }
}



# ============================================================
# 5. Helper:
#    Obtain MSigDB collection
# ============================================================

get_msig <- function(
    collection_name,
    subcollection_name = NULL
) {

    argument_names <- names(
        formals(
            msigdbr::msigdbr
        )
    )


    # --------------------------------------------------------
    # Current msigdbr API
    # --------------------------------------------------------

    if ("collection" %in% argument_names) {

        args <- list(

            species =
                "Homo sapiens",

            collection =
                collection_name
        )


        if (
            !is.null(
                subcollection_name
            )
            &&
            "subcollection"
            %in%
            argument_names
        ) {

            args$subcollection <-
                subcollection_name
        }


        result <- do.call(
            msigdbr::msigdbr,
            args
        )


    # --------------------------------------------------------
    # Compatibility with older msigdbr
    # --------------------------------------------------------

    } else {

        args <- list(

            species =
                "Homo sapiens",

            category =
                collection_name
        )


        if (
            !is.null(
                subcollection_name
            )
        ) {

            args$subcategory <-
                subcollection_name
        }


        result <- do.call(
            msigdbr::msigdbr,
            args
        )
    }


    return(
        result
    )
}



# ============================================================
# 6. Helper:
#    Read or download gene sets
# ============================================================

load_or_download_gene_sets <- function(
    cache_file,
    collection_name,
    subcollection_name = NULL,
    label
) {

    # --------------------------------------------------------
    # Use cached RDS if already available
    # --------------------------------------------------------

    if (file.exists(cache_file)) {

        cat(
            "Loading cached ",
            label,
            " gene sets:\n",
            cache_file,
            "\n\n",
            sep = ""
        )


        gene_sets <- readRDS(
            cache_file
        )


        return(
            gene_sets
        )
    }


    # --------------------------------------------------------
    # First use: download through msigdbr
    # --------------------------------------------------------

    cat(
        "Downloading ",
        label,
        " gene sets (first use only)...\n",
        sep = ""
    )


    gene_sets <- get_msig(

        collection_name =
            collection_name,

        subcollection_name =
            subcollection_name
    )


    if (
        is.null(gene_sets)
        ||
        nrow(gene_sets) == 0
    ) {

        stop(
            "No gene sets were returned for ",
            label,
            "."
        )
    }


    saveRDS(
        gene_sets,
        cache_file
    )


    cat(
        label,
        " gene sets saved to:\n",
        cache_file,
        "\n\n",
        sep = ""
    )


    return(
        gene_sets
    )
}



# ============================================================
# 7. Helper:
#    Convert ASTS signature flag safely
# ============================================================

as_logical_signature_flag <- function(x) {

    # Already logical
    if (is.logical(x)) {

        return(
            x
        )
    }


    x <- tolower(
        trimws(
            as.character(x)
        )
    )


    result <- rep(
        NA,
        length(x)
    )


    result[
        x %in%
        c(
            "true",
            "t",
            "1",
            "yes"
        )
    ] <- TRUE


    result[
        x %in%
        c(
            "false",
            "f",
            "0",
            "no"
        )
    ] <- FALSE


    return(
        result
    )
}



# ============================================================
# 8. Helper:
#    Flatten leadingEdge for TSV output
# ============================================================

prepare_fgsea_for_output <- function(result) {

    output <- copy(
        result
    )


    if (
        "leadingEdge"
        %in%
        colnames(output)
    ) {

        output[
            ,
            leadingEdge :=
                vapply(
                    leadingEdge,
                    function(x) {

                        paste(
                            x,
                            collapse = ";"
                        )
                    },
                    character(1)
                )
        ]
    }


    return(
        output
    )
}



# ============================================================
# 9. Helper:
#    Run GSEA
# ============================================================

run_gsea <- function(
    input_file,
    analysis_name,
    pathways
) {

    cat(
        "\n------------------------------------------------------------\n"
    )

    cat(
        "Running: ",
        analysis_name,
        "\n",
        sep = ""
    )


    # --------------------------------------------------------
    # Read Step14A result
    # --------------------------------------------------------

    x <- fread(
        input_file
    )


    required_columns <- c(

        "gene_symbol",

        "t_stat",

        "is_ASTS_signature"
    )


    missing_columns <- setdiff(
        required_columns,
        colnames(x)
    )


    if (
        length(
            missing_columns
        ) > 0
    ) {

        stop(
            "Required columns missing from ",
            input_file,
            ":\n",
            paste(
                missing_columns,
                collapse = ", "
            )
        )
    }


    cat(
        "Input genes: ",
        nrow(x),
        "\n",
        sep = ""
    )


    # --------------------------------------------------------
    # Convert t-statistic to numeric
    # --------------------------------------------------------

    x[
        ,
        t_stat :=
            as.numeric(
                t_stat
            )
    ]


    # --------------------------------------------------------
    # Signature flag
    # --------------------------------------------------------

    x[
        ,
        is_ASTS_signature :=
            as_logical_signature_flag(
                is_ASTS_signature
            )
    ]


    # --------------------------------------------------------
    # Remove direct ASTS signature genes
    # --------------------------------------------------------

    n_signature <- sum(
        x$is_ASTS_signature %in% TRUE,
        na.rm = TRUE
    )


    cat(
        "ASTS signature genes removed: ",
        n_signature,
        "\n",
        sep = ""
    )


    x <- x[
        is_ASTS_signature == FALSE
        |
        is.na(
            is_ASTS_signature
        )
    ]


    # --------------------------------------------------------
    # Remove invalid t statistics
    #
    # --------------------------------------------------------

    before_finite <- nrow(
        x
    )


    x <- x[
        is.finite(
            t_stat
        )
    ]


    cat(
        "Non-finite t-stat genes removed: ",
        before_finite - nrow(x),
        "\n",
        sep = ""
    )


    # --------------------------------------------------------
    # Remove empty gene symbols
    # --------------------------------------------------------

    x <- x[
        !is.na(
            gene_symbol
        )
        &
        gene_symbol != ""
    ]


    # --------------------------------------------------------
    # Duplicate gene-symbol handling
    #
    #
    # NOTE:
    # --------------------------------------------------------

    x[
        ,
        abs_t_stat :=
            abs(
                t_stat
            )
    ]


    setorder(
        x,
        -abs_t_stat
    )


    duplicated_gene_n <- sum(
        duplicated(
            x$gene_symbol
        )
    )


    x <- x[
        !duplicated(
            gene_symbol
        )
    ]


    x[
        ,
        abs_t_stat :=
            NULL
    ]


    cat(
        "Duplicated gene symbols removed: ",
        duplicated_gene_n,
        "\n",
        sep = ""
    )


    cat(
        "Genes used for GSEA: ",
        nrow(x),
        "\n",
        sep = ""
    )


    if (
        nrow(x) < 1000
    ) {

        stop(
            "Too few genes remain for GSEA: ",
            nrow(x)
        )
    }


    # --------------------------------------------------------
    # Create ranked vector
    # --------------------------------------------------------

    ranks <- x$t_stat


    names(
        ranks
    ) <- x$gene_symbol


    ranks <- sort(
        ranks,
        decreasing = TRUE
    )


    # safety check
    if (
        any(
            duplicated(
                names(ranks)
            )
        )
    ) {

        stop(
            "Duplicated gene symbols remain in ranking."
        )
    }


    if (
        any(
            !is.finite(
                ranks
            )
        )
    ) {

        stop(
            "Non-finite values remain in ranking."
        )
    }


    # --------------------------------------------------------
    # GSEA
    # --------------------------------------------------------

    result <- fgseaMultilevel(

        pathways = pathways,

        stats = ranks,

        minSize = MIN_SIZE,

        maxSize = MAX_SIZE,

        eps = 0,

        nPermSimple = N_PERM_SIMPLE
    )


        result <- as.data.table(
            result
        )
    

    # ============================================================
    # Remove pathways with non-finite GSEA statistics
    # ============================================================
    #
    #
    #
    # ============================================================

    n_result_all <- nrow(
        result
    )


    n_nonfinite <- sum(
        !is.finite(result$NES)
        |
        !is.finite(result$pval)
        |
        !is.finite(result$padj)
    )


    cat(
        "Pathways returned by fgsea: ",
        n_result_all,
        "\n",
        sep = ""
    )


    cat(
        "Pathways with non-finite statistics removed: ",
        n_nonfinite,
        "\n",
        sep = ""
    )


    result <- result[
        is.finite(NES)
        &
        is.finite(pval)
        &
        is.finite(padj)
    ]


    if (
        nrow(result) == 0
    ) {

        stop(
            "No pathways with finite GSEA statistics remain for ",
            analysis_name
        )
    }


    # ============================================================
    # Sort valid pathways
    # ============================================================

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


    # ============================================================
    # Direction
    # ============================================================

    result[
        ,
        direction :=
            fifelse(
                NES > 0,
                "ASTS_high",
                "ASTS_low"
            )
    ]


    # --------------------------------------------------------
    # Save
    # --------------------------------------------------------

    output_result <- prepare_fgsea_for_output(
        result
    )


    outfile <- file.path(

        OUT_DIR,

        paste0(
            analysis_name,
            ".tsv"
        )
    )


    fwrite(

        output_result,

        outfile,

        sep = "\t",

        quote = FALSE,

        na = "NA"
    )


    cat(
        "Pathways tested: ",
        nrow(result),
        "\n",
        sep = ""
    )


    cat(
        "FDR < 0.05: ",
        sum(
            result$padj < FDR_THRESHOLD,
            na.rm = TRUE
        ),
        "\n",
        sep = ""
    )


    return(
        result
    )
}



# ============================================================
# 10. Load Hallmark gene sets
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "Loading Hallmark gene sets\n"
)

cat(
    "============================================================\n"
)


hallmark_df <- load_or_download_gene_sets(

    cache_file =
        HALLMARK_CACHE_FILE,

    collection_name =
        "H",

    subcollection_name =
        NULL,

    label =
        "Hallmark"
)


if (
    !all(
        c(
            "gene_symbol",
            "gs_name"
        )
        %in%
        colnames(
            hallmark_df
        )
    )
) {

    stop(
        "Hallmark msigdbr output lacks gene_symbol or gs_name."
    )
}


hallmark_pathways <- split(

    hallmark_df$gene_symbol,

    hallmark_df$gs_name
)


hallmark_pathways <- lapply(
    hallmark_pathways,
    unique
)


cat(
    "Hallmark pathways loaded: ",
    length(
        hallmark_pathways
    ),
    "\n",
    sep = ""
)



# ============================================================
# 11. Load Reactome gene sets
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "Loading Reactome gene sets\n"
)

cat(
    "============================================================\n"
)


reactome_df <- load_or_download_gene_sets(

    cache_file =
        REACTOME_CACHE_FILE,

    collection_name =
        "C2",

    subcollection_name =
        "CP:REACTOME",

    label =
        "Reactome"
)


if (
    !all(
        c(
            "gene_symbol",
            "gs_name"
        )
        %in%
        colnames(
            reactome_df
        )
    )
) {

    stop(
        "Reactome msigdbr output lacks gene_symbol or gs_name."
    )
}


reactome_pathways <- split(

    reactome_df$gene_symbol,

    reactome_df$gs_name
)


reactome_pathways <- lapply(
    reactome_pathways,
    unique
)


cat(
    "Reactome pathways loaded: ",
    length(
        reactome_pathways
    ),
    "\n",
    sep = ""
)



# ============================================================
# 12. Primary GSEA
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "Primary GSEA\n"
)

cat(
    "============================================================\n"
)


hallmark_primary <- run_gsea(

    PRIMARY_FILE,

    "14b_Hallmark_primary",

    hallmark_pathways
)


reactome_primary <- run_gsea(

    PRIMARY_FILE,

    "14b_Reactome_primary",

    reactome_pathways
)



# ============================================================
# 13. Driver-adjusted GSEA
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "Driver-adjusted GSEA\n"
)

cat(
    "============================================================\n"
)


hallmark_driver <- run_gsea(

    DRIVER_FILE,

    "14b_Hallmark_driver_adjusted",

    hallmark_pathways
)


reactome_driver <- run_gsea(

    DRIVER_FILE,

    "14b_Reactome_driver_adjusted",

    reactome_pathways
)



# ============================================================
# 14. Compare Hallmark:
#     primary vs driver adjusted
# ============================================================

hallmark_compare <- merge(

    hallmark_primary[
        ,
        .(
            pathway,
            NES_primary = NES,
            padj_primary = padj
        )
    ],

    hallmark_driver[
        ,
        .(
            pathway,
            NES_driver = NES,
            padj_driver = padj
        )
    ],

    by =
        "pathway",

    all =
        FALSE
)


hallmark_compare[
    ,
    sign_concordant :=
        sign(
            NES_primary
        )
        ==
        sign(
            NES_driver
        )
]


hallmark_rho <- cor(

    hallmark_compare$NES_primary,

    hallmark_compare$NES_driver,

    method =
        "spearman",

    use =
        "complete.obs"
)


fwrite(

    hallmark_compare,

    file.path(
        OUT_DIR,
        "14b_Hallmark_primary_vs_driver.tsv"
    ),

    sep = "\t"
)



# ============================================================
# 15. Compare Reactome:
#     primary vs driver adjusted
# ============================================================

reactome_compare <- merge(

    reactome_primary[
        ,
        .(
            pathway,
            NES_primary = NES,
            padj_primary = padj
        )
    ],

    reactome_driver[
        ,
        .(
            pathway,
            NES_driver = NES,
            padj_driver = padj
        )
    ],

    by =
        "pathway",

    all =
        FALSE
)


reactome_compare[
    ,
    sign_concordant :=
        sign(
            NES_primary
        )
        ==
        sign(
            NES_driver
        )
]


reactome_rho <- cor(

    reactome_compare$NES_primary,

    reactome_compare$NES_driver,

    method =
        "spearman",

    use =
        "complete.obs"
)


fwrite(

    reactome_compare,

    file.path(
        OUT_DIR,
        "14b_Reactome_primary_vs_driver.tsv"
    ),

    sep = "\t"
)



# ============================================================
# 16. Plot helper
# ============================================================

make_gsea_plot <- function(
    result,
    title_text,
    output_file,
    top_n = 20
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
                top_n,
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

            subtitle =
                "78 ASTS signature genes excluded",

            x =
                paste0(
                    "Normalized enrichment score\n",
                    "positive = ASTS-high; ",
                    "negative = ASTS-low"
                ),

            y =
                NULL
        ) +

        theme_bw(
            base_size = 11
        )


    ggsave(

        output_file,

        p,

        width = 10,

        height = 8,

        dpi = 200
    )
}



# ============================================================
# 17. Plots
# ============================================================

make_gsea_plot(

    hallmark_primary,

    "ASTS-associated Hallmark pathways",

    file.path(
        OUT_DIR,
        "14b_Hallmark_primary.png"
    ),

    TOP_N_PLOT
)


make_gsea_plot(

    reactome_primary,

    "ASTS-associated Reactome pathways",

    file.path(
        OUT_DIR,
        "14b_Reactome_primary.png"
    ),

    TOP_N_PLOT
)


make_gsea_plot(

    hallmark_driver,

    "Driver-adjusted ASTS-associated Hallmark pathways",

    file.path(
        OUT_DIR,
        "14b_Hallmark_driver_adjusted.png"
    ),

    TOP_N_PLOT
)



# ============================================================
# 18. Summary
# ============================================================

summary_file <- file.path(
    OUT_DIR,
    "14b_summary.txt"
)


sink(
    summary_file
)


cat(
    "ASTRA-Drug Step 14B summary\n"
)

cat(
    "============================\n\n"
)


cat(
    "Primary GSEA ranking:\n"
)

cat(
    "Step14A gene-wise ASTS t-statistic\n\n"
)


cat(
    "Important:\n"
)

cat(
    "All 78 ASTS signature genes were excluded before GSEA.\n"
)

cat(
    "This reduces direct circular enrichment from genes used to calculate ASTS.\n\n"
)


cat(
    "msigdbr version: ",
    as.character(
        packageVersion("msigdbr")
    ),
    "\n\n",
    sep = ""
)


# ------------------------------------------------------------
# Hallmark
# ------------------------------------------------------------

cat(
    "Hallmark primary\n"
)

cat(
    "----------------\n"
)


cat(
    "Pathways tested: ",
    nrow(
        hallmark_primary
    ),
    "\n",
    sep = ""
)


cat(
    "FDR < 0.05: ",
    sum(
        hallmark_primary$padj
        <
        FDR_THRESHOLD,
        na.rm = TRUE
    ),
    "\n\n",
    sep = ""
)


print(

    hallmark_primary[
        seq_len(
            min(
                20,
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


# ------------------------------------------------------------
# Reactome
# ------------------------------------------------------------

cat(
    "\n\nReactome primary\n"
)

cat(
    "----------------\n"
)


cat(
    "Pathways tested: ",
    nrow(
        reactome_primary
    ),
    "\n",
    sep = ""
)


cat(
    "FDR < 0.05: ",
    sum(
        reactome_primary$padj
        <
        FDR_THRESHOLD,
        na.rm = TRUE
    ),
    "\n\n",
    sep = ""
)


print(

    reactome_primary[
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


# ------------------------------------------------------------
# Driver robustness
# ------------------------------------------------------------

cat(
    "\n\nHallmark primary vs driver-adjusted\n"
)

cat(
    "-----------------------------------\n"
)


cat(
    "Pathways compared: ",
    nrow(
        hallmark_compare
    ),
    "\n",
    sep = ""
)


cat(
    "NES Spearman r: ",
    hallmark_rho,
    "\n",
    sep = ""
)


cat(
    "Direction concordance: ",
    mean(
        hallmark_compare$sign_concordant,
        na.rm = TRUE
    ),
    "\n",
    sep = ""
)


cat(
    "\nReactome primary vs driver-adjusted\n"
)

cat(
    "-----------------------------------\n"
)


cat(
    "Pathways compared: ",
    nrow(
        reactome_compare
    ),
    "\n",
    sep = ""
)


cat(
    "NES Spearman r: ",
    reactome_rho,
    "\n",
    sep = ""
)


cat(
    "Direction concordance: ",
    mean(
        reactome_compare$sign_concordant,
        na.rm = TRUE
    ),
    "\n",
    sep = ""
)


sink()



# ============================================================
# 19. Session info
# ============================================================

sink(
    file.path(
        OUT_DIR,
        "14b_sessionInfo.txt"
    )
)


print(
    sessionInfo()
)


sink()



# ============================================================
# 20. Finished
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "Step 14B completed successfully\n"
)

cat(
    "============================================================\n\n"
)


cat(
    "Summary:\n",
    summary_file,
    "\n"
)
