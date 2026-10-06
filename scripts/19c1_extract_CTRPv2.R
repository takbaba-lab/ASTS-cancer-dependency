#!/usr/bin/env Rscript

# ============================================================
# ASTRA-Drug
# Step 19C-1
# Extract CTRPv2 drug sensitivity from PharmacoSet
# ============================================================
#
#
#
#
#
#
#   auc_recomputed
#   auc_published
#   aac_recomputed
#   aac_published
#
#
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


PSET_FILE <- file.path(
    PROJECT_DIR,
    "data/raw/ctrp/ORCESTRA/PSet_CTRPv2.rds"
)


OUT_DIR <- file.path(
    PROJECT_DIR,
    "results/19BC_mito_drug_validation/CTRPv2_extract"
)


PREFERRED_MEASURES <- c(
    "auc_recomputed",
    "auc_published",
    "aac_recomputed",
    "aac_published"
)


# ============================================================
# Packages
# ============================================================

suppressPackageStartupMessages({

    library(PharmacoGx)

})


dir.create(
    OUT_DIR,
    recursive = TRUE,
    showWarnings = FALSE
)


if (!file.exists(PSET_FILE)) {

    stop(
        "PSet file not found:\n",
        PSET_FILE
    )
}


cat(
    "Reading:\n",
    PSET_FILE,
    "\n\n"
)


# ============================================================
# Read and update object
# ============================================================

pset <- readRDS(
    PSET_FILE
)


cat(
    "Original class: ",
    paste(class(pset), collapse = ", "),
    "\n"
)


cat(
    "Updating PharmacoSet object...\n"
)


pset <- updateObject(
    pset
)


cat(
    "Updated successfully.\n\n"
)


# ============================================================
# Available sensitivity measures
# ============================================================

measures <- PharmacoGx::sensitivityMeasures(
    pset
)


cat(
    "Available sensitivity measures:\n"
)


print(
    measures
)


selected <- PREFERRED_MEASURES[
    PREFERRED_MEASURES
    %in%
    measures
]


if (length(selected) == 0) {

    stop(
        paste0(
            "No preferred AUC/AAC measure found.\n",
            "Available measures:\n",
            paste(measures, collapse = "\n")
        )
    )
}


selected <- selected[1]


cat(
    "\nSelected measure: ",
    selected,
    "\n",
    sep = ""
)


# ============================================================
# Extract metadata
# ============================================================

sample_info <- as.data.frame(
    PharmacoGx::sampleInfo(
        pset
    )
)


sample_info$`_rowname` <- rownames(
    sample_info
)


treatment_info <- as.data.frame(
    PharmacoGx::treatmentInfo(
        pset
    )
)


treatment_info$`_rowname` <- rownames(
    treatment_info
)


write.table(
    sample_info,
    file.path(
        OUT_DIR,
        "19C_CTRPv2_sample_info.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


write.table(
    treatment_info,
    file.path(
        OUT_DIR,
        "19C_CTRPv2_treatment_info.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# Curation tables
# ============================================================

cur <- PharmacoGx::curation(
    pset
)


sample_curation <- NULL

treatment_curation <- NULL


for (
    candidate
    in c(
        "sample",
        "cell",
        "cell.line"
    )
) {

    if (
        candidate
        %in%
        names(cur)
    ) {

        sample_curation <- as.data.frame(
            cur[[candidate]]
        )

        break
    }
}


for (
    candidate
    in c(
        "treatment",
        "drug"
    )
) {

    if (
        candidate
        %in%
        names(cur)
    ) {

        treatment_curation <- as.data.frame(
            cur[[candidate]]
        )

        break
    }
}


if (!is.null(sample_curation)) {

    sample_curation$`_rowname` <-
        rownames(sample_curation)

    write.table(
        sample_curation,
        file.path(
            OUT_DIR,
            "19C_CTRPv2_curation_sample.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}


if (!is.null(treatment_curation)) {

    treatment_curation$`_rowname` <-
        rownames(treatment_curation)

    write.table(
        treatment_curation,
        file.path(
            OUT_DIR,
            "19C_CTRPv2_curation_treatment.tsv"
        ),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
    )
}


# ============================================================
# Sensitivity profiles
# ============================================================

profiles <- as.data.frame(
    PharmacoGx::sensitivityProfiles(
        pset
    )
)


info <- as.data.frame(
    PharmacoGx::sensitivityInfo(
        pset
    )
)


if (
    !selected
    %in%
    colnames(profiles)
) {

    stop(
        paste0(
            "Selected measure not present ",
            "in sensitivityProfiles(): ",
            selected
        )
    )
}


# ------------------------------------------------------------
# Link profiles to sample/treatment metadata
# ------------------------------------------------------------

if (
    all(
        c(
            "sampleid",
            "treatmentid"
        )
        %in%
        colnames(profiles)
    )
) {

    response <- profiles[
        ,
        c(
            "sampleid",
            "treatmentid",
            selected
        ),
        drop = FALSE
    ]


} else if (
    nrow(profiles)
    ==
    nrow(info)
) {

    if (
        !all(
            c(
                "sampleid",
                "treatmentid"
            )
            %in%
            colnames(info)
        )
    ) {

        stop(
            paste0(
                "sensitivityInfo() does not contain ",
                "sampleid/treatmentid."
            )
        )
    }


    response <- data.frame(

        sampleid =
            info$sampleid,

        treatmentid =
            info$treatmentid,

        selected_measure =
            profiles[[selected]],

        stringsAsFactors = FALSE
    )


    colnames(response)[3] <- selected


} else {

    stop(
        paste0(
            "Could not align sensitivityProfiles ",
            "and sensitivityInfo.\n",
            "profiles rows = ",
            nrow(profiles),
            "\ninfo rows = ",
            nrow(info)
        )
    )
}


response$selected_measure <-
    as.numeric(
        response[[selected]]
    )


# ============================================================
# Harmonize direction
# ============================================================

if (
    grepl(
        "^aac",
        selected,
        ignore.case = TRUE
    )
) {

    response$response_harmonized <-
        -response$selected_measure

    sensitivity_direction <-
        "higher_original_value_more_sensitive; sign_flipped"


} else {

    response$response_harmonized <-
        response$selected_measure

    sensitivity_direction <-
        "lower_original_value_more_sensitive; unchanged"
}


response <- response[
    is.finite(
        response$response_harmonized
    ),
]


write.table(
    response[
        ,
        c(
            "sampleid",
            "treatmentid",
            "selected_measure",
            "response_harmonized"
        )
    ],
    file.path(
        OUT_DIR,
        "19C_CTRPv2_response_long.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# Measure metadata
# ============================================================

measure_info <- data.frame(

    selected_measure =
        selected,

    sensitivity_direction =
        sensitivity_direction,

    n_response_rows =
        nrow(response),

    stringsAsFactors = FALSE
)


write.table(
    measure_info,
    file.path(
        OUT_DIR,
        "19C_CTRPv2_measure_info.tsv"
    ),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)


# ============================================================
# Summary
# ============================================================

cat(
    "\n============================================================\n"
)

cat(
    "Step 19C-1 completed successfully\n"
)

cat(
    "============================================================\n"
)


cat(
    "Selected measure: ",
    selected,
    "\n",
    sep = ""
)


cat(
    "Direction: ",
    sensitivity_direction,
    "\n",
    sep = ""
)


cat(
    "Response rows: ",
    nrow(response),
    "\n",
    sep = ""
)
