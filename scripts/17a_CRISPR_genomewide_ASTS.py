#!/usr/bin/env python3

"""
ASTRA-Drug
Step 17A: Genome-wide CRISPR dependency landscape

----

    Gene Effect
        ~ ASTS
        + lineage
        + sex
        + proliferation proxy
        + MAPK driver
        + ERBB driver
        + PI3K driver
        + E2F/G2M cell-cycle score
        + MYC score


----



    ASTS beta > 0

    ASTS beta < 0


----


"""


# ============================================================
# 1. Path / Parameter
# ============================================================

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]


CRISPR_FILE = (
    PROJECT_DIR
    / "data/raw/depmap/26Q1"
    / "CRISPRGeneEffect.csv"
)


ASTS_FILE = (
    PROJECT_DIR
    / "results/08b_DepMap_ASTS_default_profile"
    / "08_DepMap_ASTS_scores.tsv"
)


DRIVER_FILE = (
    PROJECT_DIR
    / "results/13_driver_confounder"
    / "13_driver_flags.tsv"
)


CELL_CYCLE_FILE = (
    PROJECT_DIR
    / "results/15_cellcycle_confounder"
    / "15_cellcycle_scores.tsv"
)


OUT_DIR = (
    PROJECT_DIR
    / "results/17_CRISPR_genomewide"
)


# ------------------------------------------------------------
# Metadata columns
# ------------------------------------------------------------

ASTS_COLUMN = "ASTS_score"

PROLIFERATION_COLUMN = "proliferation_proxy"

LINEAGE_COLUMN = "OncotreeLineage"

SEX_COLUMN = "Sex"


# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

MIN_GLOBAL_LINEAGE_N = 5

MIN_MODELS_PER_GENE = 800

GENE_CHUNK_SIZE = 500

FDR_ALPHA = 0.05


# ============================================================
# 2. Imports
# ============================================================

import re

import numpy as np
import pandas as pd

from scipy.stats import t as t_distribution

from statsmodels.stats.multitest import multipletests


# ============================================================
# 3. Helper functions
# ============================================================

def zscore(series):
    """
    """

    x = pd.to_numeric(
        series,
        errors="coerce"
    )

    sd = x.std(
        ddof=1
    )

    if (
        not np.isfinite(sd)
        or
        sd == 0
    ):

        return pd.Series(
            np.nan,
            index=x.index
        )

    return (
        x - x.mean()
    ) / sd


def extract_gene_symbol(column_name):
    """

        EGFR (1956)
        MAP2K1 (5604)

    """

    x = str(
        column_name
    ).strip()

    match = re.match(
        r"^(.+?)\s+\(\d+\)$",
        x
    )

    if match:

        return (
            match
            .group(1)
            .strip()
        )

    return x


def detect_model_id_column(columns):
    """
    """

    candidates = [
        "ModelID",
        "DepMap_ID",
        "DepMapID",
        "model_id",
    ]

    for col in candidates:

        if col in columns:

            return col

    print(
        "WARNING: standard ModelID column not found."
    )

    print(
        "Using first column:",
        columns[0]
    )

    return columns[0]


def build_lineage_categories(metadata):
    """
    """

    counts = (
        metadata[
            LINEAGE_COLUMN
        ]
        .fillna("Unknown")
        .value_counts()
    )


    common = set(
        counts[
            counts
            >=
            MIN_GLOBAL_LINEAGE_N
        ].index
    )


    return (
        metadata[
            LINEAGE_COLUMN
        ]
        .fillna("Unknown")
        .astype(str)
        .apply(
            lambda x:
                x
                if x in common
                else "OtherRare"
        )
    )


def build_covariate_matrix(metadata):
    """

    """

    parts = []


    # Intercept
    parts.append(
        pd.DataFrame(
            {
                "Intercept":
                    np.ones(
                        len(metadata)
                    )
            },
            index=metadata.index
        )
    )


    # Continuous covariates
    parts.append(
        pd.DataFrame(
            {
                "Prolif_z":
                    zscore(
                        metadata[
                            PROLIFERATION_COLUMN
                        ]
                    ),

                "CellCycle_z":
                    zscore(
                        metadata[
                            "CellCycle_score"
                        ]
                    ),

                "MYC_z":
                    zscore(
                        metadata[
                            "MYC_V1_score"
                        ]
                    ),
            },
            index=metadata.index
        )
    )


    # Drivers
    parts.append(
        metadata[
            [
                "MAPK_driver",
                "ERBB_driver",
                "PI3K_driver",
            ]
        ].astype(float)
    )


    # Sex
    sex = (
        metadata[
            SEX_COLUMN
        ]
        .fillna("Unknown")
        .astype(str)
    )


    if sex.nunique() > 1:

        parts.append(
            pd.get_dummies(
                sex,
                prefix="Sex",
                drop_first=True,
                dtype=float
            )
        )


    # Lineage
    lineage = (
        metadata[
            "LineageModel"
        ]
        .fillna("Unknown")
        .astype(str)
    )


    if lineage.nunique() > 1:

        parts.append(
            pd.get_dummies(
                lineage,
                prefix="Lineage",
                drop_first=True,
                dtype=float
            )
        )


    X = pd.concat(
        parts,
        axis=1
    )


    return X.astype(float)


def fit_single_gene_with_missing(
    y,
    metadata
):
    """

    """

    valid = np.isfinite(
        y
    )


    if valid.sum() < MIN_MODELS_PER_GENE:

        return None


    meta = metadata.iloc[
        np.where(valid)[0]
    ].copy()


    y_valid = y[
        valid
    ]


    asts = zscore(
        meta[
            ASTS_COLUMN
        ]
    ).to_numpy(
        dtype=float
    )


    C = build_covariate_matrix(
        meta
    ).to_numpy(
        dtype=float
    )


    if (
        not np.all(
            np.isfinite(C)
        )
        or
        not np.all(
            np.isfinite(asts)
        )
    ):

        return None


    Q, _ = np.linalg.qr(
        C,
        mode="reduced"
    )


    asts_res = (
        asts
        -
        Q
        @ (
            Q.T
            @ asts
        )
    )


    asts_ss = np.sum(
        asts_res ** 2
    )


    if asts_ss <= 0:

        return None


    y_res = (
        y_valid
        -
        Q
        @ (
            Q.T
            @ y_valid
        )
    )


    beta = (
        asts_res
        @ y_res
    ) / asts_ss


    residual = (
        y_res
        -
        asts_res
        *
        beta
    )


    df_resid = (
        len(y_valid)
        -
        C.shape[1]
        -
        1
    )


    if df_resid <= 0:

        return None


    sigma2 = (
        np.sum(
            residual ** 2
        )
        /
        df_resid
    )


    se = np.sqrt(
        sigma2
        /
        asts_ss
    )


    if (
        not np.isfinite(se)
        or
        se == 0
    ):

        return None


    t_stat = (
        beta
        /
        se
    )


    p_value = (
        2.0
        *
        t_distribution.sf(
            abs(
                t_stat
            ),
            df=df_resid
        )
    )


    partial_r = (
        t_stat
        /
        np.sqrt(
            t_stat ** 2
            +
            df_resid
        )
    )


    y_sd = np.std(
        y_valid,
        ddof=1
    )


    standardized_beta = (
        beta
        /
        y_sd
        if y_sd > 0
        else np.nan
    )


    return {
        "n_models":
            len(y_valid),

        "beta":
            beta,

        "standardized_beta":
            standardized_beta,

        "SE":
            se,

        "t_stat":
            t_stat,

        "p_value":
            p_value,

        "partial_r":
            partial_r,
    }


# ============================================================
# 4. Main
# ============================================================

def main():

    print(
        "\n"
        "============================================================\n"
        "ASTRA-Drug Step 17A\n"
        "Genome-wide CRISPR dependency landscape\n"
        "============================================================"
    )


    # ========================================================
    # Input check
    # ========================================================

    for path in [
        CRISPR_FILE,
        ASTS_FILE,
        DRIVER_FILE,
        CELL_CYCLE_FILE,
    ]:

        if not path.exists():

            raise FileNotFoundError(
                f"Required file not found:\n{path}"
            )


    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )


    # ========================================================
    # Step 1: Metadata
    # ========================================================

    print(
        "\nStep 1: Preparing ASTS metadata"
    )


    asts = pd.read_csv(
        ASTS_FILE,
        sep="\t"
    )


    if asts[
        "ModelID"
    ].duplicated().any():

        raise RuntimeError(
            "Duplicated ModelID in ASTS file."
        )


    asts[
        "LineageModel"
    ] = build_lineage_categories(
        asts
    )


    drivers = pd.read_csv(
        DRIVER_FILE,
        sep="\t"
    )


    drivers = drivers[
        [
            "ModelID",
            "MAPK_driver",
            "ERBB_driver",
            "PI3K_driver",
        ]
    ].copy()


    cellcycle = pd.read_csv(
        CELL_CYCLE_FILE,
        sep="\t"
    )


    cellcycle = cellcycle[
        [
            "ModelID",
            "CellCycle_score",
            "MYC_V1_score",
        ]
    ].copy()


    metadata = (
        asts
        .merge(
            drivers,
            on="ModelID",
            how="left",
            validate="1:1"
        )
        .merge(
            cellcycle,
            on="ModelID",
            how="left",
            validate="1:1"
        )
    )


    required_metadata = [
        ASTS_COLUMN,
        PROLIFERATION_COLUMN,
        "LineageModel",
        "MAPK_driver",
        "ERBB_driver",
        "PI3K_driver",
        "CellCycle_score",
        "MYC_V1_score",
    ]


    metadata = metadata.dropna(
        subset=required_metadata
    ).copy()


    print(
        "Metadata-complete ASTS models:",
        len(
            metadata
        )
    )


    # ========================================================
    # Step 2: CRISPR header
    # ========================================================

    print(
        "\nStep 2: Reading CRISPR header"
    )


    header = pd.read_csv(
        CRISPR_FILE,
        nrows=0
    )


    all_columns = (
        header
        .columns
        .tolist()
    )


    model_id_column = (
        detect_model_id_column(
            all_columns
        )
    )


    gene_columns = [
        col
        for col in all_columns
        if col != model_id_column
    ]


    print(
        "CRISPR gene columns:",
        len(
            gene_columns
        )
    )


    # ========================================================
    # Step 3: Read full CRISPR matrix
    # ========================================================

    print(
        "\nStep 3: Reading full CRISPR matrix"
    )


    crispr = pd.read_csv(
        CRISPR_FILE,
        low_memory=False
    )


    crispr = crispr.rename(
        columns={
            model_id_column:
                "ModelID"
        }
    )


    if crispr[
        "ModelID"
    ].duplicated().any():

        raise RuntimeError(
            "Duplicated ModelID in CRISPR matrix."
        )


    data = crispr.merge(
        metadata,
        on="ModelID",
        how="inner",
        validate="1:1"
    )


    print(
        "CRISPR x metadata models:",
        len(
            data
        )
    )


    # ========================================================
    # Step 4: Fixed covariates
    # ========================================================

    print(
        "\nStep 4: Building stringent covariate model"
    )


    metadata_aligned = data[
        metadata.columns
    ].copy()


    asts_z = zscore(
        metadata_aligned[
            ASTS_COLUMN
        ]
    ).to_numpy(
        dtype=float
    )


    C = build_covariate_matrix(
        metadata_aligned
    ).to_numpy(
        dtype=float
    )


    if not np.all(
        np.isfinite(C)
    ):

        raise RuntimeError(
            "Non-finite covariate values remain."
        )


    Q, _ = np.linalg.qr(
        C,
        mode="reduced"
    )


    asts_res = (
        asts_z
        -
        Q
        @ (
            Q.T
            @ asts_z
        )
    )


    asts_ss = np.sum(
        asts_res ** 2
    )


    df_resid_complete = (
        len(data)
        -
        C.shape[1]
        -
        1
    )


    print(
        "Covariate columns:",
        C.shape[1]
    )


    print(
        "Residual df:",
        df_resid_complete
    )


    # ========================================================
    # Step 5: Genome-wide analysis
    # ========================================================

    print(
        "\nStep 5: Genome-wide gene-effect association"
    )


    results = []


    for start in range(
        0,
        len(gene_columns),
        GENE_CHUNK_SIZE
    ):

        end = min(
            start
            +
            GENE_CHUNK_SIZE,
            len(gene_columns)
        )


        cols = gene_columns[
            start:end
        ]


        print(
            f"Genes {start + 1}-{end} / "
            f"{len(gene_columns)}"
        )


        # force numeric
        Y_df = data[
            cols
        ].apply(
            pd.to_numeric,
            errors="coerce"
        )


        Y = Y_df.to_numpy(
            dtype=float
        )


        complete_mask = np.all(
            np.isfinite(Y),
            axis=0
        )


        # ----------------------------------------------------
        # Complete genes: vectorized FWL
        # ----------------------------------------------------

        if complete_mask.any():

            complete_indices = np.where(
                complete_mask
            )[0]


            Y_complete = Y[
                :,
                complete_indices
            ]


            Y_res = (
                Y_complete
                -
                Q
                @ (
                    Q.T
                    @ Y_complete
                )
            )


            beta = (
                asts_res
                @ Y_res
            ) / asts_ss


            residual = (
                Y_res
                -
                asts_res[
                    :,
                    None
                ]
                *
                beta[
                    None,
                    :
                ]
            )


            sse = np.sum(
                residual ** 2,
                axis=0
            )


            sigma2 = (
                sse
                /
                df_resid_complete
            )


            se = np.sqrt(
                sigma2
                /
                asts_ss
            )


            with np.errstate(
                divide="ignore",
                invalid="ignore"
            ):

                t_stat = (
                    beta
                    /
                    se
                )


            p_value = (
                2.0
                *
                t_distribution.sf(
                    np.abs(
                        t_stat
                    ),
                    df=df_resid_complete
                )
            )


            partial_r = (
                t_stat
                /
                np.sqrt(
                    t_stat ** 2
                    +
                    df_resid_complete
                )
            )


            y_sd = np.std(
                Y_complete,
                axis=0,
                ddof=1
            )


            with np.errstate(
                divide="ignore",
                invalid="ignore"
            ):

                standardized_beta = (
                    beta
                    /
                    y_sd
                )


            for j, local_index in enumerate(
                complete_indices
            ):

                col = cols[
                    local_index
                ]


                results.append(
                    {
                        "gene_symbol":
                            extract_gene_symbol(
                                col
                            ),

                        "CRISPR_column":
                            col,

                        "n_models":
                            len(data),

                        "beta":
                            beta[j],

                        "standardized_beta":
                            standardized_beta[j],

                        "SE":
                            se[j],

                        "t_stat":
                            t_stat[j],

                        "p_value":
                            p_value[j],

                        "partial_r":
                            partial_r[j],

                        "has_missing_gene_effect":
                            False,
                    }
                )


        # ----------------------------------------------------
        # Genes with missing values
        # ----------------------------------------------------

        incomplete_indices = np.where(
            ~complete_mask
        )[0]


        for local_index in incomplete_indices:

            col = cols[
                local_index
            ]


            y = Y[
                :,
                local_index
            ]


            result = fit_single_gene_with_missing(
                y,
                metadata_aligned
            )


            if result is None:

                continue


            results.append(
                {
                    "gene_symbol":
                        extract_gene_symbol(
                            col
                        ),

                    "CRISPR_column":
                        col,

                    **result,

                    "has_missing_gene_effect":
                        True,
                }
            )


    result_df = pd.DataFrame(
        results
    )


    # ========================================================
    # Step 6: Remove duplicate symbols conservatively
    # ========================================================

    print(
        "\nStep 6: Auditing duplicate gene symbols"
    )


    duplicate_symbols = (
        result_df[
            "gene_symbol"
        ]
        .duplicated(
            keep=False
        )
    )


    n_duplicate_rows = int(
        duplicate_symbols.sum()
    )


    print(
        "Rows with duplicated gene symbols:",
        n_duplicate_rows
    )


    # retain row with largest N, then largest |t|
    result_df[
        "abs_t"
    ] = result_df[
        "t_stat"
    ].abs()


    result_df = (
        result_df
        .sort_values(
            [
                "gene_symbol",
                "n_models",
                "abs_t",
            ],
            ascending=[
                True,
                False,
                False,
            ]
        )
        .drop_duplicates(
            "gene_symbol",
            keep="first"
        )
        .copy()
    )


    result_df = result_df.drop(
        columns=[
            "abs_t"
        ]
    )


    # ========================================================
    # Step 7: Multiple testing
    # ========================================================

    finite = np.isfinite(
        result_df[
            "p_value"
        ]
    )


    result_df[
        "FDR"
    ] = np.nan


    result_df.loc[
        finite,
        "FDR"
    ] = multipletests(
        result_df.loc[
            finite,
            "p_value"
        ],
        method="fdr_bh"
    )[1]


    result_df[
        "dependency_direction"
    ] = np.where(
        result_df[
            "beta"
        ] > 0,
        "ASTS_high_less_dependent",
        "ASTS_high_more_dependent"
    )


    # positive ranking
    result_df[
        "rank_by_t_desc"
    ] = (
        result_df[
            "t_stat"
        ]
        .rank(
            ascending=False,
            method="min"
        )
    )


    result_df = result_df.sort_values(
        "t_stat",
        ascending=False
    )


    # ========================================================
    # Step 8: Save
    # ========================================================

    output_file = (
        OUT_DIR
        / "17_CRISPR_genomewide_ASTS_association.tsv"
    )


    result_df.to_csv(
        output_file,
        sep="\t",
        index=False
    )


    # ========================================================
    # Step 9: Summary
    # ========================================================

    summary_file = (
        OUT_DIR
        / "17a_summary.txt"
    )


    focus_genes = [
        "EGFR",
        "ERBB2",
        "MAP2K1",
        "MAP2K2",
        "RAF1",
        "BRAF",
        "SOS1",
        "SHOC2",
        "PTPN11",
        "GRB2",
    ]


    with open(
        summary_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 17A summary",
            file=f
        )

        print(
            "============================",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "Primary model:",
            file=f
        )

        print(
            "CRISPR gene effect ~ ASTS + lineage + sex "
            "+ proliferation + drivers + CellCycle + MYC",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            f"CRISPR x metadata models: "
            f"{len(data)}",
            file=f
        )


        print(
            f"Genes analyzed: "
            f"{len(result_df)}",
            file=f
        )


        print(
            f"FDR < {FDR_ALPHA}: "
            f"{(result_df['FDR'] < FDR_ALPHA).sum()}",
            file=f
        )


        print(
            f"ASTS-high less dependent, FDR < {FDR_ALPHA}: "
            f"{((result_df['FDR'] < FDR_ALPHA) & (result_df['beta'] > 0)).sum()}",
            file=f
        )


        print(
            f"ASTS-high more dependent, FDR < {FDR_ALPHA}: "
            f"{((result_df['FDR'] < FDR_ALPHA) & (result_df['beta'] < 0)).sum()}",
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Top genes: ASTS-high less dependent",
            file=f
        )

        print(
            "-----------------------------------",
            file=f
        )


        print(
            result_df[
                [
                    "gene_symbol",
                    "n_models",
                    "standardized_beta",
                    "t_stat",
                    "p_value",
                    "FDR",
                ]
            ]
            .head(
                30
            )
            .to_string(
                index=False
            ),
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Top genes: ASTS-high more dependent",
            file=f
        )

        print(
            "-----------------------------------",
            file=f
        )


        print(
            result_df[
                [
                    "gene_symbol",
                    "n_models",
                    "standardized_beta",
                    "t_stat",
                    "p_value",
                    "FDR",
                ]
            ]
            .sort_values(
                "t_stat"
            )
            .head(
                30
            )
            .to_string(
                index=False
            ),
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Pre-specified signaling genes",
            file=f
        )

        print(
            "-----------------------------",
            file=f
        )


        focus = (
            result_df[
                result_df[
                    "gene_symbol"
                ].isin(
                    focus_genes
                )
            ]
            .copy()
        )


        order_map = {
            gene:
                index
            for index, gene
            in enumerate(
                focus_genes
            )
        }


        focus[
            "_order"
        ] = focus[
            "gene_symbol"
        ].map(
            order_map
        )


        focus = focus.sort_values(
            "_order"
        )


        print(
            focus[
                [
                    "gene_symbol",
                    "standardized_beta",
                    "t_stat",
                    "p_value",
                    "FDR",
                    "rank_by_t_desc",
                    "dependency_direction",
                ]
            ].to_string(
                index=False
            ),
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Interpretation:",
            file=f
        )

        print(
            "Positive beta = ASTS-high has weaker dependency.",
            file=f
        )

        print(
            "Negative beta = ASTS-high has stronger dependency.",
            file=f
        )


    print(
        "\n============================================================"
    )

    print(
        "Step 17A completed successfully"
    )

    print(
        "============================================================"
    )

    print(
        summary_file
    )


if __name__ == "__main__":
    main()
