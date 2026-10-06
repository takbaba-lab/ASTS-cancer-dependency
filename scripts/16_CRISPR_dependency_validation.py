#!/usr/bin/env python3

"""
ASTRA-Drug
Step 16: Orthogonal CRISPR dependency validation

----

    ASTS-high
        -> MEK inhibitor resistance
        -> EGFR pathway inhibitor resistance


------------------------



    ASTS coefficient > 0




Primary direct targets
----------------------
MEK:
    MAP2K1
    MAP2K2

EGFR/ERBB:
    EGFR
    ERBB2

Secondary genes
---------------
MAPK signaling:
    BRAF, RAF1, ARAF,
    KRAS, NRAS, HRAS,
    SOS1, SHOC2, PTPN11, GRB2

ERBB family:
    ERBB3, ERBB4

Exploratory non-EGFR RTK:
    MET, FGFR1-4, PDGFRA/B, KDR, KIT, AXL

Nested models
-------------
Model 1:
    GeneEffect
      ~ ASTS
      + lineage
      + sex
      + proliferation proxy

Model 2:
    Model 1
      + MAPK_driver
      + ERBB_driver
      + PI3K_driver

Model 3 (primary stringent sensitivity):
    Model 2
      + CellCycle_score
      + MYC_V1_score

----


"""


# ============================================================
# 1. Path / Parameter
# ============================================================

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]


# ------------------------------------------------------------
# DepMap 26Q1 CRISPR gene effect
# ------------------------------------------------------------

CRISPR_FILE = (
    PROJECT_DIR
    / "data/raw/depmap/26Q1"
    / "CRISPRGeneEffect.csv"
)


# ------------------------------------------------------------
# ASTS
# ------------------------------------------------------------

ASTS_FILE = (
    PROJECT_DIR
    / "results/08b_DepMap_ASTS_default_profile"
    / "08_DepMap_ASTS_scores.tsv"
)


# ------------------------------------------------------------
# Step 13 driver flags
# ------------------------------------------------------------

DRIVER_FILE = (
    PROJECT_DIR
    / "results/13_driver_confounder"
    / "13_driver_flags.tsv"
)


# ------------------------------------------------------------
# Step 15 cell-cycle scores
# ------------------------------------------------------------

CELL_CYCLE_FILE = (
    PROJECT_DIR
    / "results/15_cellcycle_confounder"
    / "15_cellcycle_scores.tsv"
)


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

OUT_DIR = (
    PROJECT_DIR
    / "results/16_CRISPR_dependency"
)

FIG_DIR = (
    OUT_DIR
    / "figures"
)


# ------------------------------------------------------------
# Metadata columns
# ------------------------------------------------------------

ASTS_COLUMN = "ASTS_score"

PROLIFERATION_COLUMN = "proliferation_proxy"

LINEAGE_COLUMN = "OncotreeLineage"

SEX_COLUMN = "Sex"


# ------------------------------------------------------------
# Analysis parameters
# ------------------------------------------------------------

MIN_GLOBAL_LINEAGE_N = 5

MIN_MODELS_PER_GENE = 500

MIN_RESIDUAL_DF = 30

FDR_ALPHA = 0.05


# ============================================================
# 2. Pre-specified target genes
# ============================================================

# Direct targets: strongest prior hypothesis
PRIMARY_TARGET_GENES = [
    "MAP2K1",
    "MAP2K2",
    "EGFR",
    "ERBB2",
]


# Secondary MAPK / ERBB signaling genes
SECONDARY_SIGNALING_GENES = [
    "BRAF",
    "RAF1",
    "ARAF",
    "KRAS",
    "NRAS",
    "HRAS",
    "SOS1",
    "SHOC2",
    "PTPN11",
    "GRB2",
    "ERBB3",
    "ERBB4",
]


# Exploratory non-EGFR RTK genes
EXPLORATORY_RTK_GENES = [
    "MET",
    "FGFR1",
    "FGFR2",
    "FGFR3",
    "FGFR4",
    "PDGFRA",
    "PDGFRB",
    "KDR",
    "KIT",
    "AXL",
]


ALL_TARGET_GENES = (
    PRIMARY_TARGET_GENES
    + SECONDARY_SIGNALING_GENES
    + EXPLORATORY_RTK_GENES
)


# ------------------------------------------------------------
# Predefined composite axes
#
#
# ------------------------------------------------------------

COMPOSITE_AXES = {

    "MEK_direct": [
        "MAP2K1",
        "MAP2K2",
    ],

    "ERBB_direct": [
        "EGFR",
        "ERBB2",
    ],

    "ERBB_family": [
        "EGFR",
        "ERBB2",
        "ERBB3",
        "ERBB4",
    ],
}


# ============================================================
# 3. Imports
# ============================================================

import re

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import statsmodels.api as sm

from statsmodels.stats.multitest import multipletests


# ============================================================
# 4. Helper functions
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
    DepMap gene column

        MAP2K1 (5604)


        MAP2K1

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

    first = columns[0]

    print(
        "WARNING: standard ModelID column not found."
    )

    print(
        "Using first column as ModelID:",
        first
    )

    return first


def build_lineage_categories(
    metadata
):
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


def build_design_matrix(
    work,
    model_type
):
    """
    """

    parts = [
        work[
            [
                "ASTS_z",
                "Prolif_z",
            ]
        ].astype(float)
    ]


    # --------------------------------------------------------
    # Driver adjustment
    # --------------------------------------------------------

    if model_type in {
        "driver_adjusted",
        "stringent",
    }:

        for col in [
            "MAPK_driver",
            "ERBB_driver",
            "PI3K_driver",
        ]:

            if (
                col in work.columns
                and
                work[
                    col
                ].nunique(
                    dropna=True
                )
                > 1
            ):

                parts.append(
                    work[
                        [col]
                    ].astype(float)
                )


    # --------------------------------------------------------
    # Cell-cycle / MYC adjustment
    # --------------------------------------------------------

    if model_type == "stringent":

        parts.append(
            work[
                [
                    "CellCycle_z",
                    "MYC_z",
                ]
            ].astype(float)
        )


    # --------------------------------------------------------
    # Sex
    # --------------------------------------------------------

    sex = (
        work[
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


    # --------------------------------------------------------
    # Lineage
    # --------------------------------------------------------

    lineage = (
        work[
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


    X = sm.add_constant(
        X,
        has_constant="add"
    )


    return X


def fit_nested_models(
    data,
    outcome_column
):
    """

    """

    required = [
        outcome_column,
        ASTS_COLUMN,
        PROLIFERATION_COLUMN,
        "LineageModel",
        "MAPK_driver",
        "ERBB_driver",
        "PI3K_driver",
        "CellCycle_score",
        "MYC_V1_score",
    ]


    work = data.dropna(
        subset=required
    ).copy()


    if len(work) < MIN_MODELS_PER_GENE:

        return None


    # --------------------------------------------------------
    # Standardized predictors
    # --------------------------------------------------------

    work[
        "ASTS_z"
    ] = zscore(
        work[
            ASTS_COLUMN
        ]
    )


    work[
        "Prolif_z"
    ] = zscore(
        work[
            PROLIFERATION_COLUMN
        ]
    )


    work[
        "CellCycle_z"
    ] = zscore(
        work[
            "CellCycle_score"
        ]
    )


    work[
        "MYC_z"
    ] = zscore(
        work[
            "MYC_V1_score"
        ]
    )


    work = work.dropna(
        subset=[
            "ASTS_z",
            "Prolif_z",
            "CellCycle_z",
            "MYC_z",
        ]
    )


    # --------------------------------------------------------
    # Outcome
    # --------------------------------------------------------

    y = pd.to_numeric(
        work[
            outcome_column
        ],
        errors="coerce"
    )


    valid_y = y.notna()

    work = work.loc[
        valid_y
    ].copy()

    y = y.loc[
        valid_y
    ]


    if len(work) < MIN_MODELS_PER_GENE:

        return None


    outcome_sd = y.std(
        ddof=1
    )


    if (
        not np.isfinite(outcome_sd)
        or
        outcome_sd == 0
    ):

        return None


    result = {
        "n_models":
            len(work),

        "gene_effect_mean":
            y.mean(),

        "gene_effect_sd":
            outcome_sd,
    }


    # --------------------------------------------------------
    # Fit nested models
    # --------------------------------------------------------

    model_types = [
        "base",
        "driver_adjusted",
        "stringent",
    ]


    for model_type in model_types:

        X = build_design_matrix(
            work,
            model_type
        )


        model = sm.OLS(
            y,
            X
        ).fit(
            cov_type="HC3"
        )


        if (
            model.df_resid
            <
            MIN_RESIDUAL_DF
        ):

            return None


        beta_raw = float(
            model.params[
                "ASTS_z"
            ]
        )


        se_raw = float(
            model.bse[
                "ASTS_z"
            ]
        )


        p_value = float(
            model.pvalues[
                "ASTS_z"
            ]
        )


        result[
            f"beta_{model_type}"
        ] = beta_raw


        result[
            f"SE_{model_type}"
        ] = se_raw


        result[
            f"p_{model_type}"
        ] = p_value


        result[
            f"standardized_beta_{model_type}"
        ] = (
            beta_raw
            /
            outcome_sd
        )


        result[
            f"standardized_SE_{model_type}"
        ] = (
            se_raw
            /
            outcome_sd
        )


    # --------------------------------------------------------
    # Attenuation
    # --------------------------------------------------------

    base_abs = abs(
        result[
            "beta_base"
        ]
    )


    if base_abs > 0:

        result[
            "attenuation_driver_percent"
        ] = (
            1
            -
            abs(
                result[
                    "beta_driver_adjusted"
                ]
            )
            /
            base_abs
        ) * 100


        result[
            "attenuation_stringent_percent"
        ] = (
            1
            -
            abs(
                result[
                    "beta_stringent"
                ]
            )
            /
            base_abs
        ) * 100


    else:

        result[
            "attenuation_driver_percent"
        ] = np.nan

        result[
            "attenuation_stringent_percent"
        ] = np.nan


    return result


# ============================================================
# 5. Main
# ============================================================

def main():

    print(
        "\n"
        "============================================================\n"
        "ASTRA-Drug Step 16\n"
        "CRISPR dependency validation\n"
        "============================================================"
    )


    # ========================================================
    # Input checks
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


    FIG_DIR.mkdir(
        parents=True,
        exist_ok=True
    )


    # ========================================================
    # Step 1:
    # Inspect CRISPR header only
    # ========================================================

    print(
        "\nStep 1: Inspecting CRISPRGeneEffect header"
    )


    header = pd.read_csv(
        CRISPR_FILE,
        nrows=0
    )


    columns = header.columns.tolist()


    model_id_column = (
        detect_model_id_column(
            columns
        )
    )


    print(
        "ModelID column:",
        model_id_column
    )


    # --------------------------------------------------------
    # Gene symbol -> CRISPR column mapping
    # --------------------------------------------------------

    symbol_to_column = {}


    for col in columns:

        if col == model_id_column:
            continue


        symbol = extract_gene_symbol(
            col
        )


        if symbol not in symbol_to_column:

            symbol_to_column[
                symbol
            ] = col


    available_genes = [
        gene
        for gene
        in ALL_TARGET_GENES
        if gene in symbol_to_column
    ]


    missing_genes = [
        gene
        for gene
        in ALL_TARGET_GENES
        if gene not in symbol_to_column
    ]


    print(
        "Requested genes:",
        len(
            ALL_TARGET_GENES
        )
    )


    print(
        "Available genes:",
        len(
            available_genes
        )
    )


    print(
        "Missing genes:",
        missing_genes
    )


    if (
        not all(
            gene
            in available_genes
            for gene
            in PRIMARY_TARGET_GENES
        )
    ):

        missing_primary = [
            gene
            for gene
            in PRIMARY_TARGET_GENES
            if gene
            not in available_genes
        ]


        raise RuntimeError(
            "Primary CRISPR target genes missing: "
            +
            ", ".join(
                missing_primary
            )
        )


    # ========================================================
    # Step 2:
    # Read only target columns
    # ========================================================

    print(
        "\nStep 2: Reading target CRISPR gene effects"
    )


    use_columns = [
        model_id_column
    ] + [
        symbol_to_column[
            gene
        ]
        for gene
        in available_genes
    ]


    crispr = pd.read_csv(
        CRISPR_FILE,
        usecols=use_columns,
        low_memory=False
    )


    crispr = crispr.rename(
        columns={
            model_id_column:
                "ModelID"
        }
    )


    # Rename gene columns to symbols
    rename_map = {
        symbol_to_column[
            gene
        ]:
            gene
        for gene
        in available_genes
    }


    crispr = crispr.rename(
        columns=rename_map
    )


    if crispr[
        "ModelID"
    ].duplicated().any():

        raise RuntimeError(
            "Duplicated ModelID in CRISPRGeneEffect.csv."
        )


    print(
        "CRISPR models:",
        len(
            crispr
        )
    )


    # ========================================================
    # Step 3:
    # Metadata
    # ========================================================

    print(
        "\nStep 3: Preparing ASTS metadata"
    )


    asts = pd.read_csv(
        ASTS_FILE,
        sep="\t"
    )


    if asts[
        "ModelID"
    ].duplicated().any():

        raise RuntimeError(
            "Duplicated ModelID in ASTS table."
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


    data = crispr.merge(
        metadata,
        on="ModelID",
        how="inner",
        validate="1:1"
    )


    print(
        "CRISPR x ASTS models:",
        len(
            data
        )
    )


    # ========================================================
    # Step 4:
    # Gene-level models
    # ========================================================

    print(
        "\nStep 4: Gene-level dependency models"
    )


    rows = []


    for gene in available_genes:

        print(
            "Analyzing:",
            gene
        )


        result = fit_nested_models(
            data,
            gene
        )


        if result is None:

            print(
                "  skipped"
            )

            continue


        if gene in PRIMARY_TARGET_GENES:

            category = (
                "primary_direct_target"
            )


        elif gene in SECONDARY_SIGNALING_GENES:

            category = (
                "secondary_signaling"
            )


        else:

            category = (
                "exploratory_RTK"
            )


        rows.append(
            {
                "gene":
                    gene,

                "category":
                    category,

                **result,
            }
        )


    gene_results = pd.DataFrame(
        rows
    )


    # --------------------------------------------------------
    # FDR across predefined gene panel
    # Primary interpretation uses stringent model.
    # --------------------------------------------------------

    finite = np.isfinite(
        gene_results[
            "p_stringent"
        ]
    )


    gene_results[
        "FDR_stringent"
    ] = np.nan


    gene_results.loc[
        finite,
        "FDR_stringent"
    ] = multipletests(
        gene_results.loc[
            finite,
            "p_stringent"
        ],
        method="fdr_bh"
    )[1]


    # --------------------------------------------------------
    # Direction label
    # --------------------------------------------------------

    gene_results[
        "ASTS_dependency_direction"
    ] = np.where(
        gene_results[
            "beta_stringent"
        ] > 0,
        "ASTS_high_less_dependent",
        "ASTS_high_more_dependent"
    )


    gene_results = gene_results.sort_values(
        "standardized_beta_stringent",
        ascending=False
    )


    gene_results.to_csv(
        OUT_DIR
        / "16_target_gene_dependency_results.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # Step 5:
    # Composite dependency axes
    # ========================================================

    print(
        "\nStep 5: Composite dependency axes"
    )


    composite_data = data.copy()


    # --------------------------------------------------------
    # Standardize each gene effect across models
    # --------------------------------------------------------

    for gene in available_genes:

        composite_data[
            f"{gene}_GE_z"
        ] = zscore(
            composite_data[
                gene
            ]
        )


    axis_rows = []


    for axis_name, genes in COMPOSITE_AXES.items():

        axis_genes = [
            gene
            for gene
            in genes
            if gene in available_genes
        ]


        if len(axis_genes) < 2:

            print(
                "Skipping axis:",
                axis_name
            )

            continue


        axis_column = (
            f"{axis_name}_score"
        )


        z_columns = [
            f"{gene}_GE_z"
            for gene
            in axis_genes
        ]


        composite_data[
            axis_column
        ] = (
            composite_data[
                z_columns
            ]
            .mean(
                axis=1,
                skipna=False
            )
        )


        result = fit_nested_models(
            composite_data,
            axis_column
        )


        if result is None:

            continue


        axis_rows.append(
            {
                "axis":
                    axis_name,

                "genes":
                    ";".join(
                        axis_genes
                    ),

                **result,
            }
        )


    axis_results = pd.DataFrame(
        axis_rows
    )


    if len(axis_results) > 0:

        finite = np.isfinite(
            axis_results[
                "p_stringent"
            ]
        )


        axis_results[
            "FDR_stringent"
        ] = np.nan


        axis_results.loc[
            finite,
            "FDR_stringent"
        ] = multipletests(
            axis_results.loc[
                finite,
                "p_stringent"
            ],
            method="fdr_bh"
        )[1]


        axis_results[
            "interpretation"
        ] = np.where(
            axis_results[
                "beta_stringent"
            ] > 0,
            "ASTS_high_less_dependent",
            "ASTS_high_more_dependent"
        )


    axis_results.to_csv(
        OUT_DIR
        / "16_composite_dependency_axes.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # Step 6:
    # Figure - primary/direct targets
    # ========================================================

    direct = gene_results[
        gene_results[
            "gene"
        ].isin(
            PRIMARY_TARGET_GENES
        )
    ].copy()


    direct = direct.sort_values(
        "standardized_beta_stringent"
    )


    if len(direct) > 0:

        fig, ax = plt.subplots(
            figsize=(7, 5)
        )


        y = np.arange(
            len(direct)
        )


        beta = direct[
            "standardized_beta_stringent"
        ].to_numpy()


        se = direct[
            "standardized_SE_stringent"
        ].to_numpy()


        ax.errorbar(
            beta,
            y,
            xerr=1.96 * se,
            fmt="o",
            capsize=3
        )


        ax.axvline(
            0,
            linestyle="--"
        )


        ax.set_yticks(
            y
        )


        ax.set_yticklabels(
            direct[
                "gene"
            ]
        )


        ax.set_xlabel(
            "ASTS beta on CRISPR gene effect "
            "(SD units)\n"
            "positive = ASTS-high less dependent"
        )


        ax.set_title(
            "Direct MEK / EGFR target dependencies"
        )


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / "16_01_direct_target_CRISPR_betas.png",
            dpi=200
        )


        plt.close(
            fig
        )


    # ========================================================
    # Step 7:
    # Figure - all target genes
    # ========================================================

    plot_data = gene_results.sort_values(
        "standardized_beta_stringent"
    )


    fig, ax = plt.subplots(
        figsize=(8, 10)
    )


    y = np.arange(
        len(plot_data)
    )


    ax.barh(
        y,
        plot_data[
            "standardized_beta_stringent"
        ]
    )


    ax.axvline(
        0,
        linestyle="--"
    )


    ax.set_yticks(
        y
    )


    ax.set_yticklabels(
        plot_data[
            "gene"
        ]
    )


    ax.set_xlabel(
        "Stringent-adjusted ASTS beta on CRISPR gene effect\n"
        "positive = ASTS-high less dependent"
    )


    ax.set_title(
        "ASTS-associated genetic dependencies"
    )


    fig.tight_layout()


    fig.savefig(
        FIG_DIR
        / "16_02_all_target_CRISPR_betas.png",
        dpi=200
    )


    plt.close(
        fig
    )


    # ========================================================
    # Step 8:
    # Nested-model comparison for direct targets
    # ========================================================

    if len(direct) > 0:

        direct_long_rows = []


        for _, row in direct.iterrows():

            for model_name in [
                "base",
                "driver_adjusted",
                "stringent",
            ]:

                direct_long_rows.append(
                    {
                        "gene":
                            row[
                                "gene"
                            ],

                        "model":
                            model_name,

                        "beta":
                            row[
                                f"standardized_beta_{model_name}"
                            ],
                    }
                )


        direct_long = pd.DataFrame(
            direct_long_rows
        )


        genes = [
            gene
            for gene
            in PRIMARY_TARGET_GENES
            if gene
            in set(
                direct_long[
                    "gene"
                ]
            )
        ]


        x = np.arange(
            len(genes)
        )


        width = 0.24


        fig, ax = plt.subplots(
            figsize=(8, 5)
        )


        for index, model_name in enumerate(
            [
                "base",
                "driver_adjusted",
                "stringent",
            ]
        ):

            values = []


            for gene in genes:

                value = direct_long.loc[
                    (
                        direct_long[
                            "gene"
                        ]
                        ==
                        gene
                    )
                    &
                    (
                        direct_long[
                            "model"
                        ]
                        ==
                        model_name
                    ),
                    "beta"
                ].iloc[0]


                values.append(
                    value
                )


            ax.bar(
                x
                +
                (
                    index - 1
                )
                *
                width,
                values,
                width,
                label=model_name
            )


        ax.axhline(
            0,
            linestyle="--"
        )


        ax.set_xticks(
            x
        )


        ax.set_xticklabels(
            genes
        )


        ax.set_ylabel(
            "Standardized ASTS beta on CRISPR gene effect\n"
            "positive = less dependency"
        )


        ax.set_title(
            "Direct target dependency across adjustment models"
        )


        ax.legend()


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / "16_03_direct_targets_nested_models.png",
            dpi=200
        )


        plt.close(
            fig
        )


    # ========================================================
    # Step 9:
    # Summary
    # ========================================================

    summary_file = (
        OUT_DIR
        / "16_summary.txt"
    )


    with open(
        summary_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 16 summary",
            file=f
        )

        print(
            "==========================",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "Purpose:",
            file=f
        )

        print(
            "Orthogonal functional validation using "
            "DepMap CRISPR gene-effect scores.",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            f"CRISPR models: "
            f"{len(crispr)}",
            file=f
        )


        print(
            f"CRISPR x ASTS models: "
            f"{len(data)}",
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
            "Positive ASTS beta = ASTS-high has higher "
            "gene-effect score = weaker genetic dependency.",
            file=f
        )

        print(
            "Negative ASTS beta = ASTS-high has stronger "
            "genetic dependency.",
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Primary direct targets",
            file=f
        )

        print(
            "----------------------",
            file=f
        )


        direct_output_columns = [
            "gene",
            "n_models",
            "standardized_beta_base",
            "p_base",
            "standardized_beta_driver_adjusted",
            "p_driver_adjusted",
            "standardized_beta_stringent",
            "p_stringent",
            "FDR_stringent",
            "ASTS_dependency_direction",
        ]


        print(
            gene_results.loc[
                gene_results[
                    "gene"
                ].isin(
                    PRIMARY_TARGET_GENES
                ),
                direct_output_columns
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
            "Composite dependency axes",
            file=f
        )

        print(
            "-------------------------",
            file=f
        )


        if len(axis_results) > 0:

            print(
                axis_results[
                    [
                        "axis",
                        "genes",
                        "n_models",
                        "standardized_beta_stringent",
                        "p_stringent",
                        "FDR_stringent",
                        "interpretation",
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
            "All predefined genes - stringent model",
            file=f
        )

        print(
            "--------------------------------------",
            file=f
        )


        print(
            gene_results[
                [
                    "gene",
                    "category",
                    "n_models",
                    "standardized_beta_stringent",
                    "p_stringent",
                    "FDR_stringent",
                    "ASTS_dependency_direction",
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
            "Interpretation note",
            file=f
        )

        print(
            "-------------------",
            file=f
        )


        print(
            "Concordant positive associations for MAP2K1/2 "
            "or EGFR/ERBB2 would support reduced genetic "
            "dependency in ASTS-high cells.",
            file=f
        )


        print(
            "A null CRISPR association does not invalidate "
            "drug-response results because gene knockout and "
            "pharmacological inhibition are not equivalent.",
            file=f
        )


        print(
            "This analysis is orthogonal validation and "
            "does not alter the ASTS signature.",
            file=f
        )


    print(
        "\n============================================================"
    )

    print(
        "Step 16 completed successfully"
    )

    print(
        "============================================================"
    )

    print(
        "\nSummary:"
    )

    print(
        summary_file
    )


if __name__ == "__main__":
    main()
