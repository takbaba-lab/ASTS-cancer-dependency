#!/usr/bin/env python3

"""
ASTRA-Drug
Step 15B: Proliferation / cell-cycle confounder analysis

----



----

    Drug response
      ~ ASTS
      + lineage
      + sex
      + original proliferation proxy
      + MAPK / ERBB / PI3K driver

    Model 1
      + independent E2F/G2M cell-cycle score

    Model 2
      + independent MYC_TARGETS_V1 score






"""


# ============================================================
# Path / Parameter
# ============================================================

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]


# ------------------------------------------------------------
# DepMap expression
# ------------------------------------------------------------

EXPRESSION_FILE = (
    PROJECT_DIR
    / "data/raw/depmap/26Q1"
    / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
)


# ------------------------------------------------------------
# ASTS / signature
# ------------------------------------------------------------

ASTS_FILE = (
    PROJECT_DIR
    / "results/08b_DepMap_ASTS_default_profile"
    / "08_DepMap_ASTS_scores.tsv"
)


SIGNATURE_FILE = (
    PROJECT_DIR
    / "results/08b_DepMap_ASTS_default_profile"
    / "08_signature_expression_mapping.tsv"
)


# ------------------------------------------------------------
# Hallmark genes from Step15A
# ------------------------------------------------------------

HALLMARK_FILE = (
    PROJECT_DIR
    / "results/15_cellcycle_confounder"
    / "15_hallmark_cellcycle_genes.tsv"
)


# ------------------------------------------------------------
# Step13 drivers
# ------------------------------------------------------------

DRIVER_FILE = (
    PROJECT_DIR
    / "results/13_driver_confounder"
    / "13_driver_flags.tsv"
)


# ------------------------------------------------------------
# PRISM / GDSC2
# ------------------------------------------------------------

PRISM_FILE = (
    PROJECT_DIR
    / "data/raw/prism/secondary"
    / "secondary-screen-dose-response-curve-parameters.csv"
)


GDSC2_FILE = (
    PROJECT_DIR
    / "data/raw/gdsc/GDSC2"
    / "GDSC2_fitted_dose_response_27Oct23.xlsx"
)


MODEL_FILE = (
    PROJECT_DIR
    / "data/raw/depmap/26Q1"
    / "Model.csv"
)


# ------------------------------------------------------------
# Step11 / Step12
# ------------------------------------------------------------

OVERLAP_FILE = (
    PROJECT_DIR
    / "results/11_cross_screen_mechanism"
    / "11_PRISM_GDSC2_overlap_annotated.tsv"
)


CLASS_FILE = (
    PROJECT_DIR
    / "results/12_drug_class_robustness"
    / "12_drug_class_membership.tsv"
)


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

OUT_DIR = (
    PROJECT_DIR
    / "results/15_cellcycle_confounder"
)

FIG_DIR = (
    OUT_DIR
    / "figures"
)


# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

DEFAULT_PROFILE_COLUMN = "IsDefaultEntryForModel"

ASTS_COLUMN = "ASTS_score"

PROLIFERATION_COLUMN = "proliferation_proxy"

LINEAGE_COLUMN = "OncotreeLineage"

SEX_COLUMN = "Sex"


MIN_GLOBAL_LINEAGE_N = 5

MIN_N_PER_DRUG = 100

MIN_RESIDUAL_DF = 30

SCORE_ROW_CHUNK = 100


# ------------------------------------------------------------
# Existing proliferation proxy genes
# ------------------------------------------------------------

PROLIFERATION_GENES = {
    "MKI67",
    "PCNA",
    "MCM2",
    "MCM3",
    "MCM4",
    "MCM5",
    "MCM6",
    "MCM7",
    "TOP2A",
    "CDK1",
}


# ============================================================
# Imports
# ============================================================

import re

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import statsmodels.api as sm

from scipy.stats import spearmanr


# ============================================================
# Helper functions
# ============================================================

def zscore(series):

    x = pd.to_numeric(
        series,
        errors="coerce"
    )

    sd = x.std(ddof=1)

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


def normalize_drug_name(x):

    if pd.isna(x):
        return ""

    return re.sub(
        r"[^a-z0-9]+",
        "",
        str(x).lower().strip()
    )


def normalize_drug_id(x):

    if pd.isna(x):
        return ""

    x = str(x).strip()

    if re.fullmatch(
        r"\d+\.0",
        x
    ):
        x = x[:-2]

    return x


def extract_gene_symbol(column_name):

    x = str(
        column_name
    )

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

    return x.strip()


def detect_signature_column(df):

    preferred = [
        "human_gene_symbol",
        "human_symbol",
        "gene_symbol",
        "GeneSymbol",
        "symbol",
        "human_gene",
        "Human_Symbol",
    ]

    for col in preferred:

        if col in df.columns:

            return col


    candidates = [
        col
        for col in df.columns
        if (
            "human"
            in col.lower()
            and
            (
                "gene"
                in col.lower()
                or
                "symbol"
                in col.lower()
            )
        )
    ]


    if len(candidates) == 1:

        return candidates[0]


    raise ValueError(
        "Could not detect human gene-symbol column.\n"
        f"Columns={df.columns.tolist()}"
    )


def build_lineage_categories(
    metadata
):

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


# ============================================================
# Expression-based rank scores
# ============================================================

def calculate_rank_scores(
    expression,
    gene_columns,
    symbol_to_column,
    score_gene_sets,
):

    """


    """

    scores = {
        name: []
        for name
        in score_gene_sets
    }


    model_ids = []


    n = len(
        expression
    )


    for start in range(
        0,
        n,
        SCORE_ROW_CHUNK
    ):

        end = min(
            start
            +
            SCORE_ROW_CHUNK,
            n
        )


        print(
            f"Scoring models "
            f"{start + 1}-{end}/{n}"
        )


        chunk = (
            expression
            .iloc[
                start:end
            ]
        )


        values = (
            chunk[
                gene_columns
            ]
        )


        ranks = values.rank(
            axis=1,
            pct=True,
            method="average"
        )


        model_ids.extend(
            chunk[
                "ModelID"
            ].tolist()
        )


        for score_name, genes in score_gene_sets.items():

            cols = [
                symbol_to_column[g]
                for g in genes
                if g in symbol_to_column
            ]


            if len(cols) == 0:

                raise RuntimeError(
                    f"No genes available for "
                    f"{score_name}"
                )


            score = (
                ranks[
                    cols
                ]
                .mean(
                    axis=1
                )
            )


            scores[
                score_name
            ].extend(
                score.tolist()
            )


    result = pd.DataFrame(
        {
            "ModelID":
                model_ids
        }
    )


    for name, values in scores.items():

        result[
            name
        ] = values


    return result


# ============================================================
# Cell-state association
# ============================================================

def fit_state_model(
    df,
    outcome_column,
):

    """
    Score ~ ASTS + lineage + sex + original proliferation + drivers

    """

    work = df.dropna(
        subset=[
            outcome_column,
            ASTS_COLUMN,
            PROLIFERATION_COLUMN,
            "LineageModel",
            "MAPK_driver",
            "ERBB_driver",
            "PI3K_driver",
        ]
    ).copy()


    work[
        "Outcome_z"
    ] = zscore(
        work[
            outcome_column
        ]
    )


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


    parts = [
        work[
            [
                "ASTS_z",
                "Prolif_z",
                "MAPK_driver",
                "ERBB_driver",
                "PI3K_driver",
            ]
        ].astype(float)
    ]


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


    lineage = (
        work[
            "LineageModel"
        ]
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


    model = sm.OLS(
        work[
            "Outcome_z"
        ],
        X
    ).fit(
        cov_type="HC3"
    )


    return {
        "outcome":
            outcome_column,

        "n":
            len(work),

        "beta_ASTS":
            model.params[
                "ASTS_z"
            ],

        "SE":
            model.bse[
                "ASTS_z"
            ],

        "p":
            model.pvalues[
                "ASTS_z"
            ],

        "R2":
            model.rsquared,
    }


# ============================================================
# Drug-response models
# ============================================================

def build_design_matrix(
    work,
    model_type,
):

    parts = [
        work[
            [
                "ASTS_z_model",
                "Prolif_z_model",
            ]
        ]
    ]


    # drivers
    for col in [
        "MAPK_driver",
        "ERBB_driver",
        "PI3K_driver",
    ]:

        if (
            work[
                col
            ].nunique()
            >
            1
        ):

            parts.append(
                work[
                    [col]
                ].astype(float)
            )


    if model_type in {
        "plus_cellcycle",
        "plus_cellcycle_myc",
    }:

        parts.append(
            work[
                [
                    "CellCycle_z_model"
                ]
            ]
        )


    if (
        model_type
        ==
        "plus_cellcycle_myc"
    ):

        parts.append(
            work[
                [
                    "MYC_z_model"
                ]
            ]
        )


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


    lineage = (
        work[
            "LineageModel"
        ]
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


    return sm.add_constant(
        X,
        has_constant="add"
    )


def fit_nested_drug_models(
    df
):

    """

    """

    required = [
        "response",
        ASTS_COLUMN,
        PROLIFERATION_COLUMN,
        "CellCycle_score",
        "MYC_V1_score",
        "LineageModel",
        "MAPK_driver",
        "ERBB_driver",
        "PI3K_driver",
    ]


    work = df.dropna(
        subset=required
    ).copy()


    if (
        len(work)
        <
        MIN_N_PER_DRUG
    ):

        return None


    work[
        "ASTS_z_model"
    ] = zscore(
        work[
            ASTS_COLUMN
        ]
    )


    work[
        "Prolif_z_model"
    ] = zscore(
        work[
            PROLIFERATION_COLUMN
        ]
    )


    work[
        "CellCycle_z_model"
    ] = zscore(
        work[
            "CellCycle_score"
        ]
    )


    work[
        "MYC_z_model"
    ] = zscore(
        work[
            "MYC_V1_score"
        ]
    )


    work = work.dropna(
        subset=[
            "ASTS_z_model",
            "Prolif_z_model",
            "CellCycle_z_model",
            "MYC_z_model",
        ]
    )


    results = {
        "n":
            len(work)
    }


    for model_type in [
        "driver_only",
        "plus_cellcycle",
        "plus_cellcycle_myc",
    ]:

        X = build_design_matrix(
            work,
            model_type
        )


        model = sm.OLS(
            work[
                "response"
            ],
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


        results[
            f"beta_{model_type}"
        ] = model.params[
            "ASTS_z_model"
        ]


        results[
            f"p_{model_type}"
        ] = model.pvalues[
            "ASTS_z_model"
        ]


    # attenuation relative to driver-only baseline
    baseline = abs(
        results[
            "beta_driver_only"
        ]
    )


    if baseline > 0:

        results[
            "attenuation_cellcycle_percent"
        ] = (
            1
            -
            abs(
                results[
                    "beta_plus_cellcycle"
                ]
            )
            /
            baseline
        ) * 100


        results[
            "attenuation_cellcycle_myc_percent"
        ] = (
            1
            -
            abs(
                results[
                    "beta_plus_cellcycle_myc"
                ]
            )
            /
            baseline
        ) * 100


    else:

        results[
            "attenuation_cellcycle_percent"
        ] = np.nan

        results[
            "attenuation_cellcycle_myc_percent"
        ] = np.nan


    return results


# ============================================================
# Main
# ============================================================

def main():

    print(
        "\n"
        "============================================================\n"
        "ASTRA-Drug Step 15\n"
        "Proliferation / cell-cycle confounder audit\n"
        "============================================================"
    )


    for path in [
        EXPRESSION_FILE,
        ASTS_FILE,
        SIGNATURE_FILE,
        HALLMARK_FILE,
        DRIVER_FILE,
        PRISM_FILE,
        GDSC2_FILE,
        MODEL_FILE,
        OVERLAP_FILE,
        CLASS_FILE,
    ]:

        if not path.exists():

            raise FileNotFoundError(
                f"File not found:\n{path}"
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
    # 1. Signature genes
    # ========================================================

    signature = pd.read_csv(
        SIGNATURE_FILE,
        sep="\t"
    )


    signature_col = (
        detect_signature_column(
            signature
        )
    )


    asts_genes = set(
        signature[
            signature_col
        ]
        .dropna()
        .astype(str)
        .str.strip()
    )


    print(
        "ASTS genes:",
        len(asts_genes)
    )


    # ========================================================
    # 2. Hallmark sets
    # ========================================================

    hallmark = pd.read_csv(
        HALLMARK_FILE,
        sep="\t"
    )


    hallmark_sets = {
        name: set(
            hallmark.loc[
                hallmark[
                    "gs_name"
                ] == name,
                "gene_symbol"
            ]
            .dropna()
            .astype(str)
        )
        for name in hallmark[
            "gs_name"
        ].unique()
    }


    exclusion_genes = (
        asts_genes
        |
        PROLIFERATION_GENES
    )


    e2f_genes = (
        hallmark_sets[
            "HALLMARK_E2F_TARGETS"
        ]
        -
        exclusion_genes
    )


    g2m_genes = (
        hallmark_sets[
            "HALLMARK_G2M_CHECKPOINT"
        ]
        -
        exclusion_genes
    )


    myc_genes = (
        hallmark_sets[
            "HALLMARK_MYC_TARGETS_V1"
        ]
        -
        exclusion_genes
    )


    cellcycle_genes = (
        e2f_genes
        |
        g2m_genes
    )


    print(
        "E2F genes after exclusions:",
        len(e2f_genes)
    )

    print(
        "G2M genes after exclusions:",
        len(g2m_genes)
    )

    print(
        "Cell-cycle union:",
        len(cellcycle_genes)
    )

    print(
        "MYC V1 genes after exclusions:",
        len(myc_genes)
    )


    # ========================================================
    # 3. Expression
    # ========================================================

    print(
        "\nReading expression..."
    )


    expression = pd.read_csv(
        EXPRESSION_FILE,
        low_memory=False
    )


    if (
        DEFAULT_PROFILE_COLUMN
        in expression.columns
    ):

        flag = (
            expression[
                DEFAULT_PROFILE_COLUMN
            ]
            .astype(str)
            .str.strip()
            .str.lower()
        )


        expression = expression[
            flag.isin(
                [
                    "yes",
                    "true",
                    "1",
                ]
            )
        ].copy()


    if (
        expression[
            "ModelID"
        ].duplicated().any()
    ):

        raise RuntimeError(
            "Duplicated ModelID remains "
            "in expression table."
        )


    # only gene columns
    gene_columns = [
        col
        for col in expression.columns
        if re.match(
            r"^.+\s+\(\d+\)$",
            str(col)
        )
    ]


    symbol_to_column = {}


    for col in gene_columns:

        symbol = (
            extract_gene_symbol(
                col
            )
        )

        if (
            symbol
            not in symbol_to_column
        ):

            symbol_to_column[
                symbol
            ] = col


    # ========================================================
    # 4. Rank scores
    # ========================================================

    score_gene_sets = {
        "E2F_score":
            e2f_genes,

        "G2M_score":
            g2m_genes,

        "CellCycle_score":
            cellcycle_genes,

        "MYC_V1_score":
            myc_genes,
    }


    scores = calculate_rank_scores(
        expression,
        gene_columns,
        symbol_to_column,
        score_gene_sets
    )


    scores.to_csv(
        OUT_DIR
        / "15_cellcycle_scores.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # 5. Metadata / drivers
    # ========================================================

    asts = pd.read_csv(
        ASTS_FILE,
        sep="\t"
    )


    if (
        asts[
            "ModelID"
        ].duplicated().any()
    ):

        raise RuntimeError(
            "Duplicated ASTS ModelID."
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


    metadata = (
        asts
        .merge(
            scores,
            on="ModelID",
            how="inner",
            validate="1:1"
        )
        .merge(
            drivers,
            on="ModelID",
            how="left",
            validate="1:1"
        )
    )


    # ========================================================
    # 6. ASTS vs proliferation / cell-cycle
    # ========================================================

    correlation_rows = []


    for variable in [
        PROLIFERATION_COLUMN,
        "E2F_score",
        "G2M_score",
        "CellCycle_score",
        "MYC_V1_score",
    ]:

        sub = metadata[
            [
                ASTS_COLUMN,
                variable,
            ]
        ].dropna()


        rho, p = spearmanr(
            sub[
                ASTS_COLUMN
            ],
            sub[
                variable
            ]
        )


        correlation_rows.append(
            {
                "variable":
                    variable,

                "n":
                    len(sub),

                "Spearman_r":
                    rho,

                "p":
                    p,
            }
        )


    correlation_df = pd.DataFrame(
        correlation_rows
    )


    correlation_df.to_csv(
        OUT_DIR
        / "15_ASTS_cellcycle_correlations.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # 7. Partial state associations
    # ========================================================

    state_rows = []


    for outcome in [
        "E2F_score",
        "G2M_score",
        "CellCycle_score",
        "MYC_V1_score",
    ]:

        state_rows.append(
            fit_state_model(
                metadata,
                outcome
            )
        )


    state_df = pd.DataFrame(
        state_rows
    )


    state_df.to_csv(
        OUT_DIR
        / "15_ASTS_cellcycle_partial_models.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # 8. Drug definitions
    # ========================================================

    overlap = pd.read_csv(
        OVERLAP_FILE,
        sep="\t"
    )


    overlap[
        "DRUG_ID_norm"
    ] = overlap[
        "DRUG_ID"
    ].apply(
        normalize_drug_id
    )


    classes = pd.read_csv(
        CLASS_FILE,
        sep="\t"
    )


    # ========================================================
    # 9. PRISM
    # ========================================================

    print(
        "\nPreparing PRISM..."
    )


    prism = pd.read_csv(
        PRISM_FILE,
        low_memory=False
    )


    prism = prism[
        [
            "depmap_id",
            "auc",
            "name",
        ]
    ].copy()


    prism[
        "auc"
    ] = pd.to_numeric(
        prism[
            "auc"
        ],
        errors="coerce"
    )


    prism = prism[
        prism[
            "auc"
        ] > 0
    ].copy()


    prism[
        "response"
    ] = np.log2(
        prism[
            "auc"
        ]
    )


    prism[
        "drug_key"
    ] = prism[
        "name"
    ].apply(
        normalize_drug_name
    )


    prism = (
        prism
        .groupby(
            [
                "drug_key",
                "depmap_id",
            ],
            as_index=False
        )
        .agg(
            response=(
                "response",
                "median"
            )
        )
    )


    prism = prism.rename(
        columns={
            "depmap_id":
                "ModelID"
        }
    )


    prism = prism.merge(
        metadata,
        on="ModelID",
        how="inner",
        validate="m:1"
    )


    # ========================================================
    # 10. GDSC2
    # ========================================================

    print(
        "Preparing GDSC2..."
    )


    gdsc = pd.read_excel(
        GDSC2_FILE,
        engine="openpyxl"
    )


    gdsc = gdsc[
        [
            "SANGER_MODEL_ID",
            "DRUG_ID",
            "LN_IC50",
        ]
    ].copy()


    gdsc[
        "response"
    ] = pd.to_numeric(
        gdsc[
            "LN_IC50"
        ],
        errors="coerce"
    )


    gdsc = gdsc[
        gdsc[
            "response"
        ].notna()
    ].copy()


    gdsc[
        "DRUG_ID_norm"
    ] = gdsc[
        "DRUG_ID"
    ].apply(
        normalize_drug_id
    )


    gdsc[
        "SANGER_MODEL_ID"
    ] = (
        gdsc[
            "SANGER_MODEL_ID"
        ]
        .astype(str)
        .str.strip()
    )


    model = pd.read_csv(
        MODEL_FILE
    )


    model_map = model[
        [
            "ModelID",
            "SangerModelID",
        ]
    ].dropna().copy()


    model_map[
        "SangerModelID"
    ] = (
        model_map[
            "SangerModelID"
        ]
        .astype(str)
        .str.strip()
    )


    counts = (
        model_map
        .groupby(
            "SangerModelID"
        )[
            "ModelID"
        ]
        .nunique()
    )


    good = set(
        counts[
            counts == 1
        ].index
    )


    model_map = (
        model_map[
            model_map[
                "SangerModelID"
            ].isin(
                good
            )
        ]
        .drop_duplicates(
            "SangerModelID"
        )
    )


    gdsc = (
        gdsc
        .groupby(
            [
                "DRUG_ID_norm",
                "SANGER_MODEL_ID",
            ],
            as_index=False
        )
        .agg(
            response=(
                "response",
                "median"
            )
        )
    )


    gdsc = gdsc.merge(
        model_map,
        left_on="SANGER_MODEL_ID",
        right_on="SangerModelID",
        how="inner",
        validate="m:1"
    )


    gdsc = gdsc.merge(
        metadata,
        on="ModelID",
        how="inner",
        validate="m:1"
    )


    # ========================================================
    # 11. Fit nested models
    # ========================================================

    rows = []


    for _, drug in overlap.iterrows():

        drug_name = drug[
            "compound_name"
        ]


        # PRISM
        p = prism[
            prism[
                "drug_key"
            ]
            ==
            drug[
                "drug_key"
            ]
        ].copy()


        result = fit_nested_drug_models(
            p
        )


        if result is not None:

            rows.append(
                {
                    "screen":
                        "PRISM",

                    "compound_name":
                        drug_name,

                    **result,
                }
            )


        # GDSC2
        g = gdsc[
            gdsc[
                "DRUG_ID_norm"
            ]
            ==
            drug[
                "DRUG_ID_norm"
            ]
        ].copy()


        result = fit_nested_drug_models(
            g
        )


        if result is not None:

            rows.append(
                {
                    "screen":
                        "GDSC2",

                    "compound_name":
                        drug_name,

                    **result,
                }
            )


    results = pd.DataFrame(
        rows
    )


    results.to_csv(
        OUT_DIR
        / "15_all_drugs_cellcycle_adjustment.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # 12. Cross-screen comparison
    # ========================================================

    prism_results = (
        results[
            results[
                "screen"
            ] == "PRISM"
        ]
        .drop(
            columns=[
                "screen"
            ]
        )
    )


    gdsc_results = (
        results[
            results[
                "screen"
            ] == "GDSC2"
        ]
        .drop(
            columns=[
                "screen"
            ]
        )
    )


    compare = prism_results.merge(
        gdsc_results,
        on="compound_name",
        suffixes=(
            "_PRISM",
            "_GDSC2"
        ),
        validate="1:1"
    )


    model_types = [
        "driver_only",
        "plus_cellcycle",
        "plus_cellcycle_myc",
    ]


    for model_type in model_types:

        pcol = (
            f"beta_{model_type}_PRISM"
        )

        gcol = (
            f"beta_{model_type}_GDSC2"
        )


        compare[
            f"PRISM_z_{model_type}"
        ] = zscore(
            compare[
                pcol
            ]
        )


        compare[
            f"GDSC2_z_{model_type}"
        ] = zscore(
            compare[
                gcol
            ]
        )


        compare[
            f"consensus_z_{model_type}"
        ] = (
            compare[
                f"PRISM_z_{model_type}"
            ]
            +
            compare[
                f"GDSC2_z_{model_type}"
            ]
        ) / 2


    compare = compare.merge(
        classes,
        on="compound_name",
        how="left",
        validate="1:1"
    )


    compare.to_csv(
        OUT_DIR
        / "15_PRISM_GDSC2_cellcycle_adjusted.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # 13. Class summary
    # ========================================================

    class_names = [
        col
        for col in [
            "MEK_target",
            "EGFR_target",
            "EGFR_pathway",
            "ERK_MAPK_pathway",
            "PI3K_MTOR_pathway",
            "RTK_pathway",
            "RTK_non_EGFR",
        ]
        if col in compare.columns
    ]


    class_rows = []


    for class_name in class_names:

        subset = compare[
            compare[
                class_name
            ] == 1
        ]


        for model_type in model_types:

            pcol = (
                f"beta_{model_type}_PRISM"
            )

            gcol = (
                f"beta_{model_type}_GDSC2"
            )


            resistant = (
                (
                    subset[
                        pcol
                    ] > 0
                )
                &
                (
                    subset[
                        gcol
                    ] > 0
                )
            )


            sensitive = (
                (
                    subset[
                        pcol
                    ] < 0
                )
                &
                (
                    subset[
                        gcol
                    ] < 0
                )
            )


            class_rows.append(
                {
                    "class":
                        class_name,

                    "model":
                        model_type,

                    "n_drugs":
                        len(subset),

                    "n_concordant_resistant":
                        int(
                            resistant.sum()
                        ),

                    "n_concordant_sensitive":
                        int(
                            sensitive.sum()
                        ),

                    "mean_consensus_z":
                        subset[
                            f"consensus_z_{model_type}"
                        ].mean(),

                    "min_consensus_z":
                        subset[
                            f"consensus_z_{model_type}"
                        ].min(),
                }
            )


    class_summary = pd.DataFrame(
        class_rows
    )


    class_summary.to_csv(
        OUT_DIR
        / "15_class_summary.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # 14. Global cross-screen summary
    # ========================================================

    global_rows = []


    for model_type in model_types:

        pcol = (
            f"beta_{model_type}_PRISM"
        )

        gcol = (
            f"beta_{model_type}_GDSC2"
        )


        rho, p = spearmanr(
            compare[
                pcol
            ],
            compare[
                gcol
            ]
        )


        sign_fraction = (
            np.sign(
                compare[
                    pcol
                ]
            )
            ==
            np.sign(
                compare[
                    gcol
                ]
            )
        ).mean()


        global_rows.append(
            {
                "model":
                    model_type,

                "n_drugs":
                    len(compare),

                "Spearman_r":
                    rho,

                "Spearman_p":
                    p,

                "sign_concordance":
                    sign_fraction,
            }
        )


    global_df = pd.DataFrame(
        global_rows
    )


    global_df.to_csv(
        OUT_DIR
        / "15_cross_screen_summary.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # 15. Figures
    # ========================================================

    for screen in [
        "PRISM",
        "GDSC2",
    ]:

        x = results[
            results[
                "screen"
            ] == screen
        ]


        fig, ax = plt.subplots(
            figsize=(6, 6)
        )


        ax.scatter(
            x[
                "beta_driver_only"
            ],
            x[
                "beta_plus_cellcycle"
            ],
            alpha=0.7
        )


        low = min(
            x[
                "beta_driver_only"
            ].min(),
            x[
                "beta_plus_cellcycle"
            ].min()
        )


        high = max(
            x[
                "beta_driver_only"
            ].max(),
            x[
                "beta_plus_cellcycle"
            ].max()
        )


        ax.plot(
            [
                low,
                high,
            ],
            [
                low,
                high,
            ],
            linestyle="--"
        )


        ax.axhline(
            0,
            linestyle=":"
        )


        ax.axvline(
            0,
            linestyle=":"
        )


        ax.set_xlabel(
            "Driver-adjusted ASTS beta"
        )


        ax.set_ylabel(
            "Driver + cell-cycle adjusted ASTS beta"
        )


        ax.set_title(
            f"{screen}: cell-cycle adjustment"
        )


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / f"15_{screen}_cellcycle_adjustment.png",
            dpi=200
        )


        plt.close(
            fig
        )


    # ========================================================
    # 16. Summary
    # ========================================================

    summary_file = (
        OUT_DIR
        / "15_summary.txt"
    )


    with open(
        summary_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 15 summary",
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
            "Test whether ASTS-drug associations are explained "
            "by proliferation/cell-cycle state.",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "Cell-cycle scores:",
            file=f
        )

        print(
            f"E2F genes after exclusions: "
            f"{len(e2f_genes)}",
            file=f
        )

        print(
            f"G2M genes after exclusions: "
            f"{len(g2m_genes)}",
            file=f
        )

        print(
            f"E2F/G2M union: "
            f"{len(cellcycle_genes)}",
            file=f
        )

        print(
            f"MYC V1 genes after exclusions: "
            f"{len(myc_genes)}",
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Raw ASTS correlations",
            file=f
        )

        print(
            "---------------------",
            file=f
        )

        print(
            correlation_df.to_string(
                index=False
            ),
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Partial cell-state models",
            file=f
        )

        print(
            "-------------------------",
            file=f
        )

        print(
            state_df.to_string(
                index=False
            ),
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Cross-screen drug effects",
            file=f
        )

        print(
            "-------------------------",
            file=f
        )

        print(
            global_df.to_string(
                index=False
            ),
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Drug-class summary",
            file=f
        )

        print(
            "------------------",
            file=f
        )

        print(
            class_summary.to_string(
                index=False
            ),
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Interpretation note:",
            file=f
        )

        print(
            "Cell-cycle adjustment is a stringent sensitivity "
            "analysis and may adjust away a biological mediator "
            "of the ASTS state.",
            file=f
        )

        print(
            "Attenuation therefore does not by itself demonstrate "
            "confounding or artifact.",
            file=f
        )


    print(
        "\n============================================================"
    )

    print(
        "Step 15 completed successfully"
    )

    print(
        "============================================================"
    )

    print(
        summary_file
    )


if __name__ == "__main__":
    main()
