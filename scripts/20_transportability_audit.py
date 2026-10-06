#!/usr/bin/env python3

"""
ASTRA-Drug
Step 20: Transportability / influence audit

----





Primary functional endpoints
----------------------------
1. CRISPR EGFR dependency
      expected ASTS beta > 0

2. CRISPR MITO_CORE_UNION dependency
      expected ASTS beta < 0

3. PRISM BAY-87-2243
      expected ASTS beta < 0

4. PRISM oligomycin A
      expected ASTS beta < 0

5. EGFR_pathway drugs defined in Step12
      expected ASTS beta > 0


Base stringent model
--------------------
Outcome
    ~ ASTS
    + Lineage
    + Sex
    + proliferation
    + MAPK_driver
    + ERBB_driver
    + PI3K_driver
    + CellCycle
    + MYC


Leave-one-lineage-out (LOLO)
----------------------------



Sex interaction
---------------
Outcome
    ~ ASTS
    + Sex
    + ASTS:Sex
    + Lineage
    + proliferation
    + drivers
    + CellCycle
    + MYC



"""


# ============================================================
# 1. Path / Parameter
# ============================================================

from pathlib import Path


PROJECT_DIR = Path(__file__).resolve().parents[1]


# ------------------------------------------------------------
# DepMap / CRISPR
# ------------------------------------------------------------

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


MITO_SET_FILE = (
    PROJECT_DIR
    / "results/18_mito_dependency_specificity"
    / "18_mitochondrial_gene_sets.tsv"
)


# ------------------------------------------------------------
# Drug-response data
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
# Step11 / Step12 drug definitions
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
    / "results/20_transportability_audit"
)


FIG_DIR = (
    OUT_DIR
    / "figures"
)


# ============================================================
# 2. Analysis parameters
# ============================================================

ASTS_COLUMN = "ASTS_score"

PROLIFERATION_COLUMN = "proliferation_proxy"

LINEAGE_COLUMN = "OncotreeLineage"

SEX_COLUMN = "Sex"


MIN_GLOBAL_LINEAGE_N = 5


MIN_LINEAGE_N_FOR_LOLO = 30


MIN_N_ENDPOINT = 100


MIN_RESIDUAL_DF = 30


MIN_MITO_GENE_FRACTION = 0.80


# ------------------------------------------------------------
# Key mitochondrial drugs
# ------------------------------------------------------------

PRISM_KEY_MITO_DRUGS = {
    "BAY-87-2243":
        "bay872243",

    "oligomycin A":
        "oligomycina",
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
# 4. Generic helper functions
# ============================================================

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


def zscore(series):

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


def normalize_sex(x):
    """

    """

    if pd.isna(x):

        return "Unknown"

    value = str(
        x
    ).strip().lower()

    if value in {
        "female",
        "f",
    }:

        return "Female"

    if value in {
        "male",
        "m",
    }:

        return "Male"

    return "Unknown"


# ============================================================
# 5. Metadata
# ============================================================

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


def prepare_metadata():

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


    asts[
        "SexModel"
    ] = asts[
        SEX_COLUMN
    ].apply(
        normalize_sex
    )


    drivers = pd.read_csv(
        DRIVER_FILE,
        sep="\t"
    )[
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
    )[
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


    return metadata


# ============================================================
# 6. Base stringent model
# ============================================================

def fit_base_model(
    df,
    outcome_column="outcome"
):

    required = [
        outcome_column,
        ASTS_COLUMN,
        PROLIFERATION_COLUMN,
        "CellCycle_score",
        "MYC_V1_score",
        "MAPK_driver",
        "ERBB_driver",
        "PI3K_driver",
        "LineageModel",
    ]


    work = df.dropna(
        subset=required
    ).copy()


    if len(
        work
    ) < MIN_N_ENDPOINT:

        return None


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


    if len(
        work
    ) < MIN_N_ENDPOINT:

        return None


    outcome_sd = y.std(
        ddof=1
    )


    if (
        not np.isfinite(
            outcome_sd
        )
        or
        outcome_sd == 0
    ):

        return None


    parts = [
        work[
            [
                "ASTS_z",
                "Prolif_z",
                "CellCycle_z",
                "MYC_z",
            ]
        ].astype(float)
    ]


    # --------------------------------------------------------
    # Drivers
    # --------------------------------------------------------

    for col in [
        "MAPK_driver",
        "ERBB_driver",
        "PI3K_driver",
    ]:

        if (
            work[
                col
            ].nunique(
                dropna=True
            )
            >
            1
        ):

            parts.append(
                work[
                    [col]
                ].astype(float)
            )


    # --------------------------------------------------------
    # Sex
    # --------------------------------------------------------

    sex = (
        work[
            "SexModel"
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


    beta = float(
        model.params[
            "ASTS_z"
        ]
    )


    se = float(
        model.bse[
            "ASTS_z"
        ]
    )


    p = float(
        model.pvalues[
            "ASTS_z"
        ]
    )


    return {
        "n":
            len(
                work
            ),

        "n_lineages":
            work[
                "LineageModel"
            ].nunique(),

        "outcome_sd":
            outcome_sd,

        "beta":
            beta,

        "SE":
            se,

        "p":
            p,

        "standardized_beta":
            beta
            /
            outcome_sd,

        "standardized_SE":
            se
            /
            outcome_sd,
    }


# ============================================================
# 7. Sex-interaction model
# ============================================================

def fit_sex_interaction(
    df,
    outcome_column="outcome"
):

    required = [
        outcome_column,
        ASTS_COLUMN,
        PROLIFERATION_COLUMN,
        "CellCycle_score",
        "MYC_V1_score",
        "MAPK_driver",
        "ERBB_driver",
        "PI3K_driver",
        "LineageModel",
        "SexModel",
    ]


    work = df.dropna(
        subset=required
    ).copy()


    work = work[
        work[
            "SexModel"
        ].isin(
            [
                "Female",
                "Male",
            ]
        )
    ].copy()


    if (
        work[
            "SexModel"
        ].nunique()
        <
        2
    ):

        return None


    if len(
        work
    ) < MIN_N_ENDPOINT:

        return None


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


    work[
        "SexMale"
    ] = (
        work[
            "SexModel"
        ]
        ==
        "Male"
    ).astype(float)


    work[
        "ASTS_x_SexMale"
    ] = (
        work[
            "ASTS_z"
        ]
        *
        work[
            "SexMale"
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


    outcome_sd = y.std(
        ddof=1
    )


    if (
        not np.isfinite(
            outcome_sd
        )
        or
        outcome_sd == 0
    ):

        return None


    parts = [
        work[
            [
                "ASTS_z",
                "SexMale",
                "ASTS_x_SexMale",
                "Prolif_z",
                "CellCycle_z",
                "MYC_z",
            ]
        ].astype(float)
    ]


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


    # --------------------------------------------------------
    # Female ASTS effect
    #
    # Female is reference:
    # beta_ASTSZ
    # --------------------------------------------------------

    female_beta = float(
        model.params[
            "ASTS_z"
        ]
    )


    female_se = float(
        model.bse[
            "ASTS_z"
        ]
    )


    female_p = float(
        model.pvalues[
            "ASTS_z"
        ]
    )


    # --------------------------------------------------------
    # Male ASTS effect
    #
    # ASTS + ASTS:Male
    # --------------------------------------------------------

    parameter_names = list(
        model.params.index
    )


    contrast = np.zeros(
        len(
            parameter_names
        )
    )


    contrast[
        parameter_names.index(
            "ASTS_z"
        )
    ] = 1.0


    contrast[
        parameter_names.index(
            "ASTS_x_SexMale"
        )
    ] = 1.0


    male_test = model.t_test(
        contrast
    )


    male_beta = float(
        np.asarray(
            male_test.effect
        ).squeeze()
    )


    male_se = float(
        np.asarray(
            male_test.sd
        ).squeeze()
    )


    male_p = float(
        np.asarray(
            male_test.pvalue
        ).squeeze()
    )


    interaction_beta = float(
        model.params[
            "ASTS_x_SexMale"
        ]
    )


    interaction_se = float(
        model.bse[
            "ASTS_x_SexMale"
        ]
    )


    interaction_p = float(
        model.pvalues[
            "ASTS_x_SexMale"
        ]
    )


    return {
        "n":
            len(
                work
            ),

        "n_female":
            int(
                (
                    work[
                        "SexModel"
                    ]
                    ==
                    "Female"
                ).sum()
            ),

        "n_male":
            int(
                (
                    work[
                        "SexModel"
                    ]
                    ==
                    "Male"
                ).sum()
            ),

        "female_standardized_beta":
            female_beta
            /
            outcome_sd,

        "female_standardized_SE":
            female_se
            /
            outcome_sd,

        "female_p":
            female_p,

        "male_standardized_beta":
            male_beta
            /
            outcome_sd,

        "male_standardized_SE":
            male_se
            /
            outcome_sd,

        "male_p":
            male_p,

        "interaction_standardized_beta":
            interaction_beta
            /
            outcome_sd,

        "interaction_standardized_SE":
            interaction_se
            /
            outcome_sd,

        "interaction_p":
            interaction_p,
    }


# ============================================================
# 8. Leave-one-lineage-out
# ============================================================

def run_lolo(
    df,
    endpoint_name,
    expected_sign
):

    lineage_counts = (
        df[
            "LineageModel"
        ]
        .value_counts()
    )


    lineages = list(
        lineage_counts[
            lineage_counts
            >=
            MIN_LINEAGE_N_FOR_LOLO
        ].index
    )


    rows = []


    full = fit_base_model(
        df
    )


    if full is None:

        return (
            pd.DataFrame(),
            None
        )


    for lineage in lineages:

        subset = df[
            df[
                "LineageModel"
            ]
            !=
            lineage
        ].copy()


        result = fit_base_model(
            subset
        )


        if result is None:

            continue


        rows.append(
            {
                "endpoint":
                    endpoint_name,

                "excluded_lineage":
                    lineage,

                "excluded_n":
                    int(
                        lineage_counts[
                            lineage
                        ]
                    ),

                "n_remaining":
                    result[
                        "n"
                    ],

                "standardized_beta":
                    result[
                        "standardized_beta"
                    ],

                "p":
                    result[
                        "p"
                    ],

                "expected_sign":
                    expected_sign,

                "support_beta":
                    (
                        result[
                            "standardized_beta"
                        ]
                        *
                        expected_sign
                    ),

                "same_sign_as_full":
                    (
                        np.sign(
                            result[
                                "standardized_beta"
                            ]
                        )
                        ==
                        np.sign(
                            full[
                                "standardized_beta"
                            ]
                        )
                    ),
            }
        )


    lolo = pd.DataFrame(
        rows
    )


    if len(
        lolo
    ) == 0:

        return (
            lolo,
            None
        )


    full_beta = full[
        "standardized_beta"
    ]


    summary = {
        "endpoint":
            endpoint_name,

        "full_standardized_beta":
            full_beta,

        "expected_sign":
            expected_sign,

        "full_support_beta":
            full_beta
            *
            expected_sign,

        "n_lineages_tested":
            len(
                lolo
            ),

        "LOLO_min_beta":
            lolo[
                "standardized_beta"
            ].min(),

        "LOLO_max_beta":
            lolo[
                "standardized_beta"
            ].max(),

        "LOLO_median_beta":
            lolo[
                "standardized_beta"
            ].median(),

        "LOLO_sign_concordance":
            lolo[
                "same_sign_as_full"
            ].mean(),

        "LOLO_min_support_beta":
            lolo[
                "support_beta"
            ].min(),

        "all_LOLO_support_expected_direction":
            bool(
                (
                    lolo[
                        "support_beta"
                    ]
                    >
                    0
                ).all()
            ),
    }


    return (
        lolo,
        summary
    )


# ============================================================
# 9. CRISPR endpoints
# ============================================================

def prepare_crispr_endpoints(
    metadata
):

    print(
        "\nPreparing CRISPR endpoints..."
    )


    mito_table = pd.read_csv(
        MITO_SET_FILE,
        sep="\t"
    )


    mito_genes = set(
        mito_table.loc[
            mito_table[
                "gs_name"
            ]
            ==
            "MITO_CORE_UNION",
            "gene_symbol"
        ]
        .dropna()
        .astype(str)
    )


    header = pd.read_csv(
        CRISPR_FILE,
        nrows=0
    )


    columns = header.columns.tolist()


    model_id_column = None


    for candidate in [
        "ModelID",
        "DepMap_ID",
        "DepMapID",
    ]:

        if candidate in columns:

            model_id_column = candidate

            break


    if model_id_column is None:

        model_id_column = columns[
            0
        ]


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


    if "EGFR" not in symbol_to_column:

        raise RuntimeError(
            "EGFR not found in CRISPR matrix."
        )


    available_mito = [
        gene
        for gene in mito_genes
        if gene in symbol_to_column
    ]


    print(
        "MITO_CORE genes available:",
        len(
            available_mito
        ),
        "/",
        len(
            mito_genes
        )
    )


    use_columns = [
        model_id_column,
        symbol_to_column[
            "EGFR"
        ],
    ] + [
        symbol_to_column[
            gene
        ]
        for gene in available_mito
    ]


    crispr = pd.read_csv(
        CRISPR_FILE,
        usecols=use_columns,
        low_memory=False
    )


    crispr = crispr.rename(
        columns={
            model_id_column:
                "ModelID",

            symbol_to_column[
                "EGFR"
            ]:
                "EGFR_GE",
        }
    )


    # --------------------------------------------------------
    # Mito gene columns -> symbols
    # --------------------------------------------------------

    rename_map = {
        symbol_to_column[
            gene
        ]:
            gene
        for gene in available_mito
    }


    crispr = crispr.rename(
        columns=rename_map
    )


    data = crispr.merge(
        metadata,
        on="ModelID",
        how="inner",
        validate="1:1"
    )


    # --------------------------------------------------------
    # Standardize each mitochondrial gene effect
    # --------------------------------------------------------

    mito_matrix = data[
        available_mito
    ].apply(
        pd.to_numeric,
        errors="coerce"
    )


    mito_z = (
        mito_matrix
        -
        mito_matrix.mean(
            axis=0
        )
    ) / mito_matrix.std(
        axis=0,
        ddof=1
    )


    coverage = (
        mito_z.notna()
        .mean(
            axis=1
        )
    )


    mito_score = mito_z.mean(
        axis=1,
        skipna=True
    )


    mito_score[
        coverage
        <
        MIN_MITO_GENE_FRACTION
    ] = np.nan


    data[
        "MITO_CORE_GE_score"
    ] = mito_score


    # --------------------------------------------------------
    # EGFR endpoint
    # --------------------------------------------------------

    egfr = data.copy()

    egfr[
        "outcome"
    ] = pd.to_numeric(
        egfr[
            "EGFR_GE"
        ],
        errors="coerce"
    )


    # --------------------------------------------------------
    # Mito endpoint
    # --------------------------------------------------------

    mito = data.copy()

    mito[
        "outcome"
    ] = mito[
        "MITO_CORE_GE_score"
    ]


    return {
        "CRISPR_EGFR":
            {
                "data":
                    egfr,

                # positive = less dependency
                "expected_sign":
                    +1,

                "family":
                    "EGFR_function",

                "primary_key":
                    True,
            },

        "CRISPR_MITO_CORE":
            {
                "data":
                    mito,

                # negative = stronger dependency
                "expected_sign":
                    -1,

                "family":
                    "Mito_function",

                "primary_key":
                    True,
            },
    }


# ============================================================
# 10. PRISM
# ============================================================

def prepare_prism(
    metadata
):

    prism = pd.read_csv(
        PRISM_FILE,
        low_memory=False
    )


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
        ]
        >
        0
    ].copy()


    prism[
        "outcome"
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
            compound_name=(
                "name",
                "first"
            ),

            outcome=(
                "outcome",
                "median"
            ),
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


    return prism


# ============================================================
# 11. GDSC2
# ============================================================

def prepare_gdsc2(
    metadata
):

    gdsc = pd.read_excel(
        GDSC2_FILE,
        engine="openpyxl"
    )


    gdsc[
        "outcome"
    ] = pd.to_numeric(
        gdsc[
            "LN_IC50"
        ],
        errors="coerce"
    )


    gdsc = gdsc[
        gdsc[
            "outcome"
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
        "SANGER_MODEL_ID_norm"
    ] = (
        gdsc[
            "SANGER_MODEL_ID"
        ]
        .astype(str)
        .str.strip()
    )


    gdsc = (
        gdsc
        .groupby(
            [
                "DRUG_ID_norm",
                "SANGER_MODEL_ID_norm",
            ],
            as_index=False
        )
        .agg(
            compound_name=(
                "DRUG_NAME",
                "first"
            ),

            outcome=(
                "outcome",
                "median"
            ),
        )
    )


    model = pd.read_csv(
        MODEL_FILE,
        low_memory=False
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


    valid_sanger = set(
        counts[
            counts
            ==
            1
        ].index
    )


    model_map = (
        model_map[
            model_map[
                "SangerModelID"
            ].isin(
                valid_sanger
            )
        ]
        .drop_duplicates(
            "SangerModelID"
        )
    )


    gdsc = gdsc.merge(
        model_map,
        left_on="SANGER_MODEL_ID_norm",
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


    return gdsc


# ============================================================
# 12. Drug endpoints
# ============================================================

def prepare_drug_endpoints(
    metadata
):

    print(
        "\nPreparing drug-response endpoints..."
    )


    prism = prepare_prism(
        metadata
    )


    gdsc = prepare_gdsc2(
        metadata
    )


    overlap = pd.read_csv(
        OVERLAP_FILE,
        sep="\t"
    )


    classes = pd.read_csv(
        CLASS_FILE,
        sep="\t"
    )


    overlap[
        "DRUG_ID_norm"
    ] = overlap[
        "DRUG_ID"
    ].apply(
        normalize_drug_id
    )


    if "drug_key" not in overlap.columns:

        overlap[
            "drug_key"
        ] = overlap[
            "compound_name"
        ].apply(
            normalize_drug_name
        )


    drug_def = overlap.merge(
        classes,
        on="compound_name",
        how="left",
        validate="1:1",
        suffixes=(
            "",
            "_class"
        )
    )


    endpoints = {}


    # --------------------------------------------------------
    # Key mitochondrial PRISM drugs
    # --------------------------------------------------------

    for name, key in (
        PRISM_KEY_MITO_DRUGS.items()
    ):

        subset = prism[
            prism[
                "drug_key"
            ]
            ==
            key
        ].copy()


        if len(
            subset
        ) >= MIN_N_ENDPOINT:

            endpoint_name = (
                "PRISM_MITO_"
                +
                normalize_drug_name(
                    name
                )
            )


            endpoints[
                endpoint_name
            ] = {
                "data":
                    subset,

                # negative = sensitivity
                "expected_sign":
                    -1,

                "family":
                    "Mito_drug",

                "primary_key":
                    True,
            }


    # --------------------------------------------------------
    # EGFR pathway drugs
    # --------------------------------------------------------

    if (
        "EGFR_pathway"
        not in drug_def.columns
    ):

        raise RuntimeError(
            "EGFR_pathway column not found "
            "in Step12 drug class table."
        )


    egfr_drugs = drug_def[
        drug_def[
            "EGFR_pathway"
        ]
        ==
        1
    ].copy()


    print(
        "EGFR pathway drugs:",
        len(
            egfr_drugs
        )
    )


    for _, drug in (
        egfr_drugs.iterrows()
    ):

        compound_name = drug[
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


        if len(
            p
        ) >= MIN_N_ENDPOINT:

            endpoint_name = (
                "PRISM_EGFR_"
                +
                normalize_drug_name(
                    compound_name
                )
            )


            endpoints[
                endpoint_name
            ] = {
                "data":
                    p,

                # positive = resistance
                "expected_sign":
                    +1,

                "family":
                    "EGFR_drug",

                "primary_key":
                    False,

                "compound_name":
                    compound_name,
            }


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


        if len(
            g
        ) >= MIN_N_ENDPOINT:

            endpoint_name = (
                "GDSC2_EGFR_"
                +
                normalize_drug_name(
                    compound_name
                )
            )


            endpoints[
                endpoint_name
            ] = {
                "data":
                    g,

                "expected_sign":
                    +1,

                "family":
                    "EGFR_drug",

                "primary_key":
                    False,

                "compound_name":
                    compound_name,
            }


    return endpoints


# ============================================================
# 13. Main
# ============================================================

def main():

    print(
        "\n"
        "============================================================\n"
        "ASTRA-Drug Step 20\n"
        "Transportability / influence audit\n"
        "============================================================"
    )


    # --------------------------------------------------------
    # Input check
    # --------------------------------------------------------

    for path in [
        CRISPR_FILE,
        ASTS_FILE,
        DRIVER_FILE,
        CELL_CYCLE_FILE,
        MITO_SET_FILE,
        PRISM_FILE,
        GDSC2_FILE,
        MODEL_FILE,
        OVERLAP_FILE,
        CLASS_FILE,
    ]:

        if not path.exists():

            raise FileNotFoundError(
                f"Required file not found:\n"
                f"{path}"
            )


    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )


    FIG_DIR.mkdir(
        parents=True,
        exist_ok=True
    )


    # --------------------------------------------------------
    # Freeze analysis plan
    # --------------------------------------------------------

    with open(
        OUT_DIR
        / "20_ANALYSIS_PLAN_FREEZE.txt",
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 20",
            file=f
        )

        print(
            "Transportability / influence audit",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "No new candidate discovery.",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Audits:",
            file=f
        )

        print(
            "1. Leave-one-lineage-out",
            file=f
        )

        print(
            "2. ASTS x Sex interaction",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Primary anchors:",
            file=f
        )

        print(
            "CRISPR EGFR",
            file=f
        )

        print(
            "CRISPR MITO_CORE_UNION",
            file=f
        )

        print(
            "PRISM BAY-87-2243",
            file=f
        )

        print(
            "PRISM oligomycin A",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Secondary:",
            file=f
        )

        print(
            "Step12 predefined EGFR_pathway drugs",
            file=f
        )


    # ========================================================
    # Metadata
    # ========================================================

    metadata = prepare_metadata()


    # ========================================================
    # Endpoints
    # ========================================================

    endpoints = {}


    endpoints.update(
        prepare_crispr_endpoints(
            metadata
        )
    )


    endpoints.update(
        prepare_drug_endpoints(
            metadata
        )
    )


    print(
        "\nNumber of endpoints:",
        len(
            endpoints
        )
    )


    # ========================================================
    # Base models / LOLO / sex interaction
    # ========================================================

    base_rows = []

    sex_rows = []

    lolo_frames = []

    lolo_summary_rows = []


    for endpoint_name, info in (
        endpoints.items()
    ):

        print(
            "\nAnalyzing:",
            endpoint_name
        )


        data = info[
            "data"
        ]


        expected_sign = info[
            "expected_sign"
        ]


        family = info[
            "family"
        ]


        primary_key = info.get(
            "primary_key",
            False
        )


        # ----------------------------------------------------
        # Base
        # ----------------------------------------------------

        base = fit_base_model(
            data
        )


        if base is None:

            print(
                "  Base model skipped."
            )

            continue


        base_rows.append(
            {
                "endpoint":
                    endpoint_name,

                "family":
                    family,

                "primary_key":
                    primary_key,

                "expected_sign":
                    expected_sign,

                **base,

                "support_beta":
                    (
                        base[
                            "standardized_beta"
                        ]
                        *
                        expected_sign
                    ),

                "supports_expected_direction":
                    (
                        base[
                            "standardized_beta"
                        ]
                        *
                        expected_sign
                        >
                        0
                    ),
            }
        )


        # ----------------------------------------------------
        # Sex interaction
        # ----------------------------------------------------

        sex_result = (
            fit_sex_interaction(
                data
            )
        )


        if sex_result is not None:

            sex_rows.append(
                {
                    "endpoint":
                        endpoint_name,

                    "family":
                        family,

                    "primary_key":
                        primary_key,

                    "expected_sign":
                        expected_sign,

                    **sex_result,

                    "female_support_beta":
                        (
                            sex_result[
                                "female_standardized_beta"
                            ]
                            *
                            expected_sign
                        ),

                    "male_support_beta":
                        (
                            sex_result[
                                "male_standardized_beta"
                            ]
                            *
                            expected_sign
                        ),
                }
            )


        # ----------------------------------------------------
        # LOLO
        # ----------------------------------------------------

        lolo, lolo_summary = (
            run_lolo(
                data,
                endpoint_name,
                expected_sign
            )
        )


        if len(
            lolo
        ) > 0:

            lolo[
                "family"
            ] = family

            lolo_frames.append(
                lolo
            )


        if (
            lolo_summary
            is not None
        ):

            lolo_summary[
                "family"
            ] = family

            lolo_summary[
                "primary_key"
            ] = primary_key

            lolo_summary_rows.append(
                lolo_summary
            )


    # ========================================================
    # Save base models
    # ========================================================

    base_df = pd.DataFrame(
        base_rows
    )


    base_df.to_csv(
        OUT_DIR
        / "20_endpoint_base_models.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # Sex interaction
    # ========================================================

    sex_df = pd.DataFrame(
        sex_rows
    )


    if len(
        sex_df
    ) > 0:

        sex_df[
            "interaction_FDR"
        ] = np.nan


        finite = np.isfinite(
            sex_df[
                "interaction_p"
            ]
        )


        sex_df.loc[
            finite,
            "interaction_FDR"
        ] = multipletests(
            sex_df.loc[
                finite,
                "interaction_p"
            ],
            method="fdr_bh"
        )[1]


    sex_df.to_csv(
        OUT_DIR
        / "20_sex_interaction.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # LOLO outputs
    # ========================================================

    if len(
        lolo_frames
    ) > 0:

        lolo_df = pd.concat(
            lolo_frames,
            ignore_index=True
        )

    else:

        lolo_df = pd.DataFrame()


    lolo_df.to_csv(
        OUT_DIR
        / "20_leave_one_lineage_out.tsv",
        sep="\t",
        index=False
    )


    lolo_summary_df = pd.DataFrame(
        lolo_summary_rows
    )


    lolo_summary_df.to_csv(
        OUT_DIR
        / "20_lineage_robustness_summary.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # EGFR drug-class transportability summary
    # ========================================================

    egfr_base = base_df[
        base_df[
            "family"
        ]
        ==
        "EGFR_drug"
    ].copy()


    egfr_lolo = (
        lolo_summary_df[
            lolo_summary_df[
                "family"
            ]
            ==
            "EGFR_drug"
        ].copy()
        if len(
            lolo_summary_df
        ) > 0
        else pd.DataFrame()
    )


    class_rows = []


    for screen_prefix in [
        "PRISM_EGFR_",
        "GDSC2_EGFR_",
    ]:

        sub = egfr_base[
            egfr_base[
                "endpoint"
            ].str.startswith(
                screen_prefix
            )
        ]


        if len(
            sub
        ) == 0:

            continue


        endpoints_screen = set(
            sub[
                "endpoint"
            ]
        )


        lolo_sub = egfr_lolo[
            egfr_lolo[
                "endpoint"
            ].isin(
                endpoints_screen
            )
        ]


        class_rows.append(
            {
                "screen":
                    (
                        "PRISM"
                        if screen_prefix.startswith(
                            "PRISM"
                        )
                        else
                        "GDSC2"
                    ),

                "n_drugs":
                    len(
                        sub
                    ),

                "n_expected_direction":
                    int(
                        sub[
                            "supports_expected_direction"
                        ].sum()
                    ),

                "mean_support_beta":
                    sub[
                        "support_beta"
                    ].mean(),

                "min_support_beta":
                    sub[
                        "support_beta"
                    ].min(),

                "n_drugs_LOLO_all_expected_direction":
                    (
                        int(
                            lolo_sub[
                                "all_LOLO_support_expected_direction"
                            ].sum()
                        )
                        if len(
                            lolo_sub
                        ) > 0
                        else
                        0
                    ),
            }
        )


    class_summary = pd.DataFrame(
        class_rows
    )


    class_summary.to_csv(
        OUT_DIR
        / "20_EGFR_class_transportability.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # Figure 1:
    # LOLO for primary key endpoints
    # ========================================================

    key_summary = lolo_summary_df[
        lolo_summary_df[
            "primary_key"
        ]
        ==
        True
    ].copy()


    if len(
        key_summary
    ) > 0:

        key_summary = key_summary.sort_values(
            "full_support_beta"
        )


        fig, ax = plt.subplots(
            figsize=(8, 5)
        )


        y = np.arange(
            len(
                key_summary
            )
        )


        full_support = (
            key_summary[
                "full_support_beta"
            ].to_numpy()
        )


        # Convert LOLO beta range into support-direction range
        support_min = []

        support_max = []


        for _, row in (
            key_summary.iterrows()
        ):

            sign = row[
                "expected_sign"
            ]


            values = np.array(
                [
                    row[
                        "LOLO_min_beta"
                    ]
                    *
                    sign,

                    row[
                        "LOLO_max_beta"
                    ]
                    *
                    sign,
                ]
            )


            support_min.append(
                values.min()
            )

            support_max.append(
                values.max()
            )


        support_min = np.array(
            support_min
        )

        support_max = np.array(
            support_max
        )


        xerr = np.vstack(
            [
                full_support
                -
                support_min,

                support_max
                -
                full_support,
            ]
        )


        xerr = np.maximum(
            xerr,
            0
        )


        ax.errorbar(
            full_support,
            y,
            xerr=xerr,
            fmt="o",
            capsize=4
        )


        ax.axvline(
            0,
            linestyle="--"
        )


        ax.set_yticks(
            y
        )


        ax.set_yticklabels(
            key_summary[
                "endpoint"
            ]
        )


        ax.set_xlabel(
            "Direction-aligned standardized ASTS beta\n"
            "positive = supports predefined phenotype"
        )


        ax.set_title(
            "Leave-one-lineage-out robustness"
        )


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / "20_01_primary_endpoints_LOLO.png",
            dpi=200
        )


        plt.close(
            fig
        )


    # ========================================================
    # Figure 2:
    # Female / Male effects
    # ========================================================

    key_sex = sex_df[
        sex_df[
            "primary_key"
        ]
        ==
        True
    ].copy()


    if len(
        key_sex
    ) > 0:

        key_sex = key_sex.reset_index(
            drop=True
        )


        fig, ax = plt.subplots(
            figsize=(8, 5)
        )


        y = np.arange(
            len(
                key_sex
            )
        )


        offset = 0.12


        female_beta = (
            key_sex[
                "female_standardized_beta"
            ]
            *
            key_sex[
                "expected_sign"
            ]
        )


        male_beta = (
            key_sex[
                "male_standardized_beta"
            ]
            *
            key_sex[
                "expected_sign"
            ]
        )


        female_se = key_sex[
            "female_standardized_SE"
        ]


        male_se = key_sex[
            "male_standardized_SE"
        ]


        ax.errorbar(
            female_beta,
            y - offset,
            xerr=1.96
            *
            female_se,
            fmt="o",
            capsize=3,
            label="Female"
        )


        ax.errorbar(
            male_beta,
            y + offset,
            xerr=1.96
            *
            male_se,
            fmt="s",
            capsize=3,
            label="Male"
        )


        ax.axvline(
            0,
            linestyle="--"
        )


        ax.set_yticks(
            y
        )


        ax.set_yticklabels(
            key_sex[
                "endpoint"
            ]
        )


        ax.set_xlabel(
            "Direction-aligned standardized ASTS beta\n"
            "positive = supports predefined phenotype"
        )


        ax.set_title(
            "ASTS effects by cell-line donor sex"
        )


        ax.legend()


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / "20_02_primary_endpoints_by_sex.png",
            dpi=200
        )


        plt.close(
            fig
        )


    # ========================================================
    # Summary
    # ========================================================

    summary_file = (
        OUT_DIR
        / "20_summary.txt"
    )


    with open(
        summary_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 20 summary",
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
            "Final transportability / influence audit "
            "for predefined ASTS functional phenotypes.",
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Base endpoint models",
            file=f
        )

        print(
            "--------------------",
            file=f
        )


        print(
            base_df[
                [
                    "endpoint",
                    "family",
                    "n",
                    "standardized_beta",
                    "p",
                    "expected_sign",
                    "support_beta",
                    "supports_expected_direction",
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
            "Leave-one-lineage-out",
            file=f
        )

        print(
            "---------------------",
            file=f
        )


        print(
            lolo_summary_df.to_string(
                index=False
            ),
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Sex interaction",
            file=f
        )

        print(
            "---------------",
            file=f
        )


        if len(
            sex_df
        ) > 0:

            print(
                sex_df[
                    [
                        "endpoint",
                        "n",
                        "n_female",
                        "n_male",
                        "female_standardized_beta",
                        "male_standardized_beta",
                        "interaction_standardized_beta",
                        "interaction_p",
                        "interaction_FDR",
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
            "EGFR pathway drug-class transportability",
            file=f
        )

        print(
            "---------------------------------------",
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
            "Interpretation guide",
            file=f
        )

        print(
            "--------------------",
            file=f
        )


        print(
            "LOLO_min_support_beta > 0 means that "
            "the predefined phenotype remains in the expected "
            "direction after removal of every sufficiently large lineage.",
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Sex-interaction FDR < 0.05 indicates evidence that "
            "the ASTS association differs between Female- and "
            "Male-derived cell lines.",
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "A non-significant interaction does not prove identical "
            "effects between sexes; it indicates that this dataset "
            "does not provide evidence for effect modification.",
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "DepMap Sex reflects donor sex and should not be "
            "interpreted as current androgen exposure.",
            file=f
        )


    print(
        "\n============================================================"
    )

    print(
        "Step 20 completed successfully"
    )

    print(
        "============================================================"
    )

    print(
        summary_file
    )


if __name__ == "__main__":
    main()
