#!/usr/bin/env python3

"""
ASTRA-Drug
Step 10: GDSC2 validation of PRISM-ASTS drug associations

----

Primary validation:

    disulfiram
    carboxyamidotriazole

Primary model:
    LN_IC50 ~ ASTS + Lineage + Sex + Proliferation


    beta_ASTS < 0
        higher ASTS -> greater sensitivity

    beta_ASTS > 0
        higher ASTS -> greater resistance


----
"""


# ============================================================
# Path / Parameter
# ============================================================

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]

ASTS_FILE = (
    PROJECT_DIR
    / "results/08b_DepMap_ASTS_default_profile"
    / "08_DepMap_ASTS_scores.tsv"
)

# DepMap 26Q1 model metadata
MODEL_FILE = (
    PROJECT_DIR
    / "data/raw/depmap/26Q1"
    / "Model.csv"
)

PRISM_RESULT_FILE = (
    PROJECT_DIR
    / "results/09b_PRISM_ASTS_default_profile"
    / "09_PRISM_ASTS_primary.tsv"
)

# GDSC2
GDSC_FILE = (
    PROJECT_DIR
    / "data/raw/gdsc/GDSC2"
    / "GDSC2_fitted_dose_response_27Oct23.xlsx"
)

OUT_DIR = (
    PROJECT_DIR
    / "results/10_GDSC2_validation"
)

FIG_DIR = OUT_DIR / "figures"


# ------------------------------------------------------------
# Analysis parameters
# ------------------------------------------------------------

PRIMARY_ASTS_COLUMN = "ASTS_score"
HIGHCONF_ASTS_COLUMN = "ASTS_high_conf_score"

LINEAGE_COLUMN = "OncotreeLineage"
SEX_COLUMN = "Sex"
PROLIFERATION_COLUMN = "proliferation_proxy"

# PRISM primary validation candidates
PRISM_PRIMARY_FDR_MAX = 0.10

EXPECTED_PRIMARY_CANDIDATE_N = 2

# Secondary validation
TOP_N_PRISM = 20

MIN_N_PER_DRUG = 100

MIN_GLOBAL_LINEAGE_N = 5

# within-lineage validation
MIN_WITHIN_LINEAGE_N = 15
MIN_LINEAGES_FOR_META = 3

MIN_RESIDUAL_DF = 30

EXPLORATORY_FDR_ALPHA = 0.05


# ============================================================
# Imports
# ============================================================

import hashlib
import math
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import statsmodels.api as sm

from scipy.stats import (
    spearmanr,
    norm,
)

from statsmodels.stats.multitest import multipletests


# ============================================================
# Helper functions
# ============================================================

def sha256_file(path, block_size=1024 * 1024):
    """Calculate the SHA256 checksum of a file."""

    h = hashlib.sha256()

    with open(path, "rb") as f:
        while True:
            block = f.read(block_size)

            if not block:
                break

            h.update(block)

    return h.hexdigest()


def normalize_drug_name(name):
    """


        YM-155 -> ym155
        YM 155 -> ym155

    """

    if pd.isna(name):
        return ""

    x = str(name).lower().strip()

    x = re.sub(
        r"[^a-z0-9]+",
        "",
        x
    )

    return x


def zscore(series):
    """Standardize values to mean 0 and standard deviation 1."""

    x = pd.to_numeric(
        series,
        errors="coerce"
    )

    sd = x.std(ddof=1)

    if (
        not np.isfinite(sd)
        or sd == 0
    ):
        return pd.Series(
            np.nan,
            index=x.index
        )

    return (
        x - x.mean()
    ) / sd


def bh_fdr(pvalues):
    """Benjamini-Hochberg FDR。"""

    p = np.asarray(
        pvalues,
        dtype=float
    )

    q = np.full(
        len(p),
        np.nan
    )

    ok = np.isfinite(p)

    if ok.sum() > 0:
        q[ok] = multipletests(
            p[ok],
            method="fdr_bh"
        )[1]

    return q


def fit_model(
    df,
    asts_column,
    include_lineage=True
):
    """

    """

    use = df.copy()

    required = [
        "LN_IC50",
        asts_column,
        PROLIFERATION_COLUMN,
    ]

    if include_lineage:
        required.append(
            "LineageModel"
        )

    use = use.dropna(
        subset=required
    ).copy()


    if len(use) < 5:
        return None


    use["ASTS_z_model"] = zscore(
        use[asts_column]
    )

    use["Prolif_z_model"] = zscore(
        use[PROLIFERATION_COLUMN]
    )


    use = use.dropna(
        subset=[
            "ASTS_z_model",
            "Prolif_z_model",
        ]
    ).copy()


    if len(use) < 5:
        return None


    X_parts = [
        use[
            [
                "ASTS_z_model",
                "Prolif_z_model",
            ]
        ]
    ]


    # Sex
    sex = (
        use[SEX_COLUMN]
        .fillna("Unknown")
        .astype(str)
    )

    if sex.nunique() > 1:

        sex_dummies = pd.get_dummies(
            sex,
            prefix="Sex",
            drop_first=True,
            dtype=float
        )

        X_parts.append(
            sex_dummies
        )


    # Lineage fixed effects
    if include_lineage:

        lineage = (
            use["LineageModel"]
            .fillna("Unknown")
            .astype(str)
        )

        if lineage.nunique() > 1:

            lineage_dummies = pd.get_dummies(
                lineage,
                prefix="Lineage",
                drop_first=True,
                dtype=float
            )

            X_parts.append(
                lineage_dummies
            )


    X = pd.concat(
        X_parts,
        axis=1
    )


    X = sm.add_constant(
        X,
        has_constant="add"
    )


    y = pd.to_numeric(
        use["LN_IC50"],
        errors="coerce"
    )


    try:

        model = sm.OLS(
            y,
            X
        ).fit(
            cov_type="HC3"
        )

    except Exception:
        return None


    if (
        "ASTS_z_model"
        not in model.params.index
    ):
        return None


    beta = float(
        model.params[
            "ASTS_z_model"
        ]
    )

    se = float(
        model.bse[
            "ASTS_z_model"
        ]
    )

    pvalue = float(
        model.pvalues[
            "ASTS_z_model"
        ]
    )

    zvalue = float(
        model.tvalues[
            "ASTS_z_model"
        ]
    )

    response_sd = y.std(ddof=1)


    if (
        np.isfinite(response_sd)
        and response_sd > 0
    ):
        standardized_beta = (
            beta / response_sd
        )
    else:
        standardized_beta = np.nan


    df_resid = float(
        model.df_resid
    )


    if df_resid > 0:
        partial_r2 = (
            zvalue ** 2
            /
            (
                zvalue ** 2
                + df_resid
            )
        )
    else:
        partial_r2 = np.nan


    return {
        "n":
            len(use),

        "n_lineages":
            use["LineageModel"].nunique()
            if include_lineage
            else 1,

        "beta_ASTS_LNIC50_per_1SD":
            beta,

        "se_ASTS":
            se,

        "z_ASTS":
            zvalue,

        "p_ASTS_two_sided":
            pvalue,

        "beta_ASTS_standardized":
            standardized_beta,

        "partial_R2_ASTS":
            partial_r2,

        "model_R2":
            float(model.rsquared),

        "model_adj_R2":
            float(model.rsquared_adj),

        "df_resid":
            df_resid,
    }


def one_sided_p_from_two_sided(
    beta,
    p_two_sided,
    expected_beta
):
    """
    """

    if (
        not np.isfinite(beta)
        or
        not np.isfinite(p_two_sided)
    ):
        return 1.0


    if (
        np.sign(beta)
        ==
        np.sign(expected_beta)
    ):
        return p_two_sided / 2.0

    return 1.0 - p_two_sided / 2.0


def fixed_effect_meta(df):
    """
    """

    x = df.copy()

    x = x[
        np.isfinite(
            x["beta_ASTS_LNIC50_per_1SD"]
        )
        &
        np.isfinite(
            x["se_ASTS"]
        )
        &
        (
            x["se_ASTS"] > 0
        )
    ]


    if len(x) < MIN_LINEAGES_FOR_META:
        return None


    beta = x[
        "beta_ASTS_LNIC50_per_1SD"
    ].to_numpy()

    se = x[
        "se_ASTS"
    ].to_numpy()


    weights = 1.0 / se ** 2


    meta_beta = (
        np.sum(weights * beta)
        /
        np.sum(weights)
    )


    meta_se = math.sqrt(
        1.0
        /
        np.sum(weights)
    )


    z = meta_beta / meta_se

    p = 2.0 * norm.sf(
        abs(z)
    )


    concordance = np.mean(
        np.sign(beta)
        ==
        np.sign(meta_beta)
    )


    return {
        "n_lineages_meta":
            len(x),

        "meta_beta":
            meta_beta,

        "meta_se":
            meta_se,

        "meta_z":
            z,

        "meta_p":
            p,

        "lineage_direction_concordance":
            concordance,
    }


# ============================================================
# Main
# ============================================================

def main():

    print(
        "\n"
        "============================================================\n"
        "ASTRA-Drug Step 10\n"
        "GDSC2 pharmacologic validation\n"
        "============================================================"
    )


    # --------------------------------------------------------
    # 0. Input check
    # --------------------------------------------------------

    for path in [
        ASTS_FILE,
        MODEL_FILE,
        PRISM_RESULT_FILE,
        GDSC_FILE,
    ]:

        if not path.exists():
            raise FileNotFoundError(
                f"File not found: {path}"
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
    # 1. Read PRISM results BEFORE reading GDSC response
    # --------------------------------------------------------

    print(
        "\nStep 1: Freezing PRISM validation targets"
    )


    prism = pd.read_csv(
        PRISM_RESULT_FILE,
        sep="\t"
    )


    prism = prism.sort_values(
        [
            "FDR_ASTS",
            "p_ASTS",
        ]
    ).reset_index(
        drop=True
    )


    prism["drug_key"] = (
        prism["compound_name"]
        .apply(
            normalize_drug_name
        )
    )


    primary_candidates = prism[
        prism["FDR_ASTS"]
        < PRISM_PRIMARY_FDR_MAX
    ].copy()


    if (
        len(primary_candidates)
        != EXPECTED_PRIMARY_CANDIDATE_N
    ):

        raise RuntimeError(
            "Expected "
            f"{EXPECTED_PRIMARY_CANDIDATE_N} "
            "primary PRISM candidates, but found "
            f"{len(primary_candidates)}. "
            "Check the Step09b input file."
        )


    secondary_candidates = (
        prism
        .head(
            TOP_N_PRISM
        )
        .copy()
    )


    print(
        "Primary candidates:"
    )

    print(
        primary_candidates[
            [
                "compound_name",
                "beta_ASTS_log2AUC_per_1SD",
                "FDR_ASTS",
            ]
        ].to_string(
            index=False
        )
    )


    # --------------------------------------------------------
    # 2. Freeze manifest BEFORE GDSC response is read
    # --------------------------------------------------------

    freeze_file = (
        OUT_DIR
        / "10_ANALYSIS_PLAN_FREEZE.txt"
    )


    with open(
        freeze_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 10 analysis plan freeze",
            file=f
        )

        print(
            "=======================================",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Primary dataset: GDSC2",
            file=f
        )

        print(
            "Primary response: LN_IC50",
            file=f
        )

        print(
            "Lower LN_IC50 = greater sensitivity",
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
            "LN_IC50 ~ ASTS + Lineage + Sex + Proliferation",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Primary PRISM candidates:",
            file=f
        )

        for _, row in primary_candidates.iterrows():

            print(
                f"  {row['compound_name']} "
                f"(PRISM beta="
                f"{row['beta_ASTS_log2AUC_per_1SD']})",
                file=f
            )

        print(
            "",
            file=f
        )

        print(
            "Primary test: one-sided test in the "
            "PRISM-predicted direction",
            file=f
        )

        print(
            "Primary multiplicity correction: Holm",
            file=f
        )

        print(
            "Unavailable/ambiguous primary candidates receive p=1",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            f"Minimum N per drug: {MIN_N_PER_DRUG}",
            file=f
        )

        print(
            f"Minimum global lineage N: "
            f"{MIN_GLOBAL_LINEAGE_N}",
            file=f
        )

        print(
            f"Minimum within-lineage N: "
            f"{MIN_WITHIN_LINEAGE_N}",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Drug matching rule:",
            file=f
        )

        print(
            "case/punctuation-insensitive canonical drug name only; "
            "no synonym substitution",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "No GDSC result may be used to redefine ASTS, "
            "candidate ranking, thresholds, or covariates.",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "ASTS checksum:",
            sha256_file(ASTS_FILE),
            file=f
        )

        print(
            "PRISM result checksum:",
            sha256_file(
                PRISM_RESULT_FILE
            ),
            file=f
        )

        print(
            "GDSC2 checksum:",
            sha256_file(
                GDSC_FILE
            ),
            file=f
        )


    # --------------------------------------------------------
    # 3. Read ASTS
    # --------------------------------------------------------

    print(
        "\nStep 2: Reading ASTS scores"
    )


    asts = pd.read_csv(
        ASTS_FILE,
        sep="\t"
    )


    if asts["ModelID"].duplicated().any():

        raise RuntimeError(
            "Duplicated ModelID found in Step08b ASTS file."
        )


    asts_columns = [
        "ModelID",
        PRIMARY_ASTS_COLUMN,
        HIGHCONF_ASTS_COLUMN,
        PROLIFERATION_COLUMN,
        LINEAGE_COLUMN,
        SEX_COLUMN,
    ]


    asts = asts[
        asts_columns
    ].copy()


    # --------------------------------------------------------
    # 4. SangerModelID -> DepMap ModelID
    # --------------------------------------------------------

    print(
        "\nStep 3: Building SangerModelID mapping"
    )


    model = pd.read_csv(
        MODEL_FILE
    )


    required_model = [
        "ModelID",
        "SangerModelID",
    ]


    for col in required_model:

        if col not in model.columns:

            raise ValueError(
                f"Missing Model.csv column: {col}"
            )


    model_map = model[
        required_model
    ].dropna().copy()


    model_map[
        "SangerModelID"
    ] = (
        model_map["SangerModelID"]
        .astype(str)
        .str.strip()
    )


    n_modelids = (
        model_map
        .groupby(
            "SangerModelID"
        )["ModelID"]
        .nunique()
    )


    unambiguous_sanger = set(
        n_modelids[
            n_modelids == 1
        ].index
    )


    model_map = model_map[
        model_map[
            "SangerModelID"
        ].isin(
            unambiguous_sanger
        )
    ].drop_duplicates(
        "SangerModelID"
    )


    print(
        "Unambiguous SangerModelIDs:",
        len(model_map)
    )


    # --------------------------------------------------------
    # 5. Read GDSC2
    # --------------------------------------------------------

    print(
        "\nStep 4: Reading GDSC2"
    )


    gdsc = pd.read_excel(
        GDSC_FILE,
        engine="openpyxl"
    )


    required_gdsc = [
        "SANGER_MODEL_ID",
        "DRUG_ID",
        "DRUG_NAME",
        "LN_IC50",
    ]


    missing = [
        x for x in required_gdsc
        if x not in gdsc.columns
    ]


    if missing:

        raise ValueError(
            "Missing GDSC2 columns: "
            + ", ".join(missing)
        )


    gdsc = gdsc[
        required_gdsc
    ].copy()


    gdsc["LN_IC50"] = pd.to_numeric(
        gdsc["LN_IC50"],
        errors="coerce"
    )


    gdsc = gdsc[
        gdsc["LN_IC50"].notna()
        &
        np.isfinite(
            gdsc["LN_IC50"]
        )
        &
        gdsc["SANGER_MODEL_ID"].notna()
        &
        gdsc["DRUG_NAME"].notna()
    ].copy()


    gdsc[
        "SANGER_MODEL_ID"
    ] = (
        gdsc["SANGER_MODEL_ID"]
        .astype(str)
        .str.strip()
    )


    gdsc["DRUG_ID"] = (
        gdsc["DRUG_ID"]
        .astype(str)
    )


    gdsc["drug_key"] = (
        gdsc["DRUG_NAME"]
        .apply(
            normalize_drug_name
        )
    )


    print(
        "Raw valid GDSC2 rows:",
        len(gdsc)
    )

    print(
        "GDSC2 drug IDs:",
        gdsc["DRUG_ID"].nunique()
    )


    # --------------------------------------------------------
    # 6. Collapse duplicate drug-model measurements
    # --------------------------------------------------------

    gdsc_collapsed = (
        gdsc
        .groupby(
            [
                "DRUG_ID",
                "DRUG_NAME",
                "drug_key",
                "SANGER_MODEL_ID",
            ],
            as_index=False
        )
        .agg(
            LN_IC50=(
                "LN_IC50",
                "median"
            ),
            n_original_rows=(
                "LN_IC50",
                "size"
            ),
        )
    )


    # --------------------------------------------------------
    # 7. Map GDSC2 -> DepMap ASTS
    # --------------------------------------------------------

    gdsc_mapped = gdsc_collapsed.merge(
        model_map,
        left_on="SANGER_MODEL_ID",
        right_on="SangerModelID",
        how="inner",
        validate="m:1"
    )


    gdsc_asts = gdsc_mapped.merge(
        asts,
        on="ModelID",
        how="inner",
        validate="m:1"
    )


    print(
        "GDSC2 × ASTS rows:",
        len(gdsc_asts)
    )

    print(
        "Models in intersection:",
        gdsc_asts[
            "ModelID"
        ].nunique()
    )


    # --------------------------------------------------------
    # 8. Define lineage categories BEFORE drug-wise fitting
    # --------------------------------------------------------

    model_metadata = (
        gdsc_asts[
            [
                "ModelID",
                LINEAGE_COLUMN,
            ]
        ]
        .drop_duplicates()
    )


    lineage_counts = (
        model_metadata[
            LINEAGE_COLUMN
        ]
        .fillna("Unknown")
        .value_counts()
    )


    common_lineages = set(
        lineage_counts[
            lineage_counts
            >= MIN_GLOBAL_LINEAGE_N
        ].index
    )


    def lineage_model(x):

        if pd.isna(x):
            return "Unknown"

        x = str(x)

        if x in common_lineages:
            return x

        return "OtherRare"


    gdsc_asts[
        "LineageModel"
    ] = (
        gdsc_asts[
            LINEAGE_COLUMN
        ]
        .apply(
            lineage_model
        )
    )


    # --------------------------------------------------------
    # 9. Analyze all GDSC2 drugs
    # --------------------------------------------------------

    print(
        "\nStep 5: Drug-wise GDSC2 models"
    )


    primary_rows = []
    highconf_rows = []


    for drug_id, df in gdsc_asts.groupby(
        "DRUG_ID"
    ):

        if (
            df["ModelID"].nunique()
            < MIN_N_PER_DRUG
        ):
            continue


        name = df[
            "DRUG_NAME"
        ].iloc[0]


        key = df[
            "drug_key"
        ].iloc[0]


        result = fit_model(
            df,
            PRIMARY_ASTS_COLUMN,
            include_lineage=True
        )


        if (
            result is not None
            and
            result["n"]
            >= MIN_N_PER_DRUG
            and
            result["df_resid"]
            >= MIN_RESIDUAL_DF
        ):

            primary_rows.append(
                {
                    "DRUG_ID":
                        drug_id,

                    "DRUG_NAME":
                        name,

                    "drug_key":
                        key,

                    **result,
                }
            )


        hc = fit_model(
            df,
            HIGHCONF_ASTS_COLUMN,
            include_lineage=True
        )


        if (
            hc is not None
            and
            hc["n"]
            >= MIN_N_PER_DRUG
            and
            hc["df_resid"]
            >= MIN_RESIDUAL_DF
        ):

            highconf_rows.append(
                {
                    "DRUG_ID":
                        drug_id,

                    "DRUG_NAME":
                        name,

                    "drug_key":
                        key,

                    **hc,
                }
            )


    results = pd.DataFrame(
        primary_rows
    )


    if len(results) == 0:

        raise RuntimeError(
            "No GDSC2 drugs passed analysis."
        )


    results[
        "FDR_exploratory"
    ] = bh_fdr(
        results[
            "p_ASTS_two_sided"
        ]
    )


    results = results.sort_values(
        [
            "FDR_exploratory",
            "p_ASTS_two_sided",
        ]
    )


    results.to_csv(
        OUT_DIR
        / "10_GDSC2_all_drugs.tsv",
        sep="\t",
        index=False
    )


    # High-confidence ASTS
    highconf = pd.DataFrame(
        highconf_rows
    )


    if len(highconf) > 0:

        highconf.to_csv(
            OUT_DIR
            / "10_GDSC2_highconf.tsv",
            sep="\t",
            index=False
        )


    # --------------------------------------------------------
    # 10. Match PRISM candidates to GDSC2
    # --------------------------------------------------------

    print(
        "\nStep 6: Primary candidate validation"
    )


    gdsc_name_to_ids = (
        gdsc_asts
        .groupby(
            "drug_key"
        )["DRUG_ID"]
        .unique()
        .to_dict()
    )


    result_by_id = (
        results
        .set_index(
            "DRUG_ID"
        )
    )


    def validate_candidates(
        candidates
    ):

        rows = []


        for _, row in candidates.iterrows():

            name = row[
                "compound_name"
            ]

            key = normalize_drug_name(
                name
            )

            prism_beta = row[
                "beta_ASTS_log2AUC_per_1SD"
            ]


            ids = gdsc_name_to_ids.get(
                key,
                []
            )


            if len(ids) == 0:

                status = "not_found_in_GDSC2"

                gdsc_id = ""
                gdsc_name = ""
                beta = np.nan
                p_two = np.nan
                p_one = 1.0
                same_direction = False
                n = 0


            elif len(ids) > 1:

                status = "ambiguous_multiple_GDSC2_drug_ids"

                gdsc_id = ";".join(
                    sorted(
                        map(str, ids)
                    )
                )

                gdsc_name = ""
                beta = np.nan
                p_two = np.nan
                p_one = 1.0
                same_direction = False
                n = 0


            else:

                gdsc_id = str(
                    ids[0]
                )


                if (
                    gdsc_id
                    not in result_by_id.index
                ):

                    status = (
                        "matched_but_failed_analysis_thresholds"
                    )

                    gdsc_name = (
                        gdsc_asts.loc[
                            gdsc_asts[
                                "DRUG_ID"
                            ] == gdsc_id,
                            "DRUG_NAME"
                        ]
                        .iloc[0]
                    )

                    beta = np.nan
                    p_two = np.nan
                    p_one = 1.0
                    same_direction = False
                    n = 0


                else:

                    status = "tested"

                    gd = result_by_id.loc[
                        gdsc_id
                    ]

                    gdsc_name = gd[
                        "DRUG_NAME"
                    ]

                    beta = gd[
                        "beta_ASTS_LNIC50_per_1SD"
                    ]

                    p_two = gd[
                        "p_ASTS_two_sided"
                    ]

                    n = int(
                        gd["n"]
                    )

                    same_direction = (
                        np.sign(beta)
                        ==
                        np.sign(prism_beta)
                    )

                    p_one = (
                        one_sided_p_from_two_sided(
                            beta,
                            p_two,
                            prism_beta
                        )
                    )


            rows.append(
                {
                    "PRISM_compound":
                        name,

                    "PRISM_FDR":
                        row["FDR_ASTS"],

                    "PRISM_beta":
                        prism_beta,

                    "drug_key":
                        key,

                    "GDSC_match_status":
                        status,

                    "GDSC_DRUG_ID":
                        gdsc_id,

                    "GDSC_DRUG_NAME":
                        gdsc_name,

                    "GDSC_n":
                        n,

                    "GDSC_beta":
                        beta,

                    "GDSC_p_two_sided":
                        p_two,

                    "direction_match":
                        same_direction,

                    "GDSC_p_one_sided_validation":
                        p_one,
                }
            )


        return pd.DataFrame(
            rows
        )


    primary_validation = (
        validate_candidates(
            primary_candidates
        )
    )


    # Holm correction across ALL predeclared primary candidates
    primary_validation[
        "Holm_p_primary"
    ] = multipletests(
        primary_validation[
            "GDSC_p_one_sided_validation"
        ],
        method="holm"
    )[1]


    primary_validation[
        "primary_validated"
    ] = (
        primary_validation[
            "direction_match"
        ]
        &
        (
            primary_validation[
                "Holm_p_primary"
            ] < 0.05
        )
    )


    primary_validation.to_csv(
        OUT_DIR
        / "10_primary_PRISM_candidates_validation.tsv",
        sep="\t",
        index=False
    )


    print(
        primary_validation.to_string(
            index=False
        )
    )


    # --------------------------------------------------------
    # 11. Secondary PRISM Top20
    # --------------------------------------------------------

    secondary = validate_candidates(
        secondary_candidates
    )


    secondary[
        "FDR_one_sided_top20"
    ] = bh_fdr(
        secondary[
            "GDSC_p_one_sided_validation"
        ]
    )


    secondary.to_csv(
        OUT_DIR
        / "10_PRISM_top20_validation.tsv",
        sep="\t",
        index=False
    )


    # --------------------------------------------------------
    # 12. Global PRISM-GDSC concordance
    # --------------------------------------------------------

    gdsc_result_name_counts = (
        results
        .groupby(
            "drug_key"
        )["DRUG_ID"]
        .nunique()
    )


    unique_keys = set(
        gdsc_result_name_counts[
            gdsc_result_name_counts == 1
        ].index
    )


    gdsc_unique = (
        results[
            results[
                "drug_key"
            ].isin(
                unique_keys
            )
        ]
        .drop_duplicates(
            "drug_key"
        )
    )


    prism_unique = (
        prism
        .drop_duplicates(
            "drug_key"
        )
    )


    concordance = (
        prism_unique.merge(
            gdsc_unique[
                [
                    "drug_key",
                    "DRUG_ID",
                    "DRUG_NAME",
                    "beta_ASTS_LNIC50_per_1SD",
                ]
            ],
            on="drug_key",
            how="inner"
        )
    )


    if len(concordance) >= 3:

        global_rho, global_p = (
            spearmanr(
                concordance[
                    "beta_ASTS_log2AUC_per_1SD"
                ],
                concordance[
                    "beta_ASTS_LNIC50_per_1SD"
                ]
            )
        )

        global_sign = np.mean(
            np.sign(
                concordance[
                    "beta_ASTS_log2AUC_per_1SD"
                ]
            )
            ==
            np.sign(
                concordance[
                    "beta_ASTS_LNIC50_per_1SD"
                ]
            )
        )

    else:

        global_rho = np.nan
        global_p = np.nan
        global_sign = np.nan


    concordance.to_csv(
        OUT_DIR
        / "10_PRISM_vs_GDSC2_overlap.tsv",
        sep="\t",
        index=False
    )


    # --------------------------------------------------------
    # 13. Within-lineage analysis for primary candidates
    # --------------------------------------------------------

    within_rows = []
    meta_rows = []


    tested_primary = primary_validation[
        primary_validation[
            "GDSC_match_status"
        ] == "tested"
    ]


    for _, candidate in tested_primary.iterrows():

        drug_id = str(
            candidate[
                "GDSC_DRUG_ID"
            ]
        )


        ddf = gdsc_asts[
            gdsc_asts[
                "DRUG_ID"
            ] == drug_id
        ].copy()


        candidate_within = []


        for lineage, ldf in ddf.groupby(
            LINEAGE_COLUMN
        ):

            if pd.isna(lineage):
                continue

            if len(ldf) < MIN_WITHIN_LINEAGE_N:
                continue


            res = fit_model(
                ldf,
                PRIMARY_ASTS_COLUMN,
                include_lineage=False
            )


            if (
                res is None
                or
                res["n"]
                < MIN_WITHIN_LINEAGE_N
            ):
                continue


            row = {
                "PRISM_compound":
                    candidate[
                        "PRISM_compound"
                    ],

                "GDSC_DRUG_ID":
                    drug_id,

                "lineage":
                    lineage,

                **res,
            }


            within_rows.append(
                row
            )

            candidate_within.append(
                row
            )


        if len(candidate_within) >= MIN_LINEAGES_FOR_META:

            meta = fixed_effect_meta(
                pd.DataFrame(
                    candidate_within
                )
            )


            if meta is not None:

                meta_rows.append(
                    {
                        "PRISM_compound":
                            candidate[
                                "PRISM_compound"
                            ],

                        "GDSC_DRUG_ID":
                            drug_id,

                        **meta,
                    }
                )


    pd.DataFrame(
        within_rows
    ).to_csv(
        OUT_DIR
        / "10_primary_candidates_within_lineage.tsv",
        sep="\t",
        index=False
    )


    pd.DataFrame(
        meta_rows
    ).to_csv(
        OUT_DIR
        / "10_primary_candidates_within_lineage_meta.tsv",
        sep="\t",
        index=False
    )


    # --------------------------------------------------------
    # 14. Primary vs high-confidence ASTS in GDSC2
    # --------------------------------------------------------

    hc_rho = np.nan
    hc_sign = np.nan


    if len(highconf) > 0:

        hc_compare = results[
            [
                "DRUG_ID",
                "beta_ASTS_LNIC50_per_1SD",
            ]
        ].merge(
            highconf[
                [
                    "DRUG_ID",
                    "beta_ASTS_LNIC50_per_1SD",
                ]
            ],
            on="DRUG_ID",
            suffixes=(
                "_primary",
                "_highconf"
            )
        )


        if len(hc_compare) >= 3:

            hc_rho, _ = spearmanr(
                hc_compare[
                    "beta_ASTS_LNIC50_per_1SD_primary"
                ],
                hc_compare[
                    "beta_ASTS_LNIC50_per_1SD_highconf"
                ]
            )


            hc_sign = np.mean(
                np.sign(
                    hc_compare[
                        "beta_ASTS_LNIC50_per_1SD_primary"
                    ]
                )
                ==
                np.sign(
                    hc_compare[
                        "beta_ASTS_LNIC50_per_1SD_highconf"
                    ]
                )
            )


    # --------------------------------------------------------
    # 15. Figures
    # --------------------------------------------------------

    # GDSC volcano
    plot_df = results[
        results[
            "p_ASTS_two_sided"
        ] > 0
    ].copy()


    fig, ax = plt.subplots(
        figsize=(8, 6)
    )


    ax.scatter(
        plot_df[
            "beta_ASTS_standardized"
        ],
        -np.log10(
            plot_df[
                "p_ASTS_two_sided"
            ]
        ),
        s=16,
        alpha=0.6
    )


    ax.axvline(
        0,
        linestyle="--"
    )


    ax.set_xlabel(
        "GDSC2 ASTS standardized beta\n"
        "negative = higher ASTS, greater sensitivity"
    )

    ax.set_ylabel(
        "-log10(p)"
    )

    ax.set_title(
        "GDSC2: ASTS drug associations"
    )


    fig.tight_layout()


    fig.savefig(
        FIG_DIR
        / "01_GDSC2_ASTS_volcano.png",
        dpi=200
    )


    plt.close(fig)


    # PRISM vs GDSC beta
    if len(concordance) >= 3:

        fig, ax = plt.subplots(
            figsize=(6, 6)
        )


        ax.scatter(
            concordance[
                "beta_ASTS_log2AUC_per_1SD"
            ],
            concordance[
                "beta_ASTS_LNIC50_per_1SD"
            ],
            s=16,
            alpha=0.6
        )


        ax.axhline(
            0,
            linestyle="--"
        )

        ax.axvline(
            0,
            linestyle="--"
        )


        ax.set_xlabel(
            "PRISM ASTS beta"
        )

        ax.set_ylabel(
            "GDSC2 ASTS beta"
        )


        ax.set_title(
            f"PRISM vs GDSC2\n"
            f"Spearman r={global_rho:.3f}"
        )


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / "02_PRISM_vs_GDSC2_beta.png",
            dpi=200
        )


        plt.close(fig)


    # --------------------------------------------------------
    # 16. Summary
    # --------------------------------------------------------

    summary_file = (
        OUT_DIR
        / "10_summary.txt"
    )


    with open(
        summary_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 10 summary",
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
            "GDSC2 response: LN_IC50",
            file=f
        )

        print(
            "Lower LN_IC50 = greater sensitivity",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            f"GDSC2 valid raw rows: "
            f"{len(gdsc)}",
            file=f
        )

        print(
            f"GDSC2 x ASTS rows: "
            f"{len(gdsc_asts)}",
            file=f
        )

        print(
            f"Models in intersection: "
            f"{gdsc_asts['ModelID'].nunique()}",
            file=f
        )

        print(
            f"Drugs tested with N >= "
            f"{MIN_N_PER_DRUG}: "
            f"{len(results)}",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "Primary PRISM candidates",
            file=f
        )

        print(
            "------------------------",
            file=f
        )

        print(
            primary_validation.to_string(
                index=False
            ),
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "PRISM-GDSC2 global concordance",
            file=f
        )

        print(
            "------------------------------",
            file=f
        )

        print(
            f"Matched drugs: "
            f"{len(concordance)}",
            file=f
        )

        print(
            f"Beta Spearman r: "
            f"{global_rho}",
            file=f
        )

        print(
            f"Beta Spearman p: "
            f"{global_p}",
            file=f
        )

        print(
            f"Sign concordance: "
            f"{global_sign}",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "78-gene vs 77-gene GDSC2 sensitivity",
            file=f
        )

        print(
            "------------------------------------",
            file=f
        )

        print(
            f"Beta Spearman r: "
            f"{hc_rho}",
            file=f
        )

        print(
            f"Sign concordance: "
            f"{hc_sign}",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "Checksums",
            file=f
        )

        print(
            "---------",
            file=f
        )

        print(
            "ASTS:",
            sha256_file(
                ASTS_FILE
            ),
            file=f
        )

        print(
            "PRISM result:",
            sha256_file(
                PRISM_RESULT_FILE
            ),
            file=f
        )

        print(
            "GDSC2:",
            sha256_file(
                GDSC_FILE
            ),
            file=f
        )


    print(
        "\n============================================================"
    )

    print(
        "Step 10 completed successfully"
    )

    print(
        "============================================================"
    )

    print(
        "\nSummary:",
        summary_file
    )


if __name__ == "__main__":
    main()
