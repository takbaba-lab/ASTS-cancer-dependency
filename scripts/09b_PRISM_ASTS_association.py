#!/usr/bin/env python3

"""
ASTRA-Drug
Step 09: PRISM drug-response association with ASTS

----

    log2(AUC) ~ ASTS + lineage + sex + proliferation


Primary:
    78-gene human ASTS v1.0

Sensitivity:
    77-gene high-confidence ASTS

--------------

    beta_ASTS < 0
        Higher ASTS -> greater sensitivity

    beta_ASTS > 0
        Higher ASTS -> greater resistance


----
"""


# ============================================================
# Path / Parameter
# ============================================================

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]


ASTS_SCORE_FILE = (
    PROJECT_DIR
    / "results/08b_DepMap_ASTS_default_profile"
    / "08_DepMap_ASTS_scores.tsv"
)


PRISM_FILE = (
    PROJECT_DIR
    / "data/raw/prism/secondary"
    / "secondary-screen-dose-response-curve-parameters.csv"
)


OUT_DIR = (
    PROJECT_DIR
    / "results/09b_PRISM_ASTS_default_profile"
)

FIG_DIR = OUT_DIR / "figures"


# ------------------------------------------------------------
# Primary analysis parameters
# ------------------------------------------------------------

PRIMARY_ASTS_COLUMN = "ASTS_score"

HIGH_CONF_ASTS_COLUMN = "ASTS_high_conf_score"

PROLIFERATION_COLUMN = "proliferation_proxy"

LINEAGE_COLUMN = "OncotreeLineage"

SEX_COLUMN = "Sex"


MIN_N_PER_DRUG = 100


MIN_GLOBAL_LINEAGE_N = 5


# within-lineage sensitivity
MIN_WITHIN_LINEAGE_N = 15


MIN_LINEAGES_FOR_META = 3


MIN_RESIDUAL_DF = 30


# Multiple testing
FDR_ALPHA = 0.05


# ============================================================
# Imports
# ============================================================

import hashlib
import math

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import statsmodels.api as sm

from scipy.stats import (
    spearmanr,
    norm,
)

from statsmodels.stats.multitest import (
    multipletests,
)


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


def zscore(series):
    """
    """

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


def combine_unique(series):
    """
    """

    x = (
        series
        .dropna()
        .astype(str)
        .str.strip()
    )

    x = x[
        x != ""
    ]

    if len(x) == 0:
        return ""

    return ";".join(
        sorted(
            set(x)
        )
    )


def bh_fdr(pvalues):
    """
    Benjamini-Hochberg FDR。
    """

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
            alpha=FDR_ALPHA,
            method="fdr_bh"
        )[1]

    return q


def fit_model(
    df,
    asts_column,
    include_lineage=True
):
    """

    Primary:
        log2_auc ~ ASTS_z + proliferation_z
                   + lineage + sex

    within-lineage:
        log2_auc ~ ASTS_z + proliferation_z
                   + sex
    """

    use = df.copy()


    required = [
        "log2_auc",
        asts_column,
        PROLIFERATION_COLUMN,
    ]

    if include_lineage:
        required.append(
            "LineageModel"
        )


    use = use.dropna(
        subset=required
    )


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
    )


    if len(use) < 5:

        return None


    # --------------------------------------------
    # Design matrix
    # --------------------------------------------

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
        use["log2_auc"],
        errors="coerce"
    )


    # --------------------------------------------
    # Fit
    # --------------------------------------------

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

    tvalue = float(
        model.tvalues[
            "ASTS_z_model"
        ]
    )


    response_sd = y.std(
        ddof=1
    )


    if (
        np.isfinite(response_sd)
        and response_sd > 0
    ):

        beta_standardized = (
            beta
            / response_sd
        )

    else:

        beta_standardized = np.nan


    df_resid = float(
        model.df_resid
    )


    if (
        np.isfinite(tvalue)
        and df_resid > 0
    ):

        partial_r2 = (
            tvalue ** 2
            /
            (
                tvalue ** 2
                + df_resid
            )
        )

    else:

        partial_r2 = np.nan


    # unadjusted correlation for diagnostic only
    rho, rho_p = spearmanr(
        use[asts_column],
        use["log2_auc"],
        nan_policy="omit"
    )


    return {
        "n":
            len(use),

        "n_lineages":
            use["LineageModel"].nunique()
            if include_lineage
            else 1,

        "beta_ASTS_log2AUC_per_1SD":
            beta,

        "se_ASTS":
            se,

        "t_ASTS":
            tvalue,

        "p_ASTS":
            pvalue,

        "beta_ASTS_standardized":
            beta_standardized,

        "partial_R2_ASTS":
            partial_r2,

        "model_R2":
            float(model.rsquared),

        "model_adj_R2":
            float(model.rsquared_adj),

        "df_resid":
            df_resid,

        "unadjusted_spearman":
            float(rho),

        "unadjusted_spearman_p":
            float(rho_p),
    }


def fixed_effect_meta(df):
    """
    """

    x = df.copy()

    x = x[
        np.isfinite(
            x["beta_ASTS_log2AUC_per_1SD"]
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
        "beta_ASTS_log2AUC_per_1SD"
    ].to_numpy()

    se = x[
        "se_ASTS"
    ].to_numpy()


    weights = (
        1.0
        / se ** 2
    )


    meta_beta = np.sum(
        weights * beta
    ) / np.sum(
        weights
    )


    meta_se = math.sqrt(
        1.0
        / np.sum(weights)
    )


    z = (
        meta_beta
        / meta_se
    )


    p = (
        2
        * norm.sf(
            abs(z)
        )
    )


    same_direction = np.mean(
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
            same_direction,
    }


# ============================================================
# Main
# ============================================================

def main():

    # --------------------------------------------------------
    # 0. Input check
    # --------------------------------------------------------

    print(
        "\n"
        "============================================================\n"
        "ASTRA-Drug Step 09\n"
        "PRISM × ASTS association\n"
        "============================================================"
    )


    for path in [
        ASTS_SCORE_FILE,
        PRISM_FILE,
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
    # 1. Freeze analysis plan BEFORE reading PRISM responses
    # --------------------------------------------------------

    freeze_file = (
        OUT_DIR
        / "09_ANALYSIS_PLAN_FREEZE.txt"
    )


    with open(
        freeze_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 09 analysis plan freeze",
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
            "Primary response:",
            "log2(PRISM Secondary AUC)",
            file=f
        )

        print(
            "Lower response = greater sensitivity",
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
            "log2(AUC) ~ ASTS + Lineage + Sex + Proliferation",
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            f"Primary ASTS column: "
            f"{PRIMARY_ASTS_COLUMN}",
            file=f
        )

        print(
            f"High-confidence sensitivity column: "
            f"{HIGH_CONF_ASTS_COLUMN}",
            file=f
        )

        print(
            f"Minimum N per drug: "
            f"{MIN_N_PER_DRUG}",
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
            f"Minimum lineages for meta-analysis: "
            f"{MIN_LINEAGES_FOR_META}",
            file=f
        )

        print(
            f"FDR alpha: "
            f"{FDR_ALPHA}",
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Primary interpretation:",
            file=f
        )

        print(
            "beta < 0: higher ASTS -> greater drug sensitivity",
            file=f
        )

        print(
            "beta > 0: higher ASTS -> greater drug resistance",
            file=f
        )


        print(
            "",
            file=f
        )

        print(
            "No drug-response result may be used to redefine ASTS, "
            "lineage inclusion, or thresholds.",
            file=f
        )


        print(
            "",
            file=f
        )

        print(
            "ASTS score checksum:",
            sha256_file(
                ASTS_SCORE_FILE
            ),
            file=f
        )

        print(
            "PRISM checksum:",
            sha256_file(
                PRISM_FILE
            ),
            file=f
        )


    print(
        "Analysis plan frozen:",
        freeze_file
    )


    # --------------------------------------------------------
    # 2. Read ASTS
    # --------------------------------------------------------

    print(
        "\nStep 1: Reading ASTS scores"
    )


    asts = pd.read_csv(
        ASTS_SCORE_FILE,
        sep="\t"
    )


    required_asts = [
        "ModelID",
        PRIMARY_ASTS_COLUMN,
        HIGH_CONF_ASTS_COLUMN,
        PROLIFERATION_COLUMN,
        LINEAGE_COLUMN,
        SEX_COLUMN,
    ]


    missing = [
        x for x in required_asts
        if x not in asts.columns
    ]


    if missing:

        raise ValueError(
            "Missing ASTS columns: "
            + ", ".join(missing)
        )


    asts = asts[
        required_asts
    ].copy()


    # --------------------------------------------------------
    # 3. Read PRISM
    # --------------------------------------------------------

    print(
        "\nStep 2: Reading PRISM Secondary"
    )


    header = pd.read_csv(
        PRISM_FILE,
        nrows=0
    )


    required_prism = [
        "broad_id",
        "depmap_id",
        "auc",
        "name",
    ]


    missing = [
        x for x in required_prism
        if x not in header.columns
    ]


    if missing:

        raise ValueError(
            "Missing PRISM columns: "
            + ", ".join(missing)
        )


    optional_columns = [
        "moa",
        "target",
        "disease.area",
        "disease_area",
        "indication",
        "phase",
        "screen_id",
        "passed_str_profiling",
    ]


    usecols = (
        required_prism
        +
        [
            x for x in optional_columns
            if x in header.columns
        ]
    )


    prism = pd.read_csv(
        PRISM_FILE,
        usecols=usecols,
        low_memory=False
    )


    print(
        "Raw PRISM rows:",
        len(prism)
    )


    prism["auc"] = pd.to_numeric(
        prism["auc"],
        errors="coerce"
    )


    prism = prism[
        prism["auc"].notna()
        &
        np.isfinite(
            prism["auc"]
        )
        &
        (
            prism["auc"] > 0
        )
    ].copy()


    prism["log2_auc"] = np.log2(
        prism["auc"]
    )


    # --------------------------------------------------------
    # 4. Collapse duplicate compound-name × cell-line entries
    # --------------------------------------------------------

    print(
        "\nStep 3: Collapsing duplicate drug/model entries"
    )


    prism["compound_name"] = (
        prism["name"]
        .fillna("")
        .astype(str)
        .str.strip()
    )


    prism = prism[
        prism["compound_name"]
        != ""
    ].copy()


    prism["compound_key"] = (
        prism["compound_name"]
        .str.upper()
    )


    annotation_columns = [
        x for x in [
            "moa",
            "target",
            "disease.area",
            "disease_area",
            "indication",
            "phase",
        ]
        if x in prism.columns
    ]


    agg_dict = {
        "log2_auc":
            "median",

        "auc":
            "median",

        "compound_name":
            "first",

        "broad_id":
            pd.Series.nunique,
    }


    for col in annotation_columns:

        agg_dict[col] = (
            combine_unique
        )


    prism_collapsed = (
        prism
        .groupby(
            [
                "compound_key",
                "depmap_id",
            ],
            as_index=False
        )
        .agg(
            agg_dict
        )
    )


    prism_collapsed = (
        prism_collapsed
        .rename(
            columns={
                "broad_id":
                    "n_broad_ids"
            }
        )
    )


    print(
        "Collapsed rows:",
        len(prism_collapsed)
    )


    print(
        "Unique compounds:",
        prism_collapsed[
            "compound_key"
        ].nunique()
    )


    # --------------------------------------------------------
    # 5. Merge ASTS
    # --------------------------------------------------------

    merged = prism_collapsed.merge(
        asts,
        left_on="depmap_id",
        right_on="ModelID",
        how="inner"
    )


    print(
        "PRISM × ASTS rows:",
        len(merged)
    )

    print(
        "Models in intersection:",
        merged[
            "ModelID"
        ].nunique()
    )


    # --------------------------------------------------------
    # 6. Predefine lineage categories
    #
    # IMPORTANT:
    # --------------------------------------------------------

    model_metadata = (
        merged[
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


    merged["LineageModel"] = (
        merged[
            LINEAGE_COLUMN
        ]
        .apply(
            lineage_model
        )
    )


    # --------------------------------------------------------
    # 7. Primary drug-wise associations
    # --------------------------------------------------------

    print(
        "\nStep 4: Primary drug-wise models"
    )


    primary_rows = []

    highconf_rows = []

    within_rows = []


    compounds = (
        merged[
            "compound_key"
        ]
        .drop_duplicates()
        .tolist()
    )


    for i, compound in enumerate(
        compounds,
        start=1
    ):

        df = merged[
            merged["compound_key"]
            == compound
        ].copy()


        if len(df) < MIN_N_PER_DRUG:
            continue


        display_name = (
            df["compound_name"]
            .iloc[0]
        )


        # --------------------------------------------
        # annotation
        # --------------------------------------------

        annotations = {
            "compound_key":
                compound,

            "compound_name":
                display_name,

            "n_total_before_model":
                len(df),

            "n_broad_ids_max":
                int(
                    df[
                        "n_broad_ids"
                    ].max()
                ),
        }


        for col in annotation_columns:

            annotations[col] = (
                combine_unique(
                    df[col]
                )
            )


        # --------------------------------------------
        # Primary
        # --------------------------------------------

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

            row = {
                **annotations,
                **result,
            }

            row["association_direction"] = (
                "higher_ASTS_more_sensitive"
                if result[
                    "beta_ASTS_log2AUC_per_1SD"
                ] < 0
                else
                "higher_ASTS_more_resistant"
            )


            primary_rows.append(
                row
            )


        # --------------------------------------------
        # High-confidence sensitivity
        # --------------------------------------------

        hc_result = fit_model(
            df,
            HIGH_CONF_ASTS_COLUMN,
            include_lineage=True
        )


        if (
            hc_result is not None
            and
            hc_result["n"]
            >= MIN_N_PER_DRUG
            and
            hc_result["df_resid"]
            >= MIN_RESIDUAL_DF
        ):

            row = {
                **annotations,
                **hc_result,
            }

            row["association_direction"] = (
                "higher_ASTS_more_sensitive"
                if hc_result[
                    "beta_ASTS_log2AUC_per_1SD"
                ] < 0
                else
                "higher_ASTS_more_resistant"
            )


            highconf_rows.append(
                row
            )


        # --------------------------------------------
        # Within-lineage sensitivity
        # --------------------------------------------

        for lineage, ldf in df.groupby(
            "LineageModel"
        ):

            if (
                lineage
                in [
                    "OtherRare",
                    "Unknown",
                ]
            ):
                continue


            if len(ldf) < MIN_WITHIN_LINEAGE_N:
                continue


            lresult = fit_model(
                ldf,
                PRIMARY_ASTS_COLUMN,
                include_lineage=False
            )


            if (
                lresult is None
                or
                lresult["n"]
                < MIN_WITHIN_LINEAGE_N
            ):
                continue


            within_rows.append(
                {
                    "compound_key":
                        compound,

                    "compound_name":
                        display_name,

                    "lineage":
                        lineage,

                    **lresult,
                }
            )


        if (
            i % 100 == 0
        ):

            print(
                f"  processed {i}/"
                f"{len(compounds)} compounds"
            )


    # --------------------------------------------------------
    # 8. Primary results + FDR
    # --------------------------------------------------------

    primary = pd.DataFrame(
        primary_rows
    )


    if len(primary) == 0:

        raise RuntimeError(
            "No drugs passed primary analysis."
        )


    primary[
        "FDR_ASTS"
    ] = bh_fdr(
        primary[
            "p_ASTS"
        ]
    )


    primary = primary.sort_values(
        [
            "FDR_ASTS",
            "p_ASTS",
        ]
    )


    primary.to_csv(
        OUT_DIR
        / "09_PRISM_ASTS_primary.tsv",
        sep="\t",
        index=False
    )


    # --------------------------------------------------------
    # 9. High-confidence sensitivity
    # --------------------------------------------------------

    highconf = pd.DataFrame(
        highconf_rows
    )


    if len(highconf) > 0:

        highconf[
            "FDR_ASTS"
        ] = bh_fdr(
            highconf[
                "p_ASTS"
            ]
        )


        highconf = (
            highconf
            .sort_values(
                [
                    "FDR_ASTS",
                    "p_ASTS",
                ]
            )
        )


        highconf.to_csv(
            OUT_DIR
            / "09_PRISM_ASTS_highconf.tsv",
            sep="\t",
            index=False
        )


    # --------------------------------------------------------
    # 10. Within-lineage analysis
    # --------------------------------------------------------

    within = pd.DataFrame(
        within_rows
    )


    if len(within) > 0:

        within[
            "FDR_within_all_tests"
        ] = bh_fdr(
            within[
                "p_ASTS"
            ]
        )


        within.to_csv(
            OUT_DIR
            / "09_PRISM_ASTS_within_lineage.tsv",
            sep="\t",
            index=False
        )


        # --------------------------------------------
        # Fixed-effect meta per drug
        # --------------------------------------------

        meta_rows = []


        for compound, x in within.groupby(
            "compound_key"
        ):

            meta = fixed_effect_meta(
                x
            )


            if meta is None:
                continue


            meta_rows.append(
                {
                    "compound_key":
                        compound,

                    "compound_name":
                        x[
                            "compound_name"
                        ].iloc[0],

                    **meta,
                }
            )


        meta_df = pd.DataFrame(
            meta_rows
        )


        if len(meta_df) > 0:

            meta_df[
                "meta_FDR"
            ] = bh_fdr(
                meta_df[
                    "meta_p"
                ]
            )


            meta_df = (
                meta_df
                .sort_values(
                    [
                        "meta_FDR",
                        "meta_p",
                    ]
                )
            )


            meta_df.to_csv(
                OUT_DIR
                / "09_PRISM_ASTS_within_lineage_meta.tsv",
                sep="\t",
                index=False
            )

    else:

        meta_df = pd.DataFrame()


    # --------------------------------------------------------
    # 11. Primary vs high-confidence concordance
    # --------------------------------------------------------

    if len(highconf) > 0:

        compare = primary[
            [
                "compound_key",
                "compound_name",
                "beta_ASTS_log2AUC_per_1SD",
                "FDR_ASTS",
            ]
        ].merge(
            highconf[
                [
                    "compound_key",
                    "beta_ASTS_log2AUC_per_1SD",
                    "FDR_ASTS",
                ]
            ],
            on="compound_key",
            suffixes=(
                "_primary",
                "_highconf"
            )
        )


        rho, rho_p = spearmanr(
            compare[
                "beta_ASTS_log2AUC_per_1SD_primary"
            ],
            compare[
                "beta_ASTS_log2AUC_per_1SD_highconf"
            ]
        )


        sign_concordance = np.mean(
            np.sign(
                compare[
                    "beta_ASTS_log2AUC_per_1SD_primary"
                ]
            )
            ==
            np.sign(
                compare[
                    "beta_ASTS_log2AUC_per_1SD_highconf"
                ]
            )
        )


        compare.to_csv(
            OUT_DIR
            / "09_primary_vs_highconf.tsv",
            sep="\t",
            index=False
        )

    else:

        rho = np.nan
        rho_p = np.nan
        sign_concordance = np.nan


    # --------------------------------------------------------
    # 12. Figures
    # --------------------------------------------------------

    print(
        "\nStep 5: Creating figures"
    )


    # Volcano
    plot_df = primary[
        primary[
            "p_ASTS"
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
                "p_ASTS"
            ]
        ),
        s=15,
        alpha=0.6
    )


    ax.axvline(
        0,
        linestyle="--"
    )


    ax.set_xlabel(
        "ASTS standardized beta\n"
        "negative = higher ASTS, greater sensitivity"
    )

    ax.set_ylabel(
        "-log10(p)"
    )

    ax.set_title(
        "PRISM Secondary: ASTS drug associations"
    )


    fig.tight_layout()

    fig.savefig(
        FIG_DIR
        / "01_PRISM_ASTS_volcano.png",
        dpi=200
    )

    plt.close(fig)


    # Primary vs high-confidence
    if len(highconf) > 0:

        fig, ax = plt.subplots(
            figsize=(6, 6)
        )


        ax.scatter(
            compare[
                "beta_ASTS_log2AUC_per_1SD_primary"
            ],
            compare[
                "beta_ASTS_log2AUC_per_1SD_highconf"
            ],
            s=12,
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
            "Primary ASTS beta"
        )

        ax.set_ylabel(
            "High-confidence ASTS beta"
        )


        ax.set_title(
            f"Primary vs high-confidence\n"
            f"Spearman r={rho:.3f}"
        )


        fig.tight_layout()

        fig.savefig(
            FIG_DIR
            / "02_primary_vs_highconf_beta.png",
            dpi=200
        )

        plt.close(fig)


    # Primary vs within-lineage meta
    if (
        "meta_df" in locals()
        and
        len(meta_df) > 0
    ):

        meta_compare = primary[
            [
                "compound_key",
                "beta_ASTS_log2AUC_per_1SD",
            ]
        ].merge(
            meta_df[
                [
                    "compound_key",
                    "meta_beta",
                ]
            ],
            on="compound_key"
        )


        fig, ax = plt.subplots(
            figsize=(6, 6)
        )


        ax.scatter(
            meta_compare[
                "beta_ASTS_log2AUC_per_1SD"
            ],
            meta_compare[
                "meta_beta"
            ],
            s=12,
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
            "Primary lineage-adjusted beta"
        )

        ax.set_ylabel(
            "Within-lineage meta beta"
        )


        ax.set_title(
            "Primary vs within-lineage meta"
        )


        fig.tight_layout()

        fig.savefig(
            FIG_DIR
            / "03_primary_vs_within_lineage_meta.png",
            dpi=200
        )

        plt.close(fig)


    # --------------------------------------------------------
    # 13. Summary
    # --------------------------------------------------------

    summary_file = (
        OUT_DIR
        / "09_summary.txt"
    )


    with open(
        summary_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 09 summary",
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
            "Primary response: log2(PRISM Secondary AUC)",
            file=f
        )

        print(
            "Lower value = greater drug sensitivity",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            f"Raw PRISM rows: "
            f"{len(prism)}",
            file=f
        )

        print(
            f"Collapsed drug/model rows: "
            f"{len(prism_collapsed)}",
            file=f
        )

        print(
            f"PRISM x ASTS rows: "
            f"{len(merged)}",
            file=f
        )

        print(
            f"Models in intersection: "
            f"{merged['ModelID'].nunique()}",
            file=f
        )

        print(
            f"Unique compounds before filtering: "
            f"{merged['compound_key'].nunique()}",
            file=f
        )

        print(
            f"Drugs tested in primary model: "
            f"{len(primary)}",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "Primary FDR results",
            file=f
        )

        print(
            "-------------------",
            file=f
        )


        print(
            f"FDR < 0.05: "
            f"{(primary['FDR_ASTS'] < 0.05).sum()}",
            file=f
        )

        print(
            f"FDR < 0.10: "
            f"{(primary['FDR_ASTS'] < 0.10).sum()}",
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Top 20 primary associations",
            file=f
        )

        print(
            "---------------------------",
            file=f
        )


        display_cols = [
            "compound_name",
            "n",
            "beta_ASTS_log2AUC_per_1SD",
            "beta_ASTS_standardized",
            "partial_R2_ASTS",
            "p_ASTS",
            "FDR_ASTS",
            "association_direction",
        ]


        print(
            primary[
                display_cols
            ]
            .head(20)
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
            "Sensitivity: 78-gene vs 77-gene",
            file=f
        )

        print(
            "------------------------------",
            file=f
        )


        print(
            f"Beta Spearman r: "
            f"{rho}",
            file=f
        )

        print(
            f"Beta Spearman p: "
            f"{rho_p}",
            file=f
        )

        print(
            f"Sign concordance: "
            f"{sign_concordance}",
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
            "ASTS scores:",
            sha256_file(
                ASTS_SCORE_FILE
            ),
            file=f
        )

        print(
            "PRISM:",
            sha256_file(
                PRISM_FILE
            ),
            file=f
        )


    print(
        "\n============================================================"
    )

    print(
        "09 completed successfully"
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
