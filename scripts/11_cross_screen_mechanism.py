#!/usr/bin/env python3

"""
ASTRA-Drug
Step 11: Cross-screen mechanistic analysis

----



----



    PUTATIVE_TARGET
    PATHWAY_NAME


"""


# ============================================================
# Path / Parameter
# ============================================================

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]


OVERLAP_FILE = (
    PROJECT_DIR
    / "results/10_GDSC2_validation"
    / "10_PRISM_vs_GDSC2_overlap.tsv"
)


# Step10: GDSC2 full association results
GDSC_RESULT_FILE = (
    PROJECT_DIR
    / "results/10_GDSC2_validation"
    / "10_GDSC2_all_drugs.tsv"
)


# Original GDSC2 workbook
GDSC2_FILE = (
    PROJECT_DIR
    / "data/raw/gdsc/GDSC2"
    / "GDSC2_fitted_dose_response_27Oct23.xlsx"
)


OUT_DIR = (
    PROJECT_DIR
    / "results/11_cross_screen_mechanism"
)

FIG_DIR = OUT_DIR / "figures"


# ------------------------------------------------------------
# Expected overlap from Step10
# ------------------------------------------------------------

EXPECTED_OVERLAP_N = 108


# ------------------------------------------------------------
# Enrichment parameters
# ------------------------------------------------------------

MIN_CATEGORY_UNIVERSE_N = 3

MIN_CATEGORY_GROUP_N = 2

FDR_ALPHA = 0.05

TOP_N_ENRICHMENT_PLOT = 15

TOP_N_DRUGS = 20


# ============================================================
# Imports
# ============================================================

import hashlib
import math
import re

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from scipy.stats import (
    spearmanr,
    fisher_exact,
    binomtest,
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


def normalize_drug_id(x):
    """
        1234
        1234.0
    """

    if pd.isna(x):
        return ""

    x = str(x).strip()

    if re.fullmatch(
        r"\d+\.0",
        x
    ):
        x = x[:-2]

    return x


def clean_annotation(x):
    """Perform basic cleaning of annotation strings."""

    if pd.isna(x):
        return ""

    x = str(x).strip()

    if x.lower() in {
        "",
        "nan",
        "none",
        "na",
        "n/a",
    }:
        return ""

    return x


def split_pathways(x):
    """

    """

    x = clean_annotation(x)

    if x == "":
        return []

    values = re.split(
        r"[;|]+",
        x
    )

    return sorted(
        {
            v.strip()
            for v in values
            if v.strip()
        }
    )


def split_targets(x):
    """

    """

    x = clean_annotation(x)

    if x == "":
        return []

    values = re.split(
        r"[,;]+",
        x
    )

    return sorted(
        {
            v.strip()
            for v in values
            if v.strip()
        }
    )


def join_unique(series):
    """Combine multiple annotations assigned to the same DRUG_ID."""

    values = []

    for value in series:

        value = clean_annotation(
            value
        )

        if value != "":
            values.append(
                value
            )

    return ";".join(
        sorted(
            set(values)
        )
    )


def classify_direction(
    prism_beta,
    gdsc_beta
):
    """Classify PRISM/GDSC2 effect quadrants."""

    if (
        not np.isfinite(prism_beta)
        or
        not np.isfinite(gdsc_beta)
    ):
        return "missing"

    if (
        prism_beta < 0
        and
        gdsc_beta < 0
    ):
        return "concordant_sensitive"

    if (
        prism_beta > 0
        and
        gdsc_beta > 0
    ):
        return "concordant_resistant"

    if (
        prism_beta < 0
        and
        gdsc_beta > 0
    ):
        return "discordant_PRISM_sensitive_GDSC_resistant"

    if (
        prism_beta > 0
        and
        gdsc_beta < 0
    ):
        return "discordant_PRISM_resistant_GDSC_sensitive"

    return "zero_effect"


def expand_annotations(
    df,
    annotation_column,
    splitter
):
    """
    """

    rows = []

    for _, row in df.iterrows():

        annotations = splitter(
            row[
                annotation_column
            ]
        )

        for annotation in annotations:

            rows.append(
                {
                    "drug_key":
                        row["drug_key"],

                    "direction_class":
                        row["direction_class"],

                    "annotation":
                        annotation,
                }
            )

    return pd.DataFrame(
        rows
    )


def enrichment_analysis(
    df,
    annotation_column,
    splitter,
    group_name,
):
    """

    Universe:

    Test:
        group vs all other overlap drugs

    Alternative:
        greater
    """

    work = df.copy()


    work[
        "_annotations"
    ] = work[
        annotation_column
    ].apply(
        splitter
    )


    universe = work[
        work[
            "_annotations"
        ].apply(
            len
        ) > 0
    ].copy()


    if len(universe) == 0:
        return pd.DataFrame()


    group_mask = (
        universe[
            "direction_class"
        ]
        ==
        group_name
    )


    n_group = int(
        group_mask.sum()
    )

    n_other = int(
        (~group_mask).sum()
    )


    if (
        n_group == 0
        or
        n_other == 0
    ):
        return pd.DataFrame()


    all_categories = sorted(
        {
            category
            for annotations
            in universe[
                "_annotations"
            ]
            for category
            in annotations
        }
    )


    rows = []


    for category in all_categories:

        has_category = (
            universe[
                "_annotations"
            ]
            .apply(
                lambda x:
                    category in x
            )
        )


        universe_hits = int(
            has_category.sum()
        )


        group_hits = int(
            (
                has_category
                &
                group_mask
            ).sum()
        )


        if (
            universe_hits
            < MIN_CATEGORY_UNIVERSE_N
        ):
            continue


        if (
            group_hits
            < MIN_CATEGORY_GROUP_N
        ):
            continue


        other_hits = (
            universe_hits
            -
            group_hits
        )


        group_nonhits = (
            n_group
            -
            group_hits
        )


        other_nonhits = (
            n_other
            -
            other_hits
        )


        table = [
            [
                group_hits,
                group_nonhits,
            ],
            [
                other_hits,
                other_nonhits,
            ],
        ]


        odds_ratio, pvalue = (
            fisher_exact(
                table,
                alternative="greater"
            )
        )


        group_fraction = (
            group_hits
            /
            n_group
        )


        background_fraction = (
            other_hits
            /
            n_other
        )


        rows.append(
            {
                "annotation":
                    category,

                "group":
                    group_name,

                "group_hits":
                    group_hits,

                "group_n":
                    n_group,

                "other_hits":
                    other_hits,

                "other_n":
                    n_other,

                "universe_hits":
                    universe_hits,

                "universe_n":
                    len(universe),

                "group_fraction":
                    group_fraction,

                "other_fraction":
                    background_fraction,

                "odds_ratio":
                    odds_ratio,

                "p_value":
                    pvalue,
            }
        )


    result = pd.DataFrame(
        rows
    )


    if len(result) == 0:
        return result


    result[
        "FDR"
    ] = multipletests(
        result[
            "p_value"
        ],
        method="fdr_bh"
    )[1]


    result = result.sort_values(
        [
            "FDR",
            "p_value",
            "odds_ratio",
        ],
        ascending=[
            True,
            True,
            False,
        ]
    )


    return result


def safe_neglog10(x):
    """Calculate -log10 while avoiding zero values."""

    return -np.log10(
        np.clip(
            x,
            1e-300,
            1
        )
    )


# ============================================================
# Main
# ============================================================

def main():

    print(
        "\n"
        "============================================================\n"
        "ASTRA-Drug Step 11\n"
        "Cross-screen mechanistic analysis\n"
        "============================================================"
    )


    # --------------------------------------------------------
    # 0. Input check
    # --------------------------------------------------------

    for path in [
        OVERLAP_FILE,
        GDSC_RESULT_FILE,
        GDSC2_FILE,
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


    # --------------------------------------------------------
    # 1. Analysis plan freeze
    # --------------------------------------------------------

    plan_file = (
        OUT_DIR
        / "11_ANALYSIS_PLAN.txt"
    )


    with open(
        plan_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 11 analysis plan",
            file=f
        )

        print(
            "================================",
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
            "Exploratory mechanistic analysis of "
            "cross-screen ASTS-drug associations.",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Universe:",
            file=f
        )

        print(
            "Drugs evaluated in both PRISM and GDSC2.",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Concordant sensitive:",
            file=f
        )

        print(
            "PRISM beta < 0 AND GDSC2 beta < 0",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Concordant resistant:",
            file=f
        )

        print(
            "PRISM beta > 0 AND GDSC2 beta > 0",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Consensus score:",
            file=f
        )

        print(
            "mean(z(PRISM beta), z(GDSC2 beta))",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Annotations:",
            file=f
        )

        print(
            "GDSC2 PUTATIVE_TARGET and PATHWAY_NAME only",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Enrichment:",
            file=f
        )

        print(
            "one-sided Fisher exact test (greater), "
            "BH correction",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            f"Minimum category universe count: "
            f"{MIN_CATEGORY_UNIVERSE_N}",
            file=f
        )

        print(
            f"Minimum category group count: "
            f"{MIN_CATEGORY_GROUP_N}",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "This analysis is exploratory and will not "
            "modify the ASTS signature.",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Overlap checksum:",
            sha256_file(
                OVERLAP_FILE
            ),
            file=f
        )

        print(
            "GDSC2 checksum:",
            sha256_file(
                GDSC2_FILE
            ),
            file=f
        )


    # --------------------------------------------------------
    # 2. Read PRISM-GDSC2 overlap
    # --------------------------------------------------------

    print(
        "\nStep 1: Reading PRISM-GDSC2 overlap"
    )


    overlap = pd.read_csv(
        OVERLAP_FILE,
        sep="\t"
    )


    required_overlap = [
        "drug_key",
        "compound_name",
        "beta_ASTS_log2AUC_per_1SD",
        "DRUG_ID",
        "DRUG_NAME",
        "beta_ASTS_LNIC50_per_1SD",
    ]


    missing = [
        x
        for x in required_overlap
        if x not in overlap.columns
    ]


    if missing:

        raise ValueError(
            "Missing overlap columns: "
            + ", ".join(missing)
        )


    if (
        len(overlap)
        != EXPECTED_OVERLAP_N
    ):

        raise RuntimeError(
            f"Expected {EXPECTED_OVERLAP_N} overlap drugs, "
            f"but found {len(overlap)}."
        )


    if (
        overlap[
            "drug_key"
        ].duplicated().any()
    ):

        raise RuntimeError(
            "Duplicated drug_key in overlap file."
        )


    overlap[
        "DRUG_ID_norm"
    ] = overlap[
        "DRUG_ID"
    ].apply(
        normalize_drug_id
    )


    # --------------------------------------------------------
    # 3. Direction classification
    # --------------------------------------------------------

    print(
        "\nStep 2: Classifying cross-screen direction"
    )


    prism_beta_col = (
        "beta_ASTS_log2AUC_per_1SD"
    )

    gdsc_beta_col = (
        "beta_ASTS_LNIC50_per_1SD"
    )


    overlap[
        "direction_class"
    ] = [
        classify_direction(
            p,
            g
        )
        for p, g
        in zip(
            overlap[
                prism_beta_col
            ],
            overlap[
                gdsc_beta_col
            ],
        )
    ]


    # --------------------------------------------------------
    # 4. Consensus standardized effect
    # --------------------------------------------------------

    overlap[
        "PRISM_beta_z"
    ] = zscore(
        overlap[
            prism_beta_col
        ]
    )


    overlap[
        "GDSC2_beta_z"
    ] = zscore(
        overlap[
            gdsc_beta_col
        ]
    )


    overlap[
        "consensus_beta_z"
    ] = (
        overlap[
            "PRISM_beta_z"
        ]
        +
        overlap[
            "GDSC2_beta_z"
        ]
    ) / 2.0


    # --------------------------------------------------------
    # 5. Cross-screen statistics
    # --------------------------------------------------------

    rho, rho_p = spearmanr(
        overlap[
            prism_beta_col
        ],
        overlap[
            gdsc_beta_col
        ]
    )


    same_direction = (
        np.sign(
            overlap[
                prism_beta_col
            ]
        )
        ==
        np.sign(
            overlap[
                gdsc_beta_col
            ]
        )
    )


    n_same = int(
        same_direction.sum()
    )

    n_total = len(
        overlap
    )


    sign_fraction = (
        n_same
        /
        n_total
    )


    sign_test = binomtest(
        k=n_same,
        n=n_total,
        p=0.5,
        alternative="greater",
    )


    # --------------------------------------------------------
    # 6. Read GDSC2 drug annotations
    # --------------------------------------------------------

    print(
        "\nStep 3: Reading GDSC2 native annotations"
    )


    gdsc_raw = pd.read_excel(
        GDSC2_FILE,
        engine="openpyxl"
    )


    required_annotation = [
        "DRUG_ID",
        "DRUG_NAME",
        "PUTATIVE_TARGET",
        "PATHWAY_NAME",
    ]


    missing = [
        x
        for x in required_annotation
        if x not in gdsc_raw.columns
    ]


    if missing:

        raise ValueError(
            "Missing GDSC2 annotation columns: "
            + ", ".join(missing)
        )


    annotation = gdsc_raw[
        required_annotation
    ].copy()


    annotation[
        "DRUG_ID_norm"
    ] = annotation[
        "DRUG_ID"
    ].apply(
        normalize_drug_id
    )


    annotation = (
        annotation
        .groupby(
            "DRUG_ID_norm",
            as_index=False
        )
        .agg(
            GDSC2_annotation_name=(
                "DRUG_NAME",
                join_unique
            ),
            PUTATIVE_TARGET=(
                "PUTATIVE_TARGET",
                join_unique
            ),
            PATHWAY_NAME=(
                "PATHWAY_NAME",
                join_unique
            ),
        )
    )


    # --------------------------------------------------------
    # 7. Merge annotation
    # --------------------------------------------------------

    overlap = overlap.merge(
        annotation,
        on="DRUG_ID_norm",
        how="left",
        validate="1:1",
    )


    # --------------------------------------------------------
    # 8. Merge GDSC2 p/FDR information
    # --------------------------------------------------------

    gdsc_result = pd.read_csv(
        GDSC_RESULT_FILE,
        sep="\t"
    )


    gdsc_result[
        "DRUG_ID_norm"
    ] = gdsc_result[
        "DRUG_ID"
    ].apply(
        normalize_drug_id
    )


    gdsc_keep = [
        "DRUG_ID_norm",
    ]


    for col in [
        "p_ASTS_two_sided",
        "FDR_exploratory",
        "beta_ASTS_standardized",
        "partial_R2_ASTS",
    ]:

        if col in gdsc_result.columns:
            gdsc_keep.append(
                col
            )


    gdsc_result_small = (
        gdsc_result[
            gdsc_keep
        ]
        .drop_duplicates(
            "DRUG_ID_norm"
        )
    )


    gdsc_result_small = (
        gdsc_result_small.rename(
            columns={
                "p_ASTS_two_sided":
                    "GDSC2_p_ASTS",

                "FDR_exploratory":
                    "GDSC2_FDR",

                "beta_ASTS_standardized":
                    "GDSC2_beta_standardized",

                "partial_R2_ASTS":
                    "GDSC2_partial_R2",
            }
        )
    )


    overlap = overlap.merge(
        gdsc_result_small,
        on="DRUG_ID_norm",
        how="left",
        validate="1:1",
    )


    # --------------------------------------------------------
    # 9. Save fully annotated overlap
    # --------------------------------------------------------

    overlap = overlap.sort_values(
        "consensus_beta_z"
    )


    overlap.to_csv(
        OUT_DIR
        / "11_PRISM_GDSC2_overlap_annotated.tsv",
        sep="\t",
        index=False,
    )


    # --------------------------------------------------------
    # 10. Concordant groups
    # --------------------------------------------------------

    sensitive = overlap[
        overlap[
            "direction_class"
        ]
        ==
        "concordant_sensitive"
    ].copy()


    resistant = overlap[
        overlap[
            "direction_class"
        ]
        ==
        "concordant_resistant"
    ].copy()


    sensitive = sensitive.sort_values(
        "consensus_beta_z"
    )


    resistant = resistant.sort_values(
        "consensus_beta_z",
        ascending=False,
    )


    sensitive.to_csv(
        OUT_DIR
        / "11_concordant_sensitive_drugs.tsv",
        sep="\t",
        index=False,
    )


    resistant.to_csv(
        OUT_DIR
        / "11_concordant_resistant_drugs.tsv",
        sep="\t",
        index=False,
    )


    # --------------------------------------------------------
    # 11. Pathway enrichment
    # --------------------------------------------------------

    print(
        "\nStep 4: Pathway enrichment"
    )


    pathway_sensitive = (
        enrichment_analysis(
            overlap,
            "PATHWAY_NAME",
            split_pathways,
            "concordant_sensitive",
        )
    )


    pathway_resistant = (
        enrichment_analysis(
            overlap,
            "PATHWAY_NAME",
            split_pathways,
            "concordant_resistant",
        )
    )


    pathway_sensitive.to_csv(
        OUT_DIR
        / "11_pathway_enrichment_sensitive.tsv",
        sep="\t",
        index=False,
    )


    pathway_resistant.to_csv(
        OUT_DIR
        / "11_pathway_enrichment_resistant.tsv",
        sep="\t",
        index=False,
    )


    # --------------------------------------------------------
    # 12. Target enrichment
    # --------------------------------------------------------

    print(
        "\nStep 5: Target enrichment"
    )


    target_sensitive = (
        enrichment_analysis(
            overlap,
            "PUTATIVE_TARGET",
            split_targets,
            "concordant_sensitive",
        )
    )


    target_resistant = (
        enrichment_analysis(
            overlap,
            "PUTATIVE_TARGET",
            split_targets,
            "concordant_resistant",
        )
    )


    target_sensitive.to_csv(
        OUT_DIR
        / "11_target_enrichment_sensitive.tsv",
        sep="\t",
        index=False,
    )


    target_resistant.to_csv(
        OUT_DIR
        / "11_target_enrichment_resistant.tsv",
        sep="\t",
        index=False,
    )


    # --------------------------------------------------------
    # 13. Figure 1
    # Cross-screen scatter
    # --------------------------------------------------------

    fig, ax = plt.subplots(
        figsize=(7, 7)
    )


    ax.scatter(
        overlap[
            prism_beta_col
        ],
        overlap[
            gdsc_beta_col
        ],
        s=28,
        alpha=0.7,
    )


    ax.axhline(
        0,
        linestyle="--",
    )


    ax.axvline(
        0,
        linestyle="--",
    )


    label_df = (
        overlap
        .assign(
            abs_consensus=lambda x:
                x[
                    "consensus_beta_z"
                ].abs()
        )
        .sort_values(
            "abs_consensus",
            ascending=False
        )
        .head(10)
    )


    for _, row in label_df.iterrows():

        ax.annotate(
            row[
                "compound_name"
            ],
            (
                row[
                    prism_beta_col
                ],
                row[
                    gdsc_beta_col
                ],
            ),
            fontsize=7,
        )


    ax.set_xlabel(
        "PRISM ASTS beta\n"
        "(negative = ASTS-high more sensitive)"
    )


    ax.set_ylabel(
        "GDSC2 ASTS beta\n"
        "(negative = ASTS-high more sensitive)"
    )


    ax.set_title(
        "Cross-screen ASTS drug effects\n"
        f"Spearman r={rho:.3f}, "
        f"sign concordance={sign_fraction:.1%}"
    )


    fig.tight_layout()


    fig.savefig(
        FIG_DIR
        / "11_01_PRISM_vs_GDSC2_quadrants.png",
        dpi=200,
    )


    plt.close(fig)


    # --------------------------------------------------------
    # 14. Figure 2
    # Top concordant consensus drugs
    # --------------------------------------------------------

    concordant = overlap[
        overlap[
            "direction_class"
        ].isin(
            [
                "concordant_sensitive",
                "concordant_resistant",
            ]
        )
    ].copy()


    concordant[
        "abs_consensus"
    ] = concordant[
        "consensus_beta_z"
    ].abs()


    top = (
        concordant
        .sort_values(
            "abs_consensus",
            ascending=False
        )
        .head(
            TOP_N_DRUGS
        )
        .sort_values(
            "consensus_beta_z"
        )
    )


    if len(top) > 0:

        fig, ax = plt.subplots(
            figsize=(8, 7)
        )


        y = np.arange(
            len(top)
        )


        ax.barh(
            y,
            top[
                "consensus_beta_z"
            ]
        )


        ax.set_yticks(
            y
        )


        ax.set_yticklabels(
            top[
                "compound_name"
            ]
        )


        ax.axvline(
            0,
            linestyle="--"
        )


        ax.set_xlabel(
            "Cross-screen consensus beta z-score\n"
            "negative = ASTS-high sensitivity"
        )


        ax.set_title(
            "Top concordant ASTS-associated drugs"
        )


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / "11_02_top_concordant_consensus_drugs.png",
            dpi=200,
        )


        plt.close(fig)


    # --------------------------------------------------------
    # 15. Enrichment plots
    # --------------------------------------------------------

    enrichment_sets = [
        (
            pathway_sensitive,
            "Sensitive pathway enrichment",
            "11_03_pathway_sensitive.png",
        ),
        (
            pathway_resistant,
            "Resistant pathway enrichment",
            "11_04_pathway_resistant.png",
        ),
        (
            target_sensitive,
            "Sensitive target enrichment",
            "11_05_target_sensitive.png",
        ),
        (
            target_resistant,
            "Resistant target enrichment",
            "11_06_target_resistant.png",
        ),
    ]


    for data, title, filename in enrichment_sets:

        if len(data) == 0:
            continue


        plot_data = (
            data
            .head(
                TOP_N_ENRICHMENT_PLOT
            )
            .copy()
        )


        plot_data[
            "minus_log10_FDR"
        ] = safe_neglog10(
            plot_data[
                "FDR"
            ]
        )


        plot_data = (
            plot_data
            .sort_values(
                "minus_log10_FDR"
            )
        )


        fig, ax = plt.subplots(
            figsize=(8, 6)
        )


        y = np.arange(
            len(plot_data)
        )


        ax.barh(
            y,
            plot_data[
                "minus_log10_FDR"
            ]
        )


        ax.set_yticks(
            y
        )


        ax.set_yticklabels(
            plot_data[
                "annotation"
            ]
        )


        ax.set_xlabel(
            "-log10(FDR)"
        )


        ax.set_title(
            title
        )


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / filename,
            dpi=200,
        )


        plt.close(fig)


    # --------------------------------------------------------
    # 16. Summary
    # --------------------------------------------------------

    counts = (
        overlap[
            "direction_class"
        ]
        .value_counts()
    )


    summary_file = (
        OUT_DIR
        / "11_summary.txt"
    )


    with open(
        summary_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 11 summary",
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
            f"PRISM-GDSC2 overlap drugs: "
            f"{len(overlap)}",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "Cross-screen effect concordance",
            file=f
        )

        print(
            "-------------------------------",
            file=f
        )

        print(
            f"Spearman r: {rho}",
            file=f
        )

        print(
            f"Spearman p: {rho_p}",
            file=f
        )

        print(
            f"Same-direction drugs: "
            f"{n_same}/{n_total}",
            file=f
        )

        print(
            f"Sign concordance: "
            f"{sign_fraction}",
            file=f
        )

        print(
            f"One-sided binomial p vs 50%: "
            f"{sign_test.pvalue}",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "Direction classes",
            file=f
        )

        print(
            "-----------------",
            file=f
        )

        for name, count in counts.items():

            print(
                f"{name}: {count}",
                file=f
            )


        print(
            "",
            file=f
        )


        print(
            "Top concordant sensitive drugs",
            file=f
        )

        print(
            "-------------------------------",
            file=f
        )

        for _, row in sensitive.head(15).iterrows():

            print(
                f"{row['compound_name']}\t"
                f"consensus_z="
                f"{row['consensus_beta_z']:.4f}\t"
                f"PRISM_beta="
                f"{row[prism_beta_col]:.6f}\t"
                f"GDSC2_beta="
                f"{row[gdsc_beta_col]:.6f}",
                file=f
            )


        print(
            "",
            file=f
        )


        print(
            "Top concordant resistant drugs",
            file=f
        )

        print(
            "-------------------------------",
            file=f
        )

        for _, row in resistant.head(15).iterrows():

            print(
                f"{row['compound_name']}\t"
                f"consensus_z="
                f"{row['consensus_beta_z']:.4f}\t"
                f"PRISM_beta="
                f"{row[prism_beta_col]:.6f}\t"
                f"GDSC2_beta="
                f"{row[gdsc_beta_col]:.6f}",
                file=f
            )


        print(
            "",
            file=f
        )


        print(
            "Pathway enrichment",
            file=f
        )

        print(
            "------------------",
            file=f
        )


        for label, data in [
            (
                "concordant_sensitive",
                pathway_sensitive,
            ),
            (
                "concordant_resistant",
                pathway_resistant,
            ),
        ]:

            print(
                "",
                label,
                sep="",
                file=f
            )

            if len(data) == 0:

                print(
                    "No categories tested.",
                    file=f
                )

            else:

                print(
                    f"Categories tested: {len(data)}",
                    file=f
                )

                print(
                    f"FDR < {FDR_ALPHA}: "
                    f"{(data['FDR'] < FDR_ALPHA).sum()}",
                    file=f
                )

                for _, row in data.head(10).iterrows():

                    print(
                        f"{row['annotation']}\t"
                        f"hits={row['group_hits']}/"
                        f"{row['group_n']}\t"
                        f"OR={row['odds_ratio']:.4f}\t"
                        f"p={row['p_value']:.6g}\t"
                        f"FDR={row['FDR']:.6g}",
                        file=f
                    )


        print(
            "",
            file=f
        )


        print(
            "Target enrichment",
            file=f
        )

        print(
            "-----------------",
            file=f
        )


        for label, data in [
            (
                "concordant_sensitive",
                target_sensitive,
            ),
            (
                "concordant_resistant",
                target_resistant,
            ),
        ]:

            print(
                "",
                label,
                sep="",
                file=f
            )

            if len(data) == 0:

                print(
                    "No categories tested.",
                    file=f
                )

            else:

                print(
                    f"Categories tested: {len(data)}",
                    file=f
                )

                print(
                    f"FDR < {FDR_ALPHA}: "
                    f"{(data['FDR'] < FDR_ALPHA).sum()}",
                    file=f
                )

                for _, row in data.head(10).iterrows():

                    print(
                        f"{row['annotation']}\t"
                        f"hits={row['group_hits']}/"
                        f"{row['group_n']}\t"
                        f"OR={row['odds_ratio']:.4f}\t"
                        f"p={row['p_value']:.6g}\t"
                        f"FDR={row['FDR']:.6g}",
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
            "Step11 is exploratory mechanistic analysis.",
            file=f
        )

        print(
            "Enrichment results should not be interpreted "
            "as causal evidence.",
            file=f
        )

        print(
            "The ASTS signature remains frozen.",
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
            "Overlap:",
            sha256_file(
                OVERLAP_FILE
            ),
            file=f
        )

        print(
            "GDSC2:",
            sha256_file(
                GDSC2_FILE
            ),
            file=f
        )


    print(
        "\n============================================================"
    )

    print(
        "Step 11 completed successfully"
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
