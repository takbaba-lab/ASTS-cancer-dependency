#!/usr/bin/env python3

"""
ASTRA-Drug
Step 12: Drug-class robustness analysis

----

    MEK / EGFR / PI3K-MTOR inhibitor resistance





    RTK signaling -> ASTS-high sensitive
    EGFR signaling -> ASTS-high resistant


    EGFR-targeting drugs
    non-EGFR RTK drugs


----

"""


# ============================================================
# Path / Parameter
# ============================================================

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]

INPUT_FILE = (
    PROJECT_DIR
    / "results/11_cross_screen_mechanism"
    / "11_PRISM_GDSC2_overlap_annotated.tsv"
)

OUT_DIR = (
    PROJECT_DIR
    / "results/12_drug_class_robustness"
)

FIG_DIR = OUT_DIR / "figures"


# ------------------------------------------------------------
# Analysis parameters
# ------------------------------------------------------------

EXPECTED_DRUG_N = 108

RANDOM_SEED = 20260926

N_PERMUTATIONS = 100000

MIN_CLASS_N = 3

FDR_ALPHA = 0.05


# ============================================================
# Imports
# ============================================================

import hashlib
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from scipy.stats import (
    fisher_exact,
    binomtest,
    spearmanr,
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


def text_contains(
    series,
    pattern
):
    """
    """

    return (
        series
        .fillna("")
        .astype(str)
        .str.contains(
            pattern,
            case=False,
            regex=True,
            na=False,
        )
    )


def classify_direction(
    prism_beta,
    gdsc_beta
):
    """Classify cross-screen effect directions."""

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


def permutation_mean_test(
    all_values,
    observed_indices,
    n_permutations,
    rng,
    alternative="greater",
):
    """

    """

    values = np.asarray(
        all_values,
        dtype=float
    )

    observed_indices = np.asarray(
        observed_indices,
        dtype=int
    )

    values = values[
        np.isfinite(values)
    ]

    class_n = len(
        observed_indices
    )

    observed_mean = np.mean(
        all_values[
            observed_indices
        ]
    )

    random_means = np.empty(
        n_permutations,
        dtype=float
    )

    n_total = len(
        all_values
    )

    for i in range(
        n_permutations
    ):

        idx = rng.choice(
            n_total,
            size=class_n,
            replace=False,
        )

        random_means[i] = np.mean(
            all_values[idx]
        )


    if alternative == "greater":

        extreme = (
            random_means
            >= observed_mean
        ).sum()

    elif alternative == "less":

        extreme = (
            random_means
            <= observed_mean
        ).sum()

    else:

        extreme = (
            np.abs(random_means)
            >= abs(observed_mean)
        ).sum()


    p = (
        extreme + 1
    ) / (
        n_permutations + 1
    )


    return {
        "observed_mean":
            observed_mean,

        "permutation_p":
            p,
    }


def analyze_class(
    df,
    class_name,
    member_mask,
    rng,
):
    """
    """

    members = df[
        member_mask
    ].copy()

    n_class = len(
        members
    )

    if n_class < MIN_CLASS_N:
        return None, None


    # --------------------------------------------------------
    # Basic counts
    # --------------------------------------------------------

    n_concordant_resistant = int(
        (
            members[
                "direction_class"
            ]
            ==
            "concordant_resistant"
        ).sum()
    )

    n_concordant_sensitive = int(
        (
            members[
                "direction_class"
            ]
            ==
            "concordant_sensitive"
        ).sum()
    )

    n_same_direction = int(
        (
            np.sign(
                members[
                    "PRISM_beta"
                ]
            )
            ==
            np.sign(
                members[
                    "GDSC2_beta"
                ]
            )
        ).sum()
    )


    # --------------------------------------------------------
    # Fisher:
    # --------------------------------------------------------

    other = df[
        ~member_mask
    ]


    a = n_concordant_resistant

    b = (
        n_class
        -
        n_concordant_resistant
    )

    c = int(
        (
            other[
                "direction_class"
            ]
            ==
            "concordant_resistant"
        ).sum()
    )

    d = len(other) - c


    fisher_or, fisher_p = fisher_exact(
        [
            [a, b],
            [c, d],
        ],
        alternative="greater",
    )


    # --------------------------------------------------------
    # Mean standardized effects
    # --------------------------------------------------------

    prism_mean = members[
        "PRISM_z"
    ].mean()

    gdsc_mean = members[
        "GDSC2_z"
    ].mean()

    consensus_mean = members[
        "consensus_z"
    ].mean()


    # --------------------------------------------------------
    # Permutation test
    # resistant hypothesis = positive effect
    # --------------------------------------------------------

    member_indices = np.flatnonzero(
        member_mask.to_numpy()
    )


    prism_perm = permutation_mean_test(
        df[
            "PRISM_z"
        ].to_numpy(),
        member_indices,
        N_PERMUTATIONS,
        rng,
        alternative="greater",
    )


    gdsc_perm = permutation_mean_test(
        df[
            "GDSC2_z"
        ].to_numpy(),
        member_indices,
        N_PERMUTATIONS,
        rng,
        alternative="greater",
    )


    consensus_perm = permutation_mean_test(
        df[
            "consensus_z"
        ].to_numpy(),
        member_indices,
        N_PERMUTATIONS,
        rng,
        alternative="greater",
    )


    # --------------------------------------------------------
    # Leave-one-drug-out
    # --------------------------------------------------------

    loo_rows = []


    for remove_idx, remove_row in members.iterrows():

        remaining = members.drop(
            index=remove_idx
        )


        if len(remaining) < 2:
            continue


        loo_rows.append(
            {
                "class":
                    class_name,

                "removed_drug":
                    remove_row[
                        "compound_name"
                    ],

                "remaining_n":
                    len(remaining),

                "mean_PRISM_z":
                    remaining[
                        "PRISM_z"
                    ].mean(),

                "mean_GDSC2_z":
                    remaining[
                        "GDSC2_z"
                    ].mean(),

                "mean_consensus_z":
                    remaining[
                        "consensus_z"
                    ].mean(),

                "concordant_resistant_fraction":
                    (
                        remaining[
                            "direction_class"
                        ]
                        ==
                        "concordant_resistant"
                    ).mean(),
            }
        )


    loo = pd.DataFrame(
        loo_rows
    )


    if len(loo) > 0:

        min_loo_consensus = loo[
            "mean_consensus_z"
        ].min()

        max_loo_consensus = loo[
            "mean_consensus_z"
        ].max()

        min_loo_resistant_fraction = loo[
            "concordant_resistant_fraction"
        ].min()

    else:

        min_loo_consensus = np.nan
        max_loo_consensus = np.nan
        min_loo_resistant_fraction = np.nan


    summary = {
        "class":
            class_name,

        "n_drugs":
            n_class,

        "n_concordant_resistant":
            n_concordant_resistant,

        "n_concordant_sensitive":
            n_concordant_sensitive,

        "n_same_direction":
            n_same_direction,

        "same_direction_fraction":
            n_same_direction / n_class,

        "concordant_resistant_fraction":
            n_concordant_resistant / n_class,

        "mean_PRISM_z":
            prism_mean,

        "mean_GDSC2_z":
            gdsc_mean,

        "mean_consensus_z":
            consensus_mean,

        "fisher_OR_resistant":
            fisher_or,

        "fisher_p_resistant":
            fisher_p,

        "perm_p_PRISM":
            prism_perm[
                "permutation_p"
            ],

        "perm_p_GDSC2":
            gdsc_perm[
                "permutation_p"
            ],

        "perm_p_consensus":
            consensus_perm[
                "permutation_p"
            ],

        "min_LOO_consensus_z":
            min_loo_consensus,

        "max_LOO_consensus_z":
            max_loo_consensus,

        "min_LOO_resistant_fraction":
            min_loo_resistant_fraction,
    }


    members = members.copy()

    members[
        "class"
    ] = class_name


    return (
        summary,
        (
            members,
            loo,
        ),
    )


# ============================================================
# Main
# ============================================================

def main():

    print(
        "\n"
        "============================================================\n"
        "ASTRA-Drug Step 12\n"
        "Drug-class robustness analysis\n"
        "============================================================"
    )


    # --------------------------------------------------------
    # Input
    # --------------------------------------------------------

    if not INPUT_FILE.exists():

        raise FileNotFoundError(
            f"Input file not found:\n{INPUT_FILE}"
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
    # Read data
    # --------------------------------------------------------

    df = pd.read_csv(
        INPUT_FILE,
        sep="\t"
    )


    if len(df) != EXPECTED_DRUG_N:

        raise RuntimeError(
            f"Expected {EXPECTED_DRUG_N} drugs, "
            f"but found {len(df)}."
        )


    # --------------------------------------------------------
    # Required columns
    # --------------------------------------------------------

    prism_col = (
        "beta_ASTS_log2AUC_per_1SD"
    )

    gdsc_col = (
        "beta_ASTS_LNIC50_per_1SD"
    )


    required = [
        "compound_name",
        prism_col,
        gdsc_col,
        "PUTATIVE_TARGET",
        "PATHWAY_NAME",
    ]


    missing = [
        x
        for x in required
        if x not in df.columns
    ]


    if missing:

        raise ValueError(
            "Missing columns: "
            + ", ".join(missing)
        )


    # --------------------------------------------------------
    # Standardized effect
    # --------------------------------------------------------

    df[
        "PRISM_beta"
    ] = pd.to_numeric(
        df[
            prism_col
        ],
        errors="coerce"
    )


    df[
        "GDSC2_beta"
    ] = pd.to_numeric(
        df[
            gdsc_col
        ],
        errors="coerce"
    )


    df[
        "PRISM_z"
    ] = zscore(
        df[
            "PRISM_beta"
        ]
    )


    df[
        "GDSC2_z"
    ] = zscore(
        df[
            "GDSC2_beta"
        ]
    )


    df[
        "consensus_z"
    ] = (
        df[
            "PRISM_z"
        ]
        +
        df[
            "GDSC2_z"
        ]
    ) / 2.0


    df[
        "direction_class"
    ] = [
        classify_direction(
            p,
            g
        )
        for p, g
        in zip(
            df[
                "PRISM_beta"
            ],
            df[
                "GDSC2_beta"
            ],
        )
    ]


    # ========================================================
    # Predefined Step11-derived classes
    # ========================================================

    # --------------------------------------------------------
    # MEK target
    # --------------------------------------------------------

    mek_mask = (
        text_contains(
            df[
                "PUTATIVE_TARGET"
            ],
            r"(^|[,;/ ])MEK1($|[,;/ ])"
        )
        |
        text_contains(
            df[
                "PUTATIVE_TARGET"
            ],
            r"(^|[,;/ ])MEK2($|[,;/ ])"
        )
    )


    # --------------------------------------------------------
    # EGFR target
    # --------------------------------------------------------

    egfr_mask = text_contains(
        df[
            "PUTATIVE_TARGET"
        ],
        r"(^|[,;/ ])EGFR($|[,;/ ])"
    )


    # --------------------------------------------------------
    # PI3K / MTOR pathway
    # --------------------------------------------------------

    pi3k_mtor_mask = text_contains(
        df[
            "PATHWAY_NAME"
        ],
        r"PI3K/MTOR signaling"
    )


    # --------------------------------------------------------
    # EGFR pathway
    # --------------------------------------------------------

    egfr_pathway_mask = text_contains(
        df[
            "PATHWAY_NAME"
        ],
        r"EGFR signaling"
    )


    # --------------------------------------------------------
    # ERK / MAPK pathway
    # --------------------------------------------------------

    erk_mask = text_contains(
        df[
            "PATHWAY_NAME"
        ],
        r"ERK MAPK signaling"
    )


    # --------------------------------------------------------
    # RTK pathway
    # --------------------------------------------------------

    rtk_mask = text_contains(
        df[
            "PATHWAY_NAME"
        ],
        r"RTK signaling"
    )


    # --------------------------------------------------------
    # RTK but not EGFR-targeting
    #
    # --------------------------------------------------------

    non_egfr_rtk_mask = (
        rtk_mask
        &
        ~egfr_mask
    )


    classes = {
        "MEK_target":
            mek_mask,

        "EGFR_target":
            egfr_mask,

        "PI3K_MTOR_pathway":
            pi3k_mtor_mask,

        "EGFR_pathway":
            egfr_pathway_mask,

        "ERK_MAPK_pathway":
            erk_mask,

        "RTK_pathway":
            rtk_mask,

        "RTK_non_EGFR":
            non_egfr_rtk_mask,
    }


    # --------------------------------------------------------
    # Save class membership matrix
    # --------------------------------------------------------

    membership = df[
        [
            "compound_name",
            "PRISM_beta",
            "GDSC2_beta",
            "PRISM_z",
            "GDSC2_z",
            "consensus_z",
            "direction_class",
            "PUTATIVE_TARGET",
            "PATHWAY_NAME",
        ]
    ].copy()


    for class_name, mask in classes.items():

        membership[
            class_name
        ] = mask.astype(int)


    membership.to_csv(
        OUT_DIR
        / "12_drug_class_membership.tsv",
        sep="\t",
        index=False,
    )


    # --------------------------------------------------------
    # Analysis
    # --------------------------------------------------------

    rng = np.random.default_rng(
        RANDOM_SEED
    )


    summary_rows = []

    member_tables = []

    loo_tables = []


    for class_name, mask in classes.items():

        print(
            f"\nAnalyzing: {class_name}"
        )

        print(
            "N:",
            int(mask.sum())
        )


        result = analyze_class(
            df,
            class_name,
            mask,
            rng,
        )


        summary, details = result


        if summary is None:

            print(
                "Skipped: class too small."
            )

            continue


        members, loo = details


        summary_rows.append(
            summary
        )

        member_tables.append(
            members
        )


        if len(loo) > 0:

            loo_tables.append(
                loo
            )


    summary_df = pd.DataFrame(
        summary_rows
    )


    # --------------------------------------------------------
    # Multiple testing correction
    #
    # --------------------------------------------------------

    if len(summary_df) > 0:

        summary_df[
            "FDR_fisher_resistant"
        ] = multipletests(
            summary_df[
                "fisher_p_resistant"
            ],
            method="fdr_bh"
        )[1]


        summary_df[
            "FDR_perm_consensus"
        ] = multipletests(
            summary_df[
                "perm_p_consensus"
            ],
            method="fdr_bh"
        )[1]


        summary_df = summary_df.sort_values(
            "mean_consensus_z",
            ascending=False,
        )


    # --------------------------------------------------------
    # Save tables
    # --------------------------------------------------------

    summary_df.to_csv(
        OUT_DIR
        / "12_class_summary.tsv",
        sep="\t",
        index=False,
    )


    if len(member_tables) > 0:

        member_df = pd.concat(
            member_tables,
            ignore_index=True,
        )

    else:

        member_df = pd.DataFrame()


    member_df.to_csv(
        OUT_DIR
        / "12_class_member_drugs.tsv",
        sep="\t",
        index=False,
    )


    if len(loo_tables) > 0:

        loo_df = pd.concat(
            loo_tables,
            ignore_index=True,
        )

    else:

        loo_df = pd.DataFrame()


    loo_df.to_csv(
        OUT_DIR
        / "12_leave_one_drug_out.tsv",
        sep="\t",
        index=False,
    )


    # ========================================================
    # Specific RTK decomposition
    # ========================================================

    rtk_subgroups = []


    for _, row in df[
        rtk_mask
    ].iterrows():

        if egfr_mask.loc[
            row.name
        ]:

            subgroup = "RTK_EGFR"

        else:

            subgroup = "RTK_non_EGFR"


        rtk_subgroups.append(
            {
                "compound_name":
                    row[
                        "compound_name"
                    ],

                "subgroup":
                    subgroup,

                "PRISM_beta":
                    row[
                        "PRISM_beta"
                    ],

                "GDSC2_beta":
                    row[
                        "GDSC2_beta"
                    ],

                "consensus_z":
                    row[
                        "consensus_z"
                    ],

                "direction_class":
                    row[
                        "direction_class"
                    ],

                "PUTATIVE_TARGET":
                    row[
                        "PUTATIVE_TARGET"
                    ],

                "PATHWAY_NAME":
                    row[
                        "PATHWAY_NAME"
                    ],
            }
        )


    rtk_df = pd.DataFrame(
        rtk_subgroups
    )


    rtk_df.to_csv(
        OUT_DIR
        / "12_RTK_decomposition.tsv",
        sep="\t",
        index=False,
    )


    # ========================================================
    # Figure 1: class mean effects
    # ========================================================

    if len(summary_df) > 0:

        plot_df = summary_df.sort_values(
            "mean_consensus_z"
        )


        fig, ax = plt.subplots(
            figsize=(8, 6)
        )


        y = np.arange(
            len(plot_df)
        )


        ax.barh(
            y,
            plot_df[
                "mean_consensus_z"
            ]
        )


        ax.set_yticks(
            y
        )


        ax.set_yticklabels(
            plot_df[
                "class"
            ]
        )


        ax.axvline(
            0,
            linestyle="--"
        )


        ax.set_xlabel(
            "Mean cross-screen consensus z-score\n"
            "positive = ASTS-high resistance"
        )


        ax.set_title(
            "Step 12: Drug-class ASTS effects"
        )


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / "12_01_class_mean_consensus.png",
            dpi=200,
        )


        plt.close(fig)


    # ========================================================
    # Figure 2: drug-level effects for selected resistance classes
    # ========================================================

    selected_classes = [
        "MEK_target",
        "EGFR_target",
        "PI3K_MTOR_pathway",
    ]


    selected = membership[
        membership[
            selected_classes
        ].sum(
            axis=1
        ) > 0
    ].copy()


    if len(selected) > 0:

        selected = selected.sort_values(
            "consensus_z"
        )


        fig, ax = plt.subplots(
            figsize=(8, 8)
        )


        y = np.arange(
            len(selected)
        )


        ax.barh(
            y,
            selected[
                "consensus_z"
            ]
        )


        ax.set_yticks(
            y
        )


        ax.set_yticklabels(
            selected[
                "compound_name"
            ]
        )


        ax.axvline(
            0,
            linestyle="--"
        )


        ax.set_xlabel(
            "Cross-screen consensus z-score\n"
            "positive = ASTS-high resistance"
        )


        ax.set_title(
            "MEK / EGFR / PI3K-MTOR class drugs"
        )


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / "12_02_resistance_class_drugs.png",
            dpi=200,
        )


        plt.close(fig)


    # ========================================================
    # Figure 3: RTK decomposition
    # ========================================================

    if len(rtk_df) > 0:

        plot_rtk = rtk_df.sort_values(
            "consensus_z"
        )


        fig, ax = plt.subplots(
            figsize=(8, 7)
        )


        y = np.arange(
            len(plot_rtk)
        )


        ax.barh(
            y,
            plot_rtk[
                "consensus_z"
            ]
        )


        ax.set_yticks(
            y
        )


        labels = [
            f"{drug} [{group}]"
            for drug, group
            in zip(
                plot_rtk[
                    "compound_name"
                ],
                plot_rtk[
                    "subgroup"
                ],
            )
        ]


        ax.set_yticklabels(
            labels
        )


        ax.axvline(
            0,
            linestyle="--"
        )


        ax.set_xlabel(
            "Cross-screen consensus z-score\n"
            "negative = sensitivity; positive = resistance"
        )


        ax.set_title(
            "RTK pathway decomposition"
        )


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / "12_03_RTK_decomposition.png",
            dpi=200,
        )


        plt.close(fig)


    # ========================================================
    # Summary
    # ========================================================

    summary_file = (
        OUT_DIR
        / "12_summary.txt"
    )


    with open(
        summary_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 12 summary",
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
            "Robustness analysis of Step11-derived "
            "drug-class signals.",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            f"Overlap drugs: {len(df)}",
            file=f
        )

        print(
            f"Permutations per class: "
            f"{N_PERMUTATIONS}",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "Class summary",
            file=f
        )

        print(
            "-------------",
            file=f
        )


        if len(summary_df) > 0:

            columns_to_print = [
                "class",
                "n_drugs",
                "n_concordant_resistant",
                "n_concordant_sensitive",
                "same_direction_fraction",
                "concordant_resistant_fraction",
                "mean_PRISM_z",
                "mean_GDSC2_z",
                "mean_consensus_z",
                "fisher_OR_resistant",
                "fisher_p_resistant",
                "FDR_fisher_resistant",
                "perm_p_consensus",
                "FDR_perm_consensus",
                "min_LOO_consensus_z",
                "min_LOO_resistant_fraction",
            ]


            print(
                summary_df[
                    columns_to_print
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
            "Interpretation guide",
            file=f
        )

        print(
            "--------------------",
            file=f
        )

        print(
            "Positive mean_consensus_z = "
            "ASTS-high resistance",
            file=f
        )

        print(
            "Negative mean_consensus_z = "
            "ASTS-high sensitivity",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "A robust resistance class should ideally show:",
            file=f
        )

        print(
            "  1. positive mean PRISM and GDSC2 z",
            file=f
        )

        print(
            "  2. high concordant-resistant fraction",
            file=f
        )

        print(
            "  3. low permutation p",
            file=f
        )

        print(
            "  4. positive min leave-one-out consensus z",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "Important:",
            file=f
        )

        print(
            "Step12 is exploratory robustness analysis, "
            "not an independent validation.",
            file=f
        )

        print(
            "The classes were motivated by Step11 results.",
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
            "Input checksum:",
            file=f
        )

        print(
            sha256_file(
                INPUT_FILE
            ),
            file=f
        )


    print(
        "\n============================================================"
    )

    print(
        "Step 12 completed successfully"
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
