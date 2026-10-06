#!/usr/bin/env python3

"""
ASTRA-Drug
Step 19: Pharmacological validation of mitochondrial vulnerability

----

    ASTS-high
        -> mitochondrial / OXPHOS genetic dependency increased




----------------------------

    drug response
      ~ ASTS
      + lineage
      + sex
      + proliferation
      + canonical drivers
      + E2F/G2M cell-cycle score
      + MYC score


negative ASTS beta

positive ASTS beta


----




Primary:
    direct / near-direct ETC-OXPHOS inhibitors

Secondary:
    Primary + broader mitochondrial perturbation drugs


----

    PRISM standardized beta
    GDSC2 standardized beta
    cross-screen consensus z



    PRISM n
    GDSC2 n



    = ASTS-high selective sensitivity

"""


# ============================================================
# 1. Path / Parameter
# ============================================================

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]


# ------------------------------------------------------------
# Step 15 cross-screen drug effects
# ------------------------------------------------------------

DRUG_EFFECT_FILE = (
    PROJECT_DIR
    / "results/15_cellcycle_confounder"
    / "15_PRISM_GDSC2_cellcycle_adjusted.tsv"
)


# ------------------------------------------------------------
# Output
# ------------------------------------------------------------

OUT_DIR = (
    PROJECT_DIR
    / "results/19_mito_drug_validation"
)


FIG_DIR = (
    OUT_DIR
    / "figures"
)


# ------------------------------------------------------------
# Primary effect columns
#
# Step15 stringent sensitivity model:
# driver + CellCycle + MYC
# ------------------------------------------------------------

PRISM_BETA_COLUMN = (
    "beta_plus_cellcycle_myc_PRISM"
)

GDSC2_BETA_COLUMN = (
    "beta_plus_cellcycle_myc_GDSC2"
)

CONSENSUS_COLUMN = (
    "consensus_z_plus_cellcycle_myc"
)

PRISM_N_COLUMN = "n_PRISM"

GDSC2_N_COLUMN = "n_GDSC2"


# ------------------------------------------------------------
# Matched permutation
# ------------------------------------------------------------

N_PERMUTATIONS = 10000

MATCH_POOL_SIZE = 50

RANDOM_SEED = 20260927

MIN_PRIMARY_DRUGS_FOR_INFERENCE = 3


# ============================================================
# 2. Mechanism-defined mitochondrial drug set
# ============================================================
#
# IMPORTANT:
#
# Tier 1:
# direct / near-direct ETC or ATP-synthesis inhibitors
#
# Tier 2:
# broader mitochondrial perturbation drugs
#
# ============================================================

MITO_DRUGS = [

    # ========================================================
    # Tier 1: direct / near-direct OXPHOS
    # ========================================================

    {
        "canonical_name": "IACS-010759",
        "tier": "tier1_direct_oxphos",
        "mechanism_group": "Complex_I",
        "aliases": [
            "IACS-010759",
            "IACS010759",
        ],
    },

    {
        "canonical_name": "BAY-87-2243",
        "tier": "tier1_direct_oxphos",
        "mechanism_group": "Complex_I",
        "aliases": [
            "BAY-87-2243",
            "BAY 87-2243",
            "BAY872243",
        ],
    },

    {
        "canonical_name": "rotenone",
        "tier": "tier1_direct_oxphos",
        "mechanism_group": "Complex_I",
        "aliases": [
            "rotenone",
        ],
    },

    {
        "canonical_name": "piericidin A",
        "tier": "tier1_direct_oxphos",
        "mechanism_group": "Complex_I",
        "aliases": [
            "piericidin A",
            "piericidin-A",
            "piericidin",
        ],
    },

    {
        "canonical_name": "antimycin A",
        "tier": "tier1_direct_oxphos",
        "mechanism_group": "Complex_III",
        "aliases": [
            "antimycin A",
            "antimycin-A",
            "antimycin",
        ],
    },

    {
        "canonical_name": "atovaquone",
        "tier": "tier1_direct_oxphos",
        "mechanism_group": "Complex_III",
        "aliases": [
            "atovaquone",
        ],
    },

    {
        "canonical_name": "oligomycin A",
        "tier": "tier1_direct_oxphos",
        "mechanism_group": "ATP_synthase",
        "aliases": [
            "oligomycin A",
            "oligomycin-A",
            "oligomycin",
        ],
    },


    # ========================================================
    # Tier 2: broader mitochondrial perturbation
    # ========================================================

    {
        "canonical_name": "phenformin",
        "tier": "tier2_broad_mito",
        "mechanism_group": "Biguanide_mito_respiration",
        "aliases": [
            "phenformin",
        ],
    },

    {
        "canonical_name": "metformin",
        "tier": "tier2_broad_mito",
        "mechanism_group": "Biguanide_mito_respiration",
        "aliases": [
            "metformin",
        ],
    },

    {
        "canonical_name": "elesclomol",
        "tier": "tier2_broad_mito",
        "mechanism_group": "Mitochondrial_copper_stress",
        "aliases": [
            "elesclomol",
        ],
    },

    {
        "canonical_name": "VLX600",
        "tier": "tier2_broad_mito",
        "mechanism_group": "Mitochondrial_respiration",
        "aliases": [
            "VLX600",
            "VLX-600",
        ],
    },

    {
        "canonical_name": "devimistat",
        "tier": "tier2_broad_mito",
        "mechanism_group": "Mitochondrial_metabolism",
        "aliases": [
            "devimistat",
            "CPI-613",
            "CPI613",
        ],
    },
]


# ============================================================
# 3. Imports
# ============================================================

import hashlib
import re

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from statsmodels.stats.multitest import multipletests


# ============================================================
# 4. Helper functions
# ============================================================

def normalize_drug_name(x):
    """

        BAY-87-2243
        BAY 87 2243
        BAY872243

    """

    if pd.isna(x):
        return ""

    return re.sub(
        r"[^a-z0-9]+",
        "",
        str(x).lower().strip()
    )


def sha256_file(path):
    """
    """

    h = hashlib.sha256()

    with open(
        path,
        "rb"
    ) as f:

        while True:

            block = f.read(
                1024 * 1024
            )

            if not block:
                break

            h.update(
                block
            )

    return h.hexdigest()


def build_candidate_table():
    """
    """

    rows = []

    for drug in MITO_DRUGS:

        rows.append(
            {
                "canonical_name":
                    drug[
                        "canonical_name"
                    ],

                "tier":
                    drug[
                        "tier"
                    ],

                "mechanism_group":
                    drug[
                        "mechanism_group"
                    ],

                "aliases":
                    ";".join(
                        drug[
                            "aliases"
                        ]
                    ),
            }
        )

    return pd.DataFrame(
        rows
    )


def build_alias_map():
    """
    normalized alias -> canonical drug name
    """

    alias_map = {}

    for drug in MITO_DRUGS:

        canonical = drug[
            "canonical_name"
        ]

        for alias in drug[
            "aliases"
        ]:

            key = normalize_drug_name(
                alias
            )

            if (
                key in alias_map
                and
                alias_map[
                    key
                ]
                !=
                canonical
            ):

                raise RuntimeError(
                    "Alias collision: "
                    f"{alias}"
                )

            alias_map[
                key
            ] = canonical

    return alias_map


def classify_direction(
    prism_beta,
    gdsc_beta
):
    """
    """

    if (
        prism_beta < 0
        and
        gdsc_beta < 0
    ):

        return (
            "concordant_sensitive"
        )

    if (
        prism_beta > 0
        and
        gdsc_beta > 0
    ):

        return (
            "concordant_resistant"
        )

    return "discordant"


def zscore_features(df):
    """
    """

    x = df.copy()

    for col in [
        PRISM_N_COLUMN,
        GDSC2_N_COLUMN,
    ]:

        log_col = (
            f"log1p_{col}"
        )

        z_col = (
            f"z_{log_col}"
        )

        x[
            log_col
        ] = np.log1p(
            pd.to_numeric(
                x[
                    col
                ],
                errors="coerce"
            )
        )


        sd = x[
            log_col
        ].std(
            ddof=1
        )


        if (
            not np.isfinite(sd)
            or
            sd == 0
        ):

            raise RuntimeError(
                "Cannot standardize matching column: "
                f"{col}"
            )


        x[
            z_col
        ] = (
            x[
                log_col
            ]
            -
            x[
                log_col
            ].mean()
        ) / sd

    return x


def matched_permutation(
    set_name,
    target_rows,
    background,
    rng
):
    """

    Statistic:
        mean consensus z

    Negative:
        ASTS-high greater sensitivity
    """

    if len(
        target_rows
    ) == 0:

        return None, None


    all_for_matching = pd.concat(
        [
            target_rows,
            background,
        ],
        ignore_index=True
    )


    all_for_matching = zscore_features(
        all_for_matching
    )


    target_keys = set(
        target_rows[
            "compound_name"
        ]
    )


    target_match = all_for_matching[
        all_for_matching[
            "compound_name"
        ].isin(
            target_keys
        )
    ].copy()


    control_match = all_for_matching[
        ~all_for_matching[
            "compound_name"
        ].isin(
            target_keys
        )
    ].copy()


    if len(
        control_match
    ) < len(
        target_match
    ):

        raise RuntimeError(
            "Too few background drugs."
        )


    feature_cols = [
        f"z_log1p_{PRISM_N_COLUMN}",
        f"z_log1p_{GDSC2_N_COLUMN}",
    ]


    control_matrix = control_match[
        feature_cols
    ].to_numpy(
        dtype=float
    )


    control_names = control_match[
        "compound_name"
    ].to_numpy(
        dtype=object
    )


    candidate_pools = {}


    # --------------------------------------------------------
    # Build nearest-neighbour pool for each mitochondrial drug
    # --------------------------------------------------------

    for _, row in target_match.iterrows():

        target_vector = row[
            feature_cols
        ].to_numpy(
            dtype=float
        )


        distance = np.sum(
            (
                control_matrix
                -
                target_vector[
                    None,
                    :
                ]
            )
            ** 2,
            axis=1
        )


        k = min(
            MATCH_POOL_SIZE,
            len(
                control_match
            )
        )


        indices = np.argpartition(
            distance,
            k - 1
        )[
            :k
        ]


        indices = indices[
            np.argsort(
                distance[
                    indices
                ]
            )
        ]


        candidate_pools[
            row[
                "compound_name"
            ]
        ] = [
            control_names[i]
            for i in indices
        ]


    background_indexed = (
        background
        .set_index(
            "compound_name"
        )
    )


    # --------------------------------------------------------
    # Observed statistics
    # --------------------------------------------------------

    observed = target_rows[
        CONSENSUS_COLUMN
    ].to_numpy(
        dtype=float
    )


    observed_mean = float(
        np.mean(
            observed
        )
    )


    observed_median = float(
        np.median(
            observed
        )
    )


    fraction_negative = float(
        np.mean(
            observed < 0
        )
    )


    # --------------------------------------------------------
    # Permutation
    # --------------------------------------------------------

    null_means = np.empty(
        N_PERMUTATIONS,
        dtype=float
    )


    target_names_original = list(
        target_rows[
            "compound_name"
        ]
    )


    all_background_names = list(
        background_indexed.index
    )


    for i in range(
        N_PERMUTATIONS
    ):

        target_order = (
            target_names_original.copy()
        )


        rng.shuffle(
            target_order
        )


        selected = set()

        sampled = []


        for target_name in target_order:

            pool = candidate_pools[
                target_name
            ]


            available = [
                x
                for x in pool
                if x not in selected
            ]


            if len(
                available
            ) > 0:

                chosen = rng.choice(
                    available
                )


            else:

                remaining = [
                    x
                    for x in all_background_names
                    if x not in selected
                ]


                chosen = rng.choice(
                    remaining
                )


            selected.add(
                chosen
            )

            sampled.append(
                chosen
            )


        null_means[
            i
        ] = (
            background_indexed
            .loc[
                sampled,
                CONSENSUS_COLUMN
            ]
            .mean()
        )


    # --------------------------------------------------------
    # Empirical test
    #
    # Primary alternative:
    # mitochondrial drugs are more negative
    # --------------------------------------------------------

    p_sensitive = (
        1
        +
        np.sum(
            null_means
            <=
            observed_mean
        )
    ) / (
        N_PERMUTATIONS
        +
        1
    )


    p_resistant = (
        1
        +
        np.sum(
            null_means
            >=
            observed_mean
        )
    ) / (
        N_PERMUTATIONS
        +
        1
    )


    p_two_sided = min(
        1.0,
        2
        *
        min(
            p_sensitive,
            p_resistant
        )
    )


    null_mean = float(
        np.mean(
            null_means
        )
    )


    null_sd = float(
        np.std(
            null_means,
            ddof=1
        )
    )


    matched_z = (
        (
            observed_mean
            -
            null_mean
        )
        /
        null_sd
        if null_sd > 0
        else np.nan
    )


    summary = {
        "drug_set":
            set_name,

        "n_drugs":
            len(
                target_rows
            ),

        "mean_consensus_z":
            observed_mean,

        "median_consensus_z":
            observed_median,

        "fraction_negative":
            fraction_negative,

        "n_concordant_sensitive":
            int(
                (
                    target_rows[
                        "cross_screen_direction"
                    ]
                    ==
                    "concordant_sensitive"
                ).sum()
            ),

        "n_concordant_resistant":
            int(
                (
                    target_rows[
                        "cross_screen_direction"
                    ]
                    ==
                    "concordant_resistant"
                ).sum()
            ),

        "matched_null_mean":
            null_mean,

        "matched_null_sd":
            null_sd,

        "matched_z":
            matched_z,

        "p_ASTShigh_sensitive":
            p_sensitive,

        "p_ASTShigh_resistant":
            p_resistant,

        "p_two_sided":
            p_two_sided,
    }


    return (
        summary,
        null_means
    )


# ============================================================
# 5. Main
# ============================================================

def main():

    print(
        "\n"
        "============================================================\n"
        "ASTRA-Drug Step 19\n"
        "Pharmacological validation of mitochondrial vulnerability\n"
        "============================================================"
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
    # Freeze mitochondrial candidate set BEFORE reading effects
    # ========================================================

    print(
        "\nStep 1: Freezing mitochondrial drug set"
    )


    candidate_table = (
        build_candidate_table()
    )


    freeze_file = (
        OUT_DIR
        / "19_MITO_DRUG_SET_FREEZE.tsv"
    )


    candidate_table.to_csv(
        freeze_file,
        sep="\t",
        index=False
    )


    plan_file = (
        OUT_DIR
        / "19_ANALYSIS_PLAN_FREEZE.txt"
    )


    with open(
        plan_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 19 analysis plan",
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
            "Primary drug set:",
            file=f
        )

        print(
            "tier1_direct_oxphos",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Secondary drug set:",
            file=f
        )

        print(
            "tier1_direct_oxphos + tier2_broad_mito",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Primary drug response model:",
            file=f
        )

        print(
            "Step15 plus_cellcycle_myc model",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            "Negative consensus z = ASTS-high sensitivity",
            file=f
        )

        print(
            "",
            file=f
        )

        print(
            f"Matched permutations: "
            f"{N_PERMUTATIONS}",
            file=f
        )

        print(
            f"Match pool size: "
            f"{MATCH_POOL_SIZE}",
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
            "This is a post-hoc mechanism-defined analysis "
            "performed after mitochondrial dependency was "
            "identified in Steps17-18.",
            file=f
        )


    print(
        "Freeze file:",
        freeze_file
    )


    # ========================================================
    # Step 2:
    # Read drug effects
    # ========================================================

    print(
        "\nStep 2: Reading Step15 cross-screen drug effects"
    )


    if not DRUG_EFFECT_FILE.exists():

        raise FileNotFoundError(
            f"Drug effect file not found:\n"
            f"{DRUG_EFFECT_FILE}"
        )


    effects = pd.read_csv(
        DRUG_EFFECT_FILE,
        sep="\t"
    )


    required_columns = [
        "compound_name",
        PRISM_BETA_COLUMN,
        GDSC2_BETA_COLUMN,
        CONSENSUS_COLUMN,
        PRISM_N_COLUMN,
        GDSC2_N_COLUMN,
    ]


    missing = [
        col
        for col in required_columns
        if col not in effects.columns
    ]


    if len(
        missing
    ) > 0:

        raise RuntimeError(
            "Required Step15 columns missing:\n"
            +
            ", ".join(
                missing
            )
            +
            "\n\nAvailable columns:\n"
            +
            "\n".join(
                effects.columns
            )
        )


    # Numeric conversion
    for col in [
        PRISM_BETA_COLUMN,
        GDSC2_BETA_COLUMN,
        CONSENSUS_COLUMN,
        PRISM_N_COLUMN,
        GDSC2_N_COLUMN,
    ]:

        effects[
            col
        ] = pd.to_numeric(
            effects[
                col
            ],
            errors="coerce"
        )


    effects = effects.dropna(
        subset=required_columns
    ).copy()


    effects[
        "drug_name_norm"
    ] = effects[
        "compound_name"
    ].apply(
        normalize_drug_name
    )


    # ========================================================
    # Step 3:
    # Match predefined candidates
    # ========================================================

    print(
        "\nStep 3: Auditing candidate availability"
    )


    alias_map = (
        build_alias_map()
    )


    effects[
        "matched_canonical"
    ] = effects[
        "drug_name_norm"
    ].map(
        alias_map
    )


    availability_rows = []

    candidate_effect_rows = []


    for drug in MITO_DRUGS:

        canonical = drug[
            "canonical_name"
        ]


        matches = effects[
            effects[
                "matched_canonical"
            ]
            ==
            canonical
        ].copy()


        if len(
            matches
        ) == 0:

            availability_rows.append(
                {
                    "canonical_name":
                        canonical,

                    "tier":
                        drug[
                            "tier"
                        ],

                    "mechanism_group":
                        drug[
                            "mechanism_group"
                        ],

                    "status":
                        "not_found",

                    "matched_name":
                        "",

                    "n_matches":
                        0,
                }
            )

            continue


        # ----------------------------------------------------
        # If multiple aliases/records exist:
        # deterministically retain highest cross-screen coverage
        # ----------------------------------------------------

        matches[
            "_coverage"
        ] = np.minimum(
            matches[
                PRISM_N_COLUMN
            ],
            matches[
                GDSC2_N_COLUMN
            ]
        )


        matches = matches.sort_values(
            [
                "_coverage",
                PRISM_N_COLUMN,
                GDSC2_N_COLUMN,
            ],
            ascending=False
        )


        selected = matches.iloc[
            0
        ]


        availability_rows.append(
            {
                "canonical_name":
                    canonical,

                "tier":
                    drug[
                        "tier"
                    ],

                "mechanism_group":
                    drug[
                        "mechanism_group"
                    ],

                "status":
                    "found",

                "matched_name":
                    selected[
                        "compound_name"
                    ],

                "n_matches":
                    len(
                        matches
                    ),
            }
        )


        candidate_effect_rows.append(
            {
                "canonical_name":
                    canonical,

                "tier":
                    drug[
                        "tier"
                    ],

                "mechanism_group":
                    drug[
                        "mechanism_group"
                    ],

                "compound_name":
                    selected[
                        "compound_name"
                    ],

                PRISM_N_COLUMN:
                    selected[
                        PRISM_N_COLUMN
                    ],

                GDSC2_N_COLUMN:
                    selected[
                        GDSC2_N_COLUMN
                    ],

                PRISM_BETA_COLUMN:
                    selected[
                        PRISM_BETA_COLUMN
                    ],

                GDSC2_BETA_COLUMN:
                    selected[
                        GDSC2_BETA_COLUMN
                    ],

                CONSENSUS_COLUMN:
                    selected[
                        CONSENSUS_COLUMN
                    ],

                "cross_screen_direction":
                    classify_direction(
                        selected[
                            PRISM_BETA_COLUMN
                        ],
                        selected[
                            GDSC2_BETA_COLUMN
                        ]
                    ),
            }
        )


    availability = pd.DataFrame(
        availability_rows
    )


    availability.to_csv(
        OUT_DIR
        / "19_mito_candidate_availability.tsv",
        sep="\t",
        index=False
    )


    candidate_effects = pd.DataFrame(
        candidate_effect_rows
    )


    candidate_effects.to_csv(
        OUT_DIR
        / "19_mito_candidate_effects.tsv",
        sep="\t",
        index=False
    )


    print(
        availability.to_string(
            index=False
        )
    )


    # ========================================================
    # Step 4:
    # Define primary / secondary sets
    # ========================================================

    primary = candidate_effects[
        candidate_effects[
            "tier"
        ]
        ==
        "tier1_direct_oxphos"
    ].copy()


    secondary = candidate_effects.copy()


    # Background excludes every mitochondrial candidate
    candidate_names_present = set(
        candidate_effects[
            "compound_name"
        ]
    )


    background = effects[
        ~effects[
            "compound_name"
        ].isin(
            candidate_names_present
        )
    ].copy()


    # ========================================================
    # Step 5:
    # Matched drug-set permutation
    # ========================================================

    print(
        "\nStep 4: Matched drug-set permutation"
    )


    rng = np.random.default_rng(
        RANDOM_SEED
    )


    permutation_rows = []

    null_distributions = {}


    analysis_sets = {
        "tier1_direct_oxphos":
            primary,

        "tier1_plus_tier2_mito":
            secondary,
    }


    for set_name, subset in analysis_sets.items():

        print(
            "\nAnalyzing:",
            set_name,
            "n =",
            len(
                subset
            )
        )


        if len(
            subset
        ) == 0:

            continue


        result, null = matched_permutation(
            set_name=set_name,
            target_rows=subset,
            background=background,
            rng=rng,
        )


        if result is None:

            continue


        permutation_rows.append(
            result
        )


        null_distributions[
            set_name
        ] = null


    permutation_df = pd.DataFrame(
        permutation_rows
    )


    if len(
        permutation_df
    ) > 0:

        permutation_df[
            "FDR_ASTShigh_sensitive"
        ] = multipletests(
            permutation_df[
                "p_ASTShigh_sensitive"
            ],
            method="fdr_bh"
        )[1]


    permutation_df.to_csv(
        OUT_DIR
        / "19_mito_drug_set_permutation.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # Step 6:
    # Candidate bar plot
    # ========================================================

    if len(
        candidate_effects
    ) > 0:

        plot_data = (
            candidate_effects
            .sort_values(
                CONSENSUS_COLUMN
            )
            .copy()
        )


        fig, ax = plt.subplots(
            figsize=(8, 6)
        )


        y = np.arange(
            len(
                plot_data
            )
        )


        ax.barh(
            y,
            plot_data[
                CONSENSUS_COLUMN
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
                "canonical_name"
            ]
        )


        ax.set_xlabel(
            "Cross-screen consensus z\n"
            "negative = ASTS-high greater sensitivity"
        )


        ax.set_title(
            "Mitochondrial-targeting drugs"
        )


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / "19_01_mito_drug_consensus.png",
            dpi=200
        )


        plt.close(
            fig
        )


    # ========================================================
    # Step 7:
    # Null-distribution plots
    # ========================================================

    for _, row in permutation_df.iterrows():

        set_name = row[
            "drug_set"
        ]


        null = null_distributions[
            set_name
        ]


        fig, ax = plt.subplots(
            figsize=(7, 5)
        )


        ax.hist(
            null,
            bins=50
        )


        ax.axvline(
            row[
                "mean_consensus_z"
            ],
            linestyle="--"
        )


        ax.set_xlabel(
            "Mean consensus z of matched random drug set"
        )


        ax.set_ylabel(
            "Permutation count"
        )


        ax.set_title(
            f"{set_name}\n"
            "observed mitochondrial set shown by dashed line"
        )


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / f"19_02_permutation_{set_name}.png",
            dpi=200
        )


        plt.close(
            fig
        )


    # ========================================================
    # Step 8:
    # Summary
    # ========================================================

    summary_file = (
        OUT_DIR
        / "19_summary.txt"
    )


    with open(
        summary_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 19 summary",
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
            "Test whether the selective mitochondrial "
            "genetic dependency observed in Steps17-18 "
            "is accompanied by pharmacological sensitivity.",
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
            "This is a post-hoc mechanism-defined "
            "pharmacological validation, not a fully "
            "independent prospective validation.",
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Drug-response model:",
            file=f
        )

        print(
            "Step15 driver + CellCycle + MYC adjusted model",
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Candidate availability",
            file=f
        )

        print(
            "----------------------",
            file=f
        )


        print(
            availability.to_string(
                index=False
            ),
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Candidate effects",
            file=f
        )

        print(
            "-----------------",
            file=f
        )


        if len(
            candidate_effects
        ) > 0:

            print(
                candidate_effects.to_string(
                    index=False
                ),
                file=f
            )


        print(
            "",
            file=f
        )


        print(
            "Matched drug-set permutation",
            file=f
        )

        print(
            "----------------------------",
            file=f
        )


        if len(
            permutation_df
        ) > 0:

            print(
                permutation_df.to_string(
                    index=False
                ),
                file=f
            )


        print(
            "",
            file=f
        )


        print(
            "Primary inference rule:",
            file=f
        )


        if len(
            primary
        ) < MIN_PRIMARY_DRUGS_FOR_INFERENCE:

            print(
                "Primary direct-OXPHOS set contains fewer "
                f"than {MIN_PRIMARY_DRUGS_FOR_INFERENCE} "
                "cross-screen drugs.",
                file=f
            )

            print(
                "Treat primary analysis as underpowered "
                "and descriptive.",
                file=f
            )

        else:

            print(
                "Primary analysis has sufficient candidate "
                "count for the pre-specified matched "
                "permutation analysis.",
                file=f
            )


        print(
            "",
            file=f
        )


        print(
            "Input checksum",
            file=f
        )

        print(
            "--------------",
            file=f
        )


        print(
            "Step15 drug effects:",
            sha256_file(
                DRUG_EFFECT_FILE
            ),
            file=f
        )


    print(
        "\n============================================================"
    )

    print(
        "Step 19 completed successfully"
    )

    print(
        "============================================================"
    )

    print(
        summary_file
    )


if __name__ == "__main__":
    main()
