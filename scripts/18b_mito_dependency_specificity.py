#!/usr/bin/env python3

"""
ASTRA-Drug
Step 18B: Mitochondrial dependency specificity audit

----

    oxidative phosphorylation
    mitochondrial translation
    respiratory electron transport








   Gene Effect
      ~ ASTS
      + lineage
      + sex
      + proliferation
      + MAPK/ERBB/PI3K drivers
      + CellCycle
      + MYC
      + global common-essential burden
      + global centered CRISPR shift


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


MITO_SET_FILE = (
    PROJECT_DIR
    / "results/18_mito_dependency_specificity"
    / "18_mitochondrial_gene_sets.tsv"
)


OUT_DIR = (
    PROJECT_DIR
    / "results/18_mito_dependency_specificity"
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

MIN_MODELS_PER_GENE = 800

GENE_CHUNK_SIZE = 500


# ------------------------------------------------------------
# Global common-essential reference
# ------------------------------------------------------------

N_REFERENCE_ESSENTIAL_GENES = 500

MIN_REFERENCE_COVERAGE = 0.95

ESSENTIAL_EFFECT_THRESHOLD = -0.5

MIN_NEGATIVE_FRACTION = 0.80


# ------------------------------------------------------------
# Matched permutation
# ------------------------------------------------------------

N_PERMUTATIONS = 10000

MATCH_POOL_SIZE = 100

RANDOM_SEED = 20260927


# ============================================================
# 2. Imports
# ============================================================

import re

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

from scipy.stats import (
    spearmanr,
    t as t_distribution,
)

from statsmodels.stats.multitest import multipletests


# ============================================================
# 3. Helper functions
# ============================================================

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


def detect_model_id_column(columns):

    candidates = [
        "ModelID",
        "DepMap_ID",
        "DepMapID",
        "model_id",
    ]

    for col in candidates:

        if col in columns:
            return col

    return columns[0]


def build_lineage_categories(metadata):

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
    Step17 stringent model
    +
    two global CRISPR burden covariates.
    """

    parts = []


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

                "GlobalCEBurden_z":
                    zscore(
                        metadata[
                            "GlobalCEBurden"
                        ]
                    ),

                "GlobalCenteredShift_z":
                    zscore(
                        metadata[
                            "GlobalCenteredShift"
                        ]
                    ),
            },
            index=metadata.index
        )
    )


    parts.append(
        metadata[
            [
                "MAPK_driver",
                "ERBB_driver",
                "PI3K_driver",
            ]
        ].astype(float)
    )


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


    return pd.concat(
        parts,
        axis=1
    ).astype(float)


def fit_single_gene_missing(
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


    asts_z = zscore(
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
            np.isfinite(asts_z)
        )
    ):

        return None


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
        se <= 0
    ):

        return None


    t_stat = (
        beta
        /
        se
    )


    p_value = (
        2
        *
        t_distribution.sf(
            abs(
                t_stat
            ),
            df=df_resid
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


    partial_r = (
        t_stat
        /
        np.sqrt(
            t_stat ** 2
            +
            df_resid
        )
    )


    return {
        "n_models":
            len(y_valid),

        "beta_burden_adjusted":
            beta,

        "standardized_beta_burden_adjusted":
            standardized_beta,

        "SE_burden_adjusted":
            se,

        "t_stat_burden_adjusted":
            t_stat,

        "p_value_burden_adjusted":
            p_value,

        "partial_r_burden_adjusted":
            partial_r,
    }


def prepare_matching_features(
    gene_stats,
    eligible_genes
):
    """
    """

    x = gene_stats.loc[
        eligible_genes
    ].copy()


    x[
        "log_sd"
    ] = np.log(
        x[
            "gene_effect_sd"
        ].clip(
            lower=1e-8
        )
    )


    mean_mu = (
        x[
            "gene_effect_median"
        ].mean()
    )

    mean_sd = (
        x[
            "gene_effect_median"
        ].std(
            ddof=1
        )
    )


    logsd_mu = (
        x[
            "log_sd"
        ].mean()
    )

    logsd_sd = (
        x[
            "log_sd"
        ].std(
            ddof=1
        )
    )


    x[
        "match_mean_z"
    ] = (
        x[
            "gene_effect_median"
        ]
        -
        mean_mu
    ) / mean_sd


    x[
        "match_logsd_z"
    ] = (
        x[
            "log_sd"
        ]
        -
        logsd_mu
    ) / logsd_sd


    return x


def matched_gene_set_permutation(
    set_name,
    target_genes,
    association,
    gene_stats,
    control_genes,
    rng,
):
    """
    """

    target_genes = [
        g
        for g in target_genes
        if (
            g in association.index
            and
            g in gene_stats.index
        )
    ]


    if len(target_genes) < 10:

        return None, None


    usable_genes = list(
        set(
            control_genes
        )
        |
        set(
            target_genes
        )
    )


    matching = prepare_matching_features(
        gene_stats,
        usable_genes
    )


    controls = [
        g
        for g in control_genes
        if g in matching.index
    ]


    control_matrix = matching.loc[
        controls,
        [
            "match_mean_z",
            "match_logsd_z",
        ]
    ].to_numpy(
        dtype=float
    )


    # --------------------------------------------------------
    # Candidate pool for every target gene
    # --------------------------------------------------------

    candidate_pools = {}


    for gene in target_genes:

        target_vector = matching.loc[
            gene,
            [
                "match_mean_z",
                "match_logsd_z",
            ]
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
            len(controls)
        )


        indices = np.argpartition(
            distance,
            k - 1
        )[
            :k
        ]


        # nearest first
        indices = indices[
            np.argsort(
                distance[
                    indices
                ]
            )
        ]


        candidate_pools[
            gene
        ] = [
            controls[i]
            for i in indices
        ]


    # --------------------------------------------------------
    # Observed target-set statistic
    # --------------------------------------------------------

    observed_values = association.loc[
        target_genes,
        "t_stat_burden_adjusted"
    ].to_numpy(
        dtype=float
    )


    observed_mean_t = float(
        np.mean(
            observed_values
        )
    )


    observed_median_t = float(
        np.median(
            observed_values
        )
    )


    fraction_negative = float(
        np.mean(
            observed_values < 0
        )
    )


    # --------------------------------------------------------
    # Matched random sets
    # --------------------------------------------------------

    null_means = np.empty(
        N_PERMUTATIONS,
        dtype=float
    )


    target_order = list(
        target_genes
    )


    all_controls = np.array(
        controls,
        dtype=object
    )


    for permutation in range(
        N_PERMUTATIONS
    ):

        rng.shuffle(
            target_order
        )


        selected = set()


        sampled_genes = []


        for gene in target_order:

            pool = candidate_pools[
                gene
            ]


            available = [
                g
                for g in pool
                if g not in selected
            ]


            if len(available) > 0:

                chosen = rng.choice(
                    available
                )


            else:

                remaining = [
                    g
                    for g in all_controls
                    if g not in selected
                ]


                chosen = rng.choice(
                    remaining
                )


            selected.add(
                chosen
            )


            sampled_genes.append(
                chosen
            )


        null_means[
            permutation
        ] = association.loc[
            sampled_genes,
            "t_stat_burden_adjusted"
        ].mean()


    # --------------------------------------------------------
    # Empirical significance
    #
    # Primary alternative:
    # target is more negative than matched random sets
    # --------------------------------------------------------

    p_more_dependent = (
        1
        +
        np.sum(
            null_means
            <=
            observed_mean_t
        )
    ) / (
        N_PERMUTATIONS
        +
        1
    )


    p_less_dependent = (
        1
        +
        np.sum(
            null_means
            >=
            observed_mean_t
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
            p_more_dependent,
            p_less_dependent
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
            observed_mean_t
            -
            null_mean
        )
        /
        null_sd
        if null_sd > 0
        else np.nan
    )


    summary = {
        "gene_set":
            set_name,

        "n_genes":
            len(
                target_genes
            ),

        "mean_t":
            observed_mean_t,

        "median_t":
            observed_median_t,

        "fraction_negative":
            fraction_negative,

        "matched_null_mean":
            null_mean,

        "matched_null_sd":
            null_sd,

        "matched_z":
            matched_z,

        "p_more_dependent":
            p_more_dependent,

        "p_less_dependent":
            p_less_dependent,

        "p_two_sided":
            p_two_sided,
    }


    return (
        summary,
        null_means,
    )


# ============================================================
# 4. Main
# ============================================================

def main():

    print(
        "\n"
        "============================================================\n"
        "ASTRA-Drug Step 18B\n"
        "Mitochondrial dependency specificity audit\n"
        "============================================================"
    )


    for path in [
        CRISPR_FILE,
        ASTS_FILE,
        DRIVER_FILE,
        CELL_CYCLE_FILE,
        MITO_SET_FILE,
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
    # Step 1: Mitochondrial target genes
    # ========================================================

    print(
        "\nStep 1: Reading mitochondrial gene sets"
    )


    mito_sets_table = pd.read_csv(
        MITO_SET_FILE,
        sep="\t"
    )


    mito_sets = {
        set_name:
            set(
                subset[
                    "gene_symbol"
                ]
                .dropna()
                .astype(str)
            )
        for set_name, subset
        in mito_sets_table.groupby(
            "gs_name"
        )
    }


    all_mito_genes = set().union(
        *mito_sets.values()
    )


    print(
        "Mitochondrial genes excluded from burden construction:",
        len(
            all_mito_genes
        )
    )


    # ========================================================
    # Step 2: Metadata
    # ========================================================

    print(
        "\nStep 2: Preparing metadata"
    )


    asts = pd.read_csv(
        ASTS_FILE,
        sep="\t"
    )


    asts[
        "LineageModel"
    ] = build_lineage_categories(
        asts
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
    ]


    cellcycle = pd.read_csv(
        CELL_CYCLE_FILE,
        sep="\t"
    )[
        [
            "ModelID",
            "CellCycle_score",
            "MYC_V1_score",
        ]
    ]


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


    # ========================================================
    # Step 3: CRISPR matrix
    # ========================================================

    print(
        "\nStep 3: Reading CRISPR matrix"
    )


    header = pd.read_csv(
        CRISPR_FILE,
        nrows=0
    )


    model_id_column = detect_model_id_column(
        header.columns.tolist()
    )


    gene_columns = [
        col
        for col in header.columns
        if col != model_id_column
    ]


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
    # Step 4: Gene matrix and baseline statistics
    # ========================================================

    print(
        "\nStep 4: Computing gene-level baseline statistics"
    )


    Y_df = data[
        gene_columns
    ].apply(
        pd.to_numeric,
        errors="coerce"
    )


    Y = Y_df.to_numpy(
        dtype=float
    )


    symbols = [
        extract_gene_symbol(
            col
        )
        for col in gene_columns
    ]


    coverage = np.mean(
        np.isfinite(
            Y
        ),
        axis=0
    )


    medians = np.nanmedian(
        Y,
        axis=0
    )


    means = np.nanmean(
        Y,
        axis=0
    )


    sds = np.nanstd(
        Y,
        axis=0,
        ddof=1
    )


    negative_fraction = np.nanmean(
        Y
        <
        ESSENTIAL_EFFECT_THRESHOLD,
        axis=0
    )


    gene_stats = pd.DataFrame(
        {
            "gene_symbol":
                symbols,

            "CRISPR_column":
                gene_columns,

            "coverage":
                coverage,

            "gene_effect_median":
                medians,

            "gene_effect_mean":
                means,

            "gene_effect_sd":
                sds,

            "fraction_below_minus0.5":
                negative_fraction,
        }
    )


    gene_stats = (
        gene_stats
        .sort_values(
            [
                "gene_symbol",
                "coverage",
            ],
            ascending=[
                True,
                False,
            ]
        )
        .drop_duplicates(
            "gene_symbol"
        )
        .set_index(
            "gene_symbol"
        )
    )


    # ========================================================
    # Step 5: Build non-mito common-essential reference
    # ========================================================

    print(
        "\nStep 5: Selecting non-mito common-essential reference genes"
    )


    candidate_stats = gene_stats[
        (
            gene_stats[
                "coverage"
            ]
            >=
            MIN_REFERENCE_COVERAGE
        )
        &
        (
            ~gene_stats.index.isin(
                all_mito_genes
            )
        )
    ].copy()


    strong_reference = candidate_stats[
        candidate_stats[
            "fraction_below_minus0.5"
        ]
        >=
        MIN_NEGATIVE_FRACTION
    ].copy()


    if (
        len(
            strong_reference
        )
        >=
        N_REFERENCE_ESSENTIAL_GENES
    ):

        reference = (
            strong_reference
            .sort_values(
                "gene_effect_median"
            )
            .head(
                N_REFERENCE_ESSENTIAL_GENES
            )
        )


        reference_rule = (
            "coverage>=0.95; "
            "fraction(effect<-0.5)>=0.80; "
            "lowest median gene effect"
        )


    else:

        reference = (
            candidate_stats
            .sort_values(
                "gene_effect_median"
            )
            .head(
                N_REFERENCE_ESSENTIAL_GENES
            )
        )


        reference_rule = (
            "fallback: coverage>=0.95; "
            "lowest median gene effect"
        )


    reference_genes = set(
        reference.index
    )


    print(
        "Reference genes:",
        len(
            reference_genes
        )
    )


    print(
        "Rule:",
        reference_rule
    )


    reference.reset_index().to_csv(
        OUT_DIR
        / "18_global_dependency_reference_genes.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # Step 6: Build global burden metrics
    # ========================================================

    print(
        "\nStep 6: Calculating global CRISPR burden"
    )


    symbol_to_column = {
        extract_gene_symbol(
            col
        ):
            col
        for col in gene_columns
    }


    reference_columns = [
        symbol_to_column[
            gene
        ]
        for gene in reference_genes
        if gene in symbol_to_column
    ]


    global_ce_burden = data[
        reference_columns
    ].median(
        axis=1,
        skipna=True
    )


    # --------------------------------------------------------
    # Second metric:
    # median gene-centered shift across non-mito genes
    # --------------------------------------------------------

    non_mito_indices = [
        i
        for i, gene
        in enumerate(
            symbols
        )
        if (
            gene
            not in all_mito_genes
        )
    ]


    Y_non_mito = Y[
        :,
        non_mito_indices
    ]


    med_non_mito = np.nanmedian(
        Y_non_mito,
        axis=0
    )


    centered = (
        Y_non_mito
        -
        med_non_mito[
            None,
            :
        ]
    )


    global_centered_shift = np.nanmedian(
        centered,
        axis=1
    )


    data[
        "GlobalCEBurden"
    ] = global_ce_burden


    data[
        "GlobalCenteredShift"
    ] = global_centered_shift


    burden_scores = data[
        [
            "ModelID",
            ASTS_COLUMN,
            "GlobalCEBurden",
            "GlobalCenteredShift",
        ]
    ].copy()


    burden_scores.to_csv(
        OUT_DIR
        / "18_global_dependency_burden_scores.tsv",
        sep="\t",
        index=False
    )


    # --------------------------------------------------------
    # Correlations
    # --------------------------------------------------------

    correlation_rows = []


    for variable in [
        "GlobalCEBurden",
        "GlobalCenteredShift",
    ]:

        rho, p = spearmanr(
            data[
                ASTS_COLUMN
            ],
            data[
                variable
            ],
            nan_policy="omit"
        )


        correlation_rows.append(
            {
                "variable":
                    variable,

                "Spearman_r_with_ASTS":
                    rho,

                "p":
                    p,
            }
        )


    rho_burdens, p_burdens = spearmanr(
        data[
            "GlobalCEBurden"
        ],
        data[
            "GlobalCenteredShift"
        ],
        nan_policy="omit"
    )


    correlation_rows.append(
        {
            "variable":
                "GlobalCEBurden_vs_GlobalCenteredShift",

            "Spearman_r_with_ASTS":
                rho_burdens,

            "p":
                p_burdens,
        }
    )


    correlation_df = pd.DataFrame(
        correlation_rows
    )


    correlation_df.to_csv(
        OUT_DIR
        / "18_global_dependency_burden_correlations.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # Step 7: Burden-adjusted genome-wide association
    # ========================================================

    print(
        "\nStep 7: Burden-adjusted genome-wide association"
    )


    metadata_columns = list(
        metadata.columns
    ) + [
        "GlobalCEBurden",
        "GlobalCenteredShift",
    ]


    metadata_aligned = data[
        metadata_columns
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
        np.isfinite(
            C
        )
    ):

        raise RuntimeError(
            "Non-finite covariates remain."
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


        print(
            f"Genes {start + 1}-{end} / "
            f"{len(gene_columns)}"
        )


        cols = gene_columns[
            start:end
        ]


        chunk = data[
            cols
        ].apply(
            pd.to_numeric,
            errors="coerce"
        ).to_numpy(
            dtype=float
        )


        complete_mask = np.all(
            np.isfinite(
                chunk
            ),
            axis=0
        )


        # ----------------------------------------------------
        # Complete genes
        # ----------------------------------------------------

        complete_indices = np.where(
            complete_mask
        )[0]


        if len(
            complete_indices
        ) > 0:

            Y_complete = chunk[
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


            sigma2 = (
                np.sum(
                    residual ** 2,
                    axis=0
                )
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
                2
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


                gene = extract_gene_symbol(
                    col
                )


                results.append(
                    {
                        "gene_symbol":
                            gene,

                        "CRISPR_column":
                            col,

                        "n_models":
                            len(data),

                        "beta_burden_adjusted":
                            beta[j],

                        "standardized_beta_burden_adjusted":
                            standardized_beta[j],

                        "SE_burden_adjusted":
                            se[j],

                        "t_stat_burden_adjusted":
                            t_stat[j],

                        "p_value_burden_adjusted":
                            p_value[j],

                        "partial_r_burden_adjusted":
                            partial_r[j],

                        "is_burden_reference_gene":
                            gene
                            in reference_genes,

                        "is_mito_target_gene":
                            gene
                            in all_mito_genes,

                        "has_missing_gene_effect":
                            False,
                    }
                )


        # ----------------------------------------------------
        # Missing genes
        # ----------------------------------------------------

        incomplete_indices = np.where(
            ~complete_mask
        )[0]


        for local_index in incomplete_indices:

            col = cols[
                local_index
            ]


            gene = extract_gene_symbol(
                col
            )


            result = fit_single_gene_missing(
                chunk[
                    :,
                    local_index
                ],
                metadata_aligned
            )


            if result is None:
                continue


            results.append(
                {
                    "gene_symbol":
                        gene,

                    "CRISPR_column":
                        col,

                    **result,

                    "is_burden_reference_gene":
                        gene
                        in reference_genes,

                    "is_mito_target_gene":
                        gene
                        in all_mito_genes,

                    "has_missing_gene_effect":
                        True,
                }
            )


    association = pd.DataFrame(
        results
    )


    association[
        "abs_t"
    ] = association[
        "t_stat_burden_adjusted"
    ].abs()


    association = (
        association
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
            "gene_symbol"
        )
        .drop(
            columns=[
                "abs_t"
            ]
        )
        .copy()
    )


    # Add baseline gene stats
    association = association.merge(
        gene_stats.reset_index()[
            [
                "gene_symbol",
                "coverage",
                "gene_effect_median",
                "gene_effect_mean",
                "gene_effect_sd",
                "fraction_below_minus0.5",
            ]
        ],
        on="gene_symbol",
        how="left",
        validate="1:1"
    )


    finite = np.isfinite(
        association[
            "p_value_burden_adjusted"
        ]
    )


    association[
        "FDR_burden_adjusted"
    ] = np.nan


    association.loc[
        finite,
        "FDR_burden_adjusted"
    ] = multipletests(
        association.loc[
            finite,
            "p_value_burden_adjusted"
        ],
        method="fdr_bh"
    )[1]


    association = association.sort_values(
        "t_stat_burden_adjusted",
        ascending=False
    )


    association.to_csv(
        OUT_DIR
        / "18_CRISPR_genomewide_burden_adjusted.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # Step 8: Matched gene-set permutation
    # ========================================================

    print(
        "\nStep 8: Matched gene-set permutation"
    )


    association_indexed = association.set_index(
        "gene_symbol"
    )


    gene_stats_indexed = gene_stats.copy()


    control_genes = [
        gene
        for gene in association_indexed.index
        if (
            gene
            not in all_mito_genes
            and
            gene
            not in reference_genes
            and
            np.isfinite(
                association_indexed.loc[
                    gene,
                    "t_stat_burden_adjusted"
                ]
            )
            and
            gene
            in gene_stats_indexed.index
            and
            np.isfinite(
                gene_stats_indexed.loc[
                    gene,
                    "gene_effect_sd"
                ]
            )
            and
            gene_stats_indexed.loc[
                gene,
                "gene_effect_sd"
            ]
            >
            0
        )
    ]


    rng = np.random.default_rng(
        RANDOM_SEED
    )


    primary_sets = [
        "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
        "REACTOME_MITOCHONDRIAL_TRANSLATION",
        "REACTOME_RESPIRATORY_ELECTRON_TRANSPORT",
        "MITO_CORE_UNION",
    ]


    permutation_rows = []

    null_distributions = {}


    for set_name in primary_sets:

        if set_name not in mito_sets:
            continue


        print(
            "Permutation:",
            set_name
        )


        summary, null_means = (
            matched_gene_set_permutation(
                set_name=set_name,
                target_genes=mito_sets[
                    set_name
                ],
                association=association_indexed,
                gene_stats=gene_stats_indexed,
                control_genes=control_genes,
                rng=rng,
            )
        )


        if summary is None:
            continue


        permutation_rows.append(
            summary
        )


        null_distributions[
            set_name
        ] = null_means


    permutation_df = pd.DataFrame(
        permutation_rows
    )


    if len(
        permutation_df
    ) > 0:

        permutation_df[
            "FDR_more_dependent"
        ] = multipletests(
            permutation_df[
                "p_more_dependent"
            ],
            method="fdr_bh"
        )[1]


    permutation_df.to_csv(
        OUT_DIR
        / "18_mito_matched_permutation.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # Step 9: Permutation figure
    # ========================================================

    if len(
        permutation_df
    ) > 0:

        fig, ax = plt.subplots(
            figsize=(9, 5)
        )


        x = np.arange(
            len(
                permutation_df
            )
        )


        null_lower = []

        null_upper = []


        for set_name in permutation_df[
            "gene_set"
        ]:

            null = null_distributions[
                set_name
            ]


            null_lower.append(
                np.quantile(
                    null,
                    0.025
                )
            )


            null_upper.append(
                np.quantile(
                    null,
                    0.975
                )
            )


        observed = permutation_df[
            "mean_t"
        ].to_numpy()


        lower_error = (
            observed
            -
            np.array(
                null_lower
            )
        )


        upper_error = (
            np.array(
                null_upper
            )
            -
            observed
        )


        # plot null 95% interval as vertical lines
        for i in range(
            len(x)
        ):

            ax.plot(
                [
                    x[i],
                    x[i],
                ],
                [
                    null_lower[i],
                    null_upper[i],
                ],
                linewidth=4,
                alpha=0.4
            )


        ax.scatter(
            x,
            observed,
            s=70
        )


        ax.axhline(
            0,
            linestyle="--"
        )


        ax.set_xticks(
            x
        )


        ax.set_xticklabels(
            permutation_df[
                "gene_set"
            ],
            rotation=30,
            ha="right"
        )


        ax.set_ylabel(
            "Mean burden-adjusted ASTS t statistic\n"
            "negative = ASTS-high more dependent"
        )


        ax.set_title(
            "Mitochondrial dependency vs matched gene sets"
        )


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / "18_mito_matched_permutation.png",
            dpi=200
        )


        plt.close(
            fig
        )


    # ========================================================
    # Step 10: Summary
    # ========================================================

    summary_file = (
        OUT_DIR
        / "18b_summary.txt"
    )


    focus_genes = [
        "EGFR",
        "MAP2K1",
        "MAP2K2",
        "SLC25A51",
        "MTIF3",
        "RMND1",
        "NDUFB4",
        "MRPS11",
        "LRPPRC",
        "NDUFA6",
        "NDUFS8",
        "COX19",
    ]


    with open(
        summary_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 18B summary",
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
            "Purpose:",
            file=f
        )

        print(
            "Test whether the ASTS-high mitochondrial "
            "dependency observed in Step17 is explained "
            "by global CRISPR dependency burden.",
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
            f"Global reference genes: "
            f"{len(reference_genes)}",
            file=f
        )


        print(
            f"Reference rule: "
            f"{reference_rule}",
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Burden correlations",
            file=f
        )

        print(
            "-------------------",
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
            "Burden-adjusted genome-wide results",
            file=f
        )

        print(
            "-----------------------------------",
            file=f
        )


        print(
            f"Genes analyzed: "
            f"{len(association)}",
            file=f
        )


        print(
            f"FDR < 0.05: "
            f"{(association['FDR_burden_adjusted'] < 0.05).sum()}",
            file=f
        )


        print(
            "",
            file=f
        )


        focus = association[
            association[
                "gene_symbol"
            ].isin(
                focus_genes
            )
        ].copy()


        order = {
            gene:
                i
            for i, gene
            in enumerate(
                focus_genes
            )
        }


        focus[
            "_order"
        ] = focus[
            "gene_symbol"
        ].map(
            order
        )


        focus = focus.sort_values(
            "_order"
        )


        print(
            "Focus genes:",
            file=f
        )


        print(
            focus[
                [
                    "gene_symbol",
                    "standardized_beta_burden_adjusted",
                    "t_stat_burden_adjusted",
                    "p_value_burden_adjusted",
                    "FDR_burden_adjusted",
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
            "Matched mitochondrial gene-set permutation",
            file=f
        )

        print(
            "-----------------------------------------",
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
            "Interpretation:",
            file=f
        )


        print(
            "If mitochondrial gene sets retain negative "
            "ASTS effects after global-burden adjustment "
            "and are significantly more negative than "
            "essentiality/variance-matched random gene sets, "
            "the Step17 signal supports a selective "
            "mitochondrial dependency state rather than a "
            "general CRISPR essentiality shift.",
            file=f
        )


    print(
        "\n============================================================"
    )

    print(
        "Step 18B completed successfully"
    )

    print(
        "============================================================"
    )

    print(
        summary_file
    )


if __name__ == "__main__":
    main()
