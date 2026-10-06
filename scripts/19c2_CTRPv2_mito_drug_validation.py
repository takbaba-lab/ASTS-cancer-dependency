#!/usr/bin/env python3

"""
ASTRA-Drug
Step 19C-2: CTRPv2 mitochondrial drug validation

----



    lower = more sensitive



    ASTS beta < 0



"""


# ============================================================
# 1. Path / Parameter
# ============================================================

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]


EXTRACT_DIR = (
    PROJECT_DIR
    / "results/19BC_mito_drug_validation/CTRPv2_extract"
)


RESPONSE_FILE = (
    EXTRACT_DIR
    / "19C_CTRPv2_response_long.tsv"
)


SAMPLE_INFO_FILE = (
    EXTRACT_DIR
    / "19C_CTRPv2_sample_info.tsv"
)


TREATMENT_INFO_FILE = (
    EXTRACT_DIR
    / "19C_CTRPv2_treatment_info.tsv"
)


SAMPLE_CURATION_FILE = (
    EXTRACT_DIR
    / "19C_CTRPv2_curation_sample.tsv"
)


TREATMENT_CURATION_FILE = (
    EXTRACT_DIR
    / "19C_CTRPv2_curation_treatment.tsv"
)


FREEZE_FILE = (
    PROJECT_DIR
    / "results/19_mito_drug_validation"
    / "19_MITO_DRUG_SET_FREEZE.tsv"
)


MODEL_FILE = (
    PROJECT_DIR
    / "data/raw/depmap/26Q1"
    / "Model.csv"
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
    / "results/19BC_mito_drug_validation"
)


FIG_DIR = OUT_DIR / "figures"


# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

ASTS_COLUMN = "ASTS_score"
PROLIFERATION_COLUMN = "proliferation_proxy"
LINEAGE_COLUMN = "OncotreeLineage"
SEX_COLUMN = "Sex"

MIN_GLOBAL_LINEAGE_N = 5
MIN_N_PER_DRUG = 100
MIN_RESIDUAL_DF = 30

N_PERMUTATIONS = 10000
MATCH_POOL_SIZE = 50
RANDOM_SEED = 20261927


# ============================================================
# 2. Imports
# ============================================================

import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm

from statsmodels.stats.multitest import multipletests


# ============================================================
# 3. Generic helpers
# ============================================================

def normalize_identifier(x):

    if pd.isna(x):
        return ""

    return re.sub(
        r"[^a-z0-9]+",
        "",
        str(x).lower().strip(),
    )


def normalize_drug_name(x):

    return normalize_identifier(x)


def zscore(series):

    x = pd.to_numeric(
        series,
        errors="coerce",
    )

    sd = x.std(ddof=1)

    if (
        not np.isfinite(sd)
        or sd == 0
    ):

        return pd.Series(
            np.nan,
            index=x.index,
        )

    return (
        x - x.mean()
    ) / sd


def load_optional_tsv(path):

    if path.exists():

        return pd.read_csv(
            path,
            sep="\t",
            low_memory=False,
        )

    return pd.DataFrame()


# ============================================================
# 4. Candidate freeze
# ============================================================

def load_candidate_freeze():

    freeze = pd.read_csv(
        FREEZE_FILE,
        sep="\t",
    )

    alias_to_canonical = {}

    for _, row in freeze.iterrows():

        aliases = str(
            row["aliases"]
        ).split(";")

        aliases.append(
            str(
                row["canonical_name"]
            )
        )

        for alias in aliases:

            key = normalize_drug_name(
                alias
            )

            if key:

                alias_to_canonical[
                    key
                ] = row[
                    "canonical_name"
                ]

    return (
        freeze,
        alias_to_canonical,
    )


def match_candidate_identifiers(
    identifiers,
    alias_map,
):
    """
    """

    normalized = {
        normalize_drug_name(x)
        for x in identifiers
        if normalize_drug_name(x)
    }

    for x in normalized:

        if x in alias_map:

            return alias_map[x]

    return None


# ============================================================
# 5. Metadata
# ============================================================

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
            counts >= MIN_GLOBAL_LINEAGE_N
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
        sep="\t",
    )

    asts["LineageModel"] = (
        build_lineage_categories(
            asts
        )
    )

    drivers = pd.read_csv(
        DRIVER_FILE,
        sep="\t",
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
        sep="\t",
    )[
        [
            "ModelID",
            "CellCycle_score",
            "MYC_V1_score",
        ]
    ]

    return (
        asts
        .merge(
            drivers,
            on="ModelID",
            how="left",
            validate="1:1",
        )
        .merge(
            cellcycle,
            on="ModelID",
            how="left",
            validate="1:1",
        )
    )


# ============================================================
# 6. CTRPv2 -> DepMap sample mapping
# ============================================================

def useful_depmap_columns(model):

    preferred = [
        "ModelID",
        "CCLEName",
        "CellLineName",
        "StrippedCellLineName",
        "CellosaurusID",
        "SangerModelID",
    ]

    return [
        x
        for x in preferred
        if x in model.columns
    ]


def build_depmap_identifier_index(model):

    index = {}

    columns = useful_depmap_columns(
        model
    )

    for _, row in model.iterrows():

        model_id = row[
            "ModelID"
        ]

        for col in columns:

            key = normalize_identifier(
                row[col]
            )

            if len(key) < 4:
                continue

            index.setdefault(
                key,
                set(),
            ).add(model_id)

    return index


def candidate_sample_columns(df):

    keep = []

    keywords = [
        "sample",
        "cell",
        "ccle",
        "depmap",
        "model",
        "name",
        "cellosaurus",
        "_rowname",
    ]

    for col in df.columns:

        lower = col.lower()

        if any(
            key in lower
            for key in keywords
        ):

            if "tissue" not in lower:

                keep.append(col)

    return keep


def build_curation_row_index(curation):

    index = {}

    if len(curation) == 0:
        return index

    for i, row in curation.iterrows():

        values = []

        for value in row:

            key = normalize_identifier(
                value
            )

            if len(key) >= 4:
                values.append(key)

        for key in values:

            index.setdefault(
                key,
                set(),
            ).add(i)

    return index


def map_ctrp_samples():

    sample_info = pd.read_csv(
        SAMPLE_INFO_FILE,
        sep="\t",
        low_memory=False,
    )

    sample_curation = (
        load_optional_tsv(
            SAMPLE_CURATION_FILE
        )
    )

    model = pd.read_csv(
        MODEL_FILE,
        low_memory=False,
    )

    depmap_index = (
        build_depmap_identifier_index(
            model
        )
    )

    curation_index = (
        build_curation_row_index(
            sample_curation
        )
    )

    sample_cols = (
        candidate_sample_columns(
            sample_info
        )
    )

    if "sampleid" not in sample_info.columns:

        raise RuntimeError(
            "sampleid not found in CTRPv2 sampleInfo."
        )

    rows = []

    for _, row in sample_info.iterrows():

        sample_id = str(
            row["sampleid"]
        )

        identifiers = set()

        for col in sample_cols:

            value = normalize_identifier(
                row[col]
            )

            if len(value) >= 4:

                identifiers.add(value)

        # Add linked curation values
        linked_curation_rows = set()

        for identifier in list(
            identifiers
        ):

            linked_curation_rows.update(
                curation_index.get(
                    identifier,
                    set(),
                )
            )

        for idx in linked_curation_rows:

            for value in (
                sample_curation
                .loc[idx]
                .tolist()
            ):

                key = normalize_identifier(
                    value
                )

                if len(key) >= 4:
                    identifiers.add(key)

        matched_models = set()

        for identifier in identifiers:

            matched_models.update(
                depmap_index.get(
                    identifier,
                    set(),
                )
            )

        if len(matched_models) == 1:

            status = "mapped"
            model_id = next(
                iter(matched_models)
            )

        elif len(matched_models) == 0:

            status = "unmapped"
            model_id = ""

        else:

            status = "ambiguous"
            model_id = ""

        rows.append(
            {
                "sampleid":
                    sample_id,

                "ModelID":
                    model_id,

                "mapping_status":
                    status,

                "n_candidate_ModelIDs":
                    len(
                        matched_models
                    ),
            }
        )

    mapping = pd.DataFrame(
        rows
    )

    mapping.to_csv(
        OUT_DIR
        / "19C_CTRPv2_sample_to_DepMap.tsv",
        sep="\t",
        index=False,
    )

    return mapping


# ============================================================
# 7. Treatment metadata
# ============================================================

def prepare_treatment_map(
    freeze,
    alias_map,
):

    treatment = pd.read_csv(
        TREATMENT_INFO_FILE,
        sep="\t",
        low_memory=False,
    )

    curation = load_optional_tsv(
        TREATMENT_CURATION_FILE
    )

    if "treatmentid" not in treatment.columns:

        raise RuntimeError(
            "treatmentid missing from treatmentInfo."
        )

    curation_index = (
        build_curation_row_index(
            curation
        )
    )

    rows = []

    for _, row in treatment.iterrows():

        treatment_id = str(
            row["treatmentid"]
        )

        identifiers = []

        # All treatment metadata values are relevant
        for value in row.tolist():

            if pd.notna(value):

                identifiers.append(
                    str(value)
                )

        normalized = {
            normalize_identifier(x)
            for x in identifiers
            if len(
                normalize_identifier(x)
            )
            >= 4
        }

        linked_rows = set()

        for key in normalized:

            linked_rows.update(
                curation_index.get(
                    key,
                    set(),
                )
            )

        for idx in linked_rows:

            identifiers.extend(
                [
                    str(x)
                    for x in (
                        curation.loc[idx]
                        .tolist()
                    )
                    if pd.notna(x)
                ]
            )

        canonical = (
            match_candidate_identifiers(
                identifiers,
                alias_map,
            )
        )

        # Display name:
        # treatmentid is kept because it is
        # guaranteed to be the PSet key.
        compound_name = treatment_id

        for preferred in [
            "drug_name",
            "treatment_name",
            "name",
            "compound_name",
        ]:

            if (
                preferred
                in treatment.columns
            ):

                value = row[
                    preferred
                ]

                if pd.notna(value):

                    compound_name = str(
                        value
                    )

                    break

        rows.append(
            {
                "treatmentid":
                    treatment_id,

                "compound_name":
                    compound_name,

                "canonical_name":
                    canonical,
            }
        )

    return pd.DataFrame(
        rows
    )


# ============================================================
# 8. Regression
# ============================================================

def fit_one_drug(df):

    required = [
        "response",
        ASTS_COLUMN,
        PROLIFERATION_COLUMN,
        "LineageModel",
        "MAPK_driver",
        "ERBB_driver",
        "PI3K_driver",
        "CellCycle_score",
        "MYC_V1_score",
    ]

    work = df.dropna(
        subset=required
    ).copy()

    if len(work) < MIN_N_PER_DRUG:
        return None

    work["ASTS_z"] = zscore(
        work[ASTS_COLUMN]
    )

    work["Prolif_z"] = zscore(
        work[PROLIFERATION_COLUMN]
    )

    work["CellCycle_z"] = zscore(
        work["CellCycle_score"]
    )

    work["MYC_z"] = zscore(
        work["MYC_V1_score"]
    )

    work = work.dropna(
        subset=[
            "ASTS_z",
            "Prolif_z",
            "CellCycle_z",
            "MYC_z",
        ]
    )

    if len(work) < MIN_N_PER_DRUG:
        return None

    response_sd = work[
        "response"
    ].std(ddof=1)

    if (
        not np.isfinite(response_sd)
        or response_sd == 0
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

    for col in [
        "MAPK_driver",
        "ERBB_driver",
        "PI3K_driver",
    ]:

        if work[col].nunique() > 1:

            parts.append(
                work[[col]].astype(float)
            )

    sex = (
        work[SEX_COLUMN]
        .fillna("Unknown")
        .astype(str)
    )

    if sex.nunique() > 1:

        parts.append(
            pd.get_dummies(
                sex,
                prefix="Sex",
                drop_first=True,
                dtype=float,
            )
        )

    lineage = (
        work["LineageModel"]
        .fillna("Unknown")
        .astype(str)
    )

    if lineage.nunique() > 1:

        parts.append(
            pd.get_dummies(
                lineage,
                prefix="Lineage",
                drop_first=True,
                dtype=float,
            )
        )

    X = pd.concat(
        parts,
        axis=1,
    )

    X = sm.add_constant(
        X,
        has_constant="add",
    )

    model = sm.OLS(
        work["response"],
        X,
    ).fit(
        cov_type="HC3"
    )

    if model.df_resid < MIN_RESIDUAL_DF:
        return None

    beta = float(
        model.params["ASTS_z"]
    )

    se = float(
        model.bse["ASTS_z"]
    )

    return {
        "n_models":
            len(work),

        "response_sd":
            response_sd,

        "beta":
            beta,

        "SE":
            se,

        "p":
            float(
                model.pvalues[
                    "ASTS_z"
                ]
            ),

        "standardized_beta":
            beta / response_sd,

        "standardized_SE":
            se / response_sd,
    }


# ============================================================
# 9. Permutation
# ============================================================

def matched_permutation(
    candidates,
    all_results,
    set_name,
):

    if len(candidates) == 0:
        return None

    rng = np.random.default_rng(
        RANDOM_SEED
    )

    candidate_ids = set(
        candidates["treatmentid"]
    )

    background = all_results[
        ~all_results[
            "treatmentid"
        ].isin(candidate_ids)
    ].copy()

    combined = pd.concat(
        [
            candidates,
            background,
        ],
        ignore_index=True,
    )

    for col in [
        "n_models",
        "response_sd",
    ]:

        v = np.log1p(
            combined[col]
        )

        combined[
            f"m_{col}"
        ] = (
            v - v.mean()
        ) / v.std(ddof=1)

    target = combined[
        combined[
            "treatmentid"
        ].isin(candidate_ids)
    ]

    controls = combined[
        ~combined[
            "treatmentid"
        ].isin(candidate_ids)
    ]

    features = [
        "m_n_models",
        "m_response_sd",
    ]

    control_matrix = controls[
        features
    ].to_numpy()

    control_ids = controls[
        "treatmentid"
    ].to_numpy(
        dtype=object
    )

    pools = {}

    for _, row in target.iterrows():

        vec = row[
            features
        ].to_numpy(
            dtype=float
        )

        d = np.sum(
            (
                control_matrix
                - vec[None, :]
            )
            ** 2,
            axis=1,
        )

        k = min(
            MATCH_POOL_SIZE,
            len(control_ids),
        )

        idx = np.argpartition(
            d,
            k - 1,
        )[:k]

        pools[
            row["treatmentid"]
        ] = [
            control_ids[i]
            for i in idx
        ]

    background_idx = (
        background
        .set_index(
            "treatmentid"
        )
    )

    observed = candidates[
        "standardized_beta"
    ].mean()

    null = np.empty(
        N_PERMUTATIONS
    )

    original = list(
        candidates[
            "treatmentid"
        ]
    )

    all_background = list(
        background_idx.index
    )

    for i in range(
        N_PERMUTATIONS
    ):

        order = original.copy()
        rng.shuffle(order)

        used = set()
        sampled = []

        for treatment_id in order:

            available = [
                x
                for x in pools[
                    treatment_id
                ]
                if x not in used
            ]

            if available:

                chosen = rng.choice(
                    available
                )

            else:

                remaining = [
                    x
                    for x in all_background
                    if x not in used
                ]

                chosen = rng.choice(
                    remaining
                )

            used.add(chosen)
            sampled.append(chosen)

        null[i] = (
            background_idx
            .loc[
                sampled,
                "standardized_beta",
            ]
            .mean()
        )

    p_sensitive = (
        1
        +
        np.sum(
            null <= observed
        )
    ) / (
        N_PERMUTATIONS + 1
    )

    null_mean = np.mean(null)
    null_sd = np.std(
        null,
        ddof=1,
    )

    return {
        "screen":
            "CTRPv2",

        "drug_set":
            set_name,

        "n_drugs":
            len(candidates),

        "mean_standardized_beta":
            observed,

        "fraction_negative":
            np.mean(
                candidates[
                    "standardized_beta"
                ]
                <
                0
            ),

        "matched_null_mean":
            null_mean,

        "matched_null_sd":
            null_sd,

        "matched_z":
            (
                observed
                -
                null_mean
            ) / null_sd,

        "p_ASTShigh_sensitive":
            p_sensitive,
    }


# ============================================================
# 10. Main
# ============================================================

def main():

    print(
        "\n"
        "============================================================\n"
        "ASTRA-Drug Step 19C-2\n"
        "CTRPv2 mitochondrial drug validation\n"
        "============================================================"
    )

    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    FIG_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    freeze, alias_map = (
        load_candidate_freeze()
    )

    metadata = prepare_metadata()

    mapping = map_ctrp_samples()

    treatment_map = (
        prepare_treatment_map(
            freeze,
            alias_map,
        )
    )

    response = pd.read_csv(
        RESPONSE_FILE,
        sep="\t",
        low_memory=False,
    )

    response["response"] = (
        pd.to_numeric(
            response[
                "response_harmonized"
            ],
            errors="coerce",
        )
    )

    response = response[
        response[
            "response"
        ].notna()
    ]

    response = response.merge(
        mapping[
            [
                "sampleid",
                "ModelID",
                "mapping_status",
            ]
        ],
        on="sampleid",
        how="left",
        validate="m:1",
    )

    response = response[
        response[
            "mapping_status"
        ]
        ==
        "mapped"
    ].copy()

    response = response.merge(
        treatment_map,
        on="treatmentid",
        how="left",
        validate="m:1",
    )

    # Collapse duplicated experiments
    response = (
        response
        .groupby(
            [
                "treatmentid",
                "ModelID",
            ],
            as_index=False,
        )
        .agg(
            compound_name=(
                "compound_name",
                "first",
            ),
            canonical_name=(
                "canonical_name",
                "first",
            ),
            response=(
                "response",
                "median",
            ),
        )
    )

    response = response.merge(
        metadata,
        on="ModelID",
        how="inner",
        validate="m:1",
    )

    print(
        "Mapped response rows:",
        len(response),
    )

    print(
        "Unique ModelIDs:",
        response[
            "ModelID"
        ].nunique(),
    )

    # --------------------------------------------------------
    # All CTRPv2 drugs
    # --------------------------------------------------------

    rows = []

    for treatment_id, subset in (
        response.groupby(
            "treatmentid"
        )
    ):

        result = fit_one_drug(
            subset
        )

        if result is None:
            continue

        rows.append(
            {
                "screen":
                    "CTRPv2",

                "treatmentid":
                    treatment_id,

                "compound_name":
                    subset[
                        "compound_name"
                    ].iloc[0],

                "canonical_name":
                    subset[
                        "canonical_name"
                    ].iloc[0],

                **result,
            }
        )

    results = pd.DataFrame(
        rows
    )

    if len(results) == 0:

        raise RuntimeError(
            "No CTRPv2 drugs passed "
            "the analysis threshold."
        )

    results["FDR_all_drugs"] = np.nan

    finite = np.isfinite(
        results["p"]
    )

    results.loc[
        finite,
        "FDR_all_drugs"
    ] = multipletests(
        results.loc[
            finite,
            "p"
        ],
        method="fdr_bh",
    )[1]

    results.to_csv(
        OUT_DIR
        / "19C_CTRPv2_all_drug_effects.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # Candidate availability
    # --------------------------------------------------------

    availability_rows = []

    candidate_results = results[
        results[
            "canonical_name"
        ].notna()
    ].copy()

    candidate_results = (
        candidate_results
        .sort_values(
            [
                "canonical_name",
                "n_models",
            ],
            ascending=[
                True,
                False,
            ],
        )
        .drop_duplicates(
            "canonical_name"
        )
    )

    if len(candidate_results) > 0:

        candidate_results[
            "FDR_candidate_set"
        ] = multipletests(
            candidate_results["p"],
            method="fdr_bh",
        )[1]

    for _, row in freeze.iterrows():

        canonical = row[
            "canonical_name"
        ]

        raw_found = (
            treatment_map[
                "canonical_name"
            ]
            ==
            canonical
        ).any()

        tested = candidate_results[
            candidate_results[
                "canonical_name"
            ]
            ==
            canonical
        ]

        if len(tested) > 0:

            status = "tested"
            matched_name = (
                tested.iloc[0][
                    "compound_name"
                ]
            )

            n_models = int(
                tested.iloc[0][
                    "n_models"
                ]
            )

        elif raw_found:

            status = (
                "found_but_not_tested"
            )

            names = (
                treatment_map.loc[
                    treatment_map[
                        "canonical_name"
                    ]
                    ==
                    canonical,
                    "compound_name",
                ]
                .astype(str)
                .unique()
            )

            matched_name = ";".join(
                names
            )

            n_models = np.nan

        else:

            status = "not_found"
            matched_name = ""
            n_models = np.nan

        availability_rows.append(
            {
                "screen":
                    "CTRPv2",

                "canonical_name":
                    canonical,

                "tier":
                    row["tier"],

                "mechanism_group":
                    row[
                        "mechanism_group"
                    ],

                "status":
                    status,

                "matched_name":
                    matched_name,

                "n_models":
                    n_models,
            }
        )

    availability = pd.DataFrame(
        availability_rows
    )

    availability.to_csv(
        OUT_DIR
        / "19C_CTRPv2_candidate_availability.tsv",
        sep="\t",
        index=False,
    )

    candidate_results = (
        candidate_results
        .merge(
            freeze[
                [
                    "canonical_name",
                    "tier",
                    "mechanism_group",
                ]
            ],
            on="canonical_name",
            how="left",
            validate="1:1",
        )
    )

    candidate_results.to_csv(
        OUT_DIR
        / "19C_CTRPv2_candidate_effects.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # Matched permutation
    # --------------------------------------------------------

    permutation_rows = []

    primary = candidate_results[
        candidate_results[
            "tier"
        ]
        ==
        "tier1_direct_oxphos"
    ]

    secondary = candidate_results.copy()

    for name, subset in {
        "tier1_direct_oxphos":
            primary,

        "tier1_plus_tier2_mito":
            secondary,
    }.items():

        result = matched_permutation(
            subset,
            results,
            name,
        )

        if result is not None:

            permutation_rows.append(
                result
            )

    permutation = pd.DataFrame(
        permutation_rows
    )

    if len(permutation) > 0:

        permutation[
            "FDR_ASTShigh_sensitive"
        ] = multipletests(
            permutation[
                "p_ASTShigh_sensitive"
            ],
            method="fdr_bh",
        )[1]

    permutation.to_csv(
        OUT_DIR
        / "19C_CTRPv2_permutation.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # Combined three-screen evidence table
    # --------------------------------------------------------

    step19b_file = (
        OUT_DIR
        / "19B_candidate_effects.tsv"
    )

    combined_parts = []

    if step19b_file.exists():

        combined_parts.append(
            pd.read_csv(
                step19b_file,
                sep="\t",
            )
        )

    combined_parts.append(
        candidate_results
    )

    combined = pd.concat(
        combined_parts,
        ignore_index=True,
        sort=False,
    )

    combined.to_csv(
        OUT_DIR
        / "19BC_candidate_evidence_all_screens.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # Figure
    # --------------------------------------------------------

    if len(combined) > 0:

        plot = combined[
            combined[
                "standardized_beta"
            ].notna()
        ].copy()

        plot["label"] = (
            plot[
                "canonical_name"
            ].astype(str)
            +
            " ["
            +
            plot[
                "screen"
            ].astype(str)
            +
            "]"
        )

        plot = plot.sort_values(
            "standardized_beta"
        )

        fig, ax = plt.subplots(
            figsize=(9, 8)
        )

        y = np.arange(
            len(plot)
        )

        ax.barh(
            y,
            plot[
                "standardized_beta"
            ],
        )

        ax.axvline(
            0,
            linestyle="--",
        )

        ax.set_yticks(y)

        ax.set_yticklabels(
            plot["label"]
        )

        ax.set_xlabel(
            "ASTS standardized beta\n"
            "negative = ASTS-high greater sensitivity"
        )

        ax.set_title(
            "Mitochondrial drug effects across screens"
        )

        fig.tight_layout()

        fig.savefig(
            FIG_DIR
            / "19BC_mito_drugs_all_screens.png",
            dpi=200,
        )

        plt.close(fig)

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    summary_file = (
        OUT_DIR
        / "19C_summary.txt"
    )

    with open(
        summary_file,
        "w",
    ) as f:

        print(
            "ASTRA-Drug Step 19C summary",
            file=f,
        )

        print(
            "============================",
            file=f,
        )

        print(
            "",
            file=f,
        )

        print(
            "CTRPv2 sample mapping",
            file=f,
        )

        print(
            "---------------------",
            file=f,
        )

        print(
            mapping[
                "mapping_status"
            ]
            .value_counts()
            .to_string(),
            file=f,
        )

        print(
            "",
            file=f,
        )

        print(
            "Candidate availability",
            file=f,
        )

        print(
            "----------------------",
            file=f,
        )

        print(
            availability.to_string(
                index=False
            ),
            file=f,
        )

        print(
            "",
            file=f,
        )

        print(
            "Candidate effects",
            file=f,
        )

        print(
            "-----------------",
            file=f,
        )

        print(
            candidate_results.to_string(
                index=False
            ),
            file=f,
        )

        print(
            "",
            file=f,
        )

        print(
            "Matched permutation",
            file=f,
        )

        print(
            "-------------------",
            file=f,
        )

        print(
            permutation.to_string(
                index=False
            ),
            file=f,
        )

    print(
        "\nStep 19C completed successfully"
    )

    print(summary_file)


if __name__ == "__main__":
    main()
