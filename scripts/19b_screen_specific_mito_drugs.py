#!/usr/bin/env python3

"""
ASTRA-Drug
Step 19B: Screen-specific mitochondrial drug validation

----




----

    19_MITO_DRUG_SET_FREEZE.tsv


------

Drug response
    ~ ASTS
    + lineage
    + sex
    + proliferation proxy
    + MAPK / ERBB / PI3K driver
    + E2F/G2M CellCycle
    + MYC

PRISM:
    response = log2(AUC)
    lower = more sensitive

GDSC2:
    response = LN_IC50
    lower = more sensitive


    ASTS beta < 0


"""


# ============================================================
# 1. Path / Parameter
# ============================================================

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]

FREEZE_FILE = (
    PROJECT_DIR
    / "results/19_mito_drug_validation"
    / "19_MITO_DRUG_SET_FREEZE.tsv"
)

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
RANDOM_SEED = 20260927

MIN_PRIMARY_DRUGS_FOR_INFERENCE = 3


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
# 3. Helper functions
# ============================================================

def normalize_drug_name(x):

    if pd.isna(x):
        return ""

    return re.sub(
        r"[^a-z0-9]+",
        "",
        str(x).lower().strip(),
    )


def normalize_id(x):

    if pd.isna(x):
        return ""

    x = str(x).strip()

    if re.fullmatch(r"\d+\.0", x):
        x = x[:-2]

    return x


def zscore(series):

    x = pd.to_numeric(
        series,
        errors="coerce",
    )

    sd = x.std(ddof=1)

    if not np.isfinite(sd) or sd == 0:

        return pd.Series(
            np.nan,
            index=x.index,
        )

    return (
        x - x.mean()
    ) / sd


def build_lineage_categories(metadata):

    counts = (
        metadata[LINEAGE_COLUMN]
        .fillna("Unknown")
        .value_counts()
    )

    common = set(
        counts[
            counts >= MIN_GLOBAL_LINEAGE_N
        ].index
    )

    return (
        metadata[LINEAGE_COLUMN]
        .fillna("Unknown")
        .astype(str)
        .apply(
            lambda x:
                x
                if x in common
                else "OtherRare"
        )
    )


def load_candidate_freeze():

    freeze = pd.read_csv(
        FREEZE_FILE,
        sep="\t",
    )

    required = [
        "canonical_name",
        "tier",
        "mechanism_group",
        "aliases",
    ]

    missing = [
        x
        for x in required
        if x not in freeze.columns
    ]

    if missing:

        raise RuntimeError(
            "Freeze file columns missing: "
            + ", ".join(missing)
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
                alias_to_canonical[key] = (
                    row["canonical_name"]
                )

    return (
        freeze,
        alias_to_canonical,
    )


def match_candidate_name(
    drug_name,
    alias_to_canonical,
):
    """

        rotenone
        dihydrorotenone
    """

    key = normalize_drug_name(
        drug_name
    )

    return alias_to_canonical.get(
        key,
        None
    ) 


# ============================================================
# 4. Metadata
# ============================================================

def prepare_metadata():

    asts = pd.read_csv(
        ASTS_FILE,
        sep="\t",
    )

    if asts["ModelID"].duplicated().any():

        raise RuntimeError(
            "Duplicated ModelID in ASTS table."
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
    ].copy()

    cellcycle = pd.read_csv(
        CELL_CYCLE_FILE,
        sep="\t",
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
            validate="1:1",
        )
        .merge(
            cellcycle,
            on="ModelID",
            how="left",
            validate="1:1",
        )
    )

    return metadata


# ============================================================
# 5. Regression
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

    p = float(
        model.pvalues["ASTS_z"]
    )

    return {
        "n_models":
            len(work),

        "n_lineages":
            work[
                "LineageModel"
            ].nunique(),

        "response_sd":
            response_sd,

        "beta":
            beta,

        "SE":
            se,

        "p":
            p,

        "standardized_beta":
            beta / response_sd,

        "standardized_SE":
            se / response_sd,
    }


# ============================================================
# 6. PRISM
# ============================================================

def prepare_prism(metadata):

    print(
        "\nPreparing PRISM..."
    )

    prism = pd.read_csv(
        PRISM_FILE,
        low_memory=False,
    )

    required = [
        "depmap_id",
        "auc",
        "name",
    ]

    for col in required:

        if col not in prism.columns:
            raise RuntimeError(
                f"PRISM column missing: {col}"
            )

    prism["auc"] = pd.to_numeric(
        prism["auc"],
        errors="coerce",
    )

    prism = prism[
        prism["auc"] > 0
    ].copy()

    prism["response"] = np.log2(
        prism["auc"]
    )

    prism["compound_name"] = (
        prism["name"]
        .astype(str)
        .str.strip()
    )

    prism["drug_key"] = (
        prism["compound_name"]
        .apply(
            normalize_drug_name
        )
    )

    prism = (
        prism
        .groupby(
            [
                "drug_key",
                "depmap_id",
            ],
            as_index=False,
        )
        .agg(
            compound_name=(
                "compound_name",
                "first",
            ),
            response=(
                "response",
                "median",
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
        validate="m:1",
    )

    return prism


# ============================================================
# 7. GDSC2
# ============================================================

def choose_first_existing(
    columns,
    candidates,
):

    for col in candidates:

        if col in columns:
            return col

    return None


def prepare_gdsc2(metadata):

    print(
        "\nPreparing GDSC2..."
    )

    gdsc = pd.read_excel(
        GDSC2_FILE,
        engine="openpyxl",
    )

    sanger_col = choose_first_existing(
        gdsc.columns,
        [
            "SANGER_MODEL_ID",
            "SangerModelID",
        ],
    )

    drug_id_col = choose_first_existing(
        gdsc.columns,
        [
            "DRUG_ID",
            "Drug ID",
        ],
    )

    drug_name_col = choose_first_existing(
        gdsc.columns,
        [
            "DRUG_NAME",
            "Drug Name",
            "drug_name",
        ],
    )

    response_col = choose_first_existing(
        gdsc.columns,
        [
            "LN_IC50",
            "ln_ic50",
        ],
    )

    if None in [
        sanger_col,
        drug_id_col,
        response_col,
    ]:

        raise RuntimeError(
            "Could not detect required "
            "GDSC2 columns."
        )

    if drug_name_col is None:

        print(
            "WARNING: GDSC2 drug-name "
            "column not found."
        )

        gdsc["__drug_name"] = (
            gdsc[drug_id_col]
            .astype(str)
        )

        drug_name_col = "__drug_name"

    gdsc["response"] = pd.to_numeric(
        gdsc[response_col],
        errors="coerce",
    )

    gdsc = gdsc[
        gdsc["response"].notna()
    ].copy()

    gdsc["drug_id"] = (
        gdsc[drug_id_col]
        .apply(normalize_id)
    )

    gdsc["compound_name"] = (
        gdsc[drug_name_col]
        .astype(str)
        .str.strip()
    )

    gdsc["SANGER_MODEL_ID_norm"] = (
        gdsc[sanger_col]
        .astype(str)
        .str.strip()
    )

    gdsc = (
        gdsc
        .groupby(
            [
                "drug_id",
                "SANGER_MODEL_ID_norm",
            ],
            as_index=False,
        )
        .agg(
            compound_name=(
                "compound_name",
                "first",
            ),
            response=(
                "response",
                "median",
            ),
        )
    )

    model = pd.read_csv(
        MODEL_FILE,
        low_memory=False,
    )

    if "SangerModelID" not in model.columns:

        raise RuntimeError(
            "SangerModelID not found "
            "in Model.csv."
        )

    model_map = model[
        [
            "ModelID",
            "SangerModelID",
        ]
    ].dropna().copy()

    model_map["SangerModelID"] = (
        model_map["SangerModelID"]
        .astype(str)
        .str.strip()
    )

    # Only one-to-one Sanger -> ModelID mapping
    counts = (
        model_map
        .groupby("SangerModelID")[
            "ModelID"
        ]
        .nunique()
    )

    unique_sanger = set(
        counts[
            counts == 1
        ].index
    )

    model_map = (
        model_map[
            model_map[
                "SangerModelID"
            ].isin(unique_sanger)
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
        validate="m:1",
    )

    gdsc = gdsc.merge(
        metadata,
        on="ModelID",
        how="inner",
        validate="m:1",
    )

    return gdsc


# ============================================================
# 8. All-drug analysis
# ============================================================

def analyze_all_drugs(
    data,
    screen_name,
):

    rows = []

    grouping_col = (
        "drug_key"
        if screen_name == "PRISM"
        else "drug_id"
    )

    for drug_key, subset in (
        data.groupby(
            grouping_col
        )
    ):

        result = fit_one_drug(
            subset
        )

        if result is None:
            continue

        compound_name = (
            subset[
                "compound_name"
            ]
            .iloc[0]
        )

        rows.append(
            {
                "screen":
                    screen_name,

                "drug_key":
                    str(drug_key),

                "compound_name":
                    compound_name,

                **result,
            }
        )

    results = pd.DataFrame(
        rows
    )

    if len(results) == 0:
        raise RuntimeError(
            f"No drugs tested in {screen_name}."
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

    return results


# ============================================================
# 9. Candidate audit
# ============================================================

def candidate_audit(
    prepared_data,
    results,
    freeze,
    alias_to_canonical,
    screen_name,
):

    prepared_names = (
        prepared_data[
            "compound_name"
        ]
        .dropna()
        .astype(str)
        .unique()
    )

    raw_candidate_names = {}

    for name in prepared_names:

        canonical = match_candidate_name(
            name,
            alias_to_canonical,
        )

        if canonical is not None:

            raw_candidate_names.setdefault(
                canonical,
                []
            ).append(name)

    tested = results.copy()

    tested["canonical_name"] = (
        tested[
            "compound_name"
        ]
        .apply(
            lambda x:
                match_candidate_name(
                    x,
                    alias_to_canonical,
                )
        )
    )

    tested_candidates = tested[
        tested[
            "canonical_name"
        ].notna()
    ].copy()

    # If multiple entries match the same drug,
    # retain the entry with the largest n.
    tested_candidates = (
        tested_candidates
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
            "canonical_name",
            keep="first",
        )
    )

    if len(tested_candidates) > 0:

        tested_candidates[
            "FDR_candidate_set"
        ] = multipletests(
            tested_candidates["p"],
            method="fdr_bh",
        )[1]

    availability_rows = []

    for _, row in freeze.iterrows():

        canonical = row[
            "canonical_name"
        ]

        found_raw = (
            canonical
            in raw_candidate_names
        )

        tested_row = (
            tested_candidates[
                tested_candidates[
                    "canonical_name"
                ]
                ==
                canonical
            ]
        )

        if len(tested_row) > 0:

            status = "tested"

            matched_name = (
                tested_row.iloc[0][
                    "compound_name"
                ]
            )

            n_models = int(
                tested_row.iloc[0][
                    "n_models"
                ]
            )

        elif found_raw:

            status = (
                "found_but_not_tested"
            )

            matched_name = ";".join(
                sorted(
                    set(
                        raw_candidate_names[
                            canonical
                        ]
                    )
                )
            )

            n_models = np.nan

        else:

            status = "not_found"
            matched_name = ""
            n_models = np.nan

        availability_rows.append(
            {
                "screen":
                    screen_name,

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

    if len(tested_candidates) > 0:

        tested_candidates = (
            tested_candidates
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

    return (
        availability,
        tested_candidates,
    )


# ============================================================
# 10. Matched drug-set permutation
# ============================================================

def matched_permutation(
    candidate_results,
    all_results,
    set_name,
    rng,
):

    if len(candidate_results) == 0:

        return None

    candidate_names = set(
        candidate_results[
            "compound_name"
        ]
    )

    background = all_results[
        ~all_results[
            "compound_name"
        ].isin(
            candidate_names
        )
    ].copy()

    match_data = pd.concat(
        [
            candidate_results,
            background,
        ],
        ignore_index=True,
    )

    for col in [
        "n_models",
        "response_sd",
    ]:

        transformed = np.log1p(
            pd.to_numeric(
                match_data[col],
                errors="coerce",
            )
        )

        sd = transformed.std(
            ddof=1
        )

        match_data[
            f"match_{col}"
        ] = (
            transformed
            -
            transformed.mean()
        ) / sd

    candidates_match = (
        match_data[
            match_data[
                "compound_name"
            ].isin(candidate_names)
        ]
        .copy()
    )

    controls_match = (
        match_data[
            ~match_data[
                "compound_name"
            ].isin(candidate_names)
        ]
        .copy()
    )

    feature_cols = [
        "match_n_models",
        "match_response_sd",
    ]

    control_matrix = (
        controls_match[
            feature_cols
        ]
        .to_numpy(
            dtype=float
        )
    )

    control_names = (
        controls_match[
            "compound_name"
        ]
        .to_numpy(
            dtype=object
        )
    )

    pools = {}

    for _, row in (
        candidates_match.iterrows()
    ):

        target = row[
            feature_cols
        ].to_numpy(
            dtype=float
        )

        distance = np.sum(
            (
                control_matrix
                -
                target[None, :]
            )
            ** 2,
            axis=1,
        )

        k = min(
            MATCH_POOL_SIZE,
            len(control_names),
        )

        idx = np.argpartition(
            distance,
            k - 1,
        )[:k]

        idx = idx[
            np.argsort(
                distance[idx]
            )
        ]

        pools[
            row["compound_name"]
        ] = [
            control_names[i]
            for i in idx
        ]

    background_index = (
        background
        .set_index(
            "compound_name"
        )
    )

    observed = float(
        candidate_results[
            "standardized_beta"
        ].mean()
    )

    null = np.empty(
        N_PERMUTATIONS,
        dtype=float,
    )

    candidate_order_original = list(
        candidate_results[
            "compound_name"
        ]
    )

    all_background = list(
        background_index.index
    )

    for i in range(
        N_PERMUTATIONS
    ):

        order = (
            candidate_order_original.copy()
        )

        rng.shuffle(order)

        selected = set()
        sampled = []

        for drug in order:

            available = [
                x
                for x in pools[drug]
                if x not in selected
            ]

            if available:

                chosen = rng.choice(
                    available
                )

            else:

                remaining = [
                    x
                    for x in all_background
                    if x not in selected
                ]

                chosen = rng.choice(
                    remaining
                )

            selected.add(chosen)
            sampled.append(chosen)

        null[i] = (
            background_index
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

    p_resistant = (
        1
        +
        np.sum(
            null >= observed
        )
    ) / (
        N_PERMUTATIONS + 1
    )

    null_mean = float(
        np.mean(null)
    )

    null_sd = float(
        np.std(
            null,
            ddof=1,
        )
    )

    matched_z = (
        (
            observed - null_mean
        ) / null_sd
        if null_sd > 0
        else np.nan
    )

    return {
        "drug_set":
            set_name,

        "n_drugs":
            len(
                candidate_results
            ),

        "mean_standardized_beta":
            observed,

        "fraction_negative":
            float(
                np.mean(
                    candidate_results[
                        "standardized_beta"
                    ]
                    <
                    0
                )
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
    }


# ============================================================
# 11. Screen summary
# ============================================================

def summarize_screen(
    screen_name,
    all_results,
    candidates,
):

    rows = []

    rng = np.random.default_rng(
        RANDOM_SEED
        +
        (
            0
            if screen_name == "PRISM"
            else 1000
        )
    )

    primary = candidates[
        candidates[
            "tier"
        ]
        ==
        "tier1_direct_oxphos"
    ].copy()

    secondary = candidates.copy()

    for set_name, subset in {
        "tier1_direct_oxphos":
            primary,

        "tier1_plus_tier2_mito":
            secondary,
    }.items():

        result = matched_permutation(
            subset,
            all_results,
            set_name,
            rng,
        )

        if result is not None:

            result["screen"] = (
                screen_name
            )

            rows.append(result)

    return pd.DataFrame(rows)


# ============================================================
# 12. Main
# ============================================================

def main():

    print(
        "\n"
        "============================================================\n"
        "ASTRA-Drug Step 19B\n"
        "Screen-specific mitochondrial drug validation\n"
        "============================================================"
    )

    required_files = [
        FREEZE_FILE,
        PRISM_FILE,
        GDSC2_FILE,
        MODEL_FILE,
        ASTS_FILE,
        DRIVER_FILE,
        CELL_CYCLE_FILE,
    ]

    for path in required_files:

        if not path.exists():

            raise FileNotFoundError(
                f"File not found:\n{path}"
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

    prism = prepare_prism(
        metadata
    )

    gdsc = prepare_gdsc2(
        metadata
    )

    prism_results = (
        analyze_all_drugs(
            prism,
            "PRISM",
        )
    )

    gdsc_results = (
        analyze_all_drugs(
            gdsc,
            "GDSC2",
        )
    )

    prism_results.to_csv(
        OUT_DIR
        / "19B_PRISM_all_drug_effects.tsv",
        sep="\t",
        index=False,
    )

    gdsc_results.to_csv(
        OUT_DIR
        / "19B_GDSC2_all_drug_effects.tsv",
        sep="\t",
        index=False,
    )

    p_avail, p_candidates = (
        candidate_audit(
            prism,
            prism_results,
            freeze,
            alias_map,
            "PRISM",
        )
    )

    g_avail, g_candidates = (
        candidate_audit(
            gdsc,
            gdsc_results,
            freeze,
            alias_map,
            "GDSC2",
        )
    )

    availability = pd.concat(
        [
            p_avail,
            g_avail,
        ],
        ignore_index=True,
    )

    availability.to_csv(
        OUT_DIR
        / "19B_candidate_availability.tsv",
        sep="\t",
        index=False,
    )

    candidate_effects = pd.concat(
        [
            p_candidates,
            g_candidates,
        ],
        ignore_index=True,
    )

    candidate_effects.to_csv(
        OUT_DIR
        / "19B_candidate_effects.tsv",
        sep="\t",
        index=False,
    )

    permutations = pd.concat(
        [
            summarize_screen(
                "PRISM",
                prism_results,
                p_candidates,
            ),
            summarize_screen(
                "GDSC2",
                gdsc_results,
                g_candidates,
            ),
        ],
        ignore_index=True,
    )

    if len(permutations) > 0:

        permutations[
            "FDR_ASTShigh_sensitive"
        ] = multipletests(
            permutations[
                "p_ASTShigh_sensitive"
            ],
            method="fdr_bh",
        )[1]

    permutations.to_csv(
        OUT_DIR
        / "19B_screen_specific_permutation.tsv",
        sep="\t",
        index=False,
    )

    # --------------------------------------------------------
    # Figure
    # --------------------------------------------------------

    if len(candidate_effects) > 0:

        plot_data = (
            candidate_effects
            .copy()
            .sort_values(
                "standardized_beta"
            )
        )

        labels = (
            plot_data[
                "canonical_name"
            ]
            +
            " ["
            +
            plot_data["screen"]
            +
            "]"
        )

        fig, ax = plt.subplots(
            figsize=(9, 7)
        )

        y = np.arange(
            len(plot_data)
        )

        ax.barh(
            y,
            plot_data[
                "standardized_beta"
            ],
        )

        ax.axvline(
            0,
            linestyle="--",
        )

        ax.set_yticks(y)
        ax.set_yticklabels(labels)

        ax.set_xlabel(
            "ASTS standardized beta\n"
            "negative = ASTS-high greater sensitivity"
        )

        ax.set_title(
            "Screen-specific mitochondrial drug effects"
        )

        fig.tight_layout()

        fig.savefig(
            FIG_DIR
            / "19B_screen_specific_candidate_effects.png",
            dpi=200,
        )

        plt.close(fig)

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    summary_file = (
        OUT_DIR
        / "19B_summary.txt"
    )

    with open(
        summary_file,
        "w",
    ) as f:

        print(
            "ASTRA-Drug Step 19B summary",
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
            candidate_effects.to_string(
                index=False
            ),
            file=f,
        )

        print(
            "",
            file=f,
        )

        print(
            "Matched permutations",
            file=f,
        )

        print(
            "--------------------",
            file=f,
        )

        print(
            permutations.to_string(
                index=False
            ),
            file=f,
        )

    print(
        "\nStep 19B completed successfully"
    )

    print(summary_file)


if __name__ == "__main__":
    main()
