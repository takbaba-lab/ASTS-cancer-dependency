#!/usr/bin/env python3

"""
ASTRA-Drug
Step 13: Oncogenic-driver confounder analysis

----

    ASTS-high
        -> MEK inhibitor resistance
        -> EGFR inhibitor resistance


    KRAS / NRAS / BRAF
    EGFR / ERBB2
    PIK3CA


Primary driver covariates
-------------------------
MAPK_driver:
    KRAS or NRAS or BRAF hotspot mutation

ERBB_driver:
    EGFR or ERBB2 hotspot mutation

PI3K_driver:
    PIK3CA hotspot mutation


Baseline:
    response ~ ASTS + Lineage + Sex + Proliferation

Driver-adjusted:
    response ~ ASTS + Lineage + Sex + Proliferation
                     + MAPK_driver
                     + ERBB_driver
                     + PI3K_driver



    ASTS ~ Lineage + Sex + Proliferation + Drivers


----

"""


# ============================================================
# Path / Parameter
# ============================================================

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]


# ------------------------------------------------------------
# DepMap ASTS
# ------------------------------------------------------------

ASTS_FILE = (
    PROJECT_DIR
    / "results/08b_DepMap_ASTS_default_profile"
    / "08_DepMap_ASTS_scores.tsv"
)


# ------------------------------------------------------------
# DepMap mutation hotspot matrix
# ------------------------------------------------------------

HOTSPOT_FILE = (
    PROJECT_DIR
    / "data/raw/depmap/26Q1"
    / "OmicsSomaticMutationsMatrixHotspot.csv"
)


# ------------------------------------------------------------
# DepMap model metadata
# ------------------------------------------------------------

MODEL_FILE = (
    PROJECT_DIR
    / "data/raw/depmap/26Q1"
    / "Model.csv"
)


# ------------------------------------------------------------
# PRISM Secondary
# ------------------------------------------------------------

PRISM_FILE = (
    PROJECT_DIR
    / "data/raw/prism/secondary"
    / "secondary-screen-dose-response-curve-parameters.csv"
)


# ------------------------------------------------------------
# GDSC2
# ------------------------------------------------------------

GDSC2_FILE = (
    PROJECT_DIR
    / "data/raw/gdsc/GDSC2"
    / "GDSC2_fitted_dose_response_27Oct23.xlsx"
)


# ------------------------------------------------------------
# Step11 overlap
# ------------------------------------------------------------

OVERLAP_FILE = (
    PROJECT_DIR
    / "results/11_cross_screen_mechanism"
    / "11_PRISM_GDSC2_overlap_annotated.tsv"
)


# ------------------------------------------------------------
# Step12 class membership
# ------------------------------------------------------------

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
    / "results/13_driver_confounder"
)

FIG_DIR = OUT_DIR / "figures"


# ============================================================
# Parameters
# ============================================================

DRIVER_GENES = [
    "KRAS",
    "NRAS",
    "BRAF",
    "EGFR",
    "ERBB2",
    "PIK3CA",
]


MIN_N_PER_DRUG = 100

MIN_GLOBAL_LINEAGE_N = 5

MIN_RESIDUAL_DF = 30

EXPECTED_OVERLAP_N = 108


PRIMARY_ASTS_COLUMN = "ASTS_score"

PROLIFERATION_COLUMN = "proliferation_proxy"

LINEAGE_COLUMN = "OncotreeLineage"

SEX_COLUMN = "Sex"


# ============================================================
# Imports
# ============================================================

import hashlib
import re

import numpy as np
import pandas as pd

import matplotlib.pyplot as plt

import statsmodels.api as sm

from scipy.stats import spearmanr


# ============================================================
# Helper functions
# ============================================================

def sha256_file(path, block_size=1024 * 1024):

    h = hashlib.sha256()

    with open(path, "rb") as f:

        while True:

            block = f.read(block_size)

            if not block:
                break

            h.update(block)

    return h.hexdigest()


def normalize_drug_name(name):

    if pd.isna(name):
        return ""

    return re.sub(
        r"[^a-z0-9]+",
        "",
        str(name).lower().strip(),
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


def detect_model_id_column(df):

    candidates = [
        "ModelID",
        "DepMap_ID",
        "DepMapID",
        "model_id",
    ]

    for col in candidates:

        if col in df.columns:
            return col

    raise ValueError(
        "Could not detect ModelID column in mutation matrix."
    )


def find_gene_column(
    columns,
    gene
):
    """

        KRAS
        KRAS (3845)

    """

    exact = [
        col
        for col in columns
        if str(col) == gene
    ]

    if len(exact) == 1:
        return exact[0]


    pattern = re.compile(
        rf"^{re.escape(gene)}(?:\s|\(|$)",
        flags=re.IGNORECASE,
    )


    hits = [
        col
        for col in columns
        if pattern.search(
            str(col)
        )
    ]


    if len(hits) != 1:

        raise ValueError(
            f"Could not uniquely identify mutation column "
            f"for {gene}. Hits={hits}"
        )


    return hits[0]


def mutation_to_binary(series):
    """
        0/1
        True/False
        Y/N
    """

    if pd.api.types.is_numeric_dtype(
        series
    ):

        x = pd.to_numeric(
            series,
            errors="coerce",
        )

        return (
            x.fillna(0) > 0
        ).astype(int)


    x = (
        series
        .astype(str)
        .str.strip()
        .str.lower()
    )


    positive = {
        "1",
        "true",
        "t",
        "yes",
        "y",
    }


    return x.isin(
        positive
    ).astype(int)


def build_lineage_model(
    asts
):
    """
    """

    metadata = (
        asts[
            [
                "ModelID",
                LINEAGE_COLUMN,
            ]
        ]
        .drop_duplicates()
    )


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


    def convert(x):

        if pd.isna(x):
            return "Unknown"

        x = str(x)

        if x in common:
            return x

        return "OtherRare"


    return asts[
        LINEAGE_COLUMN
    ].apply(
        convert
    )


def fit_response_model(
    df,
    include_drivers,
):
    """
    Drug response model。

    """

    work = df.copy()


    required = [
        "response",
        PRIMARY_ASTS_COLUMN,
        PROLIFERATION_COLUMN,
        "LineageModel",
    ]


    work = work.dropna(
        subset=required
    ).copy()


    work[
        "ASTS_z_model"
    ] = zscore(
        work[
            PRIMARY_ASTS_COLUMN
        ]
    )


    work[
        "Prolif_z_model"
    ] = zscore(
        work[
            PROLIFERATION_COLUMN
        ]
    )


    work = work.dropna(
        subset=[
            "ASTS_z_model",
            "Prolif_z_model",
        ]
    ).copy()


    if len(work) < MIN_N_PER_DRUG:
        return None


    X_parts = [
        work[
            [
                "ASTS_z_model",
                "Prolif_z_model",
            ]
        ]
    ]


    # --------------------------------------------------------
    # Sex
    # --------------------------------------------------------

    sex = (
        work[
            SEX_COLUMN
        ]
        .fillna("Unknown")
        .astype(str)
    )


    if sex.nunique() > 1:

        X_parts.append(
            pd.get_dummies(
                sex,
                prefix="Sex",
                drop_first=True,
                dtype=float,
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

        X_parts.append(
            pd.get_dummies(
                lineage,
                prefix="Lineage",
                drop_first=True,
                dtype=float,
            )
        )


    # --------------------------------------------------------
    # Driver covariates
    # --------------------------------------------------------

    used_drivers = []


    if include_drivers:

        for col in [
            "MAPK_driver",
            "ERBB_driver",
            "PI3K_driver",
        ]:

            if (
                col in work.columns
                and
                work[col].nunique(
                    dropna=True
                ) > 1
            ):

                X_parts.append(
                    work[
                        [col]
                    ].astype(float)
                )

                used_drivers.append(
                    col
                )


    X = pd.concat(
        X_parts,
        axis=1,
    )


    X = sm.add_constant(
        X,
        has_constant="add",
    )


    y = pd.to_numeric(
        work[
            "response"
        ],
        errors="coerce",
    )


    try:

        model = sm.OLS(
            y,
            X,
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


    if (
        model.df_resid
        < MIN_RESIDUAL_DF
    ):

        return None


    return {
        "n":
            len(work),

        "beta_ASTS":
            float(
                model.params[
                    "ASTS_z_model"
                ]
            ),

        "se_ASTS":
            float(
                model.bse[
                    "ASTS_z_model"
                ]
            ),

        "p_ASTS":
            float(
                model.pvalues[
                    "ASTS_z_model"
                ]
            ),

        "model_R2":
            float(
                model.rsquared
            ),

        "model_adj_R2":
            float(
                model.rsquared_adj
            ),

        "df_resid":
            float(
                model.df_resid
            ),

        "drivers_used":
            ";".join(
                used_drivers
            ),

        "n_MAPK_driver":
            int(
                work[
                    "MAPK_driver"
                ].sum()
            ),

        "n_ERBB_driver":
            int(
                work[
                    "ERBB_driver"
                ].sum()
            ),

        "n_PI3K_driver":
            int(
                work[
                    "PI3K_driver"
                ].sum()
            ),
    }


def fit_asts_driver_model(
    df
):
    """
    """

    work = df.copy()


    work = work.dropna(
        subset=[
            PRIMARY_ASTS_COLUMN,
            PROLIFERATION_COLUMN,
            "LineageModel",
        ]
    ).copy()


    work[
        "ASTS_z"
    ] = zscore(
        work[
            PRIMARY_ASTS_COLUMN
        ]
    )


    work[
        "Prolif_z"
    ] = zscore(
        work[
            PROLIFERATION_COLUMN
        ]
    )


    # --------------------------------------------------------
    # Common base design
    # --------------------------------------------------------

    X_parts = [
        work[
            [
                "Prolif_z",
            ]
        ]
    ]


    sex = (
        work[
            SEX_COLUMN
        ]
        .fillna("Unknown")
        .astype(str)
    )


    if sex.nunique() > 1:

        X_parts.append(
            pd.get_dummies(
                sex,
                prefix="Sex",
                drop_first=True,
                dtype=float,
            )
        )


    lineage = (
        work[
            "LineageModel"
        ]
        .fillna("Unknown")
        .astype(str)
    )


    if lineage.nunique() > 1:

        X_parts.append(
            pd.get_dummies(
                lineage,
                prefix="Lineage",
                drop_first=True,
                dtype=float,
            )
        )


    X_base = pd.concat(
        X_parts,
        axis=1,
    )


    X_base = sm.add_constant(
        X_base,
        has_constant="add",
    )


    y = work[
        "ASTS_z"
    ]


    base = sm.OLS(
        y,
        X_base,
    ).fit(
        cov_type="HC3"
    )


    # --------------------------------------------------------
    # Driver-adjusted model
    # --------------------------------------------------------

    driver_cols = [
        "MAPK_driver",
        "ERBB_driver",
        "PI3K_driver",
    ]


    X_driver = pd.concat(
        [
            X_base.drop(
                columns=["const"]
            ),
            work[
                driver_cols
            ].astype(float),
        ],
        axis=1,
    )


    X_driver = sm.add_constant(
        X_driver,
        has_constant="add",
    )


    model = sm.OLS(
        y,
        X_driver,
    ).fit(
        cov_type="HC3"
    )


    rows = []


    for driver in driver_cols:

        rows.append(
            {
                "driver":
                    driver,

                "n_positive":
                    int(
                        work[
                            driver
                        ].sum()
                    ),

                "beta_ASTSz":
                    float(
                        model.params[
                            driver
                        ]
                    ),

                "SE":
                    float(
                        model.bse[
                            driver
                        ]
                    ),

                "p":
                    float(
                        model.pvalues[
                            driver
                        ]
                    ),
            }
        )


    return (
        pd.DataFrame(
            rows
        ),
        {
            "n":
                len(work),

            "base_R2":
                float(
                    base.rsquared
                ),

            "driver_R2":
                float(
                    model.rsquared
                ),

            "delta_R2_drivers":
                float(
                    model.rsquared
                    -
                    base.rsquared
                ),
        },
    )


# ============================================================
# Main
# ============================================================


# ============================================================
# Step 13 mutation input repair
# ============================================================
def load_step13_driver_calls(hotspot_file, out_dir, driver_genes, conflict_policy="exclude"):
    """Audit missing and duplicated mutation calls and return one complete six-gene record per ModelID."""
    if conflict_policy not in {"exclude", "error"}:
        raise ValueError("conflict_policy must be exclude or error")
    audit_dir = out_dir / "mutation_input_audit"
    audit_dir.mkdir(parents=True, exist_ok=True)

    mut = pd.read_csv(hotspot_file, dtype="string", low_memory=False)
    model_id_col = detect_model_id_column(mut)
    gene_columns = {
        gene: find_gene_column(mut.columns, gene)
        for gene in driver_genes
    }
    print("Mutation ModelID column:", model_id_col)
    for gene, column in gene_columns.items():
        print(gene, "->", column)

    mut.insert(0, "_audit_record_number", range(2, len(mut) + 2))
    mut[model_id_col] = mut[model_id_col].str.strip()
    ids = mut[model_id_col]
    invalid_ids = ids.isna() | ids.fillna("").eq("")

    mut.loc[invalid_ids].to_csv(
        audit_dir / "invalid_model_id_rows.tsv", sep="\t", index=False
    )

    mut = mut.loc[~invalid_ids].copy()
    duplicate_rows = mut[model_id_col].duplicated(keep=False)
    mut.loc[duplicate_rows].to_csv(
        audit_dir / "duplicate_model_id_rows.tsv", sep="\t", index=False
    )

    driver = pd.DataFrame({"ModelID": mut[model_id_col]})
    unknown_rows = []
    missing_tokens = {"", "na", "nan", "n/a", "null", "none", "<na>"}
    boolean_tokens = {
        "true": 1, "t": 1, "yes": 1, "y": 1,
        "false": 0, "f": 0, "no": 0, "n": 0,
    }

    for gene, column in gene_columns.items():
        raw = mut[column]
        values = raw.str.strip().str.lower()
        missing = values.isna() | values.isin(missing_tokens)
        numeric = pd.to_numeric(values, errors="coerce")

        finite = numeric.notna() & numeric.lt(float("inf"))
        valid_numeric = (
            finite & numeric.ge(0) & numeric.mod(1).eq(0)
        ).fillna(False)
        valid_boolean = values.isin(boolean_tokens)
        unknown = ~(missing | valid_numeric | valid_boolean)

        for index in mut.index[unknown]:
            unknown_rows.append({
                "record_number": mut.at[index, "_audit_record_number"],
                "ModelID": mut.at[index, model_id_col],
                "gene": gene,
                "raw_value": raw.at[index],
            })

        calls = pd.Series(pd.NA, index=mut.index, dtype="Int64")
        calls.loc[valid_numeric] = numeric.loc[valid_numeric].gt(0).astype(int)
        calls.loc[valid_boolean] = values.loc[valid_boolean].map(boolean_tokens)
        driver[gene] = calls

    unknown_table = pd.DataFrame(
        unknown_rows, columns=["record_number", "ModelID", "gene", "raw_value"]
    )
    unknown_table.to_csv(
        audit_dir / "unrecognized_mutation_values.tsv", sep="\t", index=False
    )
    if len(unknown_table):
        raise RuntimeError(
            "Unrecognized mutation values. Check "
            + str(audit_dir / "unrecognized_mutation_values.tsv")
        )

    distinct = driver.groupby("ModelID")[driver_genes].nunique(dropna=True)
    conflicting_ids = distinct.index[distinct.gt(1).any(axis=1)]
    mut.loc[mut[model_id_col].isin(conflicting_ids)].to_csv(
        audit_dir / "conflicting_model_id_rows.tsv", sep="\t", index=False
    )

    incomplete_ids = driver.loc[
        driver[driver_genes].isna().any(axis=1), "ModelID"
    ].unique()
    mut.loc[mut[model_id_col].isin(incomplete_ids)].to_csv(
        audit_dir / "incomplete_model_id_rows.tsv", sep="\t", index=False
    )

    summary = {
        "input_rows": len(mut) + int(invalid_ids.sum()),
        "invalid_model_id_rows": int(invalid_ids.sum()),
        "valid_model_id_rows": len(mut),
        "unique_valid_model_ids": driver["ModelID"].nunique(),
        "duplicate_model_id_rows": int(duplicate_rows.sum()),
        "duplicate_model_ids": mut.loc[duplicate_rows, model_id_col].nunique(),
        "conflicting_model_ids": len(conflicting_ids),
        "incomplete_model_ids": len(incomplete_ids),
    }
    pd.Series(summary, name="count").to_csv(
        audit_dir / "input_summary.tsv", sep="\t", header=True
    )
    print("Mutation input audit:", summary)
    print("Audit directory:", audit_dir)

    excluded = pd.DataFrame({
        "ModelID": sorted(set(conflicting_ids) | set(incomplete_ids))
    })
    excluded["conflicting_calls"] = excluded["ModelID"].isin(conflicting_ids)
    excluded["incomplete_calls"] = excluded["ModelID"].isin(incomplete_ids)
    excluded.to_csv(
        audit_dir / "excluded_model_ids.tsv", sep="\t", index=False
    )
    if len(conflicting_ids) and conflict_policy == "error":
        raise RuntimeError(
            "Conflicting mutation calls across rows for the same ModelID. "
            "A documented profile-selection rule is required. Check "
            + str(audit_dir / "conflicting_model_id_rows.tsv")
        )

    if len(conflicting_ids):
        print("WARNING: Excluding entire ModelIDs with conflicting calls:", len(conflicting_ids))
        print("Analysis is restricted to models with unambiguous driver calls.")

    driver = driver.loc[~driver["ModelID"].isin(excluded["ModelID"])].copy()
    driver = driver.drop_duplicates(subset=["ModelID"])
    if driver.empty:
        raise RuntimeError("No complete mutation records remain. Check mutation_input_audit.")
    driver[driver_genes] = driver[driver_genes].astype(int)
    driver = driver.reset_index(drop=True)
    driver.to_csv(
        audit_dir / "resolved_driver_gene_calls.tsv", sep="\t", index=False
    )
    print("Mutation models retained:", len(driver))
    return driver


def main():

    print(
        "\n"
        "============================================================\n"
        "ASTRA-Drug Step 13\n"
        "Oncogenic-driver confounder analysis\n"
        "============================================================"
    )


    # ========================================================
    # 0. Input check
    # ========================================================

    required_files = [
        ASTS_FILE,
        HOTSPOT_FILE,
        MODEL_FILE,
        PRISM_FILE,
        GDSC2_FILE,
        OVERLAP_FILE,
        CLASS_FILE,
    ]


    for path in required_files:

        if not path.exists():

            raise FileNotFoundError(
                f"Required file not found:\n{path}"
            )


    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )


    FIG_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )


    # ========================================================
    # 1. Read ASTS
    # ========================================================

    print(
        "\nStep 1: Reading ASTS"
    )


    asts = pd.read_csv(
        ASTS_FILE,
        sep="\t",
    )


    if asts[
        "ModelID"
    ].duplicated().any():

        raise RuntimeError(
            "Duplicated ModelID in ASTS file."
        )


    asts[
        "LineageModel"
    ] = build_lineage_model(
        asts
    )


    # ========================================================
    # 2. Read hotspot mutations
    # ========================================================

    print(
        "\nStep 2: Reading hotspot mutation matrix"
    )


    driver = load_step13_driver_calls(
        HOTSPOT_FILE, OUT_DIR, DRIVER_GENES, conflict_policy='exclude'
    )


    # --------------------------------------------------------
    # Composite driver flags
    # --------------------------------------------------------

    driver[
        "MAPK_driver"
    ] = (
        driver[
            [
                "KRAS",
                "NRAS",
                "BRAF",
            ]
        ]
        .max(
            axis=1
        )
    )


    driver[
        "ERBB_driver"
    ] = (
        driver[
            [
                "EGFR",
                "ERBB2",
            ]
        ]
        .max(
            axis=1
        )
    )


    driver[
        "PI3K_driver"
    ] = driver[
        "PIK3CA"
    ]


    driver.to_csv(
        OUT_DIR
        / "13_driver_flags.tsv",
        sep="\t",
        index=False,
    )


    # ========================================================
    # 3. Merge ASTS + mutations
    # ========================================================

    asts_driver = asts.merge(
        driver,
        on="ModelID",
        how="inner",
        validate="1:1",
    )


    print(
        "ASTS models:",
        len(asts)
    )

    print(
        "Mutation-covered ASTS models:",
        len(asts_driver)
    )


    # --------------------------------------------------------
    # Driver prevalence
    # --------------------------------------------------------

    prevalence_rows = []


    for gene in (
        DRIVER_GENES
        +
        [
            "MAPK_driver",
            "ERBB_driver",
            "PI3K_driver",
        ]
    ):

        prevalence_rows.append(
            {
                "driver":
                    gene,

                "n_positive":
                    int(
                        asts_driver[
                            gene
                        ].sum()
                    ),

                "n_total":
                    len(
                        asts_driver
                    ),

                "fraction":
                    asts_driver[
                        gene
                    ].mean(),
            }
        )


    pd.DataFrame(
        prevalence_rows
    ).to_csv(
        OUT_DIR
        / "13_driver_prevalence.tsv",
        sep="\t",
        index=False,
    )


    # ========================================================
    # 4. Is ASTS itself predicted by driver genotype?
    # ========================================================

    print(
        "\nStep 3: ASTS ~ driver genotype"
    )


    (
        asts_driver_coef,
        asts_driver_stats,
    ) = fit_asts_driver_model(
        asts_driver
    )


    asts_driver_coef.to_csv(
        OUT_DIR
        / "13_ASTS_vs_driver_coefficients.tsv",
        sep="\t",
        index=False,
    )


    # ========================================================
    # 5. Read Step11/12 definitions
    # ========================================================

    overlap = pd.read_csv(
        OVERLAP_FILE,
        sep="\t",
    )


    if len(overlap) != EXPECTED_OVERLAP_N:

        raise RuntimeError(
            f"Expected {EXPECTED_OVERLAP_N} overlap drugs, "
            f"found {len(overlap)}."
        )


    classes = pd.read_csv(
        CLASS_FILE,
        sep="\t",
    )


    if (
        classes[
            "compound_name"
        ].duplicated().any()
    ):

        raise RuntimeError(
            "Duplicated compound_name in Step12 class file."
        )


    overlap = overlap.merge(
        classes,
        on="compound_name",
        how="left",
        suffixes=(
            "",
            "_step12",
        ),
        validate="1:1",
    )


    overlap[
        "drug_key"
    ] = overlap[
        "drug_key"
    ].astype(str)


    overlap[
        "DRUG_ID_norm"
    ] = overlap[
        "DRUG_ID"
    ].apply(
        normalize_drug_id
    )


    # ========================================================
    # 6. Prepare PRISM response
    # ========================================================

    print(
        "\nStep 4: Preparing PRISM"
    )


    prism = pd.read_csv(
        PRISM_FILE,
        low_memory=False,
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
        errors="coerce",
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
            as_index=False,
        )
        .agg(
            response=(
                "response",
                "median",
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
        asts_driver,
        on="ModelID",
        how="inner",
        validate="m:1",
    )


    # ========================================================
    # 7. Prepare GDSC2 response
    # ========================================================

    print(
        "\nStep 5: Preparing GDSC2"
    )


    gdsc = pd.read_excel(
        GDSC2_FILE,
        engine="openpyxl",
    )


    gdsc = gdsc[
        [
            "SANGER_MODEL_ID",
            "DRUG_ID",
            "DRUG_NAME",
            "LN_IC50",
        ]
    ].copy()


    gdsc[
        "response"
    ] = pd.to_numeric(
        gdsc[
            "LN_IC50"
        ],
        errors="coerce",
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


    # --------------------------------------------------------
    # Sanger -> ModelID
    # --------------------------------------------------------

    model = pd.read_csv(
        MODEL_FILE,
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


    map_counts = (
        model_map
        .groupby(
            "SangerModelID"
        )[
            "ModelID"
        ]
        .nunique()
    )


    good_sanger = set(
        map_counts[
            map_counts == 1
        ].index
    )


    model_map = (
        model_map[
            model_map[
                "SangerModelID"
            ].isin(
                good_sanger
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
            as_index=False,
        )
        .agg(
            response=(
                "response",
                "median",
            )
        )
    )


    gdsc = gdsc.merge(
        model_map,
        left_on="SANGER_MODEL_ID",
        right_on="SangerModelID",
        how="inner",
        validate="m:1",
    )


    gdsc = gdsc.merge(
        asts_driver,
        on="ModelID",
        how="inner",
        validate="m:1",
    )


    # ========================================================
    # 8. Fit all 108 overlapping drugs
    # ========================================================

    print(
        "\nStep 6: Driver-adjusted drug models"
    )


    result_rows = []


    for _, drug in overlap.iterrows():

        drug_name = drug[
            "compound_name"
        ]

        drug_key = drug[
            "drug_key"
        ]

        gdsc_id = drug[
            "DRUG_ID_norm"
        ]


        # ----------------------------------------------------
        # PRISM
        # ----------------------------------------------------

        p = prism[
            prism[
                "drug_key"
            ] == drug_key
        ].copy()


        if len(p) >= MIN_N_PER_DRUG:

            base = fit_response_model(
                p,
                include_drivers=False,
            )

            adjusted = fit_response_model(
                p,
                include_drivers=True,
            )


            if (
                base is not None
                and adjusted is not None
            ):

                result_rows.append(
                    {
                        "screen":
                            "PRISM",

                        "compound_name":
                            drug_name,

                        "drug_key":
                            drug_key,

                        "DRUG_ID_norm":
                            gdsc_id,

                        "n":
                            adjusted["n"],

                        "beta_baseline":
                            base["beta_ASTS"],

                        "p_baseline":
                            base["p_ASTS"],

                        "beta_driver_adjusted":
                            adjusted[
                                "beta_ASTS"
                            ],

                        "p_driver_adjusted":
                            adjusted[
                                "p_ASTS"
                            ],

                        "drivers_used":
                            adjusted[
                                "drivers_used"
                            ],

                        "n_MAPK_driver":
                            adjusted[
                                "n_MAPK_driver"
                            ],

                        "n_ERBB_driver":
                            adjusted[
                                "n_ERBB_driver"
                            ],

                        "n_PI3K_driver":
                            adjusted[
                                "n_PI3K_driver"
                            ],
                    }
                )


        # ----------------------------------------------------
        # GDSC2
        # ----------------------------------------------------

        g = gdsc[
            gdsc[
                "DRUG_ID_norm"
            ] == gdsc_id
        ].copy()


        if len(g) >= MIN_N_PER_DRUG:

            base = fit_response_model(
                g,
                include_drivers=False,
            )

            adjusted = fit_response_model(
                g,
                include_drivers=True,
            )


            if (
                base is not None
                and adjusted is not None
            ):

                result_rows.append(
                    {
                        "screen":
                            "GDSC2",

                        "compound_name":
                            drug_name,

                        "drug_key":
                            drug_key,

                        "DRUG_ID_norm":
                            gdsc_id,

                        "n":
                            adjusted["n"],

                        "beta_baseline":
                            base["beta_ASTS"],

                        "p_baseline":
                            base["p_ASTS"],

                        "beta_driver_adjusted":
                            adjusted[
                                "beta_ASTS"
                            ],

                        "p_driver_adjusted":
                            adjusted[
                                "p_ASTS"
                            ],

                        "drivers_used":
                            adjusted[
                                "drivers_used"
                            ],

                        "n_MAPK_driver":
                            adjusted[
                                "n_MAPK_driver"
                            ],

                        "n_ERBB_driver":
                            adjusted[
                                "n_ERBB_driver"
                            ],

                        "n_PI3K_driver":
                            adjusted[
                                "n_PI3K_driver"
                            ],
                    }
                )


    results = pd.DataFrame(
        result_rows
    )


    # --------------------------------------------------------
    # attenuation
    # --------------------------------------------------------

    results[
        "same_sign_after_adjustment"
    ] = (
        np.sign(
            results[
                "beta_baseline"
            ]
        )
        ==
        np.sign(
            results[
                "beta_driver_adjusted"
            ]
        )
    )


    results[
        "beta_abs_retention"
    ] = (
        results[
            "beta_driver_adjusted"
        ].abs()
        /
        results[
            "beta_baseline"
        ].abs()
    )


    results[
        "attenuation_percent"
    ] = (
        1
        -
        results[
            "beta_abs_retention"
        ]
    ) * 100.0


    results.to_csv(
        OUT_DIR
        / "13_all_drugs_driver_adjustment.tsv",
        sep="\t",
        index=False,
    )


    # ========================================================
    # 9. Cross-screen adjusted effects
    # ========================================================

    prism_result = (
        results[
            results[
                "screen"
            ] == "PRISM"
        ]
        .copy()
        .rename(
            columns={
                "beta_baseline":
                    "PRISM_beta_baseline",

                "beta_driver_adjusted":
                    "PRISM_beta_adjusted",

                "p_driver_adjusted":
                    "PRISM_p_adjusted",

                "attenuation_percent":
                    "PRISM_attenuation_percent",
            }
        )
    )


    gdsc_result = (
        results[
            results[
                "screen"
            ] == "GDSC2"
        ]
        .copy()
        .rename(
            columns={
                "beta_baseline":
                    "GDSC2_beta_baseline",

                "beta_driver_adjusted":
                    "GDSC2_beta_adjusted",

                "p_driver_adjusted":
                    "GDSC2_p_adjusted",

                "attenuation_percent":
                    "GDSC2_attenuation_percent",
            }
        )
    )


    compare = prism_result[
        [
            "compound_name",
            "PRISM_beta_baseline",
            "PRISM_beta_adjusted",
            "PRISM_p_adjusted",
            "PRISM_attenuation_percent",
        ]
    ].merge(
        gdsc_result[
            [
                "compound_name",
                "GDSC2_beta_baseline",
                "GDSC2_beta_adjusted",
                "GDSC2_p_adjusted",
                "GDSC2_attenuation_percent",
            ]
        ],
        on="compound_name",
        how="inner",
        validate="1:1",
    )


    # --------------------------------------------------------
    # Adjusted cross-screen z scores
    # --------------------------------------------------------

    compare[
        "PRISM_adjusted_z"
    ] = zscore(
        compare[
            "PRISM_beta_adjusted"
        ]
    )


    compare[
        "GDSC2_adjusted_z"
    ] = zscore(
        compare[
            "GDSC2_beta_adjusted"
        ]
    )


    compare[
        "adjusted_consensus_z"
    ] = (
        compare[
            "PRISM_adjusted_z"
        ]
        +
        compare[
            "GDSC2_adjusted_z"
        ]
    ) / 2.0


    compare[
        "adjusted_direction_concordant"
    ] = (
        np.sign(
            compare[
                "PRISM_beta_adjusted"
            ]
        )
        ==
        np.sign(
            compare[
                "GDSC2_beta_adjusted"
            ]
        )
    )


    # Add class membership
    class_cols = [
        col
        for col in [
            "MEK_target",
            "EGFR_target",
            "PI3K_MTOR_pathway",
            "EGFR_pathway",
            "ERK_MAPK_pathway",
            "RTK_pathway",
            "RTK_non_EGFR",
        ]
        if col in classes.columns
    ]


    compare = compare.merge(
        classes[
            [
                "compound_name",
                *class_cols,
            ]
        ],
        on="compound_name",
        how="left",
        validate="1:1",
    )


    compare.to_csv(
        OUT_DIR
        / "13_PRISM_GDSC2_driver_adjusted.tsv",
        sep="\t",
        index=False,
    )


    # ========================================================
    # 10. Class summary after adjustment
    # ========================================================

    class_rows = []


    for class_name in class_cols:

        sub = compare[
            compare[
                class_name
            ] == 1
        ].copy()


        if len(sub) == 0:
            continue


        resistant = (
            (
                sub[
                    "PRISM_beta_adjusted"
                ] > 0
            )
            &
            (
                sub[
                    "GDSC2_beta_adjusted"
                ] > 0
            )
        )


        sensitive = (
            (
                sub[
                    "PRISM_beta_adjusted"
                ] < 0
            )
            &
            (
                sub[
                    "GDSC2_beta_adjusted"
                ] < 0
            )
        )


        class_rows.append(
            {
                "class":
                    class_name,

                "n_drugs":
                    len(sub),

                "n_adjusted_concordant_resistant":
                    int(
                        resistant.sum()
                    ),

                "n_adjusted_concordant_sensitive":
                    int(
                        sensitive.sum()
                    ),

                "mean_adjusted_consensus_z":
                    sub[
                        "adjusted_consensus_z"
                    ].mean(),

                "min_adjusted_consensus_z":
                    sub[
                        "adjusted_consensus_z"
                    ].min(),

                "mean_PRISM_attenuation_percent":
                    sub[
                        "PRISM_attenuation_percent"
                    ].mean(),

                "mean_GDSC2_attenuation_percent":
                    sub[
                        "GDSC2_attenuation_percent"
                    ].mean(),
            }
        )


    class_summary = pd.DataFrame(
        class_rows
    )


    class_summary.to_csv(
        OUT_DIR
        / "13_class_summary_after_driver_adjustment.tsv",
        sep="\t",
        index=False,
    )


    # ========================================================
    # 11. Global cross-screen concordance
    # ========================================================

    adjusted_rho = np.nan
    adjusted_p = np.nan


    if len(compare) >= 3:

        (
            adjusted_rho,
            adjusted_p,
        ) = spearmanr(
            compare[
                "PRISM_beta_adjusted"
            ],
            compare[
                "GDSC2_beta_adjusted"
            ],
        )


    adjusted_sign_fraction = (
        compare[
            "adjusted_direction_concordant"
        ].mean()
    )


    # ========================================================
    # 12. Figures
    # ========================================================

    # --------------------------------------------------------
    # Figure 1:
    # baseline vs driver-adjusted beta
    # --------------------------------------------------------

    for screen in [
        "PRISM",
        "GDSC2",
    ]:

        plot = results[
            results[
                "screen"
            ] == screen
        ].copy()


        fig, ax = plt.subplots(
            figsize=(6, 6)
        )


        ax.scatter(
            plot[
                "beta_baseline"
            ],
            plot[
                "beta_driver_adjusted"
            ],
            s=24,
            alpha=0.7,
        )


        low = min(
            plot[
                "beta_baseline"
            ].min(),
            plot[
                "beta_driver_adjusted"
            ].min(),
        )


        high = max(
            plot[
                "beta_baseline"
            ].max(),
            plot[
                "beta_driver_adjusted"
            ].max(),
        )


        ax.plot(
            [low, high],
            [low, high],
            linestyle="--",
        )


        ax.axhline(
            0,
            linestyle=":",
        )


        ax.axvline(
            0,
            linestyle=":",
        )


        ax.set_xlabel(
            "Baseline ASTS beta"
        )

        ax.set_ylabel(
            "Driver-adjusted ASTS beta"
        )


        ax.set_title(
            f"{screen}: driver adjustment"
        )


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / f"13_01_{screen}_baseline_vs_adjusted.png",
            dpi=200,
        )


        plt.close(fig)


    # --------------------------------------------------------
    # Figure 2:
    # adjusted MEK/EGFR/PI3K drugs
    # --------------------------------------------------------

    selected_cols = [
        col
        for col in [
            "MEK_target",
            "EGFR_target",
            "PI3K_MTOR_pathway",
        ]
        if col in compare.columns
    ]


    if len(selected_cols) > 0:

        selected = compare[
            compare[
                selected_cols
            ].sum(
                axis=1
            ) > 0
        ].copy()


        selected = selected.sort_values(
            "adjusted_consensus_z"
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
                "adjusted_consensus_z"
            ],
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
            linestyle="--",
        )


        ax.set_xlabel(
            "Driver-adjusted cross-screen consensus z-score\n"
            "positive = ASTS-high resistance"
        )


        ax.set_title(
            "MEK / EGFR / PI3K-MTOR after driver adjustment"
        )


        fig.tight_layout()


        fig.savefig(
            FIG_DIR
            / "13_02_selected_classes_after_driver_adjustment.png",
            dpi=200,
        )


        plt.close(fig)


    # ========================================================
    # 13. Summary
    # ========================================================

    summary_file = (
        OUT_DIR
        / "13_summary.txt"
    )


    with open(
        summary_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 13 summary",
            file=f,
        )

        print(
            "==========================",
            file=f,
        )

        print(
            "",
            file=f,
        )


        print(
            "Purpose:",
            file=f,
        )

        print(
            "Test whether ASTS-drug associations are explained "
            "by canonical oncogenic driver hotspot mutations.",
            file=f,
        )

        print(
            "",
            file=f,
        )


        print(
            f"ASTS models: {len(asts)}",
            file=f,
        )

        print(
            f"Mutation-covered ASTS models: "
            f"{len(asts_driver)}",
            file=f,
        )

        print(
            "",
            file=f,
        )


        print(
            "Driver prevalence",
            file=f,
        )

        print(
            "-----------------",
            file=f,
        )


        for row in prevalence_rows:

            print(
                f"{row['driver']}\t"
                f"{row['n_positive']}/"
                f"{row['n_total']}\t"
                f"{row['fraction']:.4f}",
                file=f,
            )


        print(
            "",
            file=f,
        )


        print(
            "ASTS ~ driver genotype",
            file=f,
        )

        print(
            "----------------------",
            file=f,
        )


        print(
            f"N: {asts_driver_stats['n']}",
            file=f,
        )

        print(
            f"Base R2: "
            f"{asts_driver_stats['base_R2']}",
            file=f,
        )

        print(
            f"Driver-adjusted R2: "
            f"{asts_driver_stats['driver_R2']}",
            file=f,
        )

        print(
            f"Delta R2 from drivers: "
            f"{asts_driver_stats['delta_R2_drivers']}",
            file=f,
        )

        print(
            "",
            file=f,
        )


        print(
            asts_driver_coef.to_string(
                index=False
            ),
            file=f,
        )


        print(
            "",
            file=f,
        )


        print(
            "Drug-response models",
            file=f,
        )

        print(
            "--------------------",
            file=f,
        )

        print(
            f"PRISM drugs tested: "
            f"{(results['screen'] == 'PRISM').sum()}",
            file=f,
        )

        print(
            f"GDSC2 drugs tested: "
            f"{(results['screen'] == 'GDSC2').sum()}",
            file=f,
        )

        print(
            f"Cross-screen adjusted overlap: "
            f"{len(compare)}",
            file=f,
        )

        print(
            "",
            file=f,
        )


        print(
            "Adjusted cross-screen concordance",
            file=f,
        )

        print(
            "---------------------------------",
            file=f,
        )

        print(
            f"Spearman r: {adjusted_rho}",
            file=f,
        )

        print(
            f"Spearman p: {adjusted_p}",
            file=f,
        )

        print(
            f"Sign concordance: "
            f"{adjusted_sign_fraction}",
            file=f,
        )


        print(
            "",
            file=f,
        )


        print(
            "Class summary after driver adjustment",
            file=f,
        )

        print(
            "-------------------------------------",
            file=f,
        )


        if len(class_summary) > 0:

            print(
                class_summary.to_string(
                    index=False
                ),
                file=f,
            )


        print(
            "",
            file=f,
        )


        print(
            "Interpretation guide",
            file=f,
        )

        print(
            "--------------------",
            file=f,
        )

        print(
            "If MEK/EGFR adjusted effects remain positive "
            "in both PRISM and GDSC2, the ASTS association "
            "is not fully explained by canonical hotspot-driver "
            "mutation status.",
            file=f,
        )

        print(
            "",
            file=f,
        )

        print(
            "This does NOT prove causality.",
            file=f,
        )

        print(
            "Copy-number, fusion, epigenetic and other "
            "genomic confounders are not addressed here.",
            file=f,
        )

        print(
            "",
            file=f,
        )


        print(
            "Checksums",
            file=f,
        )

        print(
            "---------",
            file=f,
        )

        print(
            "ASTS:",
            sha256_file(
                ASTS_FILE
            ),
            file=f,
        )

        print(
            "Hotspot mutations:",
            sha256_file(
                HOTSPOT_FILE
            ),
            file=f,
        )

        print(
            "PRISM:",
            sha256_file(
                PRISM_FILE
            ),
            file=f,
        )

        print(
            "GDSC2:",
            sha256_file(
                GDSC2_FILE
            ),
            file=f,
        )


    print(
        "\n============================================================"
    )

    print(
        "Step 13 completed successfully"
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
