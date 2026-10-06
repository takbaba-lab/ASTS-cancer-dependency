#!/usr/bin/env python3

"""
ASTRA-Drug
Step 14A: ASTS-associated transcriptome

----

Primary model:
    Gene expression
        ~ ASTS
        + Lineage
        + Sex
        + Proliferation

Driver-adjusted sensitivity model:
    Gene expression
        ~ ASTS
        + Lineage
        + Sex
        + Proliferation
        + MAPK_driver
        + ERBB_driver
        + PI3K_driver

----

    is_ASTS_signature




----------

------
"""


# ============================================================
# Path / Parameter
# ============================================================

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]


# DepMap expression
EXPRESSION_FILE = (
    PROJECT_DIR
    / "data/raw/depmap/26Q1"
    / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
)


# Step08b ASTS
ASTS_FILE = (
    PROJECT_DIR
    / "results/08b_DepMap_ASTS_default_profile"
    / "08_DepMap_ASTS_scores.tsv"
)


# Step08b signature mapping
SIGNATURE_MAPPING_FILE = (
    PROJECT_DIR
    / "results/08b_DepMap_ASTS_default_profile"
    / "08_signature_expression_mapping.tsv"
)


# Step13 driver flags
DRIVER_FILE = (
    PROJECT_DIR
    / "results/13_driver_confounder"
    / "13_driver_flags.tsv"
)


OUT_DIR = (
    PROJECT_DIR
    / "results/14_ASTS_transcriptome"
)


# ------------------------------------------------------------
# Parameters
# ------------------------------------------------------------

DEFAULT_PROFILE_COLUMN = "IsDefaultEntryForModel"
DEFAULT_PROFILE_VALUE = "yes"

LINEAGE_COLUMN = "OncotreeLineage"
SEX_COLUMN = "Sex"

ASTS_COLUMN = "ASTS_score"
PROLIFERATION_COLUMN = "proliferation_proxy"

MIN_GLOBAL_LINEAGE_N = 5

GENE_CHUNK_SIZE = 1000

MIN_MODELS = 500

GENE_COLUMN_PATTERN = r"^\S[^()]*\s+\(\d+\)$"

MAX_INVALID_VALUE_EXAMPLES = 5


# ============================================================
# Imports
# ============================================================

import hashlib
import re

import numpy as np
import pandas as pd

from scipy.stats import t as t_distribution
from scipy.stats import spearmanr

from statsmodels.stats.multitest import multipletests


# ============================================================
# Helper functions
# ============================================================

def select_gene_columns(expression):
    """Return columns matching the gene-symbol (numeric-ID) naming format."""

    gene_columns = [
        col for col in expression.columns
        if re.fullmatch(GENE_COLUMN_PATTERN, str(col).strip())
    ]
    selected = set(gene_columns)
    excluded = [col for col in expression.columns if col not in selected]

    print("Excluded non-gene columns:", excluded)
    print("Detected gene expression columns:", len(gene_columns))

    if not gene_columns:
        raise ValueError(
        )

    return gene_columns


def expression_to_numeric(expression, columns):
    """Convert expression values to numeric arrays and stop with an informative error if invalid values are detected."""

    block = expression[columns]
    try:
        values = block.to_numpy(dtype=float)
    except (ValueError, TypeError) as error:
        problems = []
        for col in columns:
            numeric = pd.to_numeric(block[col], errors="coerce")
            invalid = block[col].notna() & numeric.isna()
            if invalid.any():
                examples = (
                    block.loc[invalid, col].astype(str).unique()
                    [:MAX_INVALID_VALUE_EXAMPLES].tolist()
                )
                problems.append(f"{col}: {examples}")
        raise ValueError(
            + ("\n".join(problems) or str(error))
        ) from error

    finite = np.isfinite(values)
    if not finite.all():
        problems = []
        for index, col in enumerate(columns):
            invalid = ~finite[:, index]
            if invalid.any():
                model_examples = (
                    block.index[invalid].astype(str).tolist()
                    [:MAX_INVALID_VALUE_EXAMPLES]
                )
                problems.append(
                )
        raise ValueError(
        )

    return values


def sha256_file(path, block_size=1024 * 1024):

    h = hashlib.sha256()

    with open(path, "rb") as f:

        while True:

            block = f.read(block_size)

            if not block:
                break

            h.update(block)

    return h.hexdigest()


def zscore(series):

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


def extract_gene_symbol(column_name):
    """

        KRAS (3845)


    """

    x = str(column_name)

    match = re.match(
        r"^(.+?)\s+\(\d+\)$",
        x
    )

    if match:
        return match.group(1).strip()

    return x.strip()


def detect_signature_gene_column(df):
    """
    """

    preferred = [
        "human_gene_symbol",
        "human_symbol",
        "gene_symbol",
        "GeneSymbol",
        "symbol",
        "human_gene",
        "Human_Symbol",
    ]

    for col in preferred:

        if col in df.columns:
            return col


    candidates = [
        col
        for col in df.columns
        if (
            "human" in col.lower()
            and
            (
                "symbol" in col.lower()
                or
                "gene" in col.lower()
            )
        )
    ]


    if len(candidates) == 1:
        return candidates[0]


    raise ValueError(
        "Could not uniquely detect human gene symbol column "
        f"in signature mapping.\nColumns={df.columns.tolist()}"
    )


def build_covariate_matrix(
    metadata,
    include_drivers=False,
):
    """

    """

    work = metadata.copy()


    parts = []


    # --------------------------------------------------------
    # Intercept
    # --------------------------------------------------------

    parts.append(
        pd.DataFrame(
            {
                "Intercept":
                    np.ones(
                        len(work)
                    )
            },
            index=work.index
        )
    )


    # --------------------------------------------------------
    # Proliferation
    # --------------------------------------------------------

    parts.append(
        pd.DataFrame(
            {
                "Proliferation_z":
                    zscore(
                        work[
                            PROLIFERATION_COLUMN
                        ]
                    )
            },
            index=work.index
        )
    )


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


    # --------------------------------------------------------
    # Driver flags
    # --------------------------------------------------------

    if include_drivers:

        for col in [
            "MAPK_driver",
            "ERBB_driver",
            "PI3K_driver",
        ]:

            if col not in work.columns:
                raise ValueError(
                    f"Missing driver column: {col}"
                )


            parts.append(
                pd.DataFrame(
                    {
                        col:
                            pd.to_numeric(
                                work[col],
                                errors="coerce"
                            )
                    },
                    index=work.index
                )
            )


    X = pd.concat(
        parts,
        axis=1
    )


    return X.astype(float)


def fit_gene_associations(
    expression,
    metadata,
    gene_columns,
    gene_symbols,
    signature_genes,
    include_drivers=False,
):
    """


    """

    # --------------------------------------------------------
    # Complete-case models
    # --------------------------------------------------------

    required = [
        ASTS_COLUMN,
        PROLIFERATION_COLUMN,
        "LineageModel",
    ]


    if include_drivers:

        required.extend(
            [
                "MAPK_driver",
                "ERBB_driver",
                "PI3K_driver",
            ]
        )


    valid = metadata[
        required
    ].notna().all(
        axis=1
    )


    meta = metadata.loc[
        valid
    ].copy()


    expr = expression.loc[
        valid
    ].copy()


    if len(meta) < MIN_MODELS:

        raise RuntimeError(
            f"Too few models: {len(meta)}"
        )


    # --------------------------------------------------------
    # ASTS standardization
    # --------------------------------------------------------

    asts_z = zscore(
        meta[
            ASTS_COLUMN
        ]
    ).to_numpy(
        dtype=float
    )


    # --------------------------------------------------------
    # Covariate matrix
    # --------------------------------------------------------

    C = build_covariate_matrix(
        meta,
        include_drivers=include_drivers,
    ).to_numpy(
        dtype=float
    )


    # --------------------------------------------------------
    # QR decomposition
    # --------------------------------------------------------

    Q, _ = np.linalg.qr(
        C,
        mode="reduced"
    )


    # ASTS residual after covariate removal
    asts_res = (
        asts_z
        -
        Q @ (
            Q.T @ asts_z
        )
    )


    asts_ss = np.sum(
        asts_res ** 2
    )


    if asts_ss <= 0:

        raise RuntimeError(
            "Residual ASTS variance is zero."
        )


    # full model df
    df_resid = (
        len(meta)
        -
        C.shape[1]
        -
        1
    )


    if df_resid <= 0:

        raise RuntimeError(
            "Residual degrees of freedom <= 0."
        )


    # --------------------------------------------------------
    # Gene chunks
    # --------------------------------------------------------

    rows = []


    for start in range(
        0,
        len(gene_columns),
        GENE_CHUNK_SIZE
    ):

        end = min(
            start + GENE_CHUNK_SIZE,
            len(gene_columns)
        )


        cols = gene_columns[
            start:end
        ]


        symbols = gene_symbols[
            start:end
        ]


        print(
            f"Genes {start + 1}-{end} / "
            f"{len(gene_columns)}"
        )


        Y = expression_to_numeric(expr, cols)


        # ----------------------------------------------------
        # Remove covariate effects from expression
        # ----------------------------------------------------

        Y_res = (
            Y
            -
            Q @ (
                Q.T @ Y
            )
        )


        # ----------------------------------------------------
        # ASTS slope
        # ----------------------------------------------------

        beta = (
            asts_res @ Y_res
        ) / asts_ss


        # ----------------------------------------------------
        # Full residual
        # ----------------------------------------------------

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


        sse = np.sum(
            residual ** 2,
            axis=0
        )


        sigma2 = (
            sse
            /
            df_resid
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
            2.0
            *
            t_distribution.sf(
                np.abs(
                    t_stat
                ),
                df=df_resid
            )
        )


        # partial correlation
        partial_r = (
            t_stat
            /
            np.sqrt(
                t_stat ** 2
                +
                df_resid
            )
        )


        # expression SD
        y_sd = np.nanstd(
            Y,
            axis=0,
            ddof=1
        )


        with np.errstate(
            divide="ignore",
            invalid="ignore"
        ):

            standardized_beta = (
                beta
                /
                y_sd
            )


        for i in range(
            len(cols)
        ):

            rows.append(
                {
                    "gene_symbol":
                        symbols[i],

                    "expression_column":
                        cols[i],

                    "n_models":
                        len(meta),

                    "beta_log2TPM_per_1SD_ASTS":
                        beta[i],

                    "standardized_beta":
                        standardized_beta[i],

                    "t_stat":
                        t_stat[i],

                    "p_value":
                        p_value[i],

                    "partial_r":
                        partial_r[i],

                    "is_ASTS_signature":
                        symbols[i]
                        in signature_genes,
                }
            )


    result = pd.DataFrame(
        rows
    )


    # --------------------------------------------------------
    # FDR
    # --------------------------------------------------------

    finite = np.isfinite(
        result[
            "p_value"
        ]
    )


    result[
        "FDR"
    ] = np.nan


    result.loc[
        finite,
        "FDR"
    ] = multipletests(
        result.loc[
            finite,
            "p_value"
        ],
        method="fdr_bh"
    )[1]


    result = result.sort_values(
        "t_stat",
        ascending=False
    )


    return result


# ============================================================
# Main
# ============================================================

def main():

    print(
        "\n"
        "============================================================\n"
        "ASTRA-Drug Step 14A\n"
        "ASTS-associated transcriptome\n"
        "============================================================"
    )


    # --------------------------------------------------------
    # Input checks
    # --------------------------------------------------------

    for path in [
        EXPRESSION_FILE,
        ASTS_FILE,
        SIGNATURE_MAPPING_FILE,
        DRIVER_FILE,
    ]:

        if not path.exists():

            raise FileNotFoundError(
                f"Required file not found:\n{path}"
            )


    OUT_DIR.mkdir(
        parents=True,
        exist_ok=True
    )


    # ========================================================
    # 1. Signature genes
    # ========================================================

    print(
        "\nStep 1: Reading ASTS signature"
    )


    sig = pd.read_csv(
        SIGNATURE_MAPPING_FILE,
        sep="\t"
    )


    signature_col = (
        detect_signature_gene_column(
            sig
        )
    )


    signature_genes = set(
        sig[
            signature_col
        ]
        .dropna()
        .astype(str)
        .str.strip()
    )


    print(
        "Signature column:",
        signature_col
    )

    print(
        "Signature genes:",
        len(signature_genes)
    )


    if len(signature_genes) < 70:

        raise RuntimeError(
            "Unexpectedly few signature genes. "
            "Check signature mapping."
        )


    # ========================================================
    # 2. ASTS metadata
    # ========================================================

    print(
        "\nStep 2: Reading ASTS metadata"
    )


    asts = pd.read_csv(
        ASTS_FILE,
        sep="\t"
    )


    if asts[
        "ModelID"
    ].duplicated().any():

        raise RuntimeError(
            "Duplicated ModelID in ASTS file."
        )


    # --------------------------------------------------------
    # Lineage categories
    # --------------------------------------------------------

    lineage_counts = (
        asts[
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


    asts[
        "LineageModel"
    ] = (
        asts[
            LINEAGE_COLUMN
        ]
        .fillna("Unknown")
        .astype(str)
        .apply(
            lambda x:
                x
                if x in common_lineages
                else "OtherRare"
        )
    )


    # ========================================================
    # 3. Driver flags
    # ========================================================

    print(
        "\nStep 3: Reading driver flags"
    )


    driver = pd.read_csv(
        DRIVER_FILE,
        sep="\t"
    )


    driver_keep = [
        "ModelID",
        "MAPK_driver",
        "ERBB_driver",
        "PI3K_driver",
    ]


    driver = driver[
        driver_keep
    ].copy()


    if driver[
        "ModelID"
    ].duplicated().any():

        raise RuntimeError(
            "Duplicated ModelID in driver file."
        )


    metadata = asts.merge(
        driver,
        on="ModelID",
        how="left",
        validate="1:1"
    )


    # ========================================================
    # 4. Expression
    # ========================================================

    print(
        "\nStep 4: Reading DepMap expression"
    )


    expression = pd.read_csv(
        EXPRESSION_FILE,
        low_memory=False
    )


    if DEFAULT_PROFILE_COLUMN not in expression.columns:

        raise ValueError(
            f"{DEFAULT_PROFILE_COLUMN} not found."
        )


    default_flag = (
        expression[
            DEFAULT_PROFILE_COLUMN
        ]
        .astype(str)
        .str.strip()
        .str.lower()
    )


    expression = expression[
        default_flag
        ==
        DEFAULT_PROFILE_VALUE
    ].copy()


    if (
        "ModelID"
        not in expression.columns
    ):

        raise ValueError(
            "ModelID column not found."
        )


    if expression[
        "ModelID"
    ].duplicated().any():

        duplicates = expression.loc[
            expression[
                "ModelID"
            ].duplicated(
                False
            ),
            "ModelID"
        ].unique()

        raise RuntimeError(
            "Duplicated default-profile ModelID remains: "
            + ",".join(
                duplicates[:20]
            )
        )


    print(
        "Default expression models:",
        len(expression)
    )


    # ========================================================
    # 5. Gene columns
    # ========================================================

    gene_columns_raw = select_gene_columns(expression)


    gene_symbols_raw = [
        extract_gene_symbol(
            col
        )
        for col in gene_columns_raw
    ]


    # --------------------------------------------------------
    # Avoid duplicate gene symbols
    # --------------------------------------------------------

    seen = set()

    gene_columns = []

    gene_symbols = []


    duplicated_symbols = []


    for col, symbol in zip(
        gene_columns_raw,
        gene_symbols_raw
    ):

        if symbol in seen:

            duplicated_symbols.append(
                symbol
            )

            continue


        seen.add(
            symbol
        )

        gene_columns.append(
            col
        )

        gene_symbols.append(
            symbol
        )


    print(
        "Gene columns:",
        len(gene_columns)
    )

    print(
        "Duplicate symbols removed:",
        len(
            set(
                duplicated_symbols
            )
        )
    )


    # ========================================================
    # 6. Align expression and metadata
    # ========================================================

    print(
        "\nStep 5: Aligning models"
    )


    expression = expression.set_index(
        "ModelID"
    )


    metadata = metadata.set_index(
        "ModelID"
    )


    common_models = (
        expression.index.intersection(
            metadata.index
        )
    )


    expression = expression.loc[
        common_models
    ]


    metadata = metadata.loc[
        common_models
    ]


    print(
        "Models in primary analysis:",
        len(common_models)
    )


    # ========================================================
    # 7. Primary transcriptome association
    # ========================================================

    print(
        "\nStep 6: Primary transcriptome association"
    )


    primary = fit_gene_associations(
        expression=expression,
        metadata=metadata,
        gene_columns=gene_columns,
        gene_symbols=gene_symbols,
        signature_genes=signature_genes,
        include_drivers=False,
    )


    primary.to_csv(
        OUT_DIR
        / "14_gene_association_primary.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # 8. Mutation-covered baseline
    # ========================================================

    print(
        "\nStep 7: Mutation-subset baseline"
    )


    mutation_mask = metadata[
        [
            "MAPK_driver",
            "ERBB_driver",
            "PI3K_driver",
        ]
    ].notna().all(
        axis=1
    )


    mutation_expression = expression.loc[
        mutation_mask
    ]


    mutation_metadata = metadata.loc[
        mutation_mask
    ]


    baseline_mut = fit_gene_associations(
        expression=mutation_expression,
        metadata=mutation_metadata,
        gene_columns=gene_columns,
        gene_symbols=gene_symbols,
        signature_genes=signature_genes,
        include_drivers=False,
    )


    baseline_mut.to_csv(
        OUT_DIR
        / "14_gene_association_mutation_subset_baseline.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # 9. Driver-adjusted transcriptome association
    # ========================================================

    print(
        "\nStep 8: Driver-adjusted transcriptome association"
    )


    driver_adjusted = fit_gene_associations(
        expression=mutation_expression,
        metadata=mutation_metadata,
        gene_columns=gene_columns,
        gene_symbols=gene_symbols,
        signature_genes=signature_genes,
        include_drivers=True,
    )


    driver_adjusted.to_csv(
        OUT_DIR
        / "14_gene_association_driver_adjusted.tsv",
        sep="\t",
        index=False
    )


    # ========================================================
    # 10. Compare baseline vs driver adjusted
    # ========================================================

    compare = baseline_mut[
        [
            "gene_symbol",
            "t_stat",
        ]
    ].merge(
        driver_adjusted[
            [
                "gene_symbol",
                "t_stat",
            ]
        ],
        on="gene_symbol",
        suffixes=(
            "_baseline",
            "_driver_adjusted"
        ),
        validate="1:1"
    )


    rho, rho_p = spearmanr(
        compare[
            "t_stat_baseline"
        ],
        compare[
            "t_stat_driver_adjusted"
        ]
    )


    sign_concordance = np.mean(
        np.sign(
            compare[
                "t_stat_baseline"
            ]
        )
        ==
        np.sign(
            compare[
                "t_stat_driver_adjusted"
            ]
        )
    )


    # ========================================================
    # 11. Summary
    # ========================================================

    summary_file = (
        OUT_DIR
        / "14a_summary.txt"
    )


    with open(
        summary_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 14A summary",
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
            f"Default-profile expression models: "
            f"{len(expression)}",
            file=f
        )

        print(
            f"Genes analyzed: "
            f"{len(gene_columns)}",
            file=f
        )

        print(
            f"ASTS signature genes identified: "
            f"{len(signature_genes)}",
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
            "Expression ~ ASTS + lineage + sex + proliferation",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "Primary significant genes",
            file=f
        )

        print(
            "-------------------------",
            file=f
        )

        print(
            f"FDR < 0.05: "
            f"{(primary['FDR'] < 0.05).sum()}",
            file=f
        )

        print(
            f"Positive: "
            f"{((primary['FDR'] < 0.05) & (primary['t_stat'] > 0)).sum()}",
            file=f
        )

        print(
            f"Negative: "
            f"{((primary['FDR'] < 0.05) & (primary['t_stat'] < 0)).sum()}",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "Top positive non-signature genes",
            file=f
        )

        print(
            "--------------------------------",
            file=f
        )


        top_positive = primary[
            ~primary[
                "is_ASTS_signature"
            ]
        ].head(
            20
        )


        print(
            top_positive[
                [
                    "gene_symbol",
                    "t_stat",
                    "FDR",
                    "partial_r",
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
            "Top negative non-signature genes",
            file=f
        )

        print(
            "--------------------------------",
            file=f
        )


        top_negative = (
            primary[
                ~primary[
                    "is_ASTS_signature"
                ]
            ]
            .sort_values(
                "t_stat"
            )
            .head(
                20
            )
        )


        print(
            top_negative[
                [
                    "gene_symbol",
                    "t_stat",
                    "FDR",
                    "partial_r",
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
            "Mutation-subset baseline vs driver-adjusted",
            file=f
        )

        print(
            "------------------------------------------",
            file=f
        )

        print(
            f"Spearman t-stat r: {rho}",
            file=f
        )

        print(
            f"Spearman p: {rho_p}",
            file=f
        )

        print(
            f"Sign concordance: {sign_concordance}",
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
            "Step14B primary GSEA must exclude "
            "is_ASTS_signature == True genes.",
            file=f
        )

        print(
            "This avoids direct circular enrichment from "
            "the 78 genes used to calculate ASTS.",
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
            "Expression:",
            sha256_file(
                EXPRESSION_FILE
            ),
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
            "Drivers:",
            sha256_file(
                DRIVER_FILE
            ),
            file=f
        )


    print(
        "\n============================================================"
    )

    print(
        "Step 14A completed successfully"
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
