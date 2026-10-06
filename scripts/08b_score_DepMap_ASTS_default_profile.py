#!/usr/bin/env python3

"""
ASTRA-Drug
Step 08: Score human ASTS v1.0 in DepMap expression data

----

Primary ASTS:
    DHT-relative-up   40 genes
    DHT-relative-down 38 genes
    Total             78 genes

ASTS score:

        ASTS = mean(rank of DHT-relative-up genes)
             - mean(rank of DHT-relative-down genes)



----

--------
"""


# ============================================================
# Path / Parameter
# ============================================================

from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parents[1]

DEPMAP_RELEASE = "26Q1"

SIGNATURE_FILE = (
    PROJECT_DIR
    / "results/07_ASTS_human_v1"
    / "ASTS_human_v1.0_primary.tsv"
)



EXPRESSION_FILE = (
    PROJECT_DIR
    / "data/raw/depmap"
    / "26Q1"
    / "OmicsExpressionTPMLogp1HumanProteinCodingGenes.csv"
)

MODEL_FILE = (
    PROJECT_DIR
    / "data/raw/depmap"
    / "26Q1"
    / "Model.csv"
)


OUT_DIR = (
    PROJECT_DIR
    / "results/08b_DepMap_ASTS_default_profile"
)

FIG_DIR = OUT_DIR / "figures"


CHUNK_SIZE = 100


MIN_DIRECTION_COVERAGE = 0.90


MIN_LINEAGE_N = 10


HIGH_CONFIDENCE_VALUE = 1


#
PROLIFERATION_PROXY_GENES = [
    "MKI67",
    "PCNA",
    "MCM2",
    "MCM3",
    "MCM4",
    "MCM5",
    "MCM6",
    "MCM7",
    "TOP2A",
    "CDK1",
]


DEFAULT_PROFILE_COLUMN = "IsDefaultEntryForModel"
DEFAULT_PROFILE_VALUE = "yes"


# ============================================================
# Import
# ============================================================

import re
import sys
import hashlib

import numpy as np
import pandas as pd

from scipy.stats import (
    spearmanr,
    mannwhitneyu,
    kruskal,
)

import matplotlib.pyplot as plt


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


def normalize_symbol(text):
    """

        "TP53"
        "TP53 (7157)"

    """

    text = str(text).strip()

    # "TP53 (7157)" -> "TP53"
    text = re.sub(
        r"\s+\([^()]+\)$",
        "",
        text
    )

    return text.upper()


def extract_ensembl(text):
    """
    """

    m = re.search(
        r"(ENSG\d+)",
        str(text)
    )

    if m is None:
        return None

    return m.group(1)


def detect_model_id_column(columns):
    """
    """

    candidates = [
        "ModelID",
        "model_id",
        "DepMap_ID",
        "DepMapID",
    ]

    for x in candidates:
        if x in columns:
            return x

    print(
        "WARNING: known ModelID column was not found. "
        "Using the first column:",
        columns[0]
    )

    return columns[0]


def map_signature_to_expression(
    signature,
    expression_columns
):
    """

    """

    symbol_map = {}
    ensembl_map = {}

    for col in expression_columns:

        symbol = normalize_symbol(col)

        symbol_map.setdefault(
            symbol,
            []
        ).append(col)

        ens = extract_ensembl(col)

        if ens is not None:
            ensembl_map.setdefault(
                ens,
                []
            ).append(col)


    mapping_rows = []


    for _, row in signature.iterrows():

        symbol = str(
            row["human_symbol"]
        ).upper()

        ensembl = str(
            row["human_ensembl_id"]
        )


        candidates = []


        if ensembl in ensembl_map:
            candidates.extend(
                ensembl_map[ensembl]
            )


        # symbol fallback
        if len(candidates) == 0:

            if symbol in symbol_map:
                candidates.extend(
                    symbol_map[symbol]
                )


        candidates = list(
            dict.fromkeys(candidates)
        )


        if len(candidates) == 1:

            status = "mapped"
            expression_column = candidates[0]

        elif len(candidates) == 0:

            status = "not_found"
            expression_column = None

        else:

            status = "ambiguous"
            expression_column = None


        mapping_rows.append(
            {
                "human_symbol":
                    row["human_symbol"],

                "human_ensembl_id":
                    row["human_ensembl_id"],

                "score_component":
                    row["score_component"],

                "score_weight":
                    row["score_weight"],

                "orthology_confidence":
                    row.get(
                        "orthology_confidence",
                        np.nan
                    ),

                "mapping_status":
                    status,

                "expression_column":
                    expression_column,

                "n_expression_candidates":
                    len(candidates),
            }
        )


    return pd.DataFrame(mapping_rows)


def safe_spearman(x, y):
    """
    """

    x = pd.Series(x)
    y = pd.Series(y)

    ok = (
        x.notna()
        & y.notna()
    )

    if ok.sum() < 3:
        return np.nan, np.nan

    r, p = spearmanr(
        x[ok],
        y[ok]
    )

    return r, p


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
        "ASTRA-Drug Step 08\n"
        "DepMap ASTS scoring\n"
        "============================================================"
    )

    for path in [
        SIGNATURE_FILE,
        EXPRESSION_FILE,
        MODEL_FILE,
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
    # 1. Read human ASTS
    # --------------------------------------------------------

    print("\nStep 1: Reading human ASTS")


    signature = pd.read_csv(
        SIGNATURE_FILE,
        sep="\t"
    )


    required = [
        "human_symbol",
        "human_ensembl_id",
        "score_component",
        "score_weight",
    ]


    missing = [
        x for x in required
        if x not in signature.columns
    ]


    if missing:

        raise ValueError(
            "Missing signature columns: "
            + ", ".join(missing)
        )


    if len(signature) != 78:

        print(
            "WARNING: expected 78 primary human genes, found",
            len(signature)
        )


    n_up = (
        signature["score_component"]
        == "DHT_relative_up"
    ).sum()

    n_down = (
        signature["score_component"]
        == "DHT_relative_down"
    ).sum()


    print(
        f"Primary human genes: {len(signature)}"
    )

    print(
        f"  up   : {n_up}"
    )

    print(
        f"  down : {n_down}"
    )


    # --------------------------------------------------------
    # 2. Read expression header
    # --------------------------------------------------------
    
    print("\nStep 2: Reading DepMap expression header")
    
    
    header = pd.read_csv(
        EXPRESSION_FILE,
        nrows=0
    )


    columns = list(
        header.columns
    )


    if DEFAULT_PROFILE_COLUMN not in columns:

        raise ValueError(
            f"{DEFAULT_PROFILE_COLUMN} was not found in expression file."
        )


    model_id_col = detect_model_id_column(
        columns
    )


    metadata_columns = {
        model_id_col,
        "SequencingID",
        DEFAULT_PROFILE_COLUMN,
    }


    expression_columns = [
        x for x in columns
        if x not in metadata_columns
    ]


    print(
        "Model ID column:",
        model_id_col
    )

    print(
        "Default profile column:",
        DEFAULT_PROFILE_COLUMN
    )

    print(
        "Expression features:",
        len(expression_columns)
    )


    # --------------------------------------------------------
    # 3. Map ASTS genes to DepMap columns
    # --------------------------------------------------------

    print("\nStep 3: Mapping ASTS genes")


    mapping = map_signature_to_expression(
        signature,
        expression_columns
    )


    mapping.to_csv(
        OUT_DIR
        / "08_signature_expression_mapping.tsv",
        sep="\t",
        index=False
    )


    mapped = mapping[
        mapping["mapping_status"]
        == "mapped"
    ].copy()


    mapped_up = mapped[
        mapped["score_component"]
        == "DHT_relative_up"
    ]


    mapped_down = mapped[
        mapped["score_component"]
        == "DHT_relative_down"
    ]


    print(
        f"Mapped primary genes: "
        f"{len(mapped)}/{len(signature)}"
    )

    print(
        f"  up   : {len(mapped_up)}/{n_up}"
    )

    print(
        f"  down : {len(mapped_down)}/{n_down}"
    )


    # --------------------------------------------------------
    # 4. High-confidence sensitivity set
    # --------------------------------------------------------

    if "orthology_confidence" in mapped.columns:

        high_conf = mapped[
            pd.to_numeric(
                mapped["orthology_confidence"],
                errors="coerce"
            )
            == HIGH_CONFIDENCE_VALUE
        ].copy()

    else:

        high_conf = mapped.copy()


    high_conf_up = high_conf[
        high_conf["score_component"]
        == "DHT_relative_up"
    ]


    high_conf_down = high_conf[
        high_conf["score_component"]
        == "DHT_relative_down"
    ]


    print(
        "High-confidence sensitivity genes:",
        len(high_conf)
    )


    # --------------------------------------------------------
    # 5. Proliferation proxy mapping
    # --------------------------------------------------------

    symbol_lookup = {}

    for col in expression_columns:

        symbol_lookup.setdefault(
            normalize_symbol(col),
            []
        ).append(col)


    proliferation_cols = []


    for gene in PROLIFERATION_PROXY_GENES:

        candidates = symbol_lookup.get(
            gene.upper(),
            []
        )

        if len(candidates) == 1:

            proliferation_cols.append(
                candidates[0]
            )


    print(
        "Proliferation proxy genes available:",
        len(proliferation_cols),
        "/",
        len(PROLIFERATION_PROXY_GENES)
    )


    # --------------------------------------------------------
    # 6. Score expression matrix chunk-by-chunk
    # --------------------------------------------------------

    print("\nStep 4: Calculating ASTS scores")


    primary_up_cols = list(
        mapped_up["expression_column"]
    )

    primary_down_cols = list(
        mapped_down["expression_column"]
    )


    hc_up_cols = list(
        high_conf_up["expression_column"]
    )

    hc_down_cols = list(
        high_conf_down["expression_column"]
    )


    score_rows = []


    reader = pd.read_csv(
        EXPRESSION_FILE,
        chunksize=CHUNK_SIZE
    )


    total_models = 0


    for chunk_number, chunk in enumerate(
        reader,
        start=1
    ):

        flag = (
            chunk[DEFAULT_PROFILE_COLUMN]
            .astype(str)
            .str.strip()
            .str.lower()
        )
        
        chunk = chunk[
            flag == DEFAULT_PROFILE_VALUE
        ].copy()
    
        if len(chunk) == 0:
            continue
        
        model_ids = chunk[
            model_id_col
        ].astype(str)


        expr = chunk[
            expression_columns
        ].apply(
            pd.to_numeric,
            errors="coerce"
        )


        ranks = expr.rank(
            axis=1,
            method="average",
            pct=True,
            na_option="keep"
        )


        for i in range(
            len(chunk)
        ):

            model_id = model_ids.iloc[i]


            # ----------------------------------------------
            # Primary
            # ----------------------------------------------

            up_values = ranks.iloc[i][
                primary_up_cols
            ]

            down_values = ranks.iloc[i][
                primary_down_cols
            ]


            n_up_available = int(
                up_values.notna().sum()
            )

            n_down_available = int(
                down_values.notna().sum()
            )


            up_coverage = (
                n_up_available
                / len(primary_up_cols)
                if primary_up_cols
                else 0
            )

            down_coverage = (
                n_down_available
                / len(primary_down_cols)
                if primary_down_cols
                else 0
            )


            if (
                up_coverage
                >= MIN_DIRECTION_COVERAGE
                and
                down_coverage
                >= MIN_DIRECTION_COVERAGE
            ):

                up_score = float(
                    up_values.mean()
                )

                down_score = float(
                    down_values.mean()
                )

                asts_score = (
                    up_score
                    - down_score
                )

            else:

                up_score = np.nan
                down_score = np.nan
                asts_score = np.nan


            # ----------------------------------------------
            # High-confidence sensitivity
            # ----------------------------------------------

            hc_up_values = ranks.iloc[i][
                hc_up_cols
            ]

            hc_down_values = ranks.iloc[i][
                hc_down_cols
            ]


            if (
                len(hc_up_cols) > 0
                and
                len(hc_down_cols) > 0
            ):

                hc_score = (
                    hc_up_values.mean()
                    - hc_down_values.mean()
                )

            else:

                hc_score = np.nan


            # ----------------------------------------------
            # Exploratory proliferation proxy
            # ----------------------------------------------

            if proliferation_cols:

                proliferation_proxy = float(
                    ranks.iloc[i][
                        proliferation_cols
                    ].mean()
                )

            else:

                proliferation_proxy = np.nan


            score_rows.append(
                {
                    "ModelID":
                        model_id,

                    "ASTS_score":
                        asts_score,

                    "ASTS_up_mean_rank":
                        up_score,

                    "ASTS_down_mean_rank":
                        down_score,

                    "ASTS_high_conf_score":
                        hc_score,

                    "n_up_available":
                        n_up_available,

                    "n_down_available":
                        n_down_available,

                    "up_coverage":
                        up_coverage,

                    "down_coverage":
                        down_coverage,

                    "proliferation_proxy":
                        proliferation_proxy,
                }
            )


        total_models += len(chunk)


        print(
            f"  chunk {chunk_number}: "
            f"processed {total_models} models"
        )


    scores = pd.DataFrame(
        score_rows
    )

    if scores["ModelID"].duplicated().any():
    
        duplicated_ids = (
            scores.loc[
                scores["ModelID"].duplicated(False),
                "ModelID"
            ]
            .unique()
        )
    
        raise RuntimeError(
            "Duplicated ModelID remains after default-profile filtering: "
            + ", ".join(duplicated_ids[:20])
        )

    mean_score = scores[
        "ASTS_score"
    ].mean()

    sd_score = scores[
        "ASTS_score"
    ].std(ddof=1)


    scores["ASTS_z"] = (
        scores["ASTS_score"]
        - mean_score
    ) / sd_score


    # --------------------------------------------------------
    # 7. Add Model metadata
    # --------------------------------------------------------

    print("\nStep 5: Adding Model metadata")


    model = pd.read_csv(
        MODEL_FILE
    )


    if "ModelID" not in model.columns:

        raise ValueError(
            "Model.csv does not contain ModelID."
        )


    useful_metadata = [
        "ModelID",
        "CellLineName",
        "StrippedCellLineName",
        "DepmapModelType",
        "OncotreeLineage",
        "OncotreePrimaryDisease",
        "OncotreeSubtype",
        "Sex",
        "Age",
        "AgeCategory",
        "PrimaryOrMetastasis",
    ]


    useful_metadata = [
        x for x in useful_metadata
        if x in model.columns
    ]


    model_small = model[
        useful_metadata
    ].copy()


    scores = scores.merge(
        model_small,
        on="ModelID",
        how="left"
    )


    scores.to_csv(
        OUT_DIR
        / "08_DepMap_ASTS_scores.tsv",
        sep="\t",
        index=False
    )


    # --------------------------------------------------------
    # 8. Ranked top/bottom models
    # --------------------------------------------------------

    score_valid = scores.dropna(
        subset=["ASTS_score"]
    ).copy()


    score_valid = score_valid.sort_values(
        "ASTS_score",
        ascending=False
    )


    score_valid.head(30).to_csv(
        OUT_DIR
        / "08_top30_ASTS_models.tsv",
        sep="\t",
        index=False
    )


    score_valid.tail(30).sort_values(
        "ASTS_score"
    ).to_csv(
        OUT_DIR
        / "08_bottom30_ASTS_models.tsv",
        sep="\t",
        index=False
    )


    # --------------------------------------------------------
    # 9. Lineage summary
    # --------------------------------------------------------

    if "OncotreeLineage" in scores.columns:

        lineage_summary = (
            score_valid
            .dropna(
                subset=["OncotreeLineage"]
            )
            .groupby(
                "OncotreeLineage"
            )["ASTS_score"]
            .agg(
                n="count",
                mean="mean",
                median="median",
                sd="std",
                q25=lambda x:
                    x.quantile(0.25),
                q75=lambda x:
                    x.quantile(0.75),
            )
            .reset_index()
        )


        lineage_summary = (
            lineage_summary
            .sort_values(
                "median",
                ascending=False
            )
        )


        lineage_summary.to_csv(
            OUT_DIR
            / "08_lineage_summary.tsv",
            sep="\t",
            index=False
        )

    else:

        lineage_summary = pd.DataFrame()


    # --------------------------------------------------------
    # 10. Sex summary
    # --------------------------------------------------------

    sex_test_p = np.nan


    if "Sex" in scores.columns:

        sex_summary = (
            score_valid
            .dropna(
                subset=["Sex"]
            )
            .groupby(
                "Sex"
            )["ASTS_score"]
            .agg(
                n="count",
                mean="mean",
                median="median",
                sd="std",
            )
            .reset_index()
        )


        sex_summary.to_csv(
            OUT_DIR
            / "08_sex_summary.tsv",
            sep="\t",
            index=False
        )


        female = score_valid.loc[
            score_valid["Sex"]
            .astype(str)
            .str.lower()
            == "female",
            "ASTS_score"
        ]


        male = score_valid.loc[
            score_valid["Sex"]
            .astype(str)
            .str.lower()
            == "male",
            "ASTS_score"
        ]


        if (
            len(female) >= 2
            and
            len(male) >= 2
        ):

            _, sex_test_p = mannwhitneyu(
                female,
                male,
                alternative="two-sided"
            )

    else:

        sex_summary = pd.DataFrame()


    # --------------------------------------------------------
    # 11. Correlations
    # --------------------------------------------------------

    hc_r, hc_p = safe_spearman(
        scores["ASTS_score"],
        scores["ASTS_high_conf_score"]
    )


    proliferation_r, proliferation_p = (
        safe_spearman(
            scores["ASTS_score"],
            scores[
                "proliferation_proxy"
            ]
        )
    )


    # --------------------------------------------------------
    # 12. Figures
    # --------------------------------------------------------

    print("\nStep 6: Creating figures")


    # ASTS distribution
    fig, ax = plt.subplots(
        figsize=(7, 5)
    )

    ax.hist(
        score_valid["ASTS_score"],
        bins=40
    )

    ax.set_xlabel(
        "ASTS score"
    )

    ax.set_ylabel(
        "Number of DepMap models"
    )

    ax.set_title(
        "DepMap ASTS score distribution"
    )

    fig.tight_layout()

    fig.savefig(
        FIG_DIR
        / "01_ASTS_score_distribution.png",
        dpi=200
    )

    plt.close(fig)


    # Primary vs high-confidence
    fig, ax = plt.subplots(
        figsize=(6, 6)
    )

    ax.scatter(
        scores["ASTS_score"],
        scores["ASTS_high_conf_score"],
        s=10,
        alpha=0.6
    )

    ax.set_xlabel(
        "Primary ASTS score"
    )

    ax.set_ylabel(
        "High-confidence ASTS score"
    )

    ax.set_title(
        f"Primary vs high-confidence\n"
        f"Spearman r={hc_r:.3f}"
    )

    fig.tight_layout()

    fig.savefig(
        FIG_DIR
        / "02_primary_vs_high_confidence.png",
        dpi=200
    )

    plt.close(fig)


    # Proliferation proxy
    fig, ax = plt.subplots(
        figsize=(6, 6)
    )

    ax.scatter(
        scores["proliferation_proxy"],
        scores["ASTS_score"],
        s=10,
        alpha=0.6
    )

    ax.set_xlabel(
        "Exploratory proliferation proxy"
    )

    ax.set_ylabel(
        "ASTS score"
    )

    ax.set_title(
        f"ASTS vs proliferation proxy\n"
        f"Spearman r={proliferation_r:.3f}"
    )

    fig.tight_layout()

    fig.savefig(
        FIG_DIR
        / "03_ASTS_vs_proliferation_proxy.png",
        dpi=200
    )

    plt.close(fig)


    # Lineage boxplot
    if not lineage_summary.empty:

        keep_lineages = (
            lineage_summary.loc[
                lineage_summary["n"]
                >= MIN_LINEAGE_N,
                "OncotreeLineage"
            ]
            .tolist()
        )


        line_df = score_valid[
            score_valid[
                "OncotreeLineage"
            ].isin(
                keep_lineages
            )
        ].copy()


        if len(keep_lineages) > 1:

            order = (
                line_df
                .groupby(
                    "OncotreeLineage"
                )["ASTS_score"]
                .median()
                .sort_values()
                .index
                .tolist()
            )


            box_data = [
                line_df.loc[
                    line_df[
                        "OncotreeLineage"
                    ] == x,
                    "ASTS_score"
                ].values
                for x in order
            ]


            fig, ax = plt.subplots(
                figsize=(
                    8,
                    max(
                        6,
                        0.3 * len(order)
                    )
                )
            )


            ax.boxplot(
                box_data,
                vert=False,
                labels=order,
                showfliers=False
            )


            ax.set_xlabel(
                "ASTS score"
            )

            ax.set_ylabel(
                "Oncotree lineage"
            )

            ax.set_title(
                "ASTS score by lineage"
            )

            fig.tight_layout()

            fig.savefig(
                FIG_DIR
                / "04_ASTS_by_lineage.png",
                dpi=200
            )

            plt.close(fig)


    # --------------------------------------------------------
    # 13. Summary
    # --------------------------------------------------------

    summary_file = (
        OUT_DIR
        / "08_summary.txt"
    )


    with open(
        summary_file,
        "w"
    ) as f:

        print(
            "ASTRA-Drug Step 08 summary",
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
            f"DepMap release: "
            f"{DEPMAP_RELEASE}",
            file=f
        )


        print(
            f"Expression file: "
            f"{EXPRESSION_FILE}",
            file=f
        )


        print(
            "",
            file=f
        )


        print(
            "Human ASTS primary",
            file=f
        )

        print(
            "------------------",
            file=f
        )

        print(
            f"Input genes: "
            f"{len(signature)}",
            file=f
        )

        print(
            f"Mapped genes: "
            f"{len(mapped)}",
            file=f
        )

        print(
            f"DHT-relative-up mapped: "
            f"{len(mapped_up)}/{n_up}",
            file=f
        )

        print(
            f"DHT-relative-down mapped: "
            f"{len(mapped_down)}/{n_down}",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "DepMap scoring",
            file=f
        )

        print(
            "--------------",
            file=f
        )

        print(
            f"Models processed: "
            f"{len(scores)}",
            file=f
        )

        print(
            f"Models with valid ASTS score: "
            f"{scores['ASTS_score'].notna().sum()}",
            file=f
        )

        print(
            f"Mean ASTS score: "
            f"{scores['ASTS_score'].mean():.6f}",
            file=f
        )

        print(
            f"SD ASTS score: "
            f"{scores['ASTS_score'].std():.6f}",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "Sensitivity",
            file=f
        )

        print(
            "-----------",
            file=f
        )

        print(
            f"High-confidence genes: "
            f"{len(high_conf)}",
            file=f
        )

        print(
            f"Primary vs high-confidence "
            f"Spearman r: {hc_r:.6f}",
            file=f
        )

        print(
            f"Primary vs high-confidence "
            f"p: {hc_p:.6g}",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "Exploratory proliferation check",
            file=f
        )

        print(
            "-------------------------------",
            file=f
        )

        print(
            f"Proxy genes available: "
            f"{len(proliferation_cols)}/"
            f"{len(PROLIFERATION_PROXY_GENES)}",
            file=f
        )

        print(
            f"ASTS vs proliferation proxy "
            f"Spearman r: "
            f"{proliferation_r:.6f}",
            file=f
        )

        print(
            f"p: "
            f"{proliferation_p:.6g}",
            file=f
        )

        print(
            "",
            file=f
        )


        print(
            "Sex",
            file=f
        )

        print(
            "---",
            file=f
        )

        if not sex_summary.empty:

            print(
                sex_summary.to_string(
                    index=False
                ),
                file=f
            )

            print(
                f"\nFemale vs Male "
                f"Mann-Whitney p: "
                f"{sex_test_p:.6g}",
                file=f
            )

        else:

            print(
                "Sex metadata unavailable.",
                file=f
            )


        print(
            "",
            file=f
        )


        if not lineage_summary.empty:

            print(
                "Top 10 lineages by median ASTS",
                file=f
            )

            print(
                "------------------------------",
                file=f
            )

            print(
                lineage_summary.head(
                    10
                ).to_string(
                    index=False
                ),
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
            f"Human ASTS: "
            f"{sha256_file(SIGNATURE_FILE)}",
            file=f
        )

        print(
            f"Expression: "
            f"{sha256_file(EXPRESSION_FILE)}",
            file=f
        )

        print(
            f"Model metadata: "
            f"{sha256_file(MODEL_FILE)}",
            file=f
        )


    print(
        "\n============================================================"
    )

    print(
        "08 completed successfully"
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
