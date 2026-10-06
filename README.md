# ASTRA-Drug: reproducibility repository for the ASTS study

This repository contains the analysis code and supporting metadata used to construct and evaluate the **androgen-suppressed transcriptional state (ASTS)** described in our study of transcriptional-state-associated vulnerabilities in cancer cell lines.

## Overview

The study begins with RNA-seq data from mouse adrenal zona fasciculata (zF) cells exposed to defined gonadal and hormonal conditions. We separated DHT-associated global RNA-output suppression from relative transcriptional remodeling and derived a cross-sex transcriptional signature.

The mouse signature was mapped to one-to-one human orthologs and projected into human cancer cell lines using within-cell-line expression ranks. This score is referred to as the **androgen-suppressed transcriptional state (ASTS)**.

ASTS was then integrated with pharmacogenomic and CRISPR dependency datasets to test whether this transcriptional state is associated with reproducible functional phenotypes across cancer cell lines.

### Important interpretation

ASTS represents similarity to the **relative DHT-associated transcriptional state derived from mouse adrenal zona fasciculata cells**.

ASTS should **not** be interpreted as a direct measure of:

- androgen exposure in cancer cells,
- current androgen receptor activity,
- absolute RNA-output suppression,
- NR5A1 activity in cancer cells, or
- a clinically validated treatment-response biomarker.

Associations identified using ASTS are functional associations and do not by themselves establish causality.

---

## Main analysis workflow

The analysis is organized into sequential steps:

1. Quality control of mouse zF RNA-seq data
2. Read trimming and STAR alignment
3. Gene-level quantification with featureCounts
4. Standard and ERCC-based normalization
5. Construction of the mouse DHT-associated transcriptional signature
6. Freezing of the cross-sex ASTS source signature
7. Mapping of mouse genes to one-to-one human orthologs
8. Calculation of ASTS scores in DepMap cancer cell lines
9. PRISM drug-response association analysis
10. GDSC2 validation
11. Cross-screen drug-response analysis
12. Drug-class robustness analysis
13. Adjustment for canonical oncogenic drivers
14. Transcriptome-wide ASTS association and gene-set analysis
15. Cell-cycle-related sensitivity analyses
16. Targeted CRISPR dependency analysis
17. Genome-wide CRISPR dependency analysis
18. Mitochondrial dependency specificity analyses
19. Mechanism-guided mitochondrial drug-response analyses
20. Leave-one-lineage-out and donor-sex robustness analyses

Detailed execution order and dependencies are described in [`RUN_ORDER.md`](RUN_ORDER.md).

---

## Repository structure

```text
.
 ├── README.md
 ├── RUN_ORDER.md
 ├── CITATION.cff
 ├── LICENSE
 ├── .gitignore
 │
 ├── scripts/
 │   ├── analysis scripts in Python, R, and shell
 │   └── astra_map_env.sh
 │
 ├── jobs/
 │   └── Fujitsu/PJM job scripts used on the Genkai system
 │
 ├── environment/
 │   ├── astra-map.yml
 │   ├── astra-py.yml
 │   ├── astra-r.yml
 │   └── package-lock exports
 │
 ├── metadata/
 │   └── sample_metadata.tsv
 │
 ├── data/
 │   └── README.md
 │
 ├── reference/
 │   └── README.md
 │
 └── results/
     └── README.md
```

Large raw datasets, third-party datasets, BAM files, STAR indices, and other machine-specific intermediate files are intentionally not stored in this repository.

## Data sources

### Mouse adrenal zona fasciculata RNA-seq

The mouse zF RNA-seq datasets reanalyzed in this study are publicly available from the DNA Data Bank of Japan Sequence Read Archive under:
- DRA015952
- DRA017476

The public repository does not redistribute the raw FASTQ files.

### DepMap

Human cancer cell-line expression, model annotation, and CRISPR gene-effect data were obtained from:
- DepMap Public 26Q1

Because DepMap datasets are maintained and distributed by DepMap, the original downloaded files are not redistributed in this repository.

### Pharmacogenomic datasets

Drug-response analyses used:
- PRISM Secondary screen
- GDSC2
- CTRPv2

CTRPv2 data were accessed through the ORCESTRA PharmacoSet resource used in the analysis.

Original third-party pharmacogenomic datasets are not redistributed here.

### Gene-set resources

Gene-set analyses use publicly available gene-set resources as specified in the analysis scripts and environment files. Users should obtain the relevant resources from their original providers when required.

### Mouse RNA-seq reference files

The mouse RNA-seq analysis used an mm10-based reference augmented with:
- ERCC spike-in sequences,
- EGFP, and
- mCherry.

The repository does not include large FASTA files or the STAR genome index.

By default, the public scripts expect reference resources under:
reference/
?????ucsc_mm10_ERCC_EGFP_mcherry.fa
?????mm10_ERCC_EGFP_mcherry.gtf
?????star_index/

Alternatively, external locations can be supplied where supported by environment variables such as:
export ASTRA_STAR_INDEX=/path/to/star_index
export ASTRA_GTF=/path/to/mm10_ERCC_EGFP_mcherry.gtf

See [`reference/README.md`](reference/README.md) for details.

### Computing environments

Three environment specifications are provided:
environment/astra-map.yml
environment/astra-py.yml
environment/astra-r.yml

These correspond to the principal software environments used for mapping, Python-based analyses, and R-based analyses.

Package-lock exports are also provided to document the analysis environment more completely.

The environment files were sanitized for public release and do not contain user-specific installation paths.

### Genkai job scripts

The original analyses were run in part on the Genkai supercomputer system at Kyushu University.
Genkai uses PJM/pjsub, not Slurm. Corresponding job scripts are included under:
jobs/

For example:
pjsub jobs/08b_score_DepMap_ASTS_default_profile.pjm

The public job scripts use portable project paths. When necessary, the repository root can be specified explicitly:
export ASTRA_PROJECT_DIR=/path/to/repository

The micromamba environment used by the initial QC workflow can similarly be configured through environment variables.

### Reproducing the analysis

The scripts are designed to be run sequentially according to the workflow described in [`RUN_ORDER.md`](RUN_ORDER.md).

A typical pattern is:
cd /path/to/repository
export ASTRA_PROJECT_DIR="$PWD"

followed by either direct execution of the analysis script or submission of the corresponding PJM job script.

Some analysis steps require third-party datasets that must first be downloaded from their original providers and placed in the expected local directory structure. These inputs are documented separately in [`data/README.md`](data/README.md).

The repository is intended to provide the analysis logic and reproducibility framework rather than to redistribute externally licensed or very large source datasets.

### Derived data and archived release

Processed and derived data underlying the analyses and figures are archived separately in Zenodo.

Dataset DOI: https://doi.org/10.5281/zenodo.23178149

A fixed archival release of the analysis code is also available in Zenodo.

Software DOI: https://doi.org/10.5281/zenodo.23179457

### Citation

If you use this code or the ASTS framework, please cite the associated article.

Article citation:
Miao Y, Baba T.
Manuscript citation to be added after publication.

Citation metadata for this repository will also be provided in [`CITATION.cff`](CITATION.cff).

### License

The analysis code in this repository is released under the MIT License.

Third-party datasets and resources remain subject to the licenses and terms of their original providers.

### Authors

Yingqi Miao
Kyushu University

Takashi Baba
Kyushu University

For correspondence regarding this repository or the study, please contact Takashi Baba.

### Status

This repository accompanies a manuscript currently in preparation/submission. File organization and documentation may be refined before the final archived release, but the analysis scripts corresponding to the reported study are frozen for reproducibility.
