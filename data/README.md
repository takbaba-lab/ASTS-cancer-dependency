# Data directory

This directory is reserved for local input data used by the ASTRA-Drug analysis workflow.

Large sequencing files and third-party datasets are intentionally not redistributed through this GitHub repository.

## Directory structure

The analysis scripts expect data to be organized approximately as follows:

    data/
    |-- raw/
    |   `-- fastq/
    |-- processed/
    |   |-- trimmed/
    |   |-- star/
    |   `-- counts/
    |-- external/
    |   |-- depmap/
    |   |-- prism/
    |   |-- gdsc/
    |   `-- ctrpv2/
    `-- interim/

Some directories are created automatically by the analysis scripts.

## Mouse RNA-seq data

The mouse adrenal zona fasciculata RNA-seq datasets used in this study are publicly available from the DNA Data Bank of Japan Sequence Read Archive:

- DRA015952
- DRA017476

Raw FASTQ files are not included in this repository.

For reproduction from raw sequencing data, place the FASTQ files under:

    data/raw/fastq/

The sample organization and experimental metadata are described in:

    metadata/sample_metadata.tsv

## DepMap

The study used data from DepMap Public 26Q1.

The analyses require the relevant expression, model annotation, and CRISPR gene-effect files.

Because DepMap distributes and updates these resources independently, the original files are not redistributed here.

Place locally downloaded DepMap files under:

    data/external/depmap/

The exact files used by each analysis step are specified in the corresponding scripts.

## PRISM

Drug-response analyses used the PRISM Secondary screen.

Place locally obtained PRISM input files under:

    data/external/prism/

The original PRISM files are not redistributed in this repository.

## GDSC2

GDSC2 drug-response data were used for cross-screen validation.

Place locally obtained GDSC2 files under:

    data/external/gdsc/

The original GDSC2 files are not redistributed here.

## CTRPv2

CTRPv2 data were accessed through the ORCESTRA PharmacoSet resource used in the study.

Place locally prepared CTRPv2 input files under:

    data/external/ctrpv2/

The extraction and downstream analyses are implemented in:

    scripts/19c1_extract_CTRPv2.R
    scripts/19c2_CTRPv2_mito_drug_validation.py

## Intermediate files

Temporary or derived files that are not intended for version control may be stored under:

    data/interim/

This directory is excluded from Git by `.gitignore`.

## Processed data

Selected processed and derived data required to reproduce the reported results will be deposited separately in Zenodo.

The Zenodo dataset DOI will be added before the final archived release.

Raw or licensed third-party datasets will not be redistributed through Zenodo unless redistribution is explicitly permitted by the original provider.
