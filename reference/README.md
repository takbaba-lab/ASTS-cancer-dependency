# Reference files

This directory is reserved for local reference files used for mouse RNA-seq processing.

Large reference FASTA files and STAR genome indices are intentionally not stored in this GitHub repository.

## Reference genome

The RNA-seq analysis used an mm10-based reference augmented with:

- ERCC spike-in sequences
- EGFP
- mCherry

The corresponding annotation file contains the same augmented reference components.

The default directory layout is:

    reference/
    |-- ucsc_mm10_ERCC_EGFP_mcherry.fa
    |-- mm10_ERCC_EGFP_mcherry.gtf
    `-- star_index/

## STAR index

The original analysis used STAR version 2.7.11b.

The STAR index should be generated from the same augmented genome FASTA and GTF used for the RNA-seq analysis.

The default expected location is:

    reference/star_index/

Alternatively, an external STAR index can be specified with:

    export ASTRA_STAR_INDEX=/path/to/star_index

## GTF file

The default expected GTF location is:

    reference/mm10_ERCC_EGFP_mcherry.gtf

An external annotation file can be specified with:

    export ASTRA_GTF=/path/to/mm10_ERCC_EGFP_mcherry.gtf

## FASTA file

The genome FASTA is required when rebuilding the STAR index.

The recommended local filename is:

    reference/ucsc_mm10_ERCC_EGFP_mcherry.fa

## Reproducibility

Users reproducing the analysis from raw FASTQ files should ensure that:

1. the FASTA and GTF describe the same augmented reference;
2. the STAR index is generated from those files; and
3. the STAR version and major mapping parameters are compatible with those specified in the workflow scripts.

Reference files are not redistributed in this repository because of their size and because some components originate from external resources.
