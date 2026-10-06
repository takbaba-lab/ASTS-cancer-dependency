# ASTRA-Drug analysis run order

This document describes the execution order and major dependencies of the analyses used in the ASTS study.

The workflow is organized into Steps 01-20. On the Genkai system, the corresponding PJM job scripts can be submitted with `pjsub`.

## General setup

Run commands from the repository root:

```bash
cd /path/to/repository
export ASTRA_PROJECT_DIR="$PWD"

For steps requiring STAR or featureCounts references, either place the reference files under the default reference/ directory or define:
export ASTRA_STAR_INDEX=/path/to/star_index
export ASTRA_GTF=/path/to/mm10_ERCC_EGFP_mcherry.gtf

The repository does not redistribute large raw sequencing files or third-party pharmacogenomic datasets. See [`data/README.md`](data/README.md) and [`reference/README.md`](reference/README.md).
Workflow summary
StepPurposeAnalysis script(s)PJM job
01Raw FASTQ quality controlscripts/astra_map_env.shjobs/01_fastq_qc.pjm
02Read trimming and post-trimming QCscripts/02_trimmomatic_qc.shjobs/02_trimmomatic_qc.pjm
03STAR alignmentscripts/03_STAR_mapping.shjobs/03_STAR_mapping.pjm
04Gene-level quantificationscripts/04_featureCounts.shjobs/04_featureCounts.pjm
05DESeq2 and ERCC-based normalizationscripts/05_normalization_DESeq2.Rjobs/05_normalization_DESeq2.pjm
06Build mouse DHT-associated signaturescripts/06_build_ASTS_signature.Rjobs/06_build_ASTS_signature.pjm
06bFreeze the source ASTS signaturescripts/06b_freeze_ASTS_v1.Rjobs/06b_freeze_ASTS_v1.pjm
07Map mouse signature to human orthologsscripts/07_map_ASTS_mouse_to_human.pyjobs/07_map_ASTS_mouse_to_human.pjm
08bCalculate ASTS scores in DepMapscripts/08b_score_DepMap_ASTS_default_profile.pyjobs/08b_score_DepMap_ASTS_default_profile.pjm
09bTest ASTS associations in PRISMscripts/09b_PRISM_ASTS_association.pyjobs/09b_PRISM_ASTS_association.pjm
10Validate drug associations in GDSC2scripts/10_validate_GDSC2_ASTS.pyjobs/10_validate_GDSC2_ASTS.pjm
11Integrate drug-response screensscripts/11_cross_screen_mechanism.pyjobs/11_cross_screen_mechanism.pjm
12Evaluate drug-class robustnessscripts/12_drug_class_robustness.pyjobs/12_drug_class_robustness.pjm
13Adjust for canonical oncogenic driversscripts/13_driver_confounder_analysis_fixed_v2.pyjobs/13_driver_confounder_analysis_fixed_v2.pjm
14aTranscriptome-wide ASTS associationscripts/14a_ASTS_transcriptome.pyjobs/14a_ASTS_transcriptome.pjm
14bGene-set enrichment analysisscripts/14b_ASTS_GSEA.Rjobs/14b_ASTS_GSEA.pjm
15Cell-cycle-related sensitivity analysesscripts/15a_export_cellcycle_genes.R, scripts/15b_cellcycle_confounder.py jobs/15_cellcycle_confounder.pjm
16Targeted CRISPR dependency validationscripts/16_CRISPR_dependency_validation.pyjobs/16_CRISPR_dependency_validation.pjm
17Genome-wide CRISPR dependency analysisscripts/17a_CRISPR_genomewide_ASTS.py, scripts/17b_CRISPR_dependency_GSEA.R jobs/17_CRISPR_genomewide.pjm
18Mitochondrial dependency specificityscripts/18a_export_mito_gene_sets.R, scripts/18b_mito_dependency_specificity.py, scripts/18c_burden_adjusted_GSEA.R jobs/18_mito_dependency_specificity.pjm
19Mitochondrial drug-response analysisscripts/19_mito_drug_validation.pyjobs/19_mito_drug_validation.pjm
19B-CScreen-specific and CTRPv2 mitochondrial drug analysesscripts/19b_screen_specific_mito_drugs.py, scripts/19c1_extract_CTRPv2.R, scripts/19c2_CTRPv2_mito_drug_validation.py jobs/19BC_mito_drug_validation.pjm
20Leave-one-lineage-out and donor-sex robustness analysesscripts/20_transportability_audit.pyjobs/20_transportability_audit.pjm


Detailed execution order
Step 01: Raw FASTQ quality control
Purpose
Perform initial integrity and quality checks on the mouse zF RNA-seq FASTQ files.
Required inputs
- Raw paired-end FASTQ files
- metadata/sample_metadata.tsv
- FastQC and MultiQC environment
Execution
pjsub jobs/01_fastq_qc.pjm

The job uses scripts/astra_map_env.sh to configure the micromamba environment.
Principal outputs
Quality-control reports under the project results/qc/ hierarchy.
Step 02: Read trimming and QC
Purpose
Trim sequencing reads and perform post-trimming quality control.
Depends on
Step 01.
Script
scripts/02_trimmomatic_qc.sh

Execution
pjsub jobs/02_trimmomatic_qc.pjm

Principal outputs
Trimmed paired-end FASTQ files and associated QC reports.
Step 03: STAR alignment
Purpose
Align trimmed reads to the mm10-based reference augmented with ERCC, EGFP, and mCherry sequences.
Depends on
Step 02 and a prepared STAR genome index.
Execution
pjsub jobs/03_STAR_mapping.pjm

Principal outputs
STAR alignment files and sorted BAM files.
Step 04: featureCounts
Purpose
Generate gene-level count matrices using the reference GTF.
Two counting configurations were retained to support comparison with the published analysis and strand-aware quantification.
Depends on
Step 03.
Execution
pjsub jobs/04_featureCounts.pjm

Principal outputs
Gene-level count matrices under data/processed/counts/.
Step 05: normalization and differential-expression framework
Purpose
Generate normalized expression matrices using both standard DESeq2 normalization and ERCC-based normalization.
These complementary normalization strategies are used to distinguish relative transcriptional remodeling from changes in global RNA output.
Depends on
Step 04.
Execution
pjsub jobs/05_normalization_DESeq2.pjm

Principal outputs
Normalized count matrices and differential-expression results used in downstream ASTS construction.
Step 06: construct the mouse source signature
Purpose
Identify the cross-sex DHT-associated transcriptional response in mouse zF cells and construct the source signature used for ASTS.
Depends on
Step 05.
Execution
pjsub jobs/06_build_ASTS_signature.pjm

Step 06b: freeze ASTS source signature
Purpose
Freeze the mouse source signature before downstream human pharmacogenomic analyses.
This step is important because downstream drug-response and CRISPR datasets were not used to select the source signature.
Depends on
Step 06.
Execution
pjsub jobs/06b_freeze_ASTS_v1.pjm

Principal output
Frozen mouse ASTS source-signature files used by Step 07.
Step 07: mouse-to-human ortholog mapping
Purpose
Map the frozen mouse ASTS signature to one-to-one human orthologs.
Depends on
Step 06b.
Execution
pjsub jobs/07_map_ASTS_mouse_to_human.pjm

Principal output
The human ASTS signature used for cancer cell-line scoring.
The primary mapped signature contains the human orthologs retained for the downstream ASTS calculation.
Step 08b: calculate ASTS in DepMap
Purpose
Project the human ASTS signature into DepMap cancer cell lines using within-cell-line expression ranks.
Required external input
DepMap Public 26Q1 expression and model-annotation data.
Depends on
Step 07.
Execution
pjsub jobs/08b_score_DepMap_ASTS_default_profile.pjm

Principal outputs
Cell-line-level ASTS scores and associated model annotations.
Step 09b: PRISM association analysis
Purpose
Test associations between ASTS and drug response in the PRISM Secondary screen.
Required external input
PRISM Secondary drug-response data.
Depends on
Step 08b.
Execution
pjsub jobs/09b_PRISM_ASTS_association.pjm

Step 10: GDSC2 validation
Purpose
Evaluate whether drug-response associations identified in PRISM show concordant effects in GDSC2.
Required external input
GDSC2 drug-response data.
Depends on
Step 08b.
Execution
pjsub jobs/10_validate_GDSC2_ASTS.pjm

Step 11: cross-screen integration
Purpose
Integrate PRISM and GDSC2 results and evaluate cross-screen consistency of drug-response phenotypes.
Depends on
Steps 09b and 10.
Execution
pjsub jobs/11_cross_screen_mechanism.pjm

Step 12: drug-class robustness
Purpose
Summarize and evaluate ASTS-associated drug-response patterns at the drug-class and pathway level.
Depends on
Step 11.
Execution
pjsub jobs/12_drug_class_robustness.pjm

Step 13: oncogenic-driver adjustment
Purpose
Evaluate whether the major drug-response associations remain evident after adjustment for canonical oncogenic driver variables.
Depends on
Steps 08b-12 and the required genomic annotations.
Execution
pjsub jobs/13_driver_confounder_analysis_fixed_v2.pjm

Step 14a: transcriptome-wide ASTS association
Purpose
Characterize gene-expression features associated with ASTS across cancer cell lines.
Required external input
DepMap Public 26Q1 expression data.
Depends on
Step 08b.
Execution
pjsub jobs/14a_ASTS_transcriptome.pjm

Step 14b: ASTS gene-set enrichment analysis
Purpose
Perform pathway-level interpretation of the transcriptome-wide ASTS associations.
Depends on
Step 14a.
Execution
pjsub jobs/14b_ASTS_GSEA.pjm

Step 15: cell-cycle sensitivity analyses
Purpose
Evaluate the relationship between ASTS and major cell-cycle-related transcriptional programs and assess the robustness of drug-response associations after adjustment for these variables.
Depends on
Steps 08b, 09b, 10, and 14.
Execution
pjsub jobs/15_cellcycle_confounder.pjm

The job executes:
1. 15a_export_cellcycle_genes.R
2. 15b_cellcycle_confounder.py
Step 16: targeted CRISPR dependency validation
Purpose
Test selected signaling dependencies associated with ASTS, including EGFR and downstream pathway genes.
Required external input
DepMap Public 26Q1 CRISPR gene-effect data.
Depends on
Step 08b.
Execution
pjsub jobs/16_CRISPR_dependency_validation.pjm

Step 17: genome-wide CRISPR dependency analysis
Purpose
Identify genome-wide gene dependencies associated with ASTS and perform pathway-level enrichment analysis.
Required external input
DepMap Public 26Q1 CRISPR gene-effect data.
Depends on
Step 08b.
Execution
pjsub jobs/17_CRISPR_genomewide.pjm

The job executes:
1. 17a_CRISPR_genomewide_ASTS.py
2. 17b_CRISPR_dependency_GSEA.R
Step 18: mitochondrial dependency specificity
Purpose
Test whether the genome-wide CRISPR signal is selectively enriched for mitochondrial dependencies.
The analyses include mitochondrial gene-set definition, dependency-specificity analysis, adjustment for global dependency burden, and matched-gene permutation analyses.
Depends on
Step 17.
Execution
pjsub jobs/18_mito_dependency_specificity.pjm

The job executes:
1. 18a_export_mito_gene_sets.R
2. 18b_mito_dependency_specificity.py
3. 18c_burden_adjusted_GSEA.R
Step 19: mechanism-guided mitochondrial drug analysis
Purpose
Perform post hoc, mechanism-guided pharmacologic evaluation of mitochondrial vulnerabilities suggested by the CRISPR dependency analysis.
This analysis should not be interpreted as a systematic screen for broad mitochondrial or oxidative-phosphorylation inhibitor sensitivity.
Depends on
Steps 17 and 18.
Execution
pjsub jobs/19_mito_drug_validation.pjm

Step 19B-C: screen-specific and CTRPv2 analyses
Purpose
Evaluate mitochondrial drug associations in individual pharmacogenomic screens and perform the CTRPv2 analysis.
Required external input
CTRPv2 data accessed through the ORCESTRA PharmacoSet resource.
Depends on
Steps 08b, 18, and 19.
Execution
pjsub jobs/19BC_mito_drug_validation.pjm

The job executes:
1. 19b_screen_specific_mito_drugs.py
2. 19c1_extract_CTRPv2.R
3. 19c2_CTRPv2_mito_drug_validation.py
Step 20: robustness and transportability audit
Purpose
Evaluate whether the principal ASTS-associated phenotypes are attributable to a single major cancer lineage and test for evidence of donor-sex interaction.
The leave-one-lineage-out analysis is intended to assess robustness to major lineage composition; it should not be interpreted as establishing general clinical transportability.
The donor-sex analysis tests for detectable interaction effects and should not be interpreted as demonstrating equivalence between male- and female-derived cell lines.
Depends on
The major outputs from Steps 08b-19.
Execution
pjsub jobs/20_transportability_audit.pjm

Dependency overview
The main workflow can be summarized as:
Mouse zF RNA-seq
    |
    v
Steps 01-05
    |
    v
Steps 06-06b
Mouse DHT-associated source signature
    |
    v
Step 07
Human ASTS signature
    |
    v
Step 08b
ASTS scores across cancer cell lines
    |
    +-----------------------------+
    |                             |
    v                             v
Steps 09b-13                  Steps 14-15
Drug-response analyses        Transcriptomic/state analyses
    |
    +-----------------------------+
    |
    v
Steps 16-18
CRISPR dependency analyses
    |
    v
Steps 19 / 19B-C
Mechanism-guided pharmacology
    |
    v
Step 20
Lineage and donor-sex robustness analyses

Notes on reproducibility
1. The ASTS source signature is frozen before the downstream pharmacogenomic and CRISPR analyses.
2. ASTS is calculated from relative expression ranks and should not be interpreted as a measure of absolute RNA abundance.
3. Third-party datasets are not redistributed in this repository.
4. Large sequencing intermediates such as BAM files and STAR indices are intentionally excluded.
5. The archived Zenodo dataset will contain selected processed and derived data required to reproduce the reported figures and statistical results without redistributing restricted third-party source files.
6. The public code uses portable project paths and does not require the original Genkai user directory structure.
Archived release
The final manuscript-associated code release will be archived in Zenodo.
Software DOI: to be added
Processed and derived data will be provided in a separate Zenodo dataset record.
Dataset DOI: to be added
