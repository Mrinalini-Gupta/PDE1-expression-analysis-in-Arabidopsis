PDE1 expression analysis in Arabidopsis thaliana

Reproducible transcriptomic reanalysis of PDE1 (AT1G17330) across immune-associated and developmental contexts in Arabidopsis thaliana.

This repository accompanies an MSc Bioinformatics project investigating when and where PDE1 is detected using publicly available bulk RNA-seq, single-cell/single-nucleus RNA-seq and spatial transcriptomic datasets.

This README currently documents the bulk RNA-seq immune-elicitor analysis. Single-cell and spatial scripts can be added as separate reproducible workflows.

Bulk RNA-seq dataset

The bulk analysis reuses the RNA-seq dataset from:

Bjornson M, Pimprikar P, Nürnberger T, Zipfel C. The transcriptional landscape of Arabidopsis thaliana pattern-triggered immunity. Nature Plants. 2021;7:579–586. DOI: 10.1038/s41477-021-00874-5.

Data accession: ArrayExpress E-MTAB-9694.

The analysis focuses on wild-type Col samples exposed to immune elicitors including flg22, elf18, nlp20, oligogalacturonides (OGs), Pep1, chitooctaose and 3-OH-FA.

Statistical approach

The analysis starts from raw integer RNA-seq counts and fits a new DESeq2 model specifically for the wild-type Col samples.

For formal inference, stimulus and time are combined into a single group:

group = stimulus × time
DESeq2 design = ~ group

Each immune-elicitor/time combination is compared with the matched-time mock group.

For example:

flg22_90 min vs mock_90 min

For every contrast, DESeq2 reports:

log2 fold change;

standard error;

Wald test P value;

Benjamini-Hochberg adjusted P value (FDR);

estimated fold change;

95% Wald confidence interval.

An adjusted P value below 0.05 is treated as statistical evidence for differential expression.

The expression trajectory figures use DESeq2-normalised counts for visualisation, while the treatment-effect figures use the formal DESeq2 model estimates.

Candidate genes

PDE1 is compared with PR1 (AT2G14610) as a canonical defence-responsive comparator and with eight additional candidate genes examined for HD-domain/PDEase-like relationships:

AT5G40270

AT5G40290

AT1G14520

AT2G19800

AT4G02260

AT5G56640

AT1G26160

AT2G23820

Repository structure

```text
PDE1-expression-analysis-in-Arabidopsis/
├── README.md
├── .gitignore
├── scripts/
│   └── 01_bulk_RNAseq_DESeq2_analysis.R
├── data/
│   └── README.md
└── results/
    ├── figures/
    │   ├── Fig1_PDE1_expression_with_statistics.png
    │   ├── Fig2_candidate_expression_overview.png
    │   ├── Fig3_DESeq2_log2FC_vs_mock.png
    │   ├── Fig4_DESeq2_log2FC_heatmap.png
    │   └── Fig5_ranked_max_DESeq2_effect.png
    ├── tables/
    │   ├── bulk_DESeq2_candidate_results.csv
    │   └── bulk_candidate_peak_summary.csv
    ├── analysis_parameters.txt
    └── sessionInfo.txt
```
Raw data

Raw data are not committed to this repository.

For the bulk workflow, place or retain the following files in a local data directory:

PRJEB25079_UncorrectedCounts.csv
E-MTAB-9694.sdrf.txt

The analysis script accepts the raw-data directory as a command-line argument, so the repository does not depend on a user-specific absolute path.

Example local location used during development:

C:/Users/qo25519/Documents/PDE1_bulkRNA_analysis/PDE_candidate_bulkRNA/data

Software

The workflow was developed using:

R 4.5.1

RStudio 2026.06.0+242

DESeq2

tidyverse

ggplot2

patchwork

The exact package versions used in a completed run are written automatically to:

results/sessionInfo.txt

Running the bulk analysis

1. Clone the repository

git clone https://github.com/Mrinalini-Gupta/PDE1-expression-analysis-in-Arabidopsis.git
cd PDE1-expression-analysis-in-Arabidopsis

2. Create the project folders if needed

mkdir -p scripts data results/figures results/tables

On Windows PowerShell, the folders can also be created manually.

3. Save the R script

Place:

01_bulk_RNAseq_DESeq2_analysis.R

inside:

scripts/

4. Run from the repository root

Using the local raw-data directory:

Rscript scripts/01_bulk_RNAseq_DESeq2_analysis.R "C:/Users/qo25519/Documents/PDE1_bulkRNA_analysis/PDE_candidate_bulkRNA/data"

The script will:

locate the raw count file;

identify the Arabidopsis gene-ID column;

parse genotype, stimulus, time and replicate information from sample names;

retain wild-type Col samples for formal inference;

perform minimal low-count filtering;

fit a DESeq2 stimulus-time group model;

test every available elicitor/time combination against matched-time mock;

apply Benjamini-Hochberg FDR correction within each contrast;

create normalised-expression and inferential figures;

save two compact result tables and full session information.

Main output figures

Figure 1 shows PDE1 normalised expression over time with biological sample points, mean ± SE and FDR-supported matched-time contrasts.

Figure 2 shows normalised expression trajectories for PDE1, PR1 and the additional candidate genes.

Figure 3 shows formal DESeq2 log2 fold-change estimates versus matched-time mock, with 95% Wald confidence intervals and FDR status.

Figure 4 summarises DESeq2 log2 fold changes as a heatmap; an asterisk marks adjusted P < 0.05.

Figure 5 ranks the maximum positive DESeq2 effect estimate observed for each candidate. The ranking itself is descriptive; the selected contrast retains its DESeq2 confidence interval and FDR value.

Main tables

Only two analysis tables are intended for version control:

bulk_DESeq2_all_candidate_contrasts.csv

Contains all candidate-gene treatment-vs-matched-time-mock results, including log2FC, fold change, confidence intervals, raw P values and adjusted P values.

bulk_DESeq2_candidate_peak_summary.csv

Contains the maximum estimated positive DESeq2 response for each candidate together with the stimulus, time, fold change, confidence interval and statistical support for that selected contrast.

Citation

If this workflow or repository is reused, please cite the original dataset:

Bjornson M, Pimprikar P, Nürnberger T, Zipfel C. 2021. The transcriptional landscape of Arabidopsis thaliana pattern-triggered immunity. Nature Plants. 7:579–586. https://doi.org/10.1038/s41477-021-00874-5
