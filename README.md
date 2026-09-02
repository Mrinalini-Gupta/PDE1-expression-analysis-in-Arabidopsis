# **Repository usage guide**

This README describes only how to set up, navigate and run this repository.

Folder structure
```text
PDE1-expression-analysis-in-Arabidopsis/
├── data/
│   ├── raw/          # downloaded input files
│   └── processed/    # intermediate files created by scripts
├── scripts/          # R scripts
├── results/
│   ├── figures/      # generated figures
│   └── tables/       # generated tables
├── README.md
└── .gitignore
```
Large source datasets should remain local and do not need to be uploaded to GitHub.

**Clone and configure**

git clone https://github.com/Mrinalini-Gupta/PDE1-expression-analysis-in-Arabidopsis.git
cd PDE1-expression-analysis-in-Arabidopsis

Several scripts contain BASE_DIR or base_dir. Change this to the repository location on the new computer:
```text
BASE_DIR <- "C:/Users/USERNAME/Documents/PDE1-expression-analysis-in-Arabidopsis"
```
Use forward slashes in R paths.

**R packages**

Install packages requested by the scripts if they are not already installed. Packages used across the repository include DESeq2, Seurat, SeuratObject, Matrix, tidyverse, ggplot2, dplyr, tidyr, tibble, purrr, readr, stringr, forcats, patchwork, scales, readxl, curl and gridExtra.

Example:
```text
install.packages(c("Seurat","tidyverse","patchwork","readxl","curl","gridExtra"))

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
```
## **Script 01**

scripts/01_bulk_RNAseq_DESeq2_analysis.R

Place in data/raw/:
```text
PRJEB25079_UncorrectedCounts.csv
E-MTAB-9694.sdrf.txt
```
Run:

source("scripts/01_bulk_RNAseq_DESeq2_analysis.R")

<p>Outputs are written to:<br> data/processed/<br> results/figures/<br> results/tables/<br>results/.</p>

## **Script 02**

scripts/02_final_generate_scRNA_figures.R

Place in data/raw/:
```text
GSE226826_combined_filtered.rds
```
If the local filename differs, update rds_path in the script.

Run:

source("scripts/02_final_generate_scRNA_figures.R")

<p>Outputs are written to:<br>

results/figures/<br>

results/tables/</p>

<p>Generated tables can include:<br>

Nobori_condition_time_celltype_detection.csv<br>
Nobori_candidate_cluster_statistics.csv<br>
Nobori_sample_level_detection.csv<br>
Nobori_condition_time_celltype_mean.csv<br>
Nobori_subcluster_*.csv</p>

## **Script 03**

scripts/03_final_generate_developmental_spatial_figures.R

Place the required developmental Seurat objects, Curio spatial object and supporting downloaded files under data/raw/.

Run:

source("scripts/03_final_generate_developmental_spatial_figures.R")

<p>Outputs are written to:<br>
results/figures/<br>
results/tables/</p>

<p>Generated tables can include:<br>

Lee_stage_detection_summary.csv<br>
Lee_celltype_detection_summary.csv<br>
Lee_PDE1_celltype_enrichment.csv<br>
Lee_655_subcluster_pseudobulk.csv<br>
Lee_global_PDE1_PR1_detection.csv<br>
Lee_Curio_cluster_summary.csv<br>
Lee_Curio_published_signature_summary.csv<br>
Lee_Flower_snRNA_PDE1_celltype_summary.csv<br>
Lee_Flower_PDE1_label_transfer_summary.csv</p>

## **Script 04**

scripts/04_generate_compact_supplementary_tables.R

<p>Run Scripts 01-03 first. The main inputs expected in results/tables/ are:<br>

bulk_DESeq2_candidate_results.csv<br>
Nobori_condition_time_celltype_detection.csv<br>
Nobori_candidate_cluster_statistics.csv<br>
Lee_stage_detection_summary.csv<br>
Lee_Curio_cluster_summary.csv<br>
Lee_Flower_PDE1_label_transfer_summary.csv</p>

Run:

source("scripts/04_generate_compact_supplementary_tables.R")

<p>Outputs:<br>

Supplementary_Table_S1_Bulk.csv<br>
Supplementary_Table_S1_Bulk.png<br>
Supplementary_Table_S2_Nobori_scRNA.csv<br>
Supplementary_Table_S2_Nobori_scRNA.png<br>
Supplementary_Table_S3_Lee_Developmental_Spatial.csv<br>
Supplementary_Table_S3_Lee_Developmental_Spatial.png</p>

## **Run order**

1. Clone the repository.
2. Download the required input files.
3. Place inputs in data/raw/ and required subdirectories.
4. Update BASE_DIR/base_dir in the scripts.
5. Install missing R packages.
6. Run Script 01.
7. Run Script 02.
8. Run the Flower integration workflow if its required CSV files are absent.
9. Run Script 03.
10. Run Script 04.
11. Check results/figures/ and results/tables/.

**Record the local R environment**
```text
writeLines(
  capture.output(sessionInfo()),
  file.path("results", "sessionInfo_final.txt")
)
```
**Package versions can be checked with:**
```text
packageVersion("DESeq2")
packageVersion("Seurat")
packageVersion("ggplot2")
```
**Files to keep local**

Example
```text
 .gitignore:

data/raw/*.rds
data/raw/*.rds.gz
data/raw/*.h5
data/raw/*.h5ad
data/raw/*.mtx
data/raw/*.mtx.gz
*.RData
.Rhistory
.RData
.Rproj.user/
```
Always run git status before committing files.

**Update GitHub**
```text
git status
git add README.md
git add scripts/
git add results/tables/
git commit -m "Update reproducibility files"
git push origin main
```
Avoid git add . until large local files have been excluded.

**Troubleshooting**

-File not found: check the repository path, filename, expected directory and the path specified in the script.

-Package not found: install the missing package and rerun the script.

-Object or column not found: confirm that the correct processed input object is being used.

-Flower input missing: run the Flower integration workflow and place its CSV outputs in the directory searched by Script 03.

-Supplementary-table input missing: run the preceding script that generates the required CSV before running Script 04.

-GitHub rejects a large file: remove it from Git tracking, add its path or extension to .gitignore, and keep it locally.
