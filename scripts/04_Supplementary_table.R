# =============================================================================
# 04_generate_compact_supplementary_tables.R
#
# Purpose:
# Generate three compact supplementary tables for the MSc report:
#
#   Table S1 = Bulk RNA-seq
#   Table S2 = Nobori single-cell RNA-seq
#   Table S3 = Lee developmental + spatial transcriptomics
#
# Outputs:
#   results/tables/Supplementary_Table_S1_Bulk.csv
#   results/tables/Supplementary_Table_S1_Bulk.png
#
#   results/tables/Supplementary_Table_S2_Nobori_scRNA.csv
#   results/tables/Supplementary_Table_S2_Nobori_scRNA.png
#
#   results/tables/Supplementary_Table_S3_Lee_Developmental_Spatial.csv
#   results/tables/Supplementary_Table_S3_Lee_Developmental_Spatial.png
#
# =============================================================================


rm(list = ls())
options(stringsAsFactors = FALSE)


# =============================================================================
# 1. Packages
# =============================================================================

required_pkgs <- c(
  "dplyr",
  "tidyr",
  "readr",
  "stringr",
  "gridExtra",
  "grid"
)

missing_pkgs <- setdiff(
  required_pkgs,
  rownames(installed.packages())
)

if (length(missing_pkgs) > 0) {
  install.packages(missing_pkgs, dependencies = TRUE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(gridExtra)
  library(grid)
})


# =============================================================================
# 2. Paths
# =============================================================================

BASE_DIR <-
  "C:/Users/qo25519/Documents/PDE1_bulkRNA_analysis/PDE_candidate_bulkRNA/PDE1-expression-analysis-in-Arabidopsis"

TABLE_DIR <- file.path(
  BASE_DIR,
  "results",
  "tables"
)

dir.create(
  TABLE_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)


# =============================================================================
# 3. Helper functions
# =============================================================================

read_required <- function(filename) {
  
  path <- file.path(TABLE_DIR, filename)
  
  if (!file.exists(path)) {
    stop(
      "\nRequired file not found:\n",
      path,
      "\n\nRun the relevant analysis script first."
    )
  }
  
  message("Reading: ", filename)
  
  read_csv(
    path,
    show_col_types = FALSE
  )
}


read_optional <- function(filename) {
  
  path <- file.path(TABLE_DIR, filename)
  
  if (!file.exists(path)) {
    message("Optional file not found - skipping: ", filename)
    return(NULL)
  }
  
  message("Reading: ", filename)
  
  read_csv(
    path,
    show_col_types = FALSE
  )
}


fmt_p <- function(x) {
  
  ifelse(
    is.na(x),
    "",
    ifelse(
      x < 0.001,
      format(x, scientific = TRUE, digits = 2),
      sprintf("%.3f", x)
    )
  )
}


wrap_text <- function(x, width = 25) {
  
  vapply(
    x,
    function(z) {
      
      if (is.na(z)) return("")
      
      paste(
        strwrap(
          as.character(z),
          width = width
        ),
        collapse = "\n"
      )
    },
    character(1)
  )
}


save_table_png <- function(
    df,
    filename,
    title,
    note = NULL,
    base_size = 8,
    width_px = 4200
) {
  
  df_display <- df
  
  for (j in seq_along(df_display)) {
    
    if (is.character(df_display[[j]])) {
      
      df_display[[j]] <- wrap_text(
        df_display[[j]],
        width = 28
      )
    }
  }
  
  theme_tab <- ttheme_minimal(
    base_size = base_size,
    
    core = list(
      fg_params = list(
        hjust = 0,
        x = 0.02
      ),
      padding = unit(
        c(2.5, 3),
        "mm"
      )
    ),
    
    colhead = list(
      fg_params = list(
        fontface = "bold",
        hjust = 0,
        x = 0.02
      ),
      padding = unit(
        c(3, 3),
        "mm"
      )
    )
  )
  
  tab <- tableGrob(
    df_display,
    rows = NULL,
    theme = theme_tab
  )
  
  title_grob <- textGrob(
    title,
    x = 0,
    hjust = 0,
    gp = gpar(
      fontsize = 13,
      fontface = "bold"
    )
  )
  
  if (is.null(note)) {
    
    final_grob <- arrangeGrob(
      title_grob,
      tab,
      ncol = 1,
      heights = c(0.7, 12)
    )
    
  } else {
    
    note_grob <- textGrob(
      wrap_text(note, width = 150),
      x = 0,
      hjust = 0,
      gp = gpar(
        fontsize = 7.5
      )
    )
    
    final_grob <- arrangeGrob(
      title_grob,
      tab,
      note_grob,
      ncol = 1,
      heights = c(0.7, 12, 1.1)
    )
  }
  
  
  n_rows <- nrow(df_display)
  
  height_px <- max(
    2300,
    450 + n_rows * 115
  )
  
  
  png(
    filename = file.path(
      TABLE_DIR,
      filename
    ),
    width = width_px,
    height = height_px,
    res = 300
  )
  
  grid.draw(final_grob)
  
  dev.off()
  
  message(
    "Saved PNG: ",
    file.path(TABLE_DIR, filename)
  )
}


# =============================================================================
# 4. TABLE S1 — BULK RNA-seq
# =============================================================================

bulk <- read_required(
  "bulk_DESeq2_candidate_results.csv"
)


# -----------------------------------------------------------------------------
# Keep:
# - statistically supported PDE1 results
# - statistically supported PR1 results
# - statistically supported additional candidates
#
# This gives a compact report-level summary rather than every contrast.
# -----------------------------------------------------------------------------

S1 <- bulk %>%
  
  filter(
    !is.na(padj),
    padj < 0.05
  ) %>%
  
  mutate(
    Gene = gene_label,
    
    Elicitor = stimulus,
    
    `Time (min)` = as.character(time),
    
    `log2FC` = sprintf("%.2f", log2FC),
    
    `Fold change` = sprintf("%.2f", fold_change),
    
    `95% CI` = paste0(
      sprintf("%.2f", ci95_low_log2FC),
      " to ",
      sprintf("%.2f", ci95_high_log2FC)
    ),
    
    FDR = fmt_p(padj)
  ) %>%
  
  arrange(
    factor(
      Gene,
      levels = c(
        "PDE1",
        "PR1",
        setdiff(unique(Gene), c("PDE1", "PR1"))
      )
    ),
    Elicitor,
    as.numeric(`Time (min)`)
  ) %>%
  
  select(
    Gene,
    Elicitor,
    `Time (min)`,
    `log2FC`,
    `Fold change`,
    `95% CI`,
    FDR
  )

write_csv(
  S1,
  file.path(
    TABLE_DIR,
    "Supplementary_Table_S1_Bulk.csv"
  )
)


save_table_png(
  df = S1,
  
  filename =
    "Supplementary_Table_S1_Bulk.png",
  
  title =
    "Supplementary Table S1. Statistically supported bulk RNA-seq differential-expression results",
  
  note =
    "Values are DESeq2 matched-time treatment-versus-mock contrasts. FDR represents Benjamini-Hochberg adjusted P values. Only contrasts with FDR < 0.05 are shown."
)


# =============================================================================
# 5. TABLE S2 — NOBORI SINGLE-CELL RNA-seq
# =============================================================================

nobori_detection <- read_required(
  "Nobori_condition_time_celltype_detection.csv"
)

nobori_cluster <- read_required(
  "Nobori_candidate_cluster_statistics.csv"
)


# -----------------------------------------------------------------------------
# Part A:
# Broad-cell-type detection for PDE1 and PR1.
#
# To keep compact:
# retain rows where at least one positive cell was detected.
# -----------------------------------------------------------------------------

S2_detection <- nobori_detection %>%
  
  filter(
    gene %in% c("PDE1", "PR1"),
    n_positive > 0
  ) %>%
  
  mutate(
    Section = "Broad cell type",
    
    Gene = gene,
    
    Context = paste0(
      condition,
      " | ",
      time_plot,
      " | ",
      celltype
    ),
    
    `N cells` = as.character(n_total),
    
    `Positive cells` =
      as.character(n_positive),
    
    `% detected` =
      sprintf(
        "%.2f",
        percent_detected
      ),
    
    `Statistical / descriptive result` =
      "Descriptive"
  ) %>%
  
  select(
    Section,
    Gene,
    Context,
    `N cells`,
    `Positive cells`,
    `% detected`,
    `Statistical / descriptive result`
  )


# -----------------------------------------------------------------------------
# Part B:
# FDR-supported cluster associations.
#
# Restrict to PDE1 and PR1 for a compact report table.
# -----------------------------------------------------------------------------

S2_cluster <- nobori_cluster %>%
  
  filter(
    gene_label %in% c(
      "PDE1",
      "PR1"
    ),
    !is.na(FDR),
    FDR < 0.05
  ) %>%
  
  mutate(
    Section = "Major cluster",
    
    Gene = gene_label,
    
    Context = paste0(
      "Cluster ",
      cluster,
      " | ",
      dominant_celltype
    ),
    
    `N cells` = "",
    
    `Positive cells` = "",
    
    `% detected` =
      sprintf(
        "%.2f",
        percent_expressing
      ),
    
    `Statistical / descriptive result` =
      paste0(
        "Fisher FDR = ",
        fmt_p(FDR)
      )
  ) %>%
  
  select(
    Section,
    Gene,
    Context,
    `N cells`,
    `Positive cells`,
    `% detected`,
    `Statistical / descriptive result`
  )


S2 <- bind_rows(
  S2_detection,
  S2_cluster
)


# -----------------------------------------------------------------------------
# Further compacting:
#
# Broad-cell-type rows can still be numerous.
# For each gene × condition × time, keep the two cell types
# with the largest detection percentages.
# -----------------------------------------------------------------------------

S2_broad_compact <- S2_detection %>%
  
  mutate(
    pct_numeric =
      as.numeric(`% detected`)
  ) %>%
  
  separate(
    Context,
    into = c(
      "condition_tmp",
      "time_tmp",
      "celltype_tmp"
    ),
    sep = " \\| ",
    remove = FALSE
  ) %>%
  
  group_by(
    Gene,
    condition_tmp,
    time_tmp
  ) %>%
  
  slice_max(
    order_by = pct_numeric,
    n = 2,
    with_ties = FALSE
  ) %>%
  
  ungroup() %>%
  
  select(
    -pct_numeric,
    -condition_tmp,
    -time_tmp,
    -celltype_tmp
  )


S2 <- bind_rows(
  S2_broad_compact,
  S2_cluster
) %>%
  
  arrange(
    Gene,
    Section,
    Context
  )


write_csv(
  S2,
  file.path(
    TABLE_DIR,
    "Supplementary_Table_S2_Nobori_scRNA.csv"
  )
)


save_table_png(
  df = S2,
  
  filename =
    "Supplementary_Table_S2_Nobori_scRNA.png",
  
  title =
    "Supplementary Table S2. Key PDE1 and PR1 single-cell transcript-detection results",
  
  note =
    "Broad-cell-type values are transcript-detection frequencies; for compact presentation, the two cell types with the highest detection percentage per gene, condition and time point are shown. Major-cluster rows show Fisher exact-test results with Benjamini-Hochberg correction. Non-detection should not be interpreted as absence of expression."
)


# =============================================================================
# 6. TABLE S3 — LEE DEVELOPMENTAL + SPATIAL
# =============================================================================

lee_stage <- read_required(
  "Lee_stage_detection_summary.csv"
)

curio <- read_required(
  "Lee_Curio_cluster_summary.csv"
)

flower <- read_required(
  "Lee_Flower_PDE1_label_transfer_summary.csv"
)


# -----------------------------------------------------------------------------
# Part A — all 10 developmental stages
# -----------------------------------------------------------------------------

S3_stage_PDE1 <- lee_stage %>%
  
  transmute(
    Section =
      "Development",
    
    Context =
      as.character(stage),
    
    Gene =
      "PDE1",
    
    `N total` =
      as.character(n_nuclei),
    
    `N positive` =
      as.character(PDE1_positive_n),
    
    `% detected` =
      sprintf(
        "%.2f",
        PDE1_pct_detected
      ),
    
    `Enrichment / score` =
      ""
  )


S3_stage_PR1 <- lee_stage %>%
  
  transmute(
    Section =
      "Development",
    
    Context =
      as.character(stage),
    
    Gene =
      "PR1",
    
    `N total` =
      as.character(n_nuclei),
    
    `N positive` =
      as.character(PR1_positive_n),
    
    `% detected` =
      sprintf(
        "%.2f",
        PR1_pct_detected
      ),
    
    `Enrichment / score` =
      ""
  )


# -----------------------------------------------------------------------------
# Part B — Curio spatial clusters
# -----------------------------------------------------------------------------

S3_curio <- curio %>%
  
  filter(
    PDE1_positive_n > 0
  ) %>%
  
  transmute(
    Section =
      "Curio spatial",
    
    Context =
      paste0(
        "Cluster ",
        cluster
      ),
    
    Gene =
      "PDE1",
    
    `N total` =
      as.character(n_spots),
    
    `N positive` =
      as.character(PDE1_positive_n),
    
    `% detected` =
      sprintf(
        "%.3f",
        PDE1_pct
      ),
    
    `Enrichment / score` =
      paste0(
        "log2 enrichment ",
        sprintf(
          "%.2f",
          log2_PDE1_enrichment
        )
      )
  )


# -----------------------------------------------------------------------------
# Part C — Flower spatial label transfer
#
# Keep only PDE1-associated states, and only rows with high-confidence cells.
# -----------------------------------------------------------------------------

S3_flower <- flower %>%
  
  filter(
    PDE1_associated_state,
    n_high_conf_PDE1_associated > 0
  ) %>%
  
  arrange(
    sample,
    desc(n_high_conf_PDE1_associated),
    desc(mean_transfer_score)
  ) %>%
  
  group_by(sample) %>%
  
  slice_head(
    n = 5
  ) %>%
  
  ungroup() %>%
  
  transmute(
    Section =
      "Flower spatial inference",
    
    Context =
      paste0(
        sample,
        " | ",
        predicted_flower_celltype
      ),
    
    Gene =
      "PDE1-associated state",
    
    `N total` =
      as.character(n_cells),
    
    `N positive` =
      as.character(
        n_high_conf_PDE1_associated
      ),
    
    `% detected` =
      sprintf(
        "%.2f",
        pct_high_confidence
      ),
    
    `Enrichment / score` =
      paste0(
        "Mean transfer score ",
        sprintf(
          "%.2f",
          mean_transfer_score
        )
      )
  )


S3 <- bind_rows(
  S3_stage_PDE1,
  S3_stage_PR1,
  S3_curio,
  S3_flower
)


write_csv(
  S3,
  file.path(
    TABLE_DIR,
    "Supplementary_Table_S3_Lee_Developmental_Spatial.csv"
  )
)


save_table_png(
  df = S3,
  
  filename =
    "Supplementary_Table_S3_Lee_Developmental_Spatial.png",
  
  title =
    "Supplementary Table S3. Key developmental and spatial transcriptomic results",
  
  note =
    "Developmental values represent transcript-detection frequencies across the Lee et al. life-cycle atlas. Curio values represent direct transcript detection in spatial spots. Flower values represent high-confidence snRNA-to-MERFISH transferred PDE1-associated cell states; PDE1 was not directly targeted by the Flower MERFISH panel."
)


# =============================================================================
# 7. Final summary
# =============================================================================

message("\n============================================================")
message("SUPPLEMENTARY TABLES COMPLETE")
message("============================================================")

message("\nTable S1:")
message(
  file.path(
    TABLE_DIR,
    "Supplementary_Table_S1_Bulk.csv"
  )
)
message(
  file.path(
    TABLE_DIR,
    "Supplementary_Table_S1_Bulk.png"
  )
)

message("\nTable S2:")
message(
  file.path(
    TABLE_DIR,
    "Supplementary_Table_S2_Nobori_scRNA.csv"
  )
)
message(
  file.path(
    TABLE_DIR,
    "Supplementary_Table_S2_Nobori_scRNA.png"
  )
)

message("\nTable S3:")
message(
  file.path(
    TABLE_DIR,
    "Supplementary_Table_S3_Lee_Developmental_Spatial.csv"
  )
)
message(
  file.path(
    TABLE_DIR,
    "Supplementary_Table_S3_Lee_Developmental_Spatial.png"
  )
)

message("\nRows:")
message("S1 Bulk: ", nrow(S1))
message("S2 Nobori: ", nrow(S2))
message("S3 Lee: ", nrow(S3))

message("\nDONE.")