# =============================================================================
# 01_bulk_RNAseq_DESeq2_analysis.R
# Reproducible bulk RNA-seq analysis for PDE1 in Arabidopsis thaliana
#
# Inputs:
#   data/PRJEB25079_UncorrectedCounts.csv
#   data/E-MTAB-9694.sdrf.txt
#
# Formal inference:
#   wild-type Col only
#   group = stimulus x time
#   design = ~ group
#   each stimulus/time compared with matched-time mock
#   BH-adjusted P < 0.05
#
# Outputs:
#   5 figures
#   2 main result tables
#   sessionInfo + analysis parameters
# =============================================================================

rm(list = ls())
gc()
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(ggplot2)
})

COUNTS_FILE <- "data/raw/PRJEB25079_UncorrectedCounts.csv"
SDRF_FILE   <- "data/raw/E-MTAB-9694.sdrf.txt"

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

if (!file.exists(COUNTS_FILE)) stop("Missing: ", COUNTS_FILE)
if (!file.exists(SDRF_FILE)) stop("Missing: ", SDRF_FILE)

ALPHA <- 0.05
stimulus_order <- c("mock","flg22","elf18","nlp20","OGs","Pep1","chitooct","3.OH10")
time_order <- c(0,5,10,30,90,180)

# =============================================================================
# 1. Load raw counts and SDRF
# =============================================================================

counts_raw <- read.csv(
  COUNTS_FILE,
  row.names = 1,
  check.names = FALSE
)

sdrf <- read.delim(
  SDRF_FILE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

count_test <- as.matrix(counts_raw)
storage.mode(count_test) <- "numeric"

if (anyNA(count_test)) stop("Non-numeric/NA values found in raw counts.")
if (any(count_test < 0)) stop("Negative values found in raw counts.")
if (any(abs(count_test - round(count_test)) > 1e-8)) {
  stop("Input does not contain integer raw counts.")
}

cat("\nRaw count dimensions:\n")
print(dim(counts_raw))
cat("\nRaw count range:\n")
print(range(count_test))

# =============================================================================
# 2. Rebuild metadata exactly as in the earlier working script
# =============================================================================

sample_names <- colnames(counts_raw)
sample_split <- strsplit(sample_names, "_", fixed = TRUE)

metadata <- data.frame(
  row.names = sample_names,
  sample = sample_names,
  genotype = sapply(sample_split, `[`, 1),
  stimulus = sapply(sample_split, `[`, 2),
  time_raw = sapply(sample_split, `[`, 3),
  stringsAsFactors = FALSE
)

metadata$time_raw <- as.character(metadata$time_raw)

metadata$time <- dplyr::case_when(
  metadata$time_raw %in% c("000","0001","0002","0003") ~ "0",
  metadata$time_raw %in% c("005","0051","0052","0053") ~ "5",
  metadata$time_raw %in% c("010","0101","0102","0103") ~ "10",
  metadata$time_raw %in% c("030","0301","0302","0303") ~ "30",
  metadata$time_raw %in% c("090","0901","0902","0903") ~ "90",
  metadata$time_raw %in% c("180","1801","1802","1803") ~ "180",
  TRUE ~ NA_character_
)

metadata$replicate <- dplyr::case_when(
  metadata$time_raw %in% c("000","005","010","030","090","180") ~ "rep0",
  stringr::str_ends(metadata$time_raw, "1") ~ "rep1",
  stringr::str_ends(metadata$time_raw, "2") ~ "rep2",
  stringr::str_ends(metadata$time_raw, "3") ~ "rep3",
  TRUE ~ "rep_unknown"
)

metadata$genotype <- gsub("-", ".", metadata$genotype)
metadata$genotype <- gsub("/", ".", metadata$genotype)
metadata$stimulus <- gsub("-", ".", metadata$stimulus)
metadata$stimulus <- gsub("/", ".", metadata$stimulus)

metadata$time_numeric <- as.numeric(metadata$time)

if (anyNA(metadata$time_numeric)) {
  print(metadata[is.na(metadata$time_numeric), c("sample","time_raw")])
  stop("Some sample times could not be parsed.")
}

metadata$genotype <- factor(metadata$genotype)
metadata$stimulus <- factor(metadata$stimulus)
metadata$time <- factor(metadata$time, levels = as.character(time_order))
metadata$replicate <- factor(
  metadata$replicate,
  levels = c("rep0","rep1","rep2","rep3","rep_unknown")
)

metadata_summary <- metadata |>
  dplyr::group_by(genotype, stimulus, time) |>
  dplyr::summarise(n_samples = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(genotype, stimulus, time)

print(metadata_summary, n = Inf)

write.csv(
  metadata,
  "data/processed/analysis_metadata.csv",
  row.names = FALSE
)

# =============================================================================
# 3. Collapse transcript/version suffixes to gene-level IDs
# =============================================================================

count_matrix_tx <- as.matrix(counts_raw)
storage.mode(count_matrix_tx) <- "integer"

gene_ids <- toupper(sub("\\..*$", "", rownames(count_matrix_tx)))

count_matrix <- rowsum(
  count_matrix_tx,
  group = gene_ids,
  reorder = FALSE
)

storage.mode(count_matrix) <- "integer"

cat("\nGene-level count matrix dimensions:\n")
print(dim(count_matrix))
cat("\nGene-level total counts:\n")
print(sum(count_matrix))

write.csv(
  count_matrix,
  "data/processed/processed_gene_counts_matrix.csv"
)

if (!identical(colnames(count_matrix), metadata$sample)) {
  stop("Count matrix and metadata are not aligned.")
}

# =============================================================================
# 4. Restrict formal analysis to Col
# =============================================================================

meta_col <- metadata |>
  dplyr::filter(as.character(genotype) == "Col")

col_samples <- meta_col$sample
count_col <- count_matrix[, col_samples, drop = FALSE]

if (!identical(colnames(count_col), meta_col$sample)) {
  stop("Col counts and metadata are not aligned.")
}

cat("\nCol samples:", ncol(count_col), "\n")
cat("Col total counts:", sum(count_col), "\n")

col_group_summary <- meta_col |>
  dplyr::group_by(stimulus, time, time_numeric) |>
  dplyr::summarise(n_samples = dplyr::n(), .groups = "drop") |>
  dplyr::arrange(stimulus, time_numeric)

print(col_group_summary, n = Inf)

meta_col$group <- factor(
  paste(
    as.character(meta_col$stimulus),
    meta_col$time_numeric,
    sep = "_"
  )
)

# =============================================================================
# 5. DESeq2 model
# =============================================================================

dds_col <- DESeqDataSetFromMatrix(
  countData = count_col,
  colData = as.data.frame(meta_col),
  design = ~ group
)

# ==========================================================
# Filter low-count genes while retaining pre-specified candidates
# ==========================================================

candidate_ids <- c(
  "AT1G17330",  # PDE1
  "AT2G14610",  # PR1
  "AT5G40270",  # HD-domain-containing metal-dependent phosphohydrolase
  "AT5G40290",  # HD-domain-containing metal-dependent phosphohydrolase
  "AT1G14520",  # MIOX1
  "AT2G19800",  # MIOX2
  "AT5G56640",  # MIOX5
  "AT2G23820",  # metal-dependent phosphohydrolase
  "AT1G26160",  # metal-dependent phosphohydrolase
  "AT4G02260"   # RSH1
)

# Standard low-count filter, but retain all pre-specified candidate genes
keep <- rowSums(counts(dds_col)) >= 10 
cat("\nGenes before filtering:", nrow(dds_col), "\n")
cat("Genes after filtering:", sum(keep), "\n")
cat("Genes removed:", sum(!keep), "\n")
dds_col <- dds_col[keep, ]

# Confirm all candidate genes are retained
candidate_check <- data.frame(
  gene_id = candidate_ids,
  retained = candidate_ids %in% rownames(dds_col),
  stringsAsFactors = FALSE
)

candidate_check$total_Col_counts <- sapply(
  candidate_ids,
  function(g) {
    if (g %in% rownames(count_col)) {
      sum(count_col[g, ])
    } else {
      NA_real_
    }
  }
)

print(candidate_check)

cat("\nGenes retained after count filter:", nrow(dds_col), "\n")

dds_col <- DESeq(dds_col)

saveRDS(
  dds_col,
  "data/processed/dds_Col_stimulus_time.rds"
)

normalised_counts <- counts(dds_col, normalized = TRUE)

write.csv(
  normalised_counts,
  "data/processed/deseq2_Col_normalised_gene_counts.csv"
)

# =============================================================================
# 6. Candidate genes
# =============================================================================
candidate_genes <- tibble::tribble(
  ~gene_id,      ~gene_label,   ~gene_group,
  "AT1G17330",   "PDE1",        "cGMP-stimulated phosphodiesterase",
  "AT2G14610",   "PR1",         "Pathogenesis-related protein 1",
  "AT5G40270",   "AT5G40270",   "HD-domain-containing metal-dependent phosphohydrolase",
  "AT5G40290",   "AT5G40290",   "HD-domain-containing metal-dependent phosphohydrolase",
  "AT1G14520",   "AT1G14520",   "Myo-inositol oxygenase 1 (MIOX1)",
  "AT2G19800",   "AT2G19800",   "Myo-inositol oxygenase 2 (MIOX2)",
  "AT5G56640",   "AT5G56640",   "Myo-inositol oxygenase 5 (MIOX5)",
  "AT2G23820",   "AT2G23820",   "Metal-dependent phosphohydrolase",
  "AT1G26160",   "AT1G26160",   "Metal-dependent phosphohydrolase",
  "AT4G02260",   "AT4G02260",   "RelA/SpoT homolog 1 (RSH1)"
)

candidate_found <- candidate_genes |>
  dplyr::filter(gene_id %in% rownames(dds_col))

if (nrow(candidate_found) == 0) {
  stop("No candidate genes found after filtering.")
}

gene_order <- c(
  "PDE1",
  "PR1",
  "AT5G40270",
  "AT5G40290",
  "AT1G14520",
  "AT2G19800",
  "AT5G56640",
  "AT2G23820",
  "AT1G26160",
  "AT4G02260"
)
# =============================================================================
# 7. Matched-time contrasts
# =============================================================================

available_groups <- levels(meta_col$group)
nonmock_stimuli <- setdiff(unique(as.character(meta_col$stimulus)), "mock")
available_times <- sort(unique(meta_col$time_numeric))

contrast_grid <- tidyr::expand_grid(
  stimulus = nonmock_stimuli,
  time = available_times
) |>
  dplyr::mutate(
    treatment_group = paste(stimulus, time, sep = "_"),
    mock_group = paste("mock", time, sep = "_"),
    valid = treatment_group %in% available_groups &
      mock_group %in% available_groups
  ) |>
  dplyr::filter(valid) |>
  dplyr::select(-valid)

extract_one_contrast <- function(stimulus, time, treatment_group, mock_group) {
  
  res <- DESeq2::results(
    dds_col,
    contrast = c("group", treatment_group, mock_group),
    alpha = ALPHA
  )
  
  res_df <- as.data.frame(res) |>
    tibble::rownames_to_column("gene_id")
  
  candidate_found |>
    dplyr::select(gene_id, gene_label, gene_group) |>
    dplyr::left_join(res_df, by = "gene_id") |>
    dplyr::transmute(
      gene_id,
      gene_label,
      gene_group,
      stimulus = stimulus,
      time = as.numeric(time),
      comparison = paste0(treatment_group, "_vs_", mock_group),
      baseMean,
      log2FC = log2FoldChange,
      fold_change = 2^log2FoldChange,
      lfcSE,
      stat,
      pvalue,
      padj,
      ci95_low_log2FC = log2FoldChange - 1.96 * lfcSE,
      ci95_high_log2FC = log2FoldChange + 1.96 * lfcSE,
      significant_FDR05 = !is.na(padj) & padj < ALPHA
    )
}

candidate_results <- purrr::pmap_dfr(
  contrast_grid,
  extract_one_contrast
) |>
  dplyr::mutate(
    stimulus = factor(stimulus, levels = stimulus_order),
    gene_label = factor(gene_label, levels = gene_order)
  ) |>
  dplyr::arrange(gene_label, stimulus, time)

write.csv(
  candidate_results,
  "results/tables/bulk_DESeq2_candidate_results.csv",
  row.names = FALSE
)

# =============================================================================
# 8. Peak summary
# =============================================================================

peak_summary <- candidate_results |>
  dplyr::filter(!is.na(log2FC)) |>
  dplyr::group_by(gene_id, gene_label, gene_group) |>
  dplyr::arrange(dplyr::desc(log2FC), .by_group = TRUE) |>
  dplyr::slice(1) |>
  dplyr::ungroup() |>
  dplyr::arrange(dplyr::desc(log2FC))

write.csv(
  peak_summary,
  "results/tables/bulk_candidate_peak_summary.csv",
  row.names = FALSE
)

# =============================================================================
# 9. Normalised-expression table
# =============================================================================

expression_list <- vector("list", nrow(candidate_found))

for (i in seq_len(nrow(candidate_found))) {
  gene <- candidate_found$gene_id[i]
  
  expression_list[[i]] <- data.frame(
    gene_id = gene,
    gene_label = candidate_found$gene_label[i],
    sample = colnames(normalised_counts),
    normalized_count = as.numeric(normalised_counts[gene, ]),
    stringsAsFactors = FALSE
  )
}

expression_long <- dplyr::bind_rows(expression_list) |>
  dplyr::left_join(
    meta_col |>
      dplyr::select(sample, stimulus, time_numeric, replicate),
    by = "sample"
  ) |>
  dplyr::mutate(
    log2_normalized = log2(normalized_count + 1),
    stimulus = factor(stimulus, levels = stimulus_order),
    gene_label = factor(gene_label, levels = gene_order)
  )

expression_summary <- expression_long |>
  dplyr::group_by(gene_id, gene_label, stimulus, time_numeric) |>
  dplyr::summarise(
    n_samples = dplyr::n(),
    mean_log2_normalized = mean(log2_normalized, na.rm = TRUE),
    sd_log2_normalized = sd(log2_normalized, na.rm = TRUE),
    se_log2_normalized = sd_log2_normalized / sqrt(n_samples),
    .groups = "drop"
  )

# =============================================================================
# 10. Plot theme
# =============================================================================

theme_final <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.25, colour = "grey88"),
      panel.border = element_rect(colour = "black", linewidth = 0.45),
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = base_size + 3),
      strip.background = element_rect(fill = "grey95", colour = "black"),
      strip.text = element_text(face = "bold"),
      legend.title = element_text(face = "bold")
    )
}

# =============================================================================
# 11. Figure 1: PDE1 normalised expression + FDR labels
# =============================================================================

pde1_expr <- expression_long |>
  dplyr::filter(gene_id == "AT1G17330")

pde1_summary <- expression_summary |>
  dplyr::filter(gene_id == "AT1G17330")

pde1_stats <- candidate_results |>
  dplyr::filter(gene_id == "AT1G17330") |>
  dplyr::mutate(
    sig_label = dplyr::case_when(
      is.na(padj) ~ "",
      padj < 0.001 ~ "***",
      padj < 0.01 ~ "**",
      padj < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) |>
  dplyr::filter(sig_label != "") |>
  dplyr::left_join(
    pde1_summary |>
      dplyr::select(
        stimulus,
        time_numeric,
        mean_log2_normalized,
        se_log2_normalized
      ),
    by = c("stimulus", "time" = "time_numeric")
  ) |>
  dplyr::mutate(
    y_position = mean_log2_normalized +
      dplyr::coalesce(se_log2_normalized, 0) + 0.15
  )

fig1 <- ggplot() +
  geom_point(
    data = pde1_expr,
    aes(
      x = time_numeric,
      y = log2_normalized,
      colour = stimulus,
      shape = replicate
    ),
    position = position_jitter(width = 1.6, height = 0),
    alpha = 0.55,
    size = 2
  ) +
  geom_line(
    data = pde1_summary,
    aes(
      x = time_numeric,
      y = mean_log2_normalized,
      colour = stimulus,
      group = stimulus
    ),
    linewidth = 1.1
  ) +
  geom_point(
    data = pde1_summary,
    aes(
      x = time_numeric,
      y = mean_log2_normalized,
      colour = stimulus
    ),
    size = 2.8
  ) +
  geom_errorbar(
    data = pde1_summary,
    aes(
      x = time_numeric,
      ymin = mean_log2_normalized - se_log2_normalized,
      ymax = mean_log2_normalized + se_log2_normalized,
      colour = stimulus
    ),
    width = 2.5,
    linewidth = 0.35,
    na.rm = TRUE
  ) +
  geom_text(
    data = pde1_stats,
    aes(
      x = time,
      y = y_position,
      label = sig_label,
      colour = stimulus
    ),
    show.legend = FALSE,
    fontface = "bold",
    size = 4
  ) +
  scale_x_continuous(breaks = time_order) +
  theme_final(12) +
  labs(
    title = "PDE1 expression across immune stimuli",
    subtitle = paste0(
      "Wild-type Col; mean ± SE; ",
      "* FDR<0.05, ** FDR<0.01, *** FDR<0.001 vs matched-time mock"
    ),
    x = "Time after treatment (min)",
    y = "PDE1 log2(DESeq2-normalised count + 1)",
    colour = "Stimulus",
    shape = "Replicate"
  )

ggsave(
  "results/figures/Fig1_PDE1_expression_with_statistics.png",
  fig1,
  width = 14,
  height = 9,
  dpi = 600
)

# =============================================================================
# 12. Figure 2: candidate normalised-expression overview
# =============================================================================

fig2 <- ggplot(
  expression_summary,
  aes(
    x = time_numeric,
    y = mean_log2_normalized,
    colour = gene_label,
    group = gene_label
  )
) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = mean_log2_normalized - se_log2_normalized,
      ymax = mean_log2_normalized + se_log2_normalized
    ),
    width = 2.5,
    linewidth = 0.25,
    alpha = 0.6,
    na.rm = TRUE
  ) +
  facet_wrap(~ stimulus, ncol = 4, drop = FALSE) +
  scale_x_continuous(breaks = time_order) +
  theme_final(11) +
  labs(
    title = "Expression of PDE candidate genes and PR1",
    subtitle = "Wild-type Col; mean log2(DESeq2-normalised count + 1) ± SE",
    x = "Time after treatment (min)",
    y = "Mean log2 normalised expression",
    colour = "Gene"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "results/figures/Fig2_candidate_expression_overview.png",
  fig2,
  width = 15,
  height = 9,
  dpi = 600
)

# =============================================================================
# 13. Figure 3: formal DESeq2 log2FC vs matched-time mock
# =============================================================================

fig3_data <- candidate_results |>
  dplyr::mutate(
    statistical_support = dplyr::case_when(
      is.na(padj) ~ "FDR unavailable",
      padj < ALPHA ~ "FDR < 0.05",
      TRUE ~ "FDR ≥ 0.05"
    )
  )

fig3 <- ggplot(
  fig3_data,
  aes(
    x = time,
    y = log2FC,
    colour = gene_label,
    group = gene_label
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.8, alpha = 0.85) +
  geom_errorbar(
    aes(
      ymin = ci95_low_log2FC,
      ymax = ci95_high_log2FC
    ),
    width = 2.5,
    linewidth = 0.25,
    alpha = 0.55,
    na.rm = TRUE
  ) +
  geom_point(
    aes(shape = statistical_support),
    size = 2.2
  ) +
  facet_wrap(~ stimulus, ncol = 4) +
  scale_x_continuous(breaks = time_order) +
  theme_final(11) +
  labs(
    title = "DESeq2 treatment effects for PDE candidate genes and PR1",
    subtitle = "Treatment vs matched-time mock; error bars = 95% Wald CI",
    x = "Time after treatment (min)",
    y = "DESeq2 log2 fold change vs mock",
    colour = "Gene",
    shape = "Statistical support"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "results/figures/Fig3_DESeq2_log2FC_vs_mock.png",
  fig3,
  width = 15,
  height = 9,
  dpi = 600
)

# =============================================================================
# 14. Figure 4: DESeq2 log2FC heatmap
# =============================================================================

heatmap_data <- candidate_results |>
  dplyr::mutate(
    time_factor = factor(time, levels = time_order),
    
    # Set required top-to-bottom gene order
    gene_label = factor(
      as.character(gene_label),
      levels = rev(gene_order)
    ),
    
    # Mark statistically significant results
    sig_label = ifelse(significant_FDR05, "*", "")
  )
heatmap_limit <- max(abs(heatmap_data$log2FC), na.rm = TRUE)

fig4 <- ggplot(
  heatmap_data,
  aes(
    x = time_factor,
    y = gene_label,
    fill = log2FC
  )
) +
  geom_tile(colour = "grey80", linewidth = 0.3) +
  geom_text(aes(label = sig_label), fontface = "bold", size = 3.6) +
  facet_wrap(~ stimulus, nrow = 1) +
  scale_fill_gradient2(
    low = "#8c2d2d",
    mid = "grey90",
    high = "#5e4fa2",
    midpoint = 0,
    limits = c(-heatmap_limit, heatmap_limit),
    name = "DESeq2\nlog2FC"
  ) +
  theme_final(10) +
  labs(
    title = "DESeq2 log2 fold-change heatmap of PDE candidate genes and PR1",
    subtitle = "* BH-adjusted P < 0.05; treatment vs matched-time mock",
    x = "Time after treatment (min)",
    y = "Gene"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  "results/figures/Fig4_DESeq2_log2FC_heatmap.png",
  fig4,
  width = 18,
  height = 8,
  dpi = 600
)

# =============================================================================
# 15. Figure 5: ranked maximum estimated DESeq2 response
# =============================================================================

fig5_data <- peak_summary |>
  dplyr::mutate(
    gene_label = factor(
      as.character(gene_label),
      levels = rev(as.character(gene_label))
    ),
    statistical_support = dplyr::case_when(
      is.na(padj) ~ "FDR unavailable",
      padj < ALPHA ~ "FDR < 0.05",
      TRUE ~ "FDR ≥ 0.05"
    ),
    peak_label = paste0(
      as.character(stimulus),
      ", ",
      time,
      " min\nFC=",
      sprintf("%.2f", fold_change),
      "; FDR=",
      ifelse(
        is.na(padj),
        "NA",
        formatC(padj, format = "g", digits = 2)
      )
    )
  )

fig5 <- ggplot(
  fig5_data,
  aes(
    x = gene_label,
    y = log2FC,
    fill = statistical_support
  )
) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.25) +
  geom_errorbar(
    aes(
      ymin = ci95_low_log2FC,
      ymax = ci95_high_log2FC
    ),
    width = 0.22,
    linewidth = 0.35,
    na.rm = TRUE
  ) +
  geom_text(
    aes(label = peak_label),
    hjust = -0.03,
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      "FDR < 0.05" = "#00bfc4",       # blue = significant
      "FDR ≥ 0.05" = "#f8766d"       # orange = not significant
    )
  ) +
  coord_flip() +
  theme_final(11) +
  labs(
    title = "Maximum estimated DESeq2 response of candidate genes",
    subtitle = "Descriptive ranking; each selected contrast retains its DESeq2 FDR",
    x = "Gene",
    y = "Maximum DESeq2 log2FC vs matched-time mock",
    fill = "Statistical support"
  ) +
  expand_limits(
    y = max(fig5_data$ci95_high_log2FC, na.rm = TRUE) + 1.6
  )

ggsave(
  "results/figures/Fig5_ranked_max_DESeq2_effect.png",
  fig5,
  width = 12,
  height = 8,
  dpi = 600
)

# =============================================================================
# 16. Print key PDE1 and PR1 results
# =============================================================================

cat("\nPDE1 results:\n")
print(
  candidate_results |>
    dplyr::filter(gene_id == "AT1G17330") |>
    dplyr::select(
      stimulus,
      time,
      log2FC,
      fold_change,
      pvalue,
      padj,
      significant_FDR05
    ),
  n = Inf
)

cat("\nPR1 results:\n")
print(
  candidate_results |>
    dplyr::filter(gene_id == "AT2G14610") |>
    dplyr::select(
      stimulus,
      time,
      log2FC,
      fold_change,
      pvalue,
      padj,
      significant_FDR05
    ),
  n = Inf
)

# =============================================================================
# 17. Reproducibility information
# =============================================================================

capture.output(
  sessionInfo(),
  file = "results/sessionInfo.txt"
)

writeLines(
  c(
    paste0("Raw count file: ", COUNTS_FILE),
    paste0("SDRF file: ", SDRF_FILE),
    "Transcript/version suffixes collapsed to base Arabidopsis gene IDs by summing counts.",
    "Formal inference restricted to wild-type Col samples.",
    "DESeq2 design: ~ group, where group = stimulus x time.",
    "Each non-mock group compared with matched-time mock.",
    "Low-count filter: total count >= 10 across Col samples.",
    paste0("FDR threshold: ", ALPHA),
    "Multiple testing: Benjamini-Hochberg adjustment from DESeq2.",
    "Figures 1-2: DESeq2-normalised count summaries.",
    "Figures 3-5: formal DESeq2 model estimates.",
    paste0("R version: ", R.version.string),
    paste0("DESeq2 version: ", packageVersion("DESeq2")),
    paste0("ggplot2 version: ", packageVersion("ggplot2"))
  ),
  con = "results/analysis_parameters.txt"
)

cat("\nBulk RNA-seq analysis completed successfully.\n")
