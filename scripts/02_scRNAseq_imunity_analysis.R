# ============================================================
# final_generate_scRNA_figures.R
#
# Purpose:
# Generate final scRNA figures for PDE1/PR1 analysis from
# Nobori et al. GSE226826 combined filtered Seurat object.
#
# Figures generated:
#
# 1. final_guard_response/
#    - fig2_celltypes_umap.png
#
# 2. final_umap_15_panel/
#    - fig5_PDE1_detected_cells_15panel.png
#    - fig6_PR1_detected_cells_15panel.png
#
# 3. final_scRNA_candidate_dotplots/
#    - fig3_candidate_gene_dotplot_major_clusters_celltype_yaxis.png
#
# 4. results/figures/new/nobori_publication/
#    - Figure7C_PDE1_PR1_condition_time_celltype_trend.png
#
# 5. results/figures/marker_annotation_epidermis_mesophyll/
#    - Epidermis_marker_dotplot.png
#    - Mesophyll_marker_dotplot.png
#
# 6. results/figures/subcluster_binary_8figures/
#    - PDE1_AT1G17330_Epidermis_binary_subcluster_panels.png
#    - PDE1_AT1G17330_Mesophyll_binary_subcluster_panels.png
#    - PR1_AT2G14610_Epidermis_binary_subcluster_panels.png
#    - PR1_AT2G14610_Mesophyll_binary_subcluster_panels.png
#
# Notes:
# - No tables are produced.
# - Downsampling is used for plotting only.
# - Quantitative trend figure uses sample-level summaries internally.
# - Seurat v5-compatible: uses layer = "data".
# ============================================================

rm(list = ls())

# ============================================================
# 1. Libraries
# ============================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(stringr)
library(patchwork)
library(scales)

set.seed(123)

# ============================================================
# 2. Paths
# ============================================================

base_dir <- "C:/Users/qo25519/Documents/PDE1_bulkRNA_analysis/PDE_candidate_bulkRNA/PDE1-expression-analysis-in-Arabidopsis"

raw_dir <- file.path(base_dir, "data", "raw")
rds_path <- file.path(raw_dir, "GSE226826_combined_filtered.rds")

# Single output directory for ALL figures
out_dir <- file.path(base_dir, "results", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Use same output directory for every figure category
out_fig2       <- out_dir
out_15panel    <- out_dir
out_dotplot    <- out_dir
out_trend      <- out_dir
out_marker     <- out_dir
out_subcluster <- out_dir
# ============================================================
# 3. Load data
# ============================================================

obj <- readRDS(rds_path)
DefaultAssay(obj) <- "SCT"

# ============================================================
# 4. Basic checks and helper variables
# ============================================================

if (!"celltype" %in% colnames(obj@meta.data)) {
  stop("Metadata column 'celltype' not found.")
}

if (!"sample" %in% colnames(obj@meta.data)) {
  stop("Metadata column 'sample' not found.")
}

if (!"time" %in% colnames(obj@meta.data)) {
  stop("Metadata column 'time' not found.")
}

if (!"SCT_snn_res.1" %in% colnames(obj@meta.data)) {
  stop("Metadata column 'SCT_snn_res.1' not found.")
}

reduction_use <- "umap"

message("Using reduction: ", reduction_use)

# Add parsed condition/time labels
obj$condition <- case_when(
  grepl("^00_Mock", obj$sample) ~ "Mock",
  grepl("^DC3000", obj$sample) ~ "DC3000",
  grepl("^AvrRpm1", obj$sample) ~ "AvrRpm1",
  grepl("^AvrRpt2", obj$sample) ~ "AvrRpt2",
  TRUE ~ NA_character_
)

obj$time_plot <- case_when(
  grepl("^00_Mock", obj$sample) ~ "Mock",
  grepl("04h", obj$sample) ~ "04h",
  grepl("06h", obj$sample) ~ "06h",
  grepl("09h", obj$sample) ~ "09h",
  grepl("24h", obj$sample) ~ "24h",
  TRUE ~ NA_character_
)

obj$condition <- factor(
  obj$condition,
  levels = c("Mock", "DC3000", "AvrRpm1", "AvrRpt2")
)

obj$time_plot <- factor(
  obj$time_plot,
  levels = c("Mock", "04h", "06h", "09h", "24h")
)

obj$condition_time <- case_when(
  obj$condition == "Mock" ~ "Mock",
  TRUE ~ paste(as.character(obj$condition), as.character(obj$time_plot), sep = " ")
)

obj$condition_time <- factor(
  obj$condition_time,
  levels = c(
    "Mock",
    "DC3000 04h", "DC3000 06h", "DC3000 09h", "DC3000 24h",
    "AvrRpm1 04h", "AvrRpm1 06h", "AvrRpm1 09h", "AvrRpm1 24h",
    "AvrRpt2 04h", "AvrRpt2 06h", "AvrRpt2 09h", "AvrRpt2 24h"
  )
)

# ============================================================
# 5. Helper functions
# ============================================================

get_expr_vec <- function(seu, gene_id, assay_name = "SCT") {
  DefaultAssay(seu) <- assay_name
  
  if (!(gene_id %in% rownames(seu[[assay_name]]))) {
    stop("Gene not found in assay: ", gene_id)
  }
  
  expr <- GetAssayData(seu, assay = assay_name, layer = "data")[gene_id, ]
  expr <- as.numeric(expr)
  names(expr) <- colnames(seu)
  
  return(expr)
}

get_umap_df <- function(seu, reduction_name) {
  emb <- Embeddings(seu, reduction_name) %>%
    as.data.frame() %>%
    rownames_to_column("cell")
  
  colnames(emb)[2:3] <- c("UMAP_1", "UMAP_2")
  return(emb)
}

downsample_by_group <- function(df, group_col, n_max) {
  df %>%
    group_by(.data[[group_col]]) %>%
    group_modify(~ {
      if (nrow(.x) > n_max) {
        slice_sample(.x, n = n_max)
      } else {
        .x
      }
    }) %>%
    ungroup()
}

make_15panel_df <- function(df) {
  
  mock_df <- df %>%
    filter(condition == "Mock", time_plot == "Mock")
  
  mock_repeated <- bind_rows(
    mock_df %>% mutate(panel_row = "DC3000"),
    mock_df %>% mutate(panel_row = "AvrRpm1"),
    mock_df %>% mutate(panel_row = "AvrRpt2")
  )
  
  nonmock_df <- df %>%
    filter(condition != "Mock", time_plot != "Mock") %>%
    mutate(panel_row = as.character(condition))
  
  out <- bind_rows(mock_repeated, nonmock_df)
  
  out$panel_row <- factor(out$panel_row, levels = c("DC3000", "AvrRpm1", "AvrRpt2"))
  out$time_plot <- factor(out$time_plot, levels = c("Mock", "04h", "06h", "09h", "24h"))
  
  return(out)
}

# ============================================================
# 6. Figure 2: Broad cell-type UMAP
# ============================================================

umap_df <- get_umap_df(obj, reduction_use)

meta_df <- obj@meta.data %>%
  as.data.frame() %>%
  rownames_to_column("cell")

fig2_df <- left_join(umap_df, meta_df, by = "cell")

celltype_cols <- c(
  "Epidermis" = "#F8766D",
  "Mesophyll" = "#7CAE00",
  "Unknown" = "#00BFC4",
  "Vasculature" = "#C77CFF"
)

p_fig2 <- ggplot(fig2_df, aes(x = UMAP_1, y = UMAP_2, colour = celltype)) +
  geom_point(size = 0.25, alpha = 0.85) +
  scale_colour_manual(values = celltype_cols, name = "Cell type") +
  labs(
    title = "Broad cell-type annotation from metadata",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

ggsave(
  filename = file.path(out_fig2, "Fig5_celltypes_umap.png"),
  plot = p_fig2,
  width = 7,
  height = 6,
  dpi = 300,
  bg = "white"
)

# ============================================================
# 7. Figures 5 and 6: Whole-object PDE1/PR1 15-panel UMAPs
# ============================================================

make_whole_umap_15panel <- function(seu, gene_id, gene_label, out_file) {
  
  expr_vec <- get_expr_vec(seu, gene_id, assay_name = "SCT")
  umap_df <- get_umap_df(seu, reduction_use)
  
  meta_df <- seu@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell")
  
  df <- left_join(umap_df, meta_df, by = "cell") %>%
    mutate(
      expr = expr_vec[cell],
      detected = expr > 0
    ) %>%
    filter(!is.na(condition), !is.na(time_plot))
  
  # Plotting downsample only for very large groups
  df$group_id <- ifelse(df$condition == "Mock", "Mock", as.character(df$condition_time))
  
  df_ds <- df %>%
    group_by(group_id) %>%
    group_modify(~ {
      if (.y$group_id == "Mock" && nrow(.x) > 3500) {
        slice_sample(.x, n = 3500)
      } else if (.y$group_id == "AvrRpt2 09h" && nrow(.x) > 3500) {
        slice_sample(.x, n = 3500)
      } else {
        .x
      }
    }) %>%
    ungroup()
  
  plot_df <- make_15panel_df(df_ds)
  
  p <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(
      colour = "grey88",
      size = 0.22,
      alpha = 0.65
    ) +
    geom_point(
      data = plot_df %>% filter(detected),
      colour = "#8B0000",
      size = 0.35,
      alpha = 0.95
    ) +
    facet_grid(panel_row ~ time_plot, drop = FALSE) +
    labs(
      title = paste0(gene_label, " detected cells"),
      subtitle = "Mock and pooled AvrRpt2 09h are downsampled",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_rect(fill = "grey90", colour = "black"),
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(0.45, "lines")
    )
  
  ggsave(
    filename = out_file,
    plot = p,
    width = 18,
    height = 11,
    dpi = 300,
    bg = "white"
  )
}

make_whole_umap_15panel(
  seu = obj,
  gene_id = "AT1G17330",
  gene_label = "PDE1",
  out_file = file.path(out_15panel, "Fig6_PDE1_detected_cells_15panel.png")
)

make_whole_umap_15panel(
  seu = obj,
  gene_id = "PR1",
  gene_label = "PR1",
  out_file = file.path(out_15panel, "Fig7_PR1_detected_cells_15panel.png")
)
# ============================================================
# 8. Fig. 8: Candidate gene dotplot across major clusters
# with exploratory Fisher exact test significance
# ============================================================

candidate_genes <- data.frame(
  gene_id = c(
    "AT1G17330",  # PDE1
    "PR1",        # PR1
    "MIOX1",      # AT1G14520
    "MIOX2",      # AT2G19800
    "MIOX5",      # AT5G56640
    "AT5G40270",  # 270
    "AT5G40290",  # 290
    "RSH1",       # AT4G02260
    "AT1G26160",  # 160
    "AT2G23820"   # 820
  ),
  gene_label = c(
    "PDE1",
    "PR1",
    "MIOX1",
    "MIOX2",
    "MIOX5",
    "AT5G40270",
    "AT5G40290",
    "RSH1",
    "AT1G26160",
    "AT2G23820"
  ),
  stringsAsFactors = FALSE
)

candidate_genes <- candidate_genes %>%
  filter(gene_id %in% rownames(obj[["SCT"]]))

Idents(obj) <- "SCT_snn_res.1"

# -----------------------------
# DotPlot data
# -----------------------------
dot_data <- DotPlot(
  obj,
  features = candidate_genes$gene_id,
  assay = "SCT"
)$data

dot_data <- dot_data %>%
  left_join(candidate_genes, by = c("features.plot" = "gene_id")) %>%
  mutate(
    id = as.character(id),
    id_num = as.numeric(id),
    gene_label = factor(gene_label, levels = candidate_genes$gene_label)
  )

# -----------------------------
# Fisher exact test:
# gene detected in cluster vs outside cluster
# -----------------------------
expr_mat <- GetAssayData(obj, assay = "SCT", layer = "data")

stats_list <- list()

for (i in seq_len(nrow(candidate_genes))) {
  
  gene_id <- candidate_genes$gene_id[i]
  gene_label <- candidate_genes$gene_label[i]
  
  detected <- as.numeric(expr_mat[gene_id, colnames(obj)] > 0)
  names(detected) <- colnames(obj)
  
  for (cl in sort(unique(as.character(obj$SCT_snn_res.1)))) {
    
    in_cluster <- as.character(obj$SCT_snn_res.1) == cl
    
    detected_in <- sum(detected[in_cluster] == 1, na.rm = TRUE)
    not_detected_in <- sum(detected[in_cluster] == 0, na.rm = TRUE)
    
    detected_out <- sum(detected[!in_cluster] == 1, na.rm = TRUE)
    not_detected_out <- sum(detected[!in_cluster] == 0, na.rm = TRUE)
    
    fisher_table <- matrix(
      c(
        detected_in,
        not_detected_in,
        detected_out,
        not_detected_out
      ),
      nrow = 2,
      byrow = TRUE
    )
    
    fisher_res <- fisher.test(fisher_table)
    
    stats_list[[paste(gene_label, cl, sep = "_")]] <- data.frame(
      gene_label = gene_label,
      id = cl,
      p_value = fisher_res$p.value,
      stringsAsFactors = FALSE
    )
  }
}

stats_df <- bind_rows(stats_list) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    significance = case_when(
      p_adj < 0.001 ~ "***",
      p_adj < 0.01  ~ "**",
      p_adj < 0.05  ~ "*",
      TRUE ~ ""
    ),
    id_num = as.numeric(id),
    gene_label = factor(gene_label, levels = candidate_genes$gene_label)
  )

# Add significance to dotplot data
dot_data <- dot_data %>%
  left_join(
    stats_df %>% select(gene_label, id_num, p_adj, significance),
    by = c("gene_label", "id_num")
  )

# -----------------------------
# Dominant cell type for each major cluster
# -----------------------------
dominant_celltype <- obj@meta.data %>%
  as.data.frame() %>%
  count(SCT_snn_res.1, celltype, name = "n") %>%
  group_by(SCT_snn_res.1) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    id = as.character(SCT_snn_res.1),
    id_num = as.numeric(id)
  )

cluster_levels <- sort(unique(dot_data$id_num))

dot_data$id_factor <- factor(dot_data$id_num, levels = cluster_levels)
dominant_celltype$id_factor <- factor(dominant_celltype$id_num, levels = cluster_levels)

# -----------------------------
# Left cell-type strip
# -----------------------------
p_strip <- ggplot(dominant_celltype, aes(x = 1, y = id_factor, fill = celltype)) +
  geom_tile(width = 0.7, height = 0.9) +
  scale_fill_manual(values = celltype_cols, name = "Dominant\ncell type") +
  scale_y_discrete(position = "left") +
  labs(x = NULL, y = "Major cluster") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"
  )

# -----------------------------
# Dotplot with FDR significance stars
# -----------------------------
p_dot <- ggplot(dot_data, aes(x = gene_label, y = id_factor)) +
  geom_point(aes(size = pct.exp, colour = avg.exp.scaled), alpha = 0.9) +
  geom_text(
    aes(label = significance),
    size = 3.2,
    colour = "black",
    vjust = 0.35
  ) +
  scale_colour_gradient(
    low = "grey85",
    high = "blue",
    name = "Average\nexpression"
  ) +
  scale_size(
    range = c(0.2, 8),
    name = "% expressing"
  ) +
  labs(
    title = "Candidate gene expression and cluster enrichment",
    subtitle = "Dot size = % expressing; colour = average expression; stars = Fisher exact test FDR significance",
    x = "Candidate genes",
    y = NULL
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

p_fig8 <- p_strip + p_dot + plot_layout(widths = c(0.7, 9))

ggsave(
  filename = file.path(out_dotplot, "Fig8_candidate_gene_dotplot_major_clusters_celltype_yaxis.png"),
  plot = p_fig8,
  width = 16,
  height = 11,
  dpi = 300,
  bg = "white"
)

# ============================================================
# 9. Figure 7C: Sample-level PDE1/PR1 condition-time trend
# ============================================================

make_sample_level_trend <- function(seu) {
  
  genes_for_trend <- data.frame(
    gene_id = c("AT1G17330", "PR1"),
    gene_label = c("PDE1_AT1G17330", "PR1_AT2G14610"),
    stringsAsFactors = FALSE
  )
  
  meta_df <- seu@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell") %>%
    filter(condition != "Mock") %>%
    filter(time_plot %in% c("04h", "06h", "09h", "24h")) %>%
    filter(celltype %in% c("Epidermis", "Mesophyll", "Vasculature", "Unknown"))
  
  all_gene_df <- list()
  
  for (i in seq_len(nrow(genes_for_trend))) {
    
    gene_id <- genes_for_trend$gene_id[i]
    gene_label <- genes_for_trend$gene_label[i]
    
    expr_vec <- get_expr_vec(seu, gene_id, assay_name = "SCT")
    
    gene_df <- meta_df %>%
      mutate(
        gene = gene_label,
        expr = expr_vec[cell],
        detected = expr > 0
      ) %>%
      group_by(gene, sample, condition, time_plot, celltype) %>%
      summarise(
        percent_positive = 100 * mean(detected, na.rm = TRUE),
        .groups = "drop"
      )
    
    all_gene_df[[gene_label]] <- gene_df
  }
  
  trend_df <- bind_rows(all_gene_df) %>%
    group_by(gene, condition, time_plot, celltype) %>%
    summarise(
      mean_percent_positive = mean(percent_positive, na.rm = TRUE),
      .groups = "drop"
    )
  
  trend_df$time_plot <- factor(trend_df$time_plot, levels = c("04h", "06h", "09h", "24h"))
  trend_df$condition <- factor(trend_df$condition, levels = c("DC3000", "AvrRpm1", "AvrRpt2"))
  trend_df$celltype <- factor(trend_df$celltype, levels = c("Epidermis", "Mesophyll", "Vasculature", "Unknown"))
  trend_df$gene <- factor(trend_df$gene, levels = c("PDE1_AT1G17330", "PR1_AT2G14610"))
  
  p <- ggplot(
    trend_df,
    aes(
      x = time_plot,
      y = mean_percent_positive,
      group = condition,
      colour = condition
    )
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2.4) +
    facet_grid(gene ~ celltype, scales = "free_y") +
    labs(
      title = "PDE1 and PR1 detection over infection time",
      subtitle = "Sample-level mean % positive cells by infection condition and cell type",
      x = "Time after infection",
      y = "Mean % positive cells",
      colour = "Condition"
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 18, hjust = 0),
      plot.subtitle = element_text(size = 11, hjust = 0),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.background = element_rect(fill = "grey90", colour = "black"),
      strip.text = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  return(p)
}

p_trend <- make_sample_level_trend(obj)

ggsave(
  filename = file.path(out_trend, "Fig9_PDE1_PR1_condition_time_celltype_trend.png"),
  plot = p_trend,
  width = 16,
  height = 9,
  dpi = 300,
  bg = "white"
)

# ============================================================
# 10. Marker dotplots: Epidermis and Mesophyll
# ============================================================

resolve_markers <- function(seu_obj, marker_list) {
  
  available <- rownames(seu_obj[["SCT"]])
  
  resolved <- list()
  
  for (nm in names(marker_list)) {
    candidates <- marker_list[[nm]]
    present <- candidates[candidates %in% available]
    if (length(present) > 0) {
      resolved[[nm]] <- present[1]
    }
  }
  
  data.frame(
    display_name = names(resolved),
    feature_name = unlist(resolved),
    stringsAsFactors = FALSE
  )
}

make_subcluster_object <- function(seu_obj, celltype_name) {
  
  ct_obj <- subset(seu_obj, subset = celltype == celltype_name)
  DefaultAssay(ct_obj) <- "SCT"
  
  ct_obj <- RunPCA(ct_obj, npcs = 30, verbose = FALSE)
  dims_use <- 1:min(20, ncol(Embeddings(ct_obj, "pca")))
  
  ct_obj <- FindNeighbors(ct_obj, dims = dims_use, verbose = FALSE)
  ct_obj <- FindClusters(ct_obj, resolution = 0.4, verbose = FALSE)
  ct_obj <- RunUMAP(
    ct_obj,
    dims = dims_use,
    reduction = "pca",
    reduction.name = "sub.umap",
    reduction.key = "UMAP_",
    verbose = FALSE
  )
  
  return(ct_obj)
}

epidermis_marker_list <- list(
  "PDE1 / AT1G17330" = c("AT1G17330", "PDE1"),
  "PR1 / AT2G14610" = c("PR1", "AT2G14610"),
  
  "KAT1"   = c("KAT1", "AT5G46240"),
  "SLAC1"  = c("SLAC1", "AT1G12480"),
  "ALMT12" = c("ALMT12", "QUAC1", "AT4G17970"),
  "MYB60"  = c("MYB60", "AT1G08810"),
  "FAMA"   = c("FAMA", "AT3G24140"),
  "MUTE"   = c("MUTE", "AT3G06120"),
  "SPCH"   = c("SPCH", "AT5G53210"),
  "TMM"    = c("TMM", "AT1G80080"),
  "EPF1"   = c("EPF1", "AT2G20875"),
  "EPF2"   = c("EPF2", "AT1G34245"),
  "HT1"    = c("HT1", "AT1G62400"),
  
  "ATML1"  = c("ATML1", "AT4G21750"),
  "PDF1"   = c("PDF1", "AT2G42840"),
  "GL2"    = c("GL2", "AT1G79840"),
  "LTPG1"  = c("LTPG1", "AT1G27950"),
  "LTPG2"  = c("LTPG2", "AT3G43720"),
  "FDH"    = c("FDH", "AT2G26250"),
  "LHCB1.1" = c("LHCB1.1", "AT1G29920"),
  
  "NPR1"   = c("NPR1", "AT1G64280"),
  "EDS1"   = c("EDS1", "AT3G48090"),
  "PAD4"   = c("PAD4", "AT3G52430"),
  "RBOHD"  = c("RBOHD", "AT5G47910"),
  "WRKY33" = c("WRKY33", "AT2G38470"),
  "NHL10"  = c("NHL10", "AT2G35980"),
  "CYP81F2" = c("CYP81F2", "AT5G57220")
)

mesophyll_marker_list <- list(
  "PDE1 / AT1G17330" = c("AT1G17330", "PDE1"),
  "PR1 / AT2G14610" = c("PR1", "AT2G14610"),
  
  "LHCB1.1" = c("LHCB1.1", "AT1G29920"),
  "LHCB2.1" = c("LHCB2.1"),
  "PSBO1"   = c("PSBO1"),
  "CHLH"    = c("CHLH"),
  
  "NPR1"   = c("NPR1", "AT1G64280"),
  "EDS1"   = c("EDS1", "AT3G48090"),
  "PAD4"   = c("PAD4", "AT3G52430"),
  "RBOHD"  = c("RBOHD", "AT5G47910"),
  "WRKY33" = c("WRKY33", "AT2G38470"),
  "NHL10"  = c("NHL10", "AT2G35980"),
  "CYP81F2" = c("CYP81F2", "AT5G57220"),
  "ALD1"   = c("ALD1", "AT2G13810"),
  "FMO1"   = c("FMO1", "AT1G19250"),
  
  "KAT1"   = c("KAT1", "AT5G46240"),
  "SLAC1"  = c("SLAC1", "AT1G12480"),
  "ALMT12" = c("ALMT12", "QUAC1", "AT4G17970"),
  "MYB60"  = c("MYB60", "AT1G08810"),
  
  "APL"    = c("APL", "AT1G79430"),
  "SUC2"   = c("SUC2", "AT1G22710")
)

make_marker_dotplot <- function(seu_obj, celltype_name, marker_list, out_file) {
  
  DefaultAssay(seu_obj) <- "SCT"
  
  Idents(seu_obj) <- "seurat_clusters"
  
  clust_counts <- table(Idents(seu_obj))
  keep_clusters <- names(clust_counts[clust_counts >= 30])
  seu_obj <- subset(seu_obj, idents = keep_clusters)
  
  current_ids <- as.character(Idents(seu_obj))
  numeric_ids <- suppressWarnings(as.numeric(current_ids))
  if (all(!is.na(numeric_ids))) {
    Idents(seu_obj) <- factor(current_ids, levels = as.character(sort(unique(numeric_ids))))
  }
  
  resolved_markers <- resolve_markers(seu_obj, marker_list)
  
  if (nrow(resolved_markers) < 2) {
    stop("Too few markers found for ", celltype_name)
  }
  
  p <- DotPlot(
    seu_obj,
    features = resolved_markers$feature_name,
    assay = "SCT",
    dot.scale = 5
  ) +
    RotatedAxis() +
    scale_x_discrete(labels = resolved_markers$display_name) +
    labs(
      title = paste0(celltype_name, " subcluster marker annotation"),
      subtitle = "Grouped by seurat_clusters; dot size = % detected; colour = average scaled expression",
      x = "Marker genes",
      y = paste0(celltype_name, " subclusters")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 9),
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey90"),
      panel.grid.minor = element_blank()
    )
  
  ggsave(
    filename = out_file,
    plot = p,
    width = 18,
    height = 9,
    dpi = 300,
    bg = "white"
  )
}

epi_obj <- make_subcluster_object(obj, "Epidermis")
meso_obj <- make_subcluster_object(obj, "Mesophyll")

make_marker_dotplot(
  seu_obj = epi_obj,
  celltype_name = "Epidermis",
  marker_list = epidermis_marker_list,
  out_file = file.path(out_marker, "Fig10_Epidermis_marker_dotplot.png")
)

make_marker_dotplot(
  seu_obj = meso_obj,
  celltype_name = "Mesophyll",
  marker_list = mesophyll_marker_list,
  out_file = file.path(out_marker, "Fig13_Mesophyll_marker_dotplot.png")
)

# ============================================================
# 11. Binary subcluster panels: PDE1/PR1 in Epidermis/Mesophyll
# ============================================================

make_subcluster_binary_panel <- function(ct_obj, celltype_name, gene_id, gene_label, out_file) {
  
  DefaultAssay(ct_obj) <- "SCT"
  
  expr_vec <- get_expr_vec(ct_obj, gene_id, assay_name = "SCT")
  
  umap_df <- Embeddings(ct_obj, "sub.umap") %>%
    as.data.frame() %>%
    rownames_to_column("cell")
  
  colnames(umap_df)[2:3] <- c("UMAP_1", "UMAP_2")
  
  meta_df <- ct_obj@meta.data %>%
    as.data.frame() %>%
    rownames_to_column("cell")
  
  plot_df <- left_join(meta_df, umap_df, by = "cell") %>%
    mutate(
      expr = expr_vec[cell],
      detected = expr > 0,
      condition = case_when(
        grepl("^00_Mock", sample) ~ "Mock",
        grepl("^DC3000", sample) ~ "DC3000",
        grepl("^AvrRpm1", sample) ~ "AvrRpm1",
        grepl("^AvrRpt2", sample) ~ "AvrRpt2",
        TRUE ~ NA_character_
      ),
      time_plot = case_when(
        grepl("^00_Mock", sample) ~ "Mock",
        grepl("04h", sample) ~ "04h",
        grepl("06h", sample) ~ "06h",
        grepl("09h", sample) ~ "09h",
        grepl("24h", sample) ~ "24h",
        TRUE ~ NA_character_
      ),
      group_id = ifelse(condition == "Mock", "Mock", paste(condition, time_plot, sep = "_"))
    ) %>%
    filter(!is.na(condition), !is.na(time_plot))
  
  # Downsample all condition-time groups to the smallest group size
  group_counts <- plot_df %>%
    count(group_id, name = "n")
  
  n_down <- min(group_counts$n)
  
  plot_df_ds <- plot_df %>%
    group_by(group_id) %>%
    slice_sample(n = n_down) %>%
    ungroup()
  
  # Same full background in every panel
  background_df <- plot_df %>%
    select(cell, UMAP_1, UMAP_2) %>%
    distinct()
  
  panel_grid <- expand.grid(
    panel_row = c("DC3000", "AvrRpm1", "AvrRpt2"),
    time_plot = c("Mock", "04h", "06h", "09h", "24h"),
    stringsAsFactors = FALSE
  )
  
  background_full <- tidyr::crossing(background_df, panel_grid)
  
  background_full$panel_row <- factor(
    background_full$panel_row,
    levels = c("DC3000", "AvrRpm1", "AvrRpt2")
  )
  
  background_full$time_plot <- factor(
    background_full$time_plot,
    levels = c("Mock", "04h", "06h", "09h", "24h")
  )
  
  # Detected cells
  mock_detected <- plot_df_ds %>%
    filter(group_id == "Mock", detected) %>%
    select(cell, UMAP_1, UMAP_2, time_plot) %>%
    tidyr::crossing(panel_row = c("DC3000", "AvrRpm1", "AvrRpt2"))
  
  nonmock_detected <- plot_df_ds %>%
    filter(group_id != "Mock", detected) %>%
    mutate(panel_row = as.character(condition)) %>%
    select(cell, UMAP_1, UMAP_2, panel_row, time_plot)
  
  detected_df <- bind_rows(mock_detected, nonmock_detected)
  
  detected_df$panel_row <- factor(
    detected_df$panel_row,
    levels = c("DC3000", "AvrRpm1", "AvrRpt2")
  )
  
  detected_df$time_plot <- factor(
    detected_df$time_plot,
    levels = c("Mock", "04h", "06h", "09h", "24h")
  )
  
  n_detected_total <- sum(plot_df$detected, na.rm = TRUE)
  pct_detected_total <- round(100 * mean(plot_df$detected, na.rm = TRUE), 2)
  
  p <- ggplot() +
    geom_point(
      data = background_full,
      aes(x = UMAP_1, y = UMAP_2),
      colour = "grey85",
      size = 0.28,
      alpha = 0.7
    ) +
    geom_point(
      data = detected_df,
      aes(x = UMAP_1, y = UMAP_2),
      colour = "blue3",
      size = 0.9,
      alpha = 0.95
    ) +
    facet_grid(panel_row ~ time_plot, drop = FALSE) +
    labs(
      title = paste0(gene_label, " after subclustering ", celltype_name),
      subtitle = paste0(
        "Grey = all ", celltype_name,
        " cells; blue = detected cells | pooled replicates | downsampled to n = ",
        n_down, " cells per condition-time | detected cells = ",
        n_detected_total, " (", pct_detected_total, "%)"
      ),
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      strip.background = element_rect(fill = "grey92", colour = "black"),
      strip.text = element_text(face = "bold", size = 12),
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.subtitle = element_text(size = 10.5, hjust = 0.5),
      panel.spacing = unit(0.45, "lines")
    )
  
  ggsave(
    filename = out_file,
    plot = p,
    width = 16,
    height = 11,
    dpi = 300,
    bg = "white"
  )
}

make_subcluster_binary_panel(
  ct_obj = epi_obj,
  celltype_name = "Epidermis",
  gene_id = "AT1G17330",
  gene_label = "PDE1_AT1G17330",
  out_file = file.path(out_subcluster, "Fig11_PDE1_AT1G17330_Epidermis_binary_subcluster_panels.png")
)

make_subcluster_binary_panel(
  ct_obj = meso_obj,
  celltype_name = "Mesophyll",
  gene_id = "AT1G17330",
  gene_label = "PDE1_AT1G17330",
  out_file = file.path(out_subcluster, "Fig14_PDE1_AT1G17330_Mesophyll_binary_subcluster_panels.png")
)

make_subcluster_binary_panel(
  ct_obj = epi_obj,
  celltype_name = "Epidermis",
  gene_id = "PR1",
  gene_label = "PR1_AT2G14610",
  out_file = file.path(out_subcluster, "Fig12_PR1_AT2G14610_Epidermis_binary_subcluster_panels.png")
)

make_subcluster_binary_panel(
  ct_obj = meso_obj,
  celltype_name = "Mesophyll",
  gene_id = "PR1",
  gene_label = "PR1_AT2G14610",
  out_file = file.path(out_subcluster, "Fig15_PR1_AT2G14610_Mesophyll_binary_subcluster_panels.png")
)

# ============================================================
# Done
# ============================================================

message("Done. Final figures generated only.")
message("Figure folders:")
message(out_fig2)
message(out_15panel)
message(out_dotplot)
message(out_trend)
message(out_marker)
message(out_subcluster)
