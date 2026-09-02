# =============================================================================
# final_generate_developmental_spatial_figures.R
#
# Purpose:
# Generate final report figures for Lee et al. developmental/spatial analysis
# of PDE1 / AT1G17330 and PR1 / AT2G14610.
#
# Input directory:
#   data/raw
#
# Output directory:
#   results/figures
# Notes:
#   - This script produces figures only.
#   - No CSV/table files are written.
#   - Large raw files should remain excluded from GitHub via .gitignore.
# =============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(1234)

# =============================================================================
# 0. Packages
# =============================================================================

required_pkgs <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "ggplot2",
  "patchwork",
  "dplyr",
  "tidyr",
  "tibble",
  "stringr",
  "scales",
  "forcats",
  "readxl",
  "curl"
)

missing_pkgs <- setdiff(required_pkgs, rownames(installed.packages()))

if (length(missing_pkgs) > 0) {
  stop(
    "Install missing packages first:\ninstall.packages(c(",
    paste(sprintf('"%s"', missing_pkgs), collapse = ", "),
    "))"
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(stringr)
  library(scales)
  library(forcats)
  library(readxl)
  library(curl)
})

# =============================================================================
# 1. Paths
# =============================================================================

BASE_DIR <- "C:/Users/qo25519/Documents/PDE1_bulkRNA_analysis/PDE_candidate_bulkRNA/PDE1-expression-analysis-in-Arabidopsis"

RAW_DIR <- file.path(BASE_DIR, "data", "raw")
FIG_DIR <- file.path(BASE_DIR, "results", "figures")
SUPP_DIR <- file.path(RAW_DIR, "published_supplementary_tables")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SUPP_DIR, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(RAW_DIR)) {
  stop("RAW_DIR does not exist: ", RAW_DIR)
}

# =============================================================================
# 2. Gene definitions
# =============================================================================

PDE1_ID <- "AT1G17330"
PR1_ID  <- "AT2G14610"

GENES <- c(
  PDE1 = "AT1G17330",
  PR1 = "AT2G14610",
  PDE1_like_AT5G40270 = "AT5G40270",
  PDE1_like_AT5G40290 = "AT5G40290",
  MIOX1 = "AT1G14520",
  MIOX2 = "AT2G19800",
  RSH1_candidate = "AT4G02260",
  MIOX5 = "AT5G56640",
  Candidate_AT1G26160 = "AT1G26160",
  Candidate_AT2G23820 = "AT2G23820"
)

# =============================================================================
# 3. General helper functions
# =============================================================================

find_first_file <- function(pattern, search_dir = RAW_DIR, recursive = TRUE) {
  hits <- list.files(
    search_dir,
    pattern = pattern,
    recursive = recursive,
    full.names = TRUE,
    ignore.case = TRUE
  )
  
  if (length(hits) == 0) {
    return(NA_character_)
  }
  
  hits[1]
}

save_png <- function(plot, filename, width, height, dpi = 500) {
  ggsave(
    filename = file.path(FIG_DIR, filename),
    plot = plot,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  message("Saved: ", file.path(FIG_DIR, filename))
}

read_rds_any <- function(path, max_layers = 4) {
  
  if (!file.exists(path)) {
    stop("File not found: ", path)
  }
  
  obj <- tryCatch(
    readRDS(path),
    error = function(e) NULL
  )
  
  if (!is.null(obj)) {
    return(obj)
  }
  
  current_file <- path
  temp_files <- character(0)
  
  on.exit({
    if (length(temp_files) > 0) {
      unlink(temp_files, force = TRUE)
    }
  }, add = TRUE)
  
  for (layer in seq_len(max_layers)) {
    
    con <- file(current_file, "rb")
    magic <- readBin(con, what = "raw", n = 2)
    close(con)
    
    is_gzip <- (
      length(magic) == 2 &&
        as.integer(magic[1]) == 0x1f &&
        as.integer(magic[2]) == 0x8b
    )
    
    if (!is_gzip) {
      obj <- readRDS(current_file)
      return(obj)
    }
    
    tmp <- tempfile(pattern = "unzipped_", fileext = ".rds")
    temp_files <- c(temp_files, tmp)
    
    in_con <- gzfile(current_file, "rb")
    out_con <- file(tmp, "wb")
    
    repeat {
      buffer <- readBin(in_con, what = "raw", n = 8 * 1024 * 1024)
      if (length(buffer) == 0) break
      writeBin(buffer, out_con)
    }
    
    close(in_con)
    close(out_con)
    
    current_file <- tmp
    
    obj <- tryCatch(
      readRDS(current_file),
      error = function(e) NULL
    )
    
    if (!is.null(obj)) {
      return(obj)
    }
  }
  
  stop("Could not read RDS after removing gzip layers: ", path)
}

get_assay_names <- function(obj) {
  out <- tryCatch(Assays(obj), error = function(e) NULL)
  if (is.null(out)) out <- names(obj@assays)
  as.character(out)
}

pick_reduction <- function(obj) {
  reds <- tryCatch(Reductions(obj), error = function(e) character(0))
  
  if ("umap" %in% reds) return("umap")
  um <- reds[grepl("umap", reds, ignore.case = TRUE)]
  if (length(um) > 0) return(um[1])
  if ("tsne" %in% reds) return("tsne")
  
  stop("No UMAP/tSNE reduction found. Available: ", paste(reds, collapse = ", "))
}

pick_global_reduction <- function(obj) {
  reds <- tryCatch(Reductions(obj), error = function(e) character(0))
  
  if ("tsne" %in% reds) return("tsne")
  if ("umap" %in% reds) return("umap")
  um <- reds[grepl("umap", reds, ignore.case = TRUE)]
  if (length(um) > 0) return(um[1])
  
  stop("No global tSNE/UMAP reduction found.")
}

pick_annotation_col <- function(obj) {
  nms <- colnames(obj@meta.data)
  
  candidates <- c(
    "celltype", "cell_type", "CellType", "cell.type", "Cell.Type",
    "celltype_annotation", "cell_type_annotation", "annotation",
    "cluster_annotation", "celltype_manual", "predicted.celltype",
    "predicted_celltype"
  )
  
  hit <- candidates[candidates %in% nms]
  if (length(hit) > 0) return(hit[1])
  
  regex_hit <- grep("cell.?type|annot", nms, ignore.case = TRUE, value = TRUE)
  if (length(regex_hit) > 0) return(regex_hit[1])
  
  if ("seurat_clusters" %in% nms) return("seurat_clusters")
  
  NA_character_
}

resolve_feature <- function(obj, gene_id) {
  assays <- get_assay_names(obj)
  preferred <- unique(c("RNA", "SCT", "MERFISH", "Spatial", assays))
  preferred <- preferred[preferred %in% assays]
  
  for (assay in preferred) {
    rn <- tryCatch(rownames(obj[[assay]]), error = function(e) character(0))
    if (length(rn) == 0) next
    
    exact <- which(toupper(rn) == toupper(gene_id))
    if (length(exact) > 0) {
      return(list(assay = assay, feature = rn[exact[1]]))
    }
    
    stripped <- sub("\\.[0-9]+$", "", rn)
    hit <- which(toupper(stripped) == toupper(gene_id))
    if (length(hit) > 0) {
      return(list(assay = assay, feature = rn[hit[1]]))
    }
  }
  
  NULL
}

fetch_vec <- function(obj, assay, feature, what = c("data", "counts")) {
  what <- match.arg(what)
  
  old_assay <- DefaultAssay(obj)
  on.exit({
    try(DefaultAssay(obj) <- old_assay, silent = TRUE)
  }, add = TRUE)
  
  DefaultAssay(obj) <- assay
  
  x <- tryCatch({
    z <- FetchData(obj, vars = feature, layer = what)
    setNames(as.numeric(z[[1]]), rownames(z))
  }, error = function(e) NULL)
  
  if (!is.null(x)) return(x)
  
  x <- tryCatch({
    z <- FetchData(obj, vars = feature, slot = what)
    setNames(as.numeric(z[[1]]), rownames(z))
  }, error = function(e) NULL)
  
  x
}

get_gene_values <- function(obj, gene_id) {
  hit <- resolve_feature(obj, gene_id)
  
  if (is.null(hit)) {
    return(list(
      present = FALSE,
      assay = NA_character_,
      feature = NA_character_,
      expr = NULL,
      counts = NULL,
      detected = NULL
    ))
  }
  
  expr <- fetch_vec(obj, hit$assay, hit$feature, "data")
  counts <- fetch_vec(obj, hit$assay, hit$feature, "counts")
  
  if (is.null(expr) && !is.null(counts)) {
    expr <- log1p(counts)
    names(expr) <- names(counts)
  }
  
  if (!is.null(counts)) {
    detected <- counts > 0
  } else if (!is.null(expr)) {
    detected <- expr > 0
  } else {
    detected <- NULL
  }
  
  list(
    present = TRUE,
    assay = hit$assay,
    feature = hit$feature,
    expr = expr,
    counts = counts,
    detected = detected
  )
}

align_numeric <- function(v, cells) {
  if (is.null(v)) return(rep(NA_real_, length(cells)))
  as.numeric(v[cells])
}

align_logical <- function(v, cells) {
  if (is.null(v)) return(rep(FALSE, length(cells)))
  out <- v[cells]
  out[is.na(out)] <- FALSE
  as.logical(out)
}

safe_q99 <- function(x) {
  z <- x[is.finite(x) & x > 0]
  if (length(z) == 0) return(1)
  
  q <- as.numeric(quantile(z, probs = 0.99, na.rm = TRUE, names = FALSE))
  if (!is.finite(q) || q <= 0) q <- max(z, na.rm = TRUE)
  if (!is.finite(q) || q <= 0) q <- 1
  
  q
}

make_embedding_df <- function(obj, reduction) {
  emb <- Embeddings(obj, reduction = reduction)
  
  data.frame(
    cell = rownames(emb),
    x = emb[, 1],
    y = emb[, 2],
    stringsAsFactors = FALSE
  )
}

pick_cluster_vector <- function(obj, cells) {
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    return(as.character(obj@meta.data[cells, "seurat_clusters", drop = TRUE]))
  }
  
  id <- as.character(Idents(obj))
  names(id) <- colnames(obj)
  
  id[cells]
}

umap_base_theme <- theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11, hjust = 0),
    plot.subtitle = element_text(size = 8.5, colour = "grey30"),
    plot.caption = element_text(size = 7.5, colour = "grey35"),
    plot.margin = margin(4, 4, 4, 4),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

plot_feature_umap <- function(df, value_col, gene_label, gene_id, stage) {
  v <- df[[value_col]]
  cap <- safe_q99(v)
  
  d <- df
  d$plot_value <- pmin(pmax(v, 0), cap)
  
  d0 <- d[!is.finite(d$plot_value) | d$plot_value <= 0, , drop = FALSE]
  d1 <- d[is.finite(d$plot_value) & d$plot_value > 0, , drop = FALSE]
  
  if (nrow(d1) > 0) {
    d1 <- d1[order(d1$plot_value), , drop = FALSE]
  }
  
  p <- ggplot() +
    geom_point(
      data = d0,
      aes(x = x, y = y),
      colour = "#C8C8C8",
      size = 0.24,
      alpha = 0.55,
      stroke = 0
    )
  
  if (nrow(d1) > 0) {
    p <- p +
      geom_point(
        data = d1,
        aes(x = x, y = y, colour = plot_value),
        size = 0.48,
        alpha = 1,
        stroke = 0
      ) +
      scale_colour_gradientn(
        colours = c("#56B4E9", "#2C7FB8", "#0868AC", "#084081", "#00204D"),
        limits = c(0, cap),
        oob = scales::squish,
        name = "Expression"
      )
  }
  
  p +
    labs(
      title = paste0(gene_label, "  ", gene_id),
      subtitle = paste0(
        stage,
        " | detected: ",
        sprintf("%.2f", mean(v > 0, na.rm = TRUE) * 100),
        "%"
      ),
      caption = "Grey = not detected; blue = detected; colour capped at positive-cell q99"
    ) +
    umap_base_theme
}

# =============================================================================
# 4. Lee developmental atlas: 10 lifecycle UMAPs + summaries for S13
# =============================================================================

atlas <- tribble(
  ~stage_key,        ~stage,             ~file,
  "01_seed_0d",      "Seed 0 d",          "GSE226097_seed_0d_230221.rds.gz",
  "02_seed_125d",    "Seed 1.25 d",       "GSE226097_seed_125d_250805.rds",
  "03_seedling_3d",  "Seedling 3 d",      "GSE226097_seedling_3d_250805.rds",
  "04_seedling_6d",  "Seedling 6 d",      "GSE226097_seedling_6d_230221.rds.gz",
  "05_seedling_12d", "Seedling 12 d",     "GSE226097_seedling_12d_230221.rds.gz",
  "06_rosette_21d",  "Rosette 21 d",      "GSE226097_rosette_21d_230221.rds.gz",
  "07_rosette_30d",  "Rosette 30 d",      "GSE226097_rosette_30d_230221.rds.gz",
  "08_stem",         "Stem",              "GSE226097_stem_250805.rds",
  "09_flower",       "Flower",            "GSE226097_flower_250805.rds",
  "10_silique",      "Silique",           "GSE226097_silique_230221.rds.gz"
) %>%
  mutate(path = file.path(RAW_DIR, file))

STAGE_LEVELS <- atlas$stage

missing_stage <- atlas %>% filter(!file.exists(path))

if (nrow(missing_stage) > 0) {
  stop(
    "Missing developmental RDS files in data/raw:\n",
    paste(missing_stage$file, collapse = "\n")
  )
}

pde1_umap_plots <- list()
pr1_umap_plots <- list()
stage_summary_list <- list()
group_summary_list <- list()

message("\n============================================================")
message("Generating Fig_S1, Fig_S2 and Fig_S13 inputs")
message("============================================================")

for (i in seq_len(nrow(atlas))) {
  
  stage <- atlas$stage[i]
  path <- atlas$path[i]
  
  message("\nStage: ", stage)
  
  dev_obj <- read_rds_any(path)
  reduction <- pick_reduction(dev_obj)
  ann_col <- pick_annotation_col(dev_obj)
  
  pde <- get_gene_values(dev_obj, PDE1_ID)
  pr1 <- get_gene_values(dev_obj, PR1_ID)
  
  d <- make_embedding_df(dev_obj, reduction)
  d$cluster <- pick_cluster_vector(dev_obj, d$cell)
  
  if (!is.na(ann_col)) {
    d$annotation <- as.character(dev_obj@meta.data[d$cell, ann_col, drop = TRUE])
  } else {
    d$annotation <- d$cluster
  }
  
  d$annotation[is.na(d$annotation) | d$annotation == ""] <- "Unknown"
  
  d$PDE1_expr <- if (pde$present) align_numeric(pde$expr, d$cell) else 0
  d$PR1_expr  <- if (pr1$present) align_numeric(pr1$expr, d$cell) else 0
  d$PDE1_det  <- if (pde$present) align_logical(pde$detected, d$cell) else FALSE
  d$PR1_det   <- if (pr1$present) align_logical(pr1$detected, d$cell) else FALSE
  
  pde1_umap_plots[[stage]] <-
    plot_feature_umap(d, "PDE1_expr", "PDE1", PDE1_ID, stage) +
    theme(legend.position = "none")
  
  pr1_umap_plots[[stage]] <-
    plot_feature_umap(d, "PR1_expr", "PR1", PR1_ID, stage) +
    theme(legend.position = "none")
  
  stage_summary_list[[stage]] <- tibble(
    stage = stage,
    n_nuclei = nrow(d),
    pde1_pct_detected = mean(d$PDE1_det, na.rm = TRUE) * 100
  )
  
  group_summary_list[[stage]] <- d %>%
    group_by(annotation) %>%
    summarise(
      n = n(),
      PDE1_pct = mean(PDE1_det, na.rm = TRUE) * 100,
      PDE1_mean = mean(PDE1_expr, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(stage = stage) %>%
    filter(n >= 25)
  
  rm(dev_obj, d, pde, pr1)
  invisible(gc())
}

stage_summary <- bind_rows(stage_summary_list) %>%
  mutate(stage = factor(stage, levels = STAGE_LEVELS))

group_summary <- bind_rows(group_summary_list) %>%
  mutate(stage = factor(stage, levels = STAGE_LEVELS))

fig_s1 <- wrap_plots(pde1_umap_plots, ncol = 5) +
  plot_annotation(
    title = "PDE1 (AT1G17330) across the Arabidopsis life cycle",
    subtitle = "Published stage-specific embeddings; colour scales capped at stage q99"
  )

save_png(
  fig_s1,
  "fig17_PDE1_lifecycle_10_UMAPs.png",
  width = 15,
  height = 7.2,
  dpi = 500
)

fig_s2 <- wrap_plots(pr1_umap_plots, ncol = 5) +
  plot_annotation(
    title = "PR1 (AT2G14610) across the Arabidopsis life cycle",
    subtitle = "Developmental baseline comparison for defence-associated expression"
  )

save_png(
  fig_s2,
  "fig18_PR1_lifecycle_10_UMAPs.png",
  width = 15,
  height = 7.2,
  dpi = 500
)

# =============================================================================
# 5. Fig_S13: PDE1 cell-type enrichment heatmap
# =============================================================================

message("\n============================================================")
message("Generating Fig_S13_PDE1_celltype_enrichment_heatmap")
message("============================================================")

stage_baseline <- stage_summary %>%
  transmute(
    stage = factor(as.character(stage), levels = STAGE_LEVELS),
    stage_PDE1_pct = as.numeric(pde1_pct_detected),
    stage_PDE1_fraction = as.numeric(pde1_pct_detected) / 100
  )

enrich_dat <- group_summary %>%
  mutate(
    stage = factor(as.character(stage), levels = STAGE_LEVELS),
    annotation = as.character(annotation),
    n = as.numeric(n),
    PDE1_pct = as.numeric(PDE1_pct)
  ) %>%
  filter(
    !is.na(stage),
    !is.na(annotation),
    annotation != "",
    annotation != "Unknown",
    is.finite(n),
    is.finite(PDE1_pct),
    n >= 50
  ) %>%
  left_join(stage_baseline, by = "stage") %>%
  mutate(
    observed_PDE1_positive = n * PDE1_pct / 100,
    expected_PDE1_positive = n * stage_PDE1_fraction,
    enrichment_ratio = (observed_PDE1_positive + 0.5) /
      (expected_PDE1_positive + 0.5),
    log2_enrichment = log2(enrichment_ratio)
  )

annotation_enrichment_stats <- enrich_dat %>%
  group_by(annotation) %>%
  summarise(
    max_log2_enrichment = max(log2_enrichment, na.rm = TRUE),
    max_absolute_enrichment = max(abs(log2_enrichment), na.rm = TRUE),
    max_PDE1_pct = max(PDE1_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(max_PDE1_pct >= 1 | max_absolute_enrichment >= 1) %>%
  arrange(desc(max_absolute_enrichment))

MAX_ANNOTATIONS_ENRICH <- 35

keep_annotations <- annotation_enrichment_stats %>%
  slice_head(n = MAX_ANNOTATIONS_ENRICH) %>%
  pull(annotation)

annotation_order <- annotation_enrichment_stats %>%
  filter(annotation %in% keep_annotations) %>%
  arrange(max_log2_enrichment) %>%
  pull(annotation)

plot_enrich <- enrich_dat %>%
  filter(annotation %in% keep_annotations) %>%
  mutate(
    stage = factor(stage, levels = STAGE_LEVELS),
    annotation = factor(annotation, levels = annotation_order)
  ) %>%
  select(stage, annotation, log2_enrichment) %>%
  tidyr::complete(annotation, stage) %>%
  mutate(plot_log2_enrichment = pmax(pmin(log2_enrichment, 3), -3))

fig_s13 <- ggplot(
  plot_enrich,
  aes(x = stage, y = annotation, fill = plot_log2_enrichment)
) +
  geom_tile(colour = "white", linewidth = 0.35) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-3, 3),
    breaks = c(-3, -2, -1, 0, 1, 2, 3),
    labels = c("≤ -3", "-2", "-1", "0", "+1", "+2", "≥ +3"),
    na.value = "grey85",
    name = expression(log[2] * " enrichment")
  ) +
  labs(
    title = "Cell-type enrichment of PDE1 detection across Arabidopsis development",
    subtitle = "Red = detected more often than stage expectation; blue = depleted; white = expected frequency",
    x = NULL,
    y = NULL,
    caption = "Observed/expected enrichment calculated using stage-wide PDE1 detection; 0.5-cell pseudocount. Grey = annotation unavailable."
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9.5),
    axis.text.y = element_text(size = 8.5),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10, colour = "grey30"),
    plot.caption = element_text(size = 8, colour = "grey40"),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )

save_png(
  fig_s13,
  "fig19_PDE1_celltype_enrichment_heatmap.png",
  width = 11.5,
  height = max(6, 0.27 * length(annotation_order) + 2.8),
  dpi = 500
)

# =============================================================================
# 6. Fig_S19H: PDE1/PR1 across 655 published subclusters
# =============================================================================

message("\n============================================================")
message("Generating Fig_S19H_PDE1_PR1_655_subclusters")
message("============================================================")

SUPP6_FILE <- file.path(SUPP_DIR, "Lee2025_Supplementary_Table6_655_subclusters.xlsx")

SUPP6_URL <- paste0(
  "https://media.springernature.com/original/",
  "springer-static/esm/art%3A10.1038%2Fs41477-025-02072-z/",
  "MediaObjects/41477_2025_2072_MOESM8_ESM.xlsx"
)

if (!file.exists(SUPP6_FILE)) {
  message("Downloading Lee Supplementary Table 6...")
  curl::curl_download(SUPP6_URL, destfile = SUPP6_FILE, mode = "wb", quiet = FALSE)
}

supp6 <- readxl::read_excel(
  SUPP6_FILE,
  sheet = "sub_cluster_pseudobulk_rawcount",
  skip = 2
)

names(supp6)[1] <- "gene_id"
supp6$gene_id <- as.character(supp6$gene_id)

count_matrix <- as.matrix(supp6[, -1, drop = FALSE])
storage.mode(count_matrix) <- "numeric"

libsize <- colSums(count_matrix, na.rm = TRUE)
gene_clean <- toupper(sub("\\.[0-9]+$", "", supp6$gene_id))

pde_idx <- which(gene_clean == PDE1_ID)
pr1_idx <- which(gene_clean == PR1_ID)

if (length(pde_idx) == 0) stop("PDE1 not found in Supplementary Table 6.")
if (length(pr1_idx) == 0) stop("PR1 not found in Supplementary Table 6.")

h_s19 <- tibble(
  subcluster = colnames(count_matrix),
  original_order = seq_len(ncol(count_matrix)),
  PDE1_CPM = as.numeric(count_matrix[pde_idx[1], ]) / pmax(libsize, 1) * 1e6,
  PR1_CPM  = as.numeric(count_matrix[pr1_idx[1], ]) / pmax(libsize, 1) * 1e6
)

rm(supp6, count_matrix)
invisible(gc())

h_s19 <- h_s19 %>%
  mutate(
    stage = case_when(
      grepl("^seed_3d|Seed 0", subcluster, ignore.case = TRUE) ~ "Seed 0 d",
      grepl("^seed_425d|Seed 1.25", subcluster, ignore.case = TRUE) ~ "Seed 1.25 d",
      grepl("^seedling_6d|Seedling 3", subcluster, ignore.case = TRUE) ~ "Seedling 3 d",
      grepl("^seedling_9d|Seedling 6", subcluster, ignore.case = TRUE) ~ "Seedling 6 d",
      grepl("^seedling_15d|Seedling 12", subcluster, ignore.case = TRUE) ~ "Seedling 12 d",
      grepl("^rosette_21d|Rosette 21", subcluster, ignore.case = TRUE) ~ "Rosette 21 d",
      grepl("^rosette_30d|Rosette 30", subcluster, ignore.case = TRUE) ~ "Rosette 30 d",
      grepl("^stem", subcluster, ignore.case = TRUE) ~ "Stem",
      grepl("^flower", subcluster, ignore.case = TRUE) ~ "Flower",
      grepl("^silique", subcluster, ignore.case = TRUE) ~ "Silique",
      TRUE ~ "Other"
    ),
    stage = factor(stage, levels = c(STAGE_LEVELS, "Other"))
  ) %>%
  arrange(stage, original_order) %>%
  mutate(subcluster_index = row_number())

h_long <- h_s19 %>%
  select(subcluster, subcluster_index, stage, PDE1_CPM, PR1_CPM) %>%
  pivot_longer(
    cols = c(PDE1_CPM, PR1_CPM),
    names_to = "gene",
    values_to = "CPM"
  ) %>%
  mutate(
    gene = recode(gene, PDE1_CPM = "PDE1", PR1_CPM = "PR1"),
    log2_CPM = log2(CPM + 1)
  )

boundaries <- h_s19 %>%
  group_by(stage) %>%
  summarise(end_x = max(subcluster_index), .groups = "drop")

boundary_x <- if (nrow(boundaries) > 1) {
  boundaries$end_x[-nrow(boundaries)] + 0.5
} else {
  numeric(0)
}

fig_s19h <- ggplot(
  h_long,
  aes(x = subcluster_index, y = log2_CPM, colour = stage)
) +
  geom_vline(xintercept = boundary_x, colour = "grey85", linewidth = 0.35) +
  geom_point(size = 1.25, alpha = 0.85) +
  facet_grid(gene ~ ., scales = "fixed") +
  scale_colour_viridis_d(option = "turbo", end = 0.92, name = "Developmental stage") +
  labs(
    title = "PDE1 and PR1 expression across 655 developmental subclusters",
    subtitle = "Published Lee et al. subcluster pseudobulk counts normalised to CPM",
    x = "655 published developmental subclusters",
    y = expression(log[2] * "(CPM + 1)"),
    caption = "Each point is one published subcluster; points are coloured by developmental dataset."
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text.y = element_text(face = "bold", angle = 0),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 9.5, colour = "grey30"),
    legend.position = "bottom"
  )

save_png(
  fig_s19h,
  "fig20_PDE1_PR1_655_subclusters.png",
  width = 14,
  height = 7.5,
  dpi = 500
)

# =============================================================================
# 7. Fig_S19A: Global PDE1/PR1 atlas
# =============================================================================

message("\n============================================================")
message("Generating Fig_S19A_global_PDE1_PR1")
message("============================================================")

GLOBAL_FILE <- find_first_file(
  "^GSE226097_global_integration_221009.*\\.rds$",
  search_dir = RAW_DIR,
  recursive = TRUE
)

if (is.na(GLOBAL_FILE)) {
  stop(
    "Global integration object not found. Put this file in data/raw:\n",
    "GSE226097_global_integration_221009.rds"
  )
}

global_obj <- readRDS(GLOBAL_FILE)

if (!inherits(global_obj, "Seurat")) {
  stop("Global file is not a Seurat object: ", GLOBAL_FILE)
}

global_reduction <- pick_global_reduction(global_obj)

pde_global <- get_gene_values(global_obj, PDE1_ID)
pr1_global <- get_gene_values(global_obj, PR1_ID)

global_emb <- make_embedding_df(global_obj, global_reduction)

global_emb$PDE1_det <- if (pde_global$present) {
  align_logical(pde_global$detected, global_emb$cell)
} else {
  FALSE
}

global_emb$PR1_det <- if (pr1_global$present) {
  align_logical(pr1_global$detected, global_emb$cell)
} else {
  FALSE
}

global_emb$state <- case_when(
  global_emb$PDE1_det & global_emb$PR1_det ~ "PDE1+ / PR1+",
  global_emb$PDE1_det & !global_emb$PR1_det ~ "PDE1+ only",
  !global_emb$PDE1_det & global_emb$PR1_det ~ "PR1+ only",
  TRUE ~ "Neither"
)

global_emb$state <- factor(
  global_emb$state,
  levels = c("Neither", "PDE1+ only", "PR1+ only", "PDE1+ / PR1+")
)

global_background <- global_emb %>% filter(state == "Neither")
global_positive <- global_emb %>% filter(state != "Neither")

fig_s19a <- ggplot() +
  geom_point(
    data = global_background,
    aes(x = x, y = y),
    colour = "grey80",
    size = 0.06,
    alpha = 0.4
  ) +
  geom_point(
    data = global_positive,
    aes(x = x, y = y, colour = state),
    size = 0.35,
    alpha = 0.95
  ) +
  scale_colour_manual(
    values = c(
      "PDE1+ only" = "#0072B2",
      "PR1+ only" = "#CC79A7",
      "PDE1+ / PR1+" = "#111111"
    ),
    name = "Detection state"
  ) +
  coord_equal() +
  labs(
    title = "PDE1 and PR1 across the globally integrated Arabidopsis atlas",
    subtitle = paste0("Published global ", toupper(global_reduction), " coordinates"),
    x = NULL,
    y = NULL,
    caption = "Grey = neither gene detected; blue = PDE1 only; magenta = PR1 only; black = co-detection."
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 9.5, colour = "grey30"),
    legend.position = "right"
  )

save_png(
  fig_s19a,
  "fig16_global_PDE1_PR1.png",
  width = 10,
  height = 8,
  dpi = 500
)

rm(global_obj, global_emb, global_background, global_positive, pde_global, pr1_global)
invisible(gc())

# =============================================================================
# 8. Curio spatial helpers for Fig_SP13 and Fig_SP15
# =============================================================================

message("\n============================================================")
message("Generating Curio spatial figures SP13 and SP15")
message("============================================================")

CURIO_FILE <- find_first_file(
  "GSE226097_221121_curio_seed_seurat.*\\.rds(\\.gz)?$",
  search_dir = RAW_DIR,
  recursive = TRUE
)

if (is.na(CURIO_FILE)) {
  stop(
    "Curio seed Seurat object not found. Put this file in data/raw:\n",
    "GSE226097_221121_curio_seed_seurat.rds.gz"
  )
}

curio <- read_rds_any(CURIO_FILE)

if (!inherits(curio, "Seurat")) {
  stop("Curio file is not a Seurat object.")
}

get_physical_coordinates <- function(obj) {
  
  cells <- colnames(obj)
  
  imgs <- tryCatch(Images(obj), error = function(e) character(0))
  
  if (length(imgs) > 0) {
    for (img in imgs) {
      cc <- tryCatch(
        SeuratObject::GetTissueCoordinates(obj[[img]]) %>% as.data.frame(),
        error = function(e) NULL
      )
      
      if (!is.null(cc) && nrow(cc) > 0) {
        cc$cell <- rownames(cc)
        
        coordinate_pairs <- list(
          c("x", "y"),
          c("imagecol", "imagerow"),
          c("col", "row"),
          c("pxl_col_in_fullres", "pxl_row_in_fullres")
        )
        
        for (pair in coordinate_pairs) {
          if (all(pair %in% colnames(cc))) {
            out <- cc %>%
              transmute(
                cell = cell,
                x = as.numeric(.data[[pair[1]]]),
                y = as.numeric(.data[[pair[2]]])
              ) %>%
              filter(cell %in% cells)
            
            if (nrow(out) > 0) {
              attr(out, "coordinate_source") <- paste0(
                "Seurat image/FOV ",
                img,
                ": ",
                pair[1],
                "/",
                pair[2]
              )
              return(out)
            }
          }
        }
      }
    }
  }
  
  md <- obj@meta.data
  md$cell <- rownames(md)
  low <- tolower(colnames(md))
  names(low) <- colnames(md)
  
  find_ci <- function(name) {
    hit <- names(low)[low == tolower(name)]
    if (length(hit) == 0) NA_character_ else hit[1]
  }
  
  coordinate_pairs <- list(
    c("x", "y"),
    c("global_x", "global_y"),
    c("x_global", "y_global"),
    c("x_coord", "y_coord"),
    c("center_x", "center_y"),
    c("centroid_x", "centroid_y"),
    c("imagecol", "imagerow"),
    c("col", "row")
  )
  
  for (pair in coordinate_pairs) {
    xcol <- find_ci(pair[1])
    ycol <- find_ci(pair[2])
    
    if (!is.na(xcol) && !is.na(ycol)) {
      if (is.numeric(md[[xcol]]) && is.numeric(md[[ycol]])) {
        out <- md %>%
          transmute(
            cell = cell,
            x = as.numeric(.data[[xcol]]),
            y = as.numeric(.data[[ycol]])
          )
        
        attr(out, "coordinate_source") <- paste0("metadata: ", xcol, "/", ycol)
        return(out)
      }
    }
  }
  
  stop("Could not identify physical spatial coordinates in Curio object.")
}

pick_spatial_cluster_vector <- function(obj) {
  
  nms <- colnames(obj@meta.data)
  
  preferred <- c(
    "banksy_cluster", "banksy_clusters",
    "BANKSY_cluster", "BANKSY_clusters",
    "spatial_cluster", "spatial_clusters",
    "seurat_clusters", "cluster", "clusters"
  )
  
  hit <- preferred[preferred %in% nms]
  
  if (length(hit) > 0) {
    x <- as.character(obj@meta.data[, hit[1], drop = TRUE])
    names(x) <- rownames(obj@meta.data)
    return(list(column = hit[1], value = x))
  }
  
  regex_hit <- grep(
    "banksy.*cluster|spatial.*cluster|seurat.*cluster|cluster",
    nms,
    ignore.case = TRUE,
    value = TRUE
  )
  
  if (length(regex_hit) > 0) {
    x <- as.character(obj@meta.data[, regex_hit[1], drop = TRUE])
    names(x) <- rownames(obj@meta.data)
    return(list(column = regex_hit[1], value = x))
  }
  
  x <- as.character(Idents(obj))
  names(x) <- colnames(obj)
  
  list(column = "active Seurat identities", value = x)
}

coords <- get_physical_coordinates(curio)

cluster_info <- pick_spatial_cluster_vector(curio)
coords$cluster <- cluster_info$value[coords$cell]
coords$cluster[is.na(coords$cluster) | coords$cluster == ""] <- "Unknown"

pde_curio <- get_gene_values(curio, PDE1_ID)
pr1_curio <- get_gene_values(curio, PR1_ID)

if (!pde_curio$present) {
  stop("PDE1 / AT1G17330 was not found in the Curio object.")
}

coords$PDE1_expr <- align_numeric(pde_curio$expr, coords$cell)
coords$PDE1_det <- align_logical(pde_curio$detected, coords$cell)

if (pr1_curio$present) {
  coords$PR1_expr <- align_numeric(pr1_curio$expr, coords$cell)
  coords$PR1_det <- align_logical(pr1_curio$detected, coords$cell)
} else {
  coords$PR1_expr <- NA_real_
  coords$PR1_det <- FALSE
}

MIN_CLUSTER_SPOTS <- 25

cluster_summary <- coords %>%
  group_by(cluster) %>%
  summarise(
    n_spots = n(),
    PDE1_positive_n = sum(PDE1_det, na.rm = TRUE),
    PDE1_pct = mean(PDE1_det, na.rm = TRUE) * 100,
    PR1_positive_n = sum(PR1_det, na.rm = TRUE),
    PR1_pct = mean(PR1_det, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  filter(n_spots >= MIN_CLUSTER_SPOTS)

GLOBAL_PDE1_PCT <- mean(coords$PDE1_det, na.rm = TRUE) * 100

cluster_summary <- cluster_summary %>%
  mutate(
    PDE1_enrichment = (PDE1_pct + 0.05) / (GLOBAL_PDE1_PCT + 0.05),
    log2_PDE1_enrichment = log2(PDE1_enrichment)
  )

# =============================================================================
# 9. Map Curio clusters to published Curio signatures, then Fig_SP13
# =============================================================================

SUPP2_FILE <- file.path(SUPP_DIR, "Lee2025_Supplementary_Table2.xlsx")

SUPP2_URL <- paste0(
  "https://media.springernature.com/original/",
  "springer-static/esm/art%3A10.1038%2Fs41477-025-02072-z/",
  "MediaObjects/41477_2025_2072_MOESM4_ESM.xlsx"
)

if (!file.exists(SUPP2_FILE)) {
  message("Downloading Lee Supplementary Table 2...")
  curl::curl_download(SUPP2_URL, destfile = SUPP2_FILE, mode = "wb", quiet = FALSE)
}

pub_markers <- read_excel(
  SUPP2_FILE,
  sheet = "curio_seed_cluster_markers",
  skip = 2
)

pub_markers <- pub_markers %>%
  filter(!is.na(cluster)) %>%
  mutate(cluster = as.integer(cluster))

rna_features <- rownames(curio[["RNA"]])
rna_features_base <- sub("\\.[0-9]+$", "", rna_features)

resolve_pub_marker <- function(gene, agi) {
  
  candidates <- unique(na.omit(c(as.character(agi), as.character(gene))))
  
  for (candidate in candidates) {
    if (is.na(candidate) || candidate == "") next
    
    hit <- which(toupper(rna_features) == toupper(candidate))
    if (length(hit) > 0) return(rna_features[hit[1]])
    
    hit <- which(toupper(rna_features_base) == toupper(candidate))
    if (length(hit) > 0) return(rna_features[hit[1]])
  }
  
  NA_character_
}

pub_markers$matched_feature <- mapply(
  resolve_pub_marker,
  pub_markers$gene,
  pub_markers$AGI,
  USE.NAMES = FALSE
)

marker_sets <- split(pub_markers$matched_feature, pub_markers$cluster)
marker_sets <- lapply(marker_sets, function(x) unique(x[!is.na(x)]))

rna_counts <- tryCatch(
  LayerData(curio[["RNA"]], layer = "counts"),
  error = function(e) NULL
)

if (is.null(rna_counts)) {
  rna_counts <- tryCatch(
    GetAssayData(curio, assay = "RNA", slot = "counts"),
    error = function(e) NULL
  )
}

if (is.null(rna_counts)) {
  stop("Could not obtain RNA counts from Curio object.")
}

current_cluster <- as.character(curio$seurat_clusters)

cluster_levels <- sort(unique(as.integer(current_cluster)))
cluster_levels <- as.character(cluster_levels)

current_cluster <- factor(current_cluster, levels = cluster_levels)

design <- Matrix::sparse.model.matrix(~ 0 + current_cluster)
colnames(design) <- cluster_levels

pb_counts <- rna_counts %*% design

marker_union <- unique(unlist(marker_sets))
marker_union <- intersect(marker_union, rownames(pb_counts))

pb_small <- as.matrix(pb_counts[marker_union, , drop = FALSE])
library_size <- Matrix::colSums(pb_counts)
library_size[library_size == 0] <- 1

logcpm <- log1p(t(t(pb_small) / library_size * 1e6))

marker_z <- t(scale(t(logcpm)))
marker_z[!is.finite(marker_z)] <- 0

score_matrix <- sapply(
  names(marker_sets),
  function(pub_cluster) {
    features <- intersect(marker_sets[[pub_cluster]], rownames(marker_z))
    
    if (length(features) == 0) {
      return(rep(NA_real_, ncol(marker_z)))
    }
    
    colMeans(marker_z[features, , drop = FALSE], na.rm = TRUE)
  }
)

rownames(score_matrix) <- colnames(marker_z)
colnames(score_matrix) <- paste0("Published_", colnames(score_matrix))

score_long <- as.data.frame(score_matrix) %>%
  rownames_to_column("our_cluster") %>%
  pivot_longer(
    cols = starts_with("Published_"),
    names_to = "published_cluster",
    values_to = "signature_score"
  ) %>%
  mutate(
    published_number = sub("Published_", "", published_cluster)
  )

published_labels <- c(
  "0" = "Published cluster 0",
  "1" = "Published cluster 1",
  "2" = "Published cluster 2",
  "3" = "Published cluster 3",
  "4" = "Published cluster 4 — Stele",
  "5" = "Published cluster 5 — Cotyledon",
  "6" = "Published cluster 6 — Unannotated"
)

score_long <- score_long %>%
  mutate(published_label = published_labels[published_number])

best_mapping <- score_long %>%
  group_by(our_cluster) %>%
  slice_max(signature_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(as.integer(our_cluster)) %>%
  select(our_cluster, published_number, published_label, signature_score)

sp13 <- cluster_summary %>%
  left_join(
    best_mapping %>%
      select(our_cluster, published_number, published_label, signature_score),
    by = c("cluster" = "our_cluster")
  ) %>%
  group_by(published_number, published_label) %>%
  summarise(
    mean_signature_score = weighted.mean(signature_score, w = n_spots, na.rm = TRUE),
    n_spots = sum(n_spots),
    PDE1_positive_n = sum(PDE1_positive_n),
    PR1_positive_n = sum(PR1_positive_n),
    .groups = "drop"
  ) %>%
  mutate(
    PDE1_pct = 100 * PDE1_positive_n / n_spots,
    PR1_pct = 100 * PR1_positive_n / n_spots
  )

global_pde1 <- 100 * sum(sp13$PDE1_positive_n) / sum(sp13$n_spots)

sp13 <- sp13 %>%
  mutate(
    fold_enrichment = (PDE1_pct + 0.001) / (global_pde1 + 0.001),
    log2_enrichment = log2(fold_enrichment),
    label = case_when(
      published_number == "4" ~ "Cluster 4 — Stele",
      published_number == "5" ~ "Cluster 5 — Cotyledon",
      published_number == "6" ~ "Cluster 6 — Unannotated",
      TRUE ~ paste0("Published cluster ", published_number)
    )
  ) %>%
  arrange(desc(PDE1_pct))

fig_sp13 <- sp13 %>%
  mutate(label = factor(label, levels = rev(label))) %>%
  ggplot(aes(x = PDE1_pct, y = label)) +
  geom_vline(
    xintercept = global_pde1,
    linetype = "dashed",
    colour = "grey45",
    linewidth = 0.7
  ) +
  geom_col(aes(fill = log2_enrichment), width = 0.72) +
  geom_text(
    aes(label = paste0(PDE1_positive_n, "/", n_spots)),
    hjust = -0.1,
    size = 3.5
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "grey90",
    high = "#B2182B",
    midpoint = 0,
    name = expression(log[2] * " enrichment")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "PDE1 across published Curio spatial signatures",
    subtitle = paste0(
      "Dashed line = whole-section PDE1 detection (",
      sprintf("%.3f", global_pde1),
      "%)"
    ),
    x = "% spatial spots detecting PDE1",
    y = NULL,
    caption = "Numbers indicate PDE1-positive spots / total spots assigned to each published marker signature."
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.y = element_text(size = 10)
  )

save_png(
  fig_sp13,
  "fig21_PDE1_published_Curio_signatures.png",
  width = 10,
  height = 6.5,
  dpi = 500
)

# =============================================================================
# 10. Fig_SP15: PDE1 stele spatial overlay
# =============================================================================

HIGH_CONF_STELE_CLUSTER <- "2"

sp15_dat <- coords %>%
  mutate(
    stele_highconf = cluster == HIGH_CONF_STELE_CLUSTER,
    display_group = case_when(
      PDE1_det & stele_highconf ~ "PDE1+ in stele-associated domain",
      PDE1_det & !stele_highconf ~ "PDE1+ outside stele-associated domain",
      !PDE1_det & stele_highconf ~ "Stele-associated domain",
      TRUE ~ "Other spatial spots"
    )
  )

sp15_other <- sp15_dat %>% filter(display_group == "Other spatial spots")
sp15_stele <- sp15_dat %>% filter(display_group == "Stele-associated domain")
sp15_pde_out <- sp15_dat %>% filter(display_group == "PDE1+ outside stele-associated domain")
sp15_pde_stele <- sp15_dat %>% filter(display_group == "PDE1+ in stele-associated domain")

fig_sp15 <- ggplot() +
  geom_point(
    data = sp15_other,
    aes(x = x, y = y),
    colour = "grey88",
    size = 0.25,
    alpha = 0.45
  ) +
  geom_point(
    data = sp15_stele,
    aes(x = x, y = y),
    colour = "#A6CEE3",
    size = 0.38,
    alpha = 0.75
  ) +
  geom_point(
    data = sp15_pde_out,
    aes(x = x, y = y),
    colour = "#F28E2B",
    size = 1.15,
    alpha = 1
  ) +
  geom_point(
    data = sp15_pde_stele,
    aes(x = x, y = y),
    colour = "#08306B",
    size = 1.5,
    alpha = 1
  ) +
  coord_equal() +
  labs(
    title = "PDE1 localization relative to the stele-associated spatial domain",
    subtitle = paste0(
      "High-confidence Curio cluster 2: ",
      sum(sp15_dat$PDE1_det & sp15_dat$stele_highconf),
      "/",
      sum(sp15_dat$stele_highconf),
      " spots detect PDE1"
    ),
    x = NULL,
    y = NULL,
    caption = "Light blue = stele-associated domain; dark blue = PDE1+ within this domain; orange = PDE1+ elsewhere."
  ) +
  theme_void() +
  theme(
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 10, colour = "grey25"),
    plot.caption = element_text(size = 8, colour = "grey40")
  )

save_png(
  fig_sp15,
  "fig22_PDE1_stele_spatial_overlay.png",
  width = 10,
  height = 8,
  dpi = 500
)

rm(curio, coords, pde_curio, pr1_curio)
invisible(gc())

# =============================================================================
# 11. Flower label-transfer figures: Fl10A/B/C
# =============================================================================
# =============================================================================
# Fig_Fl10A/B/C — FINAL FLOWER snRNA–MERFISH SPATIAL INFERENCE FIGURES
#
# LEFT  = all MERFISH cells labelled by transferred Flower snRNA cell state
# RIGHT = high-confidence PDE1-associated Flower cell states only
#
# IMPORTANT:
# PDE1 was NOT directly targeted by Flower MERFISH.
# These are inferred PDE1-associated spatial cell states from label transfer.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

# -------------------------------------------------------------------------
# Output directory for report/GitHub project
# -------------------------------------------------------------------------

BASE_DIR <- "C:/Users/qo25519/Documents/PDE1_bulkRNA_analysis/PDE_candidate_bulkRNA/PDE1-expression-analysis-in-Arabidopsis"

FIG_OUT <- file.path(
  BASE_DIR,
  "results",
  "figures"
)

dir.create(
  FIG_OUT,
  recursive = TRUE,
  showWarnings = FALSE
)
# =============================================================================
# Reload Flower label-transfer outputs for Fig_Fl10A/B/C
#
# These files are created by 04_Lee2025_Flower_spatial_integration.R
# They contain the already transferred Flower cell-type labels and coordinates.
# =============================================================================

FLOWER_TRANSFER_DIR <- file.path(
  RAW_DIR,
  "Flower_spatial_integration_results"
)

# If your previous Flower results are still in Downloads, use this fallback
FLOWER_TRANSFER_DIR_FALLBACK <- "C:/Users/qo25519/Downloads/sc_dev/Flower_spatial_integration_results"

COORDS_FILE <- file.path(
  FLOWER_TRANSFER_DIR,
  "02_tables",
  "Flower_spatial_predictions_with_coordinates.csv"
)

REF_SUMMARY_FILE <- file.path(
  FLOWER_TRANSFER_DIR,
  "02_tables",
  "Flower_snRNA_PDE1_PR1_by_celltype.csv"
)

if (!file.exists(COORDS_FILE)) {
  COORDS_FILE <- file.path(
    FLOWER_TRANSFER_DIR_FALLBACK,
    "02_tables",
    "Flower_spatial_predictions_with_coordinates.csv"
  )
}

if (!file.exists(REF_SUMMARY_FILE)) {
  REF_SUMMARY_FILE <- file.path(
    FLOWER_TRANSFER_DIR_FALLBACK,
    "02_tables",
    "Flower_snRNA_PDE1_PR1_by_celltype.csv"
  )
}

if (!file.exists(COORDS_FILE)) {
  stop(
    "Cannot find Flower_spatial_predictions_with_coordinates.csv.\n",
    "Run 04_Lee2025_Flower_spatial_integration.R first, or copy the Flower_spatial_integration_results folder into data/raw."
  )
}

if (!file.exists(REF_SUMMARY_FILE)) {
  stop(
    "Cannot find Flower_snRNA_PDE1_PR1_by_celltype.csv.\n",
    "Run 04_Lee2025_Flower_spatial_integration.R first, or copy the Flower_spatial_integration_results folder into data/raw."
  )
}

coords_plot <- readr::read_csv(
  COORDS_FILE,
  show_col_types = FALSE
)

ref_summary <- readr::read_csv(
  REF_SUMMARY_FILE,
  show_col_types = FALSE
)

top_pde_types <- ref_summary %>%
  filter(PDE1_n > 0) %>%
  arrange(
    desc(PDE1_pct),
    desc(PDE1_n)
  ) %>%
  slice_head(n = 5) %>%
  pull(annotation)

message("Loaded Flower transferred coordinates from:")
message(COORDS_FILE)

message("Loaded Flower snRNA PDE1 summary from:")
message(REF_SUMMARY_FILE)

message("Top PDE1-associated Flower cell states:")
print(top_pde_types)
# -------------------------------------------------------------------------
# Safety checks
# -------------------------------------------------------------------------

if (!exists("coords_plot")) {
  stop("coords_plot not found. Run the Flower label-transfer section first.")
}

if (!exists("top_pde_types")) {
  stop("top_pde_types not found. Run the Flower snRNA PDE1 cell-type ranking first.")
}

required_cols <- c(
  "sample",
  "x",
  "y",
  "predicted_flower_celltype",
  "transfer_score"
)

missing_cols <- setdiff(
  required_cols,
  colnames(coords_plot)
)

if (length(missing_cols) > 0) {
  stop(
    "coords_plot is missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# -------------------------------------------------------------------------
# Clean plotting data
# -------------------------------------------------------------------------

coords_plot <- coords_plot %>%
  mutate(
    sample = as.character(sample),
    predicted_flower_celltype = as.character(predicted_flower_celltype),
    transfer_score = as.numeric(transfer_score),
    x = as.numeric(x),
    y = as.numeric(y)
  ) %>%
  filter(
    !is.na(sample),
    !is.na(predicted_flower_celltype),
    is.finite(x),
    is.finite(y)
  )

# -------------------------------------------------------------------------
# Fixed cell-type colours across all three panels
# -------------------------------------------------------------------------

all_celltypes <- sort(
  unique(
    coords_plot$predicted_flower_celltype
  )
)

all_celltypes <- all_celltypes[
  !is.na(all_celltypes) &
    all_celltypes != ""
]

celltype_colours <- scales::hue_pal()(
  length(all_celltypes)
)

names(celltype_colours) <- all_celltypes

# -------------------------------------------------------------------------
# Figure labels and output names
# -------------------------------------------------------------------------

flower_samples <- c(
  "region1",
  "region2",
  "region2_long"
)

panel_letters <- c(
  region1 = "A",
  region2 = "B",
  region2_long = "C"
)

sample_labels <- c(
  region1 = "Flower region 1",
  region2 = "Flower region 2",
  region2_long = "Flower region 2 - longitudinal section"
)

output_names <- c(
  region1 = "fig23_Flower_region1_reference_vs_PDE1",
  region2 = "fig24_Flower_region2_reference_vs_PDE1",
  region2_long = "fig25_Flower_region2_long_reference_vs_PDE1"
)

# -------------------------------------------------------------------------
# Create one figure per Flower region
# -------------------------------------------------------------------------

for (nm in flower_samples) {
  
  message("\n============================================================")
  message("Creating ", output_names[[nm]])
  message("============================================================")
  
  d <- coords_plot %>%
    filter(sample == nm)
  
  if (nrow(d) == 0) {
    stop("No rows found in coords_plot for sample: ", nm)
  }
  
  d <- d %>%
    mutate(
      predicted_flower_celltype = factor(
        predicted_flower_celltype,
        levels = all_celltypes
      )
    )
  
  # -----------------------------------------------------------------------
  # Sample-specific coordinate limits
  # This prevents large empty white space.
  # -----------------------------------------------------------------------
  
  xr <- range(d$x, na.rm = TRUE)
  yr <- range(d$y, na.rm = TRUE)
  
  xpad <- 0.04 * diff(xr)
  ypad <- 0.04 * diff(yr)
  
  if (!is.finite(xpad) || xpad == 0) xpad <- 1
  if (!is.finite(ypad) || ypad == 0) ypad <- 1
  
  # -----------------------------------------------------------------------
  # LEFT PANEL: all inferred Flower cell states
  # -----------------------------------------------------------------------
  
  p_left <- ggplot(
    d,
    aes(
      x = x,
      y = y,
      colour = predicted_flower_celltype
    )
  ) +
    geom_point(
      size = 0.75,
      alpha = 0.90,
      stroke = 0
    ) +
    scale_colour_manual(
      values = celltype_colours,
      limits = all_celltypes,
      drop = FALSE,
      name = "Predicted\ncell type"
    ) +
    coord_equal(
      xlim = c(xr[1] - xpad, xr[2] + xpad),
      ylim = c(yr[1] - ypad, yr[2] + ypad),
      expand = FALSE
    ) +
    labs(
      title = "All inferred Flower cell states",
      subtitle = "Flower snRNA \u2192 MERFISH label transfer",
      x = NULL,
      y = NULL
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 13
      ),
      plot.subtitle = element_text(
        size = 9
      ),
      legend.title = element_text(
        size = 8
      ),
      legend.text = element_text(
        size = 7
      ),
      legend.key.height = unit(
        0.38,
        "cm"
      ),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  # -----------------------------------------------------------------------
  # RIGHT PANEL: high-confidence PDE1-associated states only
  # -----------------------------------------------------------------------
  
  foreground <- d %>%
    filter(
      transfer_score >= 0.5,
      predicted_flower_celltype %in% top_pde_types
    )
  
  background <- d %>%
    filter(
      !(
        transfer_score >= 0.5 &
          predicted_flower_celltype %in% top_pde_types
      )
    )
  
  foreground <- foreground %>%
    mutate(
      predicted_flower_celltype = factor(
        predicted_flower_celltype,
        levels = all_celltypes
      )
    )
  
  message("Total cells: ", nrow(d))
  message("High-confidence PDE1-associated cells: ", nrow(foreground))
  
  p_right <- ggplot() +
    geom_point(
      data = background,
      aes(
        x = x,
        y = y
      ),
      colour = "grey88",
      size = 0.65,
      alpha = 0.45,
      stroke = 0
    ) +
    geom_point(
      data = foreground,
      aes(
        x = x,
        y = y,
        colour = predicted_flower_celltype
      ),
      size = 1.05,
      alpha = 0.98,
      stroke = 0
    ) +
    scale_colour_manual(
      values = celltype_colours,
      limits = all_celltypes,
      drop = TRUE,
      name = "PDE1-associated\ncell state"
    ) +
    coord_equal(
      xlim = c(xr[1] - xpad, xr[2] + xpad),
      ylim = c(yr[1] - ypad, yr[2] + ypad),
      expand = FALSE
    ) +
    labs(
      title = "PDE1-associated cell states",
      subtitle = "Only assignments with prediction score \u2265 0.5",
      x = NULL,
      y = NULL
    ) +
    theme_void(base_size = 11) +
    theme(
      plot.title = element_text(
        face = "bold",
        size = 13
      ),
      plot.subtitle = element_text(
        size = 9
      ),
      legend.title = element_text(
        size = 8
      ),
      legend.text = element_text(
        size = 7
      ),
      legend.key.height = unit(
        0.38,
        "cm"
      ),
      plot.margin = margin(5, 5, 5, 5)
    )
  
  # -----------------------------------------------------------------------
  # Combine panels
  # -----------------------------------------------------------------------
  
  final_flower <- (
    p_left | p_right
  ) +
    plot_layout(
      widths = c(1, 1)
    ) +
    plot_annotation(
      title = paste0(
        panel_letters[[nm]],
        " | ",
        sample_labels[[nm]],
        ": Flower snRNA\u2013MERFISH spatial inference"
      ),
      subtitle = "Left: all transferred Flower cell states. Right: high-confidence PDE1-associated cell states.",
      caption = "PDE1 was not directly targeted by Flower MERFISH. Spatial localization is inferred from Flower snRNA-to-MERFISH label transfer using shared measured genes.",
      theme = theme(
        plot.title = element_text(
          face = "bold",
          size = 16
        ),
        plot.subtitle = element_text(
          size = 10
        ),
        plot.caption = element_text(
          size = 8,
          colour = "grey35"
        )
      )
    )
  
  print(final_flower)
  
  # -----------------------------------------------------------------------
  # Save PNG + PDF
  # -----------------------------------------------------------------------
  
  ggsave(
    filename = file.path(
      FIG_OUT,
      paste0(output_names[[nm]], ".png")
    ),
    plot = final_flower,
    width = 15,
    height = 7.5,
    dpi = 500,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(
      FIG_OUT,
      paste0(output_names[[nm]], ".pdf")
    ),
    plot = final_flower,
    width = 15,
    height = 7.5,
    bg = "white"
  )
  
  message(
    "Saved: ",
    file.path(
      FIG_OUT,
      paste0(output_names[[nm]], ".png")
    )
  )
}

message("\nFinal Flower figures completed.")
# =============================================================================
# Done
# =============================================================================

message("\n============================================================")
message("DONE: final developmental/spatial report figures generated.")
message("Output directory:")
message(FIG_DIR)
message("============================================================")