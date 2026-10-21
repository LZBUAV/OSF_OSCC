# Shared paths and helpers, sourced by Step 2 onwards
suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2); library(dplyr)
  library(patchwork); library(pheatmap); library(ggrepel)
})
options(future.globals.maxSize = 50 * 1024^3)
set.seed(1234)

ROOT <- Sys.getenv("OSF_OSCC_ROOT", normalizePath(".."))
SC   <- file.path(ROOT, "data/scRNA")
ST   <- file.path(ROOT, "data/spatial")
BULK <- file.path(ROOT, "data/bulk")
METH <- file.path(ROOT, "data/methylation")
OUT  <- file.path(ROOT, "output")
FIG  <- file.path(OUT, "figures"); TAB <- file.path(OUT, "tables"); RDS <- file.path(OUT, "rds")

fwrite <- function(x, file, ...) tryCatch(data.table::fwrite(x, file, ...), error = function(e) {
  alt <- sub("\\.csv$", "_new.csv", file); warning("File locked, writing ", alt); data.table::fwrite(x, alt, ...)
})

save_plot <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h, limitsize = FALSE)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300, limitsize = FALSE)
}

# pheatmap centres the title over the whole canvas and clips it when the canvas is too narrow:
# widen the canvas so that the full title always fits (~0.09 inch per character at default font size)
save_pheatmap <- function(name, w, h, ...) {
  main <- list(...)$main
  if (!is.null(main)) w <- max(w, 0.14 * max(nchar(strsplit(main, "\n")[[1]])) + 1.5)
  for (ext in c("pdf", "png")) pheatmap(..., filename = file.path(FIG, paste0(name, ".", ext)), width = w, height = h)
}

pal_stage <- c(Normal = "#E64B35", OSF = "#00A087", OSCC = "#3C5488")
pal_group <- c(Normal = "#E64B35", OSF = "#00A087", OSF_OSCC = "#F39B7F", deNovo_OSCC = "#8491B4", OSCC_NOS = "#B09C85")
pal_many  <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000",
               "#7E6148", "#B09C85", "#FFDC91", "#6F99AD", "#EE4C97", "#20854E")
umap_theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold"))

umap_plot <- function(obj, red = "umap", by, cols = NULL, title = NULL, label = FALSE, legend = TRUE, size = 0.1) {
  d <- data.frame(Embeddings(obj, red)[, 1:2], g = relab(obj@meta.data[[by]]))
  if (!is.null(cols)) names(cols) <- relab(names(cols))
  colnames(d)[1:2] <- c("UMAP_1", "UMAP_2"); d <- d[sample(nrow(d)), ]
  p <- ggplot(d, aes(UMAP_1, UMAP_2, colour = g)) + geom_point(size = size, stroke = 0, shape = 16) +
    guides(colour = guide_legend(override.aes = list(size = 3))) + labs(colour = NULL, title = title) +
    theme_classic() + umap_theme
  p <- p + if (is.null(cols)) scale_colour_hue() else scale_colour_manual(values = cols)
  if (label) {
    cen <- aggregate(d[, 1:2], list(g = d$g), median)
    p <- p + geom_text_repel(data = cen, aes(label = g), colour = "black", size = 3.5, fontface = "bold")
  }
  if (!legend) p <- p + theme(legend.position = "none")
  p
}

# ---- display names used in all figures (code names stay unchanged in objects and tables) ----
ct_label <- c(Epithelial = "Epithelial cell", Fibroblast = "Fibroblast", Endothelial = "Endothelial cell",
              Pericyte_SMC = "Pericyte/SMC", Myeloid = "Myeloid cell", T_NK = "T/NK cell", B_cell = "B cell",
              Plasma = "Plasma cell", Mast = "Mast cell")
group_label <- c(Normal = "Normal", OSF = "OSF", OSF_OSCC = "OSF-associated OSCC", deNovo_OSCC = "De novo OSCC",
                 OSCC_NOS = "OSCC (OSF status unknown)", OSCC = "OSCC", OSF_dysplasia = "OSF with dysplasia")
cohort_label <- c(GSE125866_CV = "GSE125866 (CV)", TCGA_oral = "TCGA oral cavity", Mean = "Mean")
relab <- function(x, map = c(ct_label, group_label, cohort_label)) {
  y <- as.character(x); hit <- y %in% names(map); y[hit] <- map[y[hit]]
  if (is.factor(x)) factor(y, levels = unique(relab(levels(x), map))) else y
}
relab_dimnames <- function(m) { rownames(m) <- relab(rownames(m)); colnames(m) <- relab(colnames(m)); m }
pal_group_lab <- setNames(pal_group, relab(names(pal_group)))
pal_ct_lab <- setNames(c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148"),
                       unname(ct_label))

fmt_p <- function(p) if (is.na(p)) "P = NA" else if (p < 2.2e-16) "P < 2.2e-16" else sprintf("P = %.2g", p)

calc_roe <- function(meta, ct, grp) {
  tab <- table(meta[[ct]], meta[[grp]]); tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  exp <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  list(roe = unclass(tab / exp), p = suppressWarnings(chisq.test(tab)$p.value), counts = unclass(tab))
}
roe_symbol <- function(x) ifelse(x >= 1.5, "+++", ifelse(x >= 1.2, "++", ifelse(x > 1, "+",
                          ifelse(x > 0.8, "+/-", ifelse(x > 0.2, "-", "---")))))
plot_roe <- function(r, name, title = "Ro/e", w = 4.5, h = 5) {
  r$roe <- relab_dimnames(r$roe)
  lab <- matrix(paste0(sprintf("%.2f", r$roe), "\n", roe_symbol(r$roe)), nrow = nrow(r$roe))
  save_pheatmap(name, w, h, r$roe, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = lab,
                number_color = "black", fontsize_number = 9,
                color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
                breaks = seq(0, 2, length.out = 101), border_color = "grey60",
                main = sprintf("%s\nchi-square %s", title, fmt_p(r$p)), angle_col = 0)
}

# Cluster-level module-score annotation; returns the label per cluster
annotate_clusters <- function(obj, markers, cluster_col = "seurat_clusters", min_margin = 0.05) {
  markers <- lapply(markers, intersect, rownames(obj)); markers <- markers[lengths(markers) > 0]
  obj <- AddModuleScore(obj, features = markers, name = "AMS_")
  sc <- paste0("AMS_", seq_along(markers))
  s <- aggregate(obj@meta.data[, sc, drop = FALSE], list(cluster = obj@meta.data[[cluster_col]]), mean)
  colnames(s)[-1] <- names(markers)
  best <- apply(s[, -1, drop = FALSE], 1, function(x) names(markers)[which.max(x)])
  mar  <- apply(s[, -1, drop = FALSE], 1, function(x) { y <- sort(x, TRUE); y[1] - y[2] })
  data.frame(cluster = s$cluster, label = best, margin = round(mar, 3), ambiguous = mar < min_margin, round(s[, -1], 3))
}

volcano_plot <- function(de, title, fc = 0.5, p = 0.05, n_lab = 10) {
  de$sig <- ifelse(de$p_val_adj < p & de$avg_log2FC > fc, "Up", ifelse(de$p_val_adj < p & de$avg_log2FC < -fc, "Down", "NS"))
  de$y <- -log10(pmax(de$p_val_adj, 1e-300))
  ambient <- "^(HB[ABDGMQZ]|IG[HKL]|JCHAIN|MT-|RP[SL]|GNLY|NKG7)"   # ambient RNA, not labelled
  lab <- bind_rows(de |> filter(sig == "Up", !grepl(ambient, gene)) |> arrange(desc(avg_log2FC)) |> head(n_lab),
                   de |> filter(sig == "Down", !grepl(ambient, gene)) |> arrange(avg_log2FC) |> head(n_lab))
  ggplot(de, aes(avg_log2FC, y, colour = sig)) + geom_point(size = 0.8, alpha = 0.7) +
    scale_colour_manual(values = c(Up = "#E64B35", Down = "#3C5488", NS = "grey75")) +
    geom_vline(xintercept = c(-fc, fc), linetype = 2, colour = "grey50") +
    geom_hline(yintercept = -log10(p), linetype = 2, colour = "grey50") +
    geom_text_repel(data = lab, aes(label = gene), size = 3, max.overlaps = 30, colour = "black") +
    labs(x = "log2 fold change", y = "-log10 adjusted P", colour = NULL, title = title) + theme_classic() + umap_theme
}

kegg_sets <- function() {
  k <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
  data.frame(term = gsub("_", " ", sub("^KEGG_", "", k$gs_name)), gene = k$gene_symbol)
}
hallmark_sets <- function() {
  k <- msigdbr::msigdbr(species = "Homo sapiens", collection = "H")
  split(k$gene_symbol, gsub("_", " ", sub("^HALLMARK_", "", k$gs_name)))
}

load_atlas <- function() {
  obj <- readRDS(file.path(RDS, "01_sc_integrated.rds"))
  man <- file.path(TAB, "Step1_celltype_manual.csv")
  if (file.exists(man)) { m <- data.table::fread(man); obj$celltype <- factor(m$celltype[match(obj$seurat_clusters, m$cluster)]) }
  obj
}

# ---- bulk loaders + cached TCGA oral-cavity cohort (used by 09_machine_learning.R and 11_methylation.R) ----
source(file.path(ROOT, "scripts", "bulk_loaders.R"))
get_TCGA_oral <- function() {
  f <- file.path(RDS, "bulk_TCGA_oral.rds")
  if (!file.exists(f)) saveRDS(load_TCGA_oral(), f)
  readRDS(f)
}
for (d in c(FIG, TAB, RDS, file.path(OUT, "logs"))) dir.create(d, showWarnings = FALSE, recursive = TRUE)