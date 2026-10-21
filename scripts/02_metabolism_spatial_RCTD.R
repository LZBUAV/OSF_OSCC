# =============================================================================
# Step 1b: metabolic reprogramming + spatial deconvolution (RCTD)
# Reproduces: Fig 2H (metabolic pathway GSEA per cell type), Fig 2I (spatial co-localisation),
#             Fig S2 (RCTD maps of all spatial samples)
# Needs: analysis/rds/01_sc_integrated.rds (from Step 1)
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(Matrix); library(data.table); library(ggplot2); library(dplyr)
  library(patchwork); library(pheatmap); library(fgsea); library(msigdbr); library(spacexr); library(png)
})
set.seed(1234)

ROOT <- Sys.getenv("OSF_OSCC_ROOT", normalizePath(".."))
ST   <- file.path(ROOT, "data/spatial")
OUT  <- file.path(ROOT, "output")
FIG  <- file.path(OUT, "figures"); TAB <- file.path(OUT, "tables"); RDS <- file.path(OUT, "rds")
fwrite <- function(x, file, ...) tryCatch(data.table::fwrite(x, file, ...), error = function(e) {
  alt <- sub("\\.csv$", "_new.csv", file); warning("File locked, writing ", alt); data.table::fwrite(x, alt, ...)
})
save_plot <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

obj <- readRDS(file.path(RDS, "01_sc_integrated.rds"))
ann_file <- file.path(TAB, "Step1_celltype_manual.csv")   # optional manual corrections (cluster,celltype)
if (file.exists(ann_file)) {
  man <- fread(ann_file); obj$celltype <- factor(man$celltype[match(obj$seurat_clusters, man$cluster)])
}

# -----------------------------------------------------------------------------
# Fig 2H: metabolic pathway GSEA per cell type, OSCC vs OSF
#   (original: GSEA of metabolism pathways per cell type, MASH vs MASLD)
# -----------------------------------------------------------------------------
kegg <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
met_kw <- "METABOLISM|GLYCOLYSIS|OXIDATIVE_PHOSPHORYLATION|CITRATE_CYCLE|BIOSYNTHESIS|DEGRADATION|PENTOSE|FATTY_ACID|GLUTATHIONE"
kegg <- kegg[grepl(met_kw, kegg$gs_name), ]
met_sets <- split(kegg$gene_symbol, kegg$gs_name)
met_sets <- met_sets[lengths(met_sets) >= 10]

met_gsea <- function(col, g1, g2, tag, title) {
  res <- list()
  for (ct in levels(obj$celltype)) {
    cells <- colnames(obj)[obj$celltype == ct & obj@meta.data[[col]] %in% c(g1, g2)]
    n <- table(factor(obj@meta.data[cells, col], c(g1, g2)))
    if (min(n) < 20) { message(tag, ": skip ", ct, " (", paste(n, collapse = "/"), " cells)"); next }
    sub <- subset(obj, cells = cells); Idents(sub) <- col
    fc <- FoldChange(sub, ident.1 = g1, ident.2 = g2)
    rk <- sort(setNames(fc$avg_log2FC, rownames(fc)), decreasing = TRUE)
    g  <- fgsea(met_sets, rk, minSize = 10, maxSize = 500)
    res[[ct]] <- data.frame(celltype = ct, pathway = g$pathway, NES = g$NES, padj = g$padj)
  }
  met <- bind_rows(res); met$pathway <- gsub("_", " ", sub("^KEGG_", "", met$pathway))
  fwrite(met, file.path(TAB, paste0("Step1b_metabolic_GSEA_", tag, ".csv")))
  top_path <- met |> group_by(pathway) |> summarise(m = min(padj, na.rm = TRUE)) |>
    arrange(m) |> slice_head(n = 30) |> pull(pathway)
  nes  <- reshape2::acast(subset(met, pathway %in% top_path), pathway ~ celltype, value.var = "NES")
  padj <- reshape2::acast(subset(met, pathway %in% top_path), pathway ~ celltype, value.var = "padj")
  nes[is.na(nes)] <- 0; padj[is.na(padj)] <- 1
  star <- ifelse(padj < 0.001, "***", ifelse(padj < 0.01, "**", ifelse(padj < 0.05, "*", "")))
  for (ext in c("pdf", "png"))
    pheatmap(nes, display_numbers = star, number_color = "black", cluster_cols = FALSE,
             color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
             breaks = seq(-3, 3, length.out = 101), border_color = "grey70", fontsize_row = 8,
             main = title, angle_col = 45,
             filename = file.path(FIG, paste0("Fig2H_Heatmap_metabolic_pathways_", tag, ".", ext)), width = 9.5, height = 8)
}
if (!file.exists(file.path(FIG, "Fig2H_Heatmap_metabolic_pathways_OSF_OSCC_vs_deNovo.png")) || Sys.getenv("REDO_2H") == "1") {
  met_gsea("stage", "OSCC", "OSF", "OSCC_vs_OSF", "Metabolic pathways: OSCC vs OSF (NES)")
  met_gsea("group", "OSF_OSCC", "deNovo_OSCC", "OSF_OSCC_vs_deNovo", "Metabolic pathways: OSF-OSCC vs de novo OSCC (NES)")
}

# -----------------------------------------------------------------------------
# RCTD reference from the scRNA atlas (max 3000 cells per cell type)
# -----------------------------------------------------------------------------
cells <- unlist(lapply(split(colnames(obj), obj$celltype), function(x) sample(x, min(3000, length(x)))))
ref_counts <- LayerData(obj, assay = "RNA", layer = "counts")[, cells]
ct <- droplevels(obj$celltype[cells]); names(ct) <- cells
reference <- Reference(ref_counts, ct, nUMI = Matrix::colSums(ref_counts))
rm(ref_counts); invisible(gc())

# -----------------------------------------------------------------------------
# GSE220978 Visium (4 OSF-associated OSCC patients)
# -----------------------------------------------------------------------------
read_visium <- function(pdir, gsm, pt) {
  pre <- file.path(pdir, paste0(gsm, "_", pt, "_"))
  m <- ReadMtx(mtx = paste0(pre, "matrix.mtx.gz"), cells = paste0(pre, "barcodes.tsv.gz"),
               features = paste0(pre, "features.tsv.gz"))
  img_dir <- file.path(tempdir(), pt); dir.create(img_dir, showWarnings = FALSE)
  untar(paste0(pre, "image.tar.gz"), exdir = img_dir)
  idir <- file.path(img_dir, paste0(pt, "_image"))
  pos <- fread(file.path(idir, "tissue_positions_list.csv"), header = FALSE,
               col.names = c("barcode", "in_tissue", "row", "col", "pxl_row", "pxl_col"))
  sf  <- jsonlite::fromJSON(file.path(idir, "scalefactors_json.json"))
  pos <- pos[in_tissue == 1 & barcode %in% colnames(m)]
  list(counts = m[, pos$barcode], pos = pos, img = readPNG(file.path(idir, "tissue_hires_image.png")),
       scalef = sf$tissue_hires_scalef)
}

spatial_plot <- function(v, value, title, limits = NULL, pal = c("grey95", "#FDE725", "#E64B35", "#7E0000")) {
  d <- data.frame(x = v$pos$pxl_col * v$scalef, y = v$pos$pxl_row * v$scalef, value = value)
  h <- dim(v$img)[1]; w <- dim(v$img)[2]
  ggplot(d, aes(x, y, colour = value)) +
    annotation_raster(v$img, xmin = 0, xmax = w, ymin = -h, ymax = 0) +
    geom_point(aes(y = -y), size = 0.9, alpha = 0.9) +
    scale_colour_gradientn(colours = pal, limits = limits) +
    coord_fixed(xlim = range(d$x) + c(-20, 20), ylim = -rev(range(d$y)) + c(-20, 20)) +
    theme_void() + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, size = 10))
}

pts <- data.frame(pt = paste0("Patient", 1:4), gsm = paste0("GSM68334", 84:87))
key_ct <- intersect(c("Myeloid", "Fibroblast", "Epithelial"), levels(ct))
rctd_w <- list(); coloc <- list(); maps <- list()
for (i in seq_len(nrow(pts))) {
  pt <- pts$pt[i]; message("RCTD ", pt)
  v <- read_visium(file.path(ST, "GSE220978", pt), pts$gsm[i], pt)
  coords <- data.frame(x = v$pos$pxl_col, y = v$pos$pxl_row, row.names = v$pos$barcode)
  puck <- SpatialRNA(coords, v$counts, nUMI = Matrix::colSums(v$counts))
  rctd <- create.RCTD(puck, reference, max_cores = 8)
  rctd <- run.RCTD(rctd, doublet_mode = "full")
  w <- as.matrix(normalize_weights(rctd@results$weights))
  w <- w[v$pos$barcode[v$pos$barcode %in% rownames(w)], ]
  v$pos <- v$pos[match(rownames(w), v$pos$barcode)]
  rctd_w[[pt]] <- data.frame(sample = pt, barcode = rownames(w), w, check.names = FALSE)

  # Fig 2I / S2: weight maps of key cell types
  maps[[pt]] <- wrap_plots(lapply(key_ct, function(k) spatial_plot(v, w[, k], paste(pt, "-", k), c(0, max(w[, k])))), nrow = 1)

  # Fig S2: dominant cell type per spot
  dom <- colnames(w)[max.col(w)]
  d <- data.frame(x = v$pos$pxl_col * v$scalef, y = v$pos$pxl_row * v$scalef, celltype = dom)
  pS2 <- ggplot(d, aes(x, -y, colour = celltype)) +
    annotation_raster(v$img, xmin = 0, xmax = dim(v$img)[2], ymin = -dim(v$img)[1], ymax = 0) +
    geom_point(size = 0.9) + coord_fixed() + theme_void() + ggtitle(paste(pt, "- RCTD dominant cell type")) +
    guides(colour = guide_legend(override.aes = list(size = 3)))
  save_plot(pS2, paste0("FigS2_RCTD_dominant_celltype_", pt), 6.5, 5.5)

  ctst <- cor.test(w[, "Myeloid"], w[, "Fibroblast"], method = "spearman")
  coloc[[pt]] <- data.frame(sample = pt, rho_Myeloid_Fibroblast = ctst$estimate, p = ctst$p.value)
}
save_plot(wrap_plots(maps, ncol = 1), "Fig2I_Spatial_RCTD_Myeloid_Fibroblast_Epithelial", 12, 14)
save_plot(maps[[1]], "Fig2I_Spatial_RCTD_Patient1", 12, 4)
fwrite(bind_rows(rctd_w), file.path(TAB, "Step1b_RCTD_weights_GSE220978.csv"))
fwrite(bind_rows(coloc), file.path(TAB, "Step1b_spatial_colocalisation_Myeloid_Fibroblast.csv"))

# Co-localisation matrix of RCTD weights across all spots (supports Fig 2I)
allw <- bind_rows(rctd_w)[, -(1:2)]
cm <- cor(allw, method = "spearman")
for (ext in c("pdf", "png"))
  pheatmap(cm, display_numbers = TRUE, number_format = "%.2f",
           color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100), breaks = seq(-1, 1, length.out = 101),
           main = "Spatial co-localisation (Spearman, RCTD weights)",
           filename = file.path(FIG, paste0("FigS2_Colocalisation_matrix_GSE220978.", ext)), width = 6.5, height = 6)

saveRDS(list(reference = reference, weights = rctd_w), file.path(RDS, "02_RCTD_GSE220978.rds"))
cat("Step 1b done\n")
