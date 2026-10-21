# =============================================================================
# Step 1: scRNA-seq integration, QC, clustering, major cell-type annotation
# Reproduces: Fig 2A-G, Supplementary Table 1 (samples/QC), Supplementary Table 3 (proportions)
# Run:  "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" 01_scRNA_integration.R
# Output: analysis/figures/Fig2*, analysis/tables/*, analysis/rds/01_sc_integrated.rds
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(Matrix); library(data.table)
  library(ggplot2); library(dplyr); library(patchwork); library(pheatmap)
})
options(future.globals.maxSize = 50 * 1024^3)
set.seed(1234)

ROOT <- Sys.getenv("OSF_OSCC_ROOT", normalizePath(".."))
SC   <- file.path(ROOT, "data/scRNA")
ST   <- file.path(ROOT, "data/spatial")
OUT  <- file.path(ROOT, "output")
FIG  <- file.path(OUT, "figures"); TAB <- file.path(OUT, "tables"); RDS <- file.path(OUT, "rds")
for (d in c(FIG, TAB, RDS)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# fwrite that does not crash when the target file is open in Excel (writes *_new.csv instead)
fwrite <- function(x, file, ...) tryCatch(data.table::fwrite(x, file, ...), error = function(e) {
  alt <- sub("\\.csv$", "_new.csv", file); warning("File locked, writing ", alt); data.table::fwrite(x, alt, ...)
})

save_plot <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

# -----------------------------------------------------------------------------
# 1. Sample sheet
#    stage  : Normal / OSF / OSCC                         (Fig 2E, Ro/e across stages)
#    group  : Normal / OSF / OSF_OSCC / deNovo_OSCC / OSCC_NOS
#             OSF_OSCC  = OSCC with concomitant OSMF (GSE215403 P8, P10, P12; Kurkalang 2023)
#             deNovo_OSCC = the other 9 GSE215403 tumours (OSMF absent per the paper)
#             OSCC_NOS  = OSCC with unknown OSF status (GSE195832 pre-treatment, GSE310797)
# -----------------------------------------------------------------------------
d262 <- file.path(SC, "GSE262446"); d215 <- file.path(SC, "GSE215403", "OSCC")
d195 <- file.path(SC, "GSE195832", "scRNA_Pre"); d310 <- file.path(ST, "GSE310797", "scRNA")

ss <- rbind(
  data.frame(dataset = "GSE262446", sample = "Normal_LB", gsm = "GSM8169899", patient = "LB",
             group = "Normal", dir = file.path(d262, "Control"), prefix = "GSM8169899_LB_", type = "10x"),
  data.frame(dataset = "GSE262446", sample = "OSF_691nm", gsm = "GSM8169898", patient = "6-9-1nm",
             group = "OSF", dir = file.path(d262, "OSF"), prefix = "GSM8169898_6-9-1nm_", type = "10x"),
  data.frame(dataset = "GSE215403",
             sample  = paste0("GB_P", 1:12),
             gsm     = paste0("GSM66348", 69:80),
             patient = paste0("P", 1:12),
             group   = ifelse(1:12 %in% c(8, 10, 12), "OSF_OSCC", "deNovo_OSCC"),
             dir     = d215,
             prefix  = paste0("GSM66348", 69:80, "_",
                              c("OSCC", "scB1", "scB2", "scB5", "scB7", "scB8", "scB9",
                                "scB10", "scB12", "scB13", "scB14", "scB15"), "_"),
             type = "10x"),
  data.frame(dataset = "GSE195832", sample = paste0("VU_Pt", 1:4, "_Pre"),
             gsm = c("GSM5851565", "GSM5851567", "GSM5851569", "GSM5851571"),
             patient = paste0("Pt", 1:4), group = "OSCC_NOS", dir = d195,
             prefix = c("GSM5851565_3219-SR-1", "GSM5851567_3219-SR-3",
                        "GSM5851569_3219-SR-5", "GSM5851571_3521-SR-1"), type = "10x"),
  data.frame(dataset = "GSE310797", sample = paste0("CS_PT", 1:16),
             gsm = paste0("GSM93092", 14:29), patient = paste0("PT", 1:16), group = "OSCC_NOS",
             dir = d310, prefix = paste0("GSM93092", 14:29, "_PT", 1:16, "_"), type = "dense")
)
ss$stage <- ifelse(ss$group %in% c("Normal", "OSF"), ss$group, "OSCC")
# GSM9309226_PT13_count.txt.gz is truncated on GEO itself (only 9,577 of 51,795 genes) -> excluded
ss <- ss[ss$sample != "CS_PT13", ]

# -----------------------------------------------------------------------------
# 2. Readers. 10x files are read by Ensembl ID and renamed with the GRCh38-2020-A
#    symbols of GSE215403 so that old (GSE262446/GSE195832) and new annotations agree.
# -----------------------------------------------------------------------------
ref_feat <- fread(file.path(d215, "GSM6634869_OSCC_features.tsv.gz"), header = FALSE)
ens2sym  <- setNames(make.unique(ref_feat$V2), ref_feat$V1)

read_10x <- function(dir, prefix) {
  m <- ReadMtx(mtx      = file.path(dir, paste0(prefix, "matrix.mtx.gz")),
               cells    = file.path(dir, paste0(prefix, "barcodes.tsv.gz")),
               features = file.path(dir, paste0(prefix, "features.tsv.gz")),
               feature.column = 1)
  m <- m[rownames(m) %in% names(ens2sym), ]
  rownames(m) <- ens2sym[rownames(m)]
  m
}

read_dense <- function(dir, prefix, chunk = 2000) {
  dt <- fread(file.path(dir, paste0(prefix, "count.txt.gz")))
  genes <- make.unique(dt[[1]]); dt[, 1 := NULL]
  idx <- split(seq_len(ncol(dt)), ceiling(seq_len(ncol(dt)) / chunk))
  m <- do.call(cbind, lapply(idx, function(j) as(as.matrix(dt[, ..j]), "dgCMatrix")))
  rownames(m) <- genes
  m
}

mad_out <- function(x, n = 3, side = c("both", "higher")) {
  side <- match.arg(side); med <- median(x); md <- mad(x)
  if (side == "higher") x > med + n * md else x < med - n * md | x > med + n * md
}

# -----------------------------------------------------------------------------
# 3. Load every sample, per-sample MAD-based QC (3 MADs, as in the original paper)
# -----------------------------------------------------------------------------
f_clu <- file.path(RDS, "01_sc_clustered.rds")   # checkpoint: delete it to recompute from raw data
if (!file.exists(f_clu)) {
mats <- list(); qc <- list()
for (i in seq_len(nrow(ss))) {
  s <- ss[i, ]; message("Reading ", s$sample)
  m <- if (s$type == "10x") read_10x(s$dir, s$prefix) else read_dense(s$dir, s$prefix)
  colnames(m) <- paste0(s$sample, "_", colnames(m))
  nC <- Matrix::colSums(m); nF <- Matrix::colSums(m > 0)
  mt <- Matrix::colSums(m[grepl("^MT-", rownames(m)), , drop = FALSE]) / pmax(nC, 1) * 100
  keep <- nF >= 200 & !(mad_out(log10(nC + 1)) | mad_out(log10(nF + 1)) | mad_out(mt, side = "higher"))
  qc[[s$sample]] <- data.frame(sample = s$sample, cells_raw = ncol(m), cells_pass = sum(keep),
                               median_nCount = median(nC[keep]), median_nFeature = median(nF[keep]),
                               median_pct_mt = round(median(mt[keep]), 2))
  mats[[s$sample]] <- m[, keep]
}

common_genes <- Reduce(intersect, lapply(mats, rownames))
message("Genes shared by all datasets: ", length(common_genes))

objs <- lapply(names(mats), function(n) {
  o <- CreateSeuratObject(mats[[n]][common_genes, ], project = n)
  meta <- ss[ss$sample == n, c("dataset", "sample", "gsm", "patient", "group", "stage")]
  for (k in names(meta)) o[[k]] <- meta[[k]]
  o
})
rm(mats); invisible(gc())

obj <- merge(objs[[1]], objs[-1]); rm(objs); invisible(gc())
obj <- JoinLayers(obj)
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
obj$stage <- factor(obj$stage, levels = c("Normal", "OSF", "OSCC"))
obj$group <- factor(obj$group, levels = c("Normal", "OSF", "OSF_OSCC", "deNovo_OSCC", "OSCC_NOS"))

# Supplementary Table 1: sample information + QC
supp1 <- merge(ss[, c("dataset", "gsm", "sample", "patient", "group", "stage")],
               do.call(rbind, qc), by = "sample")
fwrite(supp1, file.path(TAB, "SuppTable1_scRNA_samples_QC.csv"))

# -----------------------------------------------------------------------------
# 4. Normalisation, HVGs (3000), PCA, UMAP before / after Harmony
# -----------------------------------------------------------------------------
obj <- NormalizeData(obj) |> FindVariableFeatures(nfeatures = 3000) |> ScaleData() |>
  RunPCA(npcs = 50, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, reduction.name = "umap.unintegrated", verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = "sample", reduction.use = "pca", dims.use = 1:30, verbose = FALSE)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30) |> FindClusters(resolution = 0.5)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap", verbose = FALSE)
saveRDS(obj, f_clu)
} else obj <- readRDS(f_clu)

# -----------------------------------------------------------------------------
# 5. Major cell-type annotation by canonical markers (cluster-level module scores)
# -----------------------------------------------------------------------------
markers <- list(
  Epithelial  = c("KRT5", "KRT14", "KRT17", "KRT6A", "SFN"),
  Fibroblast  = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "PLVAP"),
  Pericyte_SMC = c("RGS5", "ACTA2", "MYH11", "PDGFRB", "TAGLN"),
  Myeloid     = c("CD68", "CD14", "LYZ", "CD163", "C1QA"),
  T_NK        = c("CD3D", "CD3E", "CD2", "NKG7", "GNLY"),
  B_cell      = c("MS4A1", "CD79A", "CD19", "BANK1"),
  Plasma      = c("MZB1", "JCHAIN", "IGHG1", "SDC1"),
  Mast        = c("TPSAB1", "CPA3", "KIT", "MS4A2")
)
markers <- lapply(markers, intersect, rownames(obj))
obj <- AddModuleScore(obj, features = markers, name = "MS_")
sc_cols <- paste0("MS_", seq_along(markers))
cl_score <- aggregate(obj@meta.data[, sc_cols], list(cluster = obj$seurat_clusters), mean)
colnames(cl_score)[-1] <- names(markers)
best <- apply(cl_score[, -1], 1, function(x) names(markers)[which.max(x)])
margin <- apply(cl_score[, -1], 1, function(x) { y <- sort(x, decreasing = TRUE); y[1] - y[2] })
anno <- data.frame(cluster = cl_score$cluster, celltype = best, score_margin = round(margin, 3),
                   ambiguous = margin < 0.1)
fwrite(cbind(anno, round(cl_score[, -1], 3)), file.path(TAB, "Step1_cluster_annotation_scores.csv"))
obj$celltype <- factor(anno$celltype[match(obj$seurat_clusters, anno$cluster)], levels = names(markers))
obj@meta.data[, sc_cols] <- NULL

if (!file.exists(file.path(TAB, "Step1_cluster_markers.csv"))) {
  cl_markers <- FindAllMarkers(obj, only.pos = TRUE, max.cells.per.ident = 500,
                               logfc.threshold = 0.5, min.pct = 0.25, verbose = FALSE)
  fwrite(cl_markers, file.path(TAB, "Step1_cluster_markers.csv"))
}

# Manual curation after reviewing cluster markers (Step1_celltype_manual.csv); "Excluded" = doublets / low quality
man_file <- file.path(TAB, "Step1_celltype_manual.csv")
if (file.exists(man_file)) {
  man <- fread(man_file)
  obj$celltype <- man$celltype[match(as.character(obj$seurat_clusters), as.character(man$cluster))]
  obj <- subset(obj, celltype != "Excluded")
  obj$celltype <- factor(obj$celltype, levels = names(markers))
}

# -----------------------------------------------------------------------------
# 6. Figures
# -----------------------------------------------------------------------------
pal_ds    <- c(GSE262446 = "#E64B35", GSE215403 = "#4DBBD5", GSE195832 = "#00A087", GSE310797 = "#3C5488")
pal_stage <- c(Normal = "#E64B35", OSF = "#00A087", OSCC = "#3C5488")
pal_group <- c(Normal = "#E64B35", OSF = "#00A087", OSF_OSCC = "#F39B7F", deNovo_OSCC = "#8491B4", OSCC_NOS = "#B09C85")
pal_ct    <- setNames(c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
                        "#91D1C2", "#DC0000", "#7E6148"), names(markers))

umap_theme <- theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# UMAP drawn directly with ggplot2 (Seurat::DimPlot legends break under ggplot2 >= 4.0)
umap_plot <- function(obj, red = "umap", by, cols = NULL, title = NULL, label = FALSE, legend = TRUE) {
  d <- data.frame(Embeddings(obj, red)[, 1:2], g = obj@meta.data[[by]])
  colnames(d)[1:2] <- c("UMAP_1", "UMAP_2"); d <- d[sample(nrow(d)), ]
  p <- ggplot(d, aes(UMAP_1, UMAP_2, colour = g)) + geom_point(size = 0.05, stroke = 0, shape = 16) +
    guides(colour = guide_legend(override.aes = list(size = 3))) + labs(colour = NULL, title = title) +
    theme_classic() + umap_theme
  p <- p + if (is.null(cols)) scale_colour_hue() else scale_colour_manual(values = cols)
  if (label) {
    cen <- aggregate(d[, 1:2], list(g = d$g), median)
    p <- p + ggrepel::geom_text_repel(data = cen, aes(label = g), colour = "black", size = 4, fontface = "bold")
  }
  if (!legend) p <- p + theme(legend.position = "none")
  p
}

# Fig 2A / 2B: before and after batch correction
save_plot(umap_plot(obj, "umap.unintegrated", "dataset", pal_ds, "Before batch correction"), "Fig2A_UMAP_before_Harmony", 6.5, 5)
save_plot(umap_plot(obj, "umap", "dataset", pal_ds, "After Harmony integration"), "Fig2B_UMAP_after_Harmony", 6.5, 5)

# Fig 2C: disease status
save_plot(umap_plot(obj, "umap", "stage", pal_stage, "Disease status"), "Fig2C_UMAP_disease_status", 6.5, 5)
save_plot(umap_plot(obj, "umap", "group", pal_group, "Sample group"), "Fig2C_supp_UMAP_sample_group", 7, 5)

# Fig 2D: major cell types
save_plot(umap_plot(obj, "umap", "celltype", pal_ct, "Major cell types", label = TRUE), "Fig2D_UMAP_major_celltypes", 7, 5.5)
save_plot(umap_plot(obj, "umap", "seurat_clusters", NULL, "Clusters", label = TRUE, legend = FALSE), "Step1_check_UMAP_clusters", 6, 5.5)

# Ro/e (observed / expected), same definition as the original paper
calc_roe <- function(meta, ct = "celltype", grp = "stage") {
  tab <- table(meta[[ct]], meta[[grp]]); tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  exp <- outer(rowSums(tab), colSums(tab)) / sum(tab)
  list(roe = unclass(tab / exp), p = suppressWarnings(chisq.test(tab)$p.value), counts = unclass(tab))
}
roe_symbol <- function(x) ifelse(x >= 1.5, "+++", ifelse(x >= 1.2, "++", ifelse(x > 1, "+",
                          ifelse(x > 0.8, "+/-", ifelse(x > 0.2, "-", "---")))))
plot_roe <- function(r, name, title) {
  lab <- matrix(paste0(sprintf("%.2f", r$roe), "\n", roe_symbol(r$roe)), nrow = nrow(r$roe))
  args <- list(r$roe, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = lab,
               number_color = "black", fontsize_number = 9,
               color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
               breaks = seq(0, 2, length.out = 101), border_color = "grey60",
               main = sprintf("%s (chi-square p = %.2g)", title, r$p), angle_col = 0)
  for (ext in c("pdf", "png")) do.call(pheatmap, c(args, list(filename = file.path(FIG, paste0(name, ".", ext)),
                                                              width = 4.5, height = 5)))
}

# Fig 2E: Ro/e across Normal / OSF / OSCC
r2e <- calc_roe(obj@meta.data, "celltype", "stage")
plot_roe(r2e, "Fig2E_RoE_major_celltypes_by_stage", "Ro/e")
# Normal/OSF samples of GSE262446 are ~97% epithelial (epithelium-focused dissociation), so the
# stromal/immune compartment is also compared after removing epithelial cells
mne <- subset(obj@meta.data, celltype != "Epithelial"); mne$celltype <- droplevels(mne$celltype)
plot_roe(calc_roe(mne, "celltype", "stage"), "Fig2E_supp_RoE_nonepithelial_by_stage", "Ro/e (non-epithelial)")

# Fig 2F: Ro/e OSF-associated OSCC vs de novo OSCC (both from GSE215403 -> no batch confound)
meta215 <- subset(obj@meta.data, dataset == "GSE215403")
meta215$group <- droplevels(meta215$group)
r2f <- calc_roe(meta215, "celltype", "group")
plot_roe(r2f, "Fig2F_RoE_OSF_OSCC_vs_deNovo_OSCC", "Ro/e")

# Fig 2G: dot plot of canonical markers
top_mk <- unlist(lapply(markers, head, 3))
p2g <- DotPlot(obj, features = unique(top_mk), group.by = "celltype", cols = c("lightgrey", "#E64B35")) +
  RotatedAxis() + theme(axis.title = element_blank()) + ggtitle("Canonical marker genes") + umap_theme
save_plot(p2g, "Fig2G_DotPlot_canonical_markers", 11, 4.5)

# Supplementary Table 3: cell-type proportions per stage, group and sample
prop_tab <- function(by) {
  t <- as.data.frame.matrix(table(obj$celltype, obj@meta.data[[by]]))
  p <- sweep(t, 2, colSums(t), "/") * 100
  cbind(celltype = rownames(t), setNames(t, paste0("n_", colnames(t))), setNames(round(p, 2), paste0("pct_", colnames(p))))
}
fwrite(prop_tab("stage"),  file.path(TAB, "SuppTable3_celltype_proportion_by_stage.csv"))
fwrite(prop_tab("group"),  file.path(TAB, "SuppTable3_celltype_proportion_by_group.csv"))
fwrite(prop_tab("sample"), file.path(TAB, "SuppTable3_celltype_proportion_by_sample.csv"))
fwrite(data.frame(r2e$roe, check.names = FALSE), file.path(TAB, "Step1_RoE_by_stage.csv"), row.names = TRUE)
fwrite(data.frame(r2f$roe, check.names = FALSE), file.path(TAB, "Step1_RoE_OSF_OSCC_vs_deNovo.csv"), row.names = TRUE)

# -----------------------------------------------------------------------------
# 7. Save
# -----------------------------------------------------------------------------
saveRDS(obj, file.path(RDS, "01_sc_integrated.rds"))
writeLines(capture.output(sessionInfo()), file.path(OUT, "logs", "01_sessionInfo.txt"))
cat("\nDone. Cells:", ncol(obj), " Genes:", nrow(obj), "\n")
print(table(obj$celltype, obj$stage))
