# =============================================================================
# Step 6: ligand-high vs ligand-low key macrophages, and spatial co-localisation with target fibroblasts
# Reproduces: Fig 7A-E, Fig S9 (other spatial samples), Fig S10 (ligand-receptor spatial co-expression),
#             Supplementary Table 12
# Needs: 04_myeloid.rds, 03_fibroblast.rds, tables/Step5_axis_selected.csv, Visium GSE220978
# =============================================================================
source("utils.R")
suppressPackageStartupMessages({ library(fgsea); library(png) })

key  <- readLines(file.path(TAB, "Step3_key_subset.txt"))[1]
axis <- fread(file.path(TAB, "Step5_axis_selected.csv"))
lig <- strsplit(axis$ligand, "_")[[1]][1]; rec <- strsplit(axis$receptor, "_")[[1]][1]; target <- axis$target
message("Axis: ", axis$TF, " -> ", lig, " -> ", rec, " | ", key, " -> ", target)

mye <- readRDS(file.path(RDS, "04_myeloid.rds"))
km <- subset(mye, subtype == key & stage == "OSCC"); rm(mye); invisible(gc())

# Fig 7A: samples split at the 75th percentile of mean ligand expression in key macrophages
smean <- tapply(LayerData(km, layer = "data")[lig, ], km$sample, mean)
ncell <- table(km$sample); smean <- smean[names(ncell)[ncell >= 10]]
cut75 <- quantile(smean, 0.75)
km <- subset(km, sample %in% names(smean)); km$lig_group <- ifelse(smean[km$sample] >= cut75, "High", "Low")
fwrite(data.frame(sample = names(smean), mean_ligand = smean, group = ifelse(smean >= cut75, "High", "Low")),
       file.path(TAB, "Step6_samples_ligand_high_low.csv"))
Idents(km) <- "lig_group"
de <- FindMarkers(km, ident.1 = "High", ident.2 = "Low", logfc.threshold = 0, min.pct = 0.1, verbose = FALSE)
de$gene <- rownames(de); fwrite(de, file.path(TAB, "SuppTable12_DEG_keymacro_ligand_high_vs_low.csv"))
save_plot(volcano_plot(de, paste0(key, ": ", lig, "-high vs ", lig, "-low samples")), "Fig7A_Volcano_ligand_high_vs_low", 6.5, 5.5)

# Fig 7B / 7C: GSEA (KEGG)
k2g <- kegg_sets(); ks <- split(k2g$gene, k2g$term)
rk <- sort(setNames(de$avg_log2FC, de$gene), decreasing = TRUE)
gs <- fgsea(ks, rk, minSize = 10, maxSize = 500); gs <- gs[order(gs$padj), ]
fwrite(gs[, c("pathway", "NES", "pval", "padj", "size")], file.path(TAB, "Step6_GSEA_KEGG_ligand_high_vs_low.csv"))
top <- rbind(head(gs[gs$NES > 0, ], 10), head(gs[gs$NES < 0, ], 10)); top$pathway <- factor(top$pathway, top$pathway[order(top$NES)])
p7b <- ggplot(top, aes(NES, pathway, fill = -log10(padj))) + geom_col() + scale_fill_gradient(low = "#91D1C2", high = "#E64B35") +
  labs(y = NULL, title = paste("GSEA (KEGG):", lig, "-high vs low")) + theme_classic() + umap_theme
save_plot(p7b, "Fig7B_GSEA_KEGG_bar", 8, 6)
p_top <- as.character(gs$pathway[gs$NES > 0][1])
p7c <- plotEnrichment(ks[[p_top]], rk) + labs(title = sprintf("%s\nNES = %.2f, padj = %.2g", p_top, gs$NES[gs$pathway == p_top], gs$padj[gs$pathway == p_top])) +
  theme_classic() + umap_theme
save_plot(p7c, "Fig7C_GSEA_enrichment_plot", 6, 4.5)
cc_path <- grep("COMPLEMENT", gs$pathway, value = TRUE)[1]
if (!is.na(cc_path)) save_plot(plotEnrichment(ks[[cc_path]], rk) + ggtitle(sprintf("%s (NES %.2f, padj %.2g)", cc_path, gs$NES[gs$pathway == cc_path], gs$padj[gs$pathway == cc_path])) +
                                 theme_classic() + umap_theme, "Fig7C_supp_GSEA_complement_coagulation", 6, 4.5)

# -----------------------------------------------------------------------------
# Spatial: key-macrophage vs target-fibroblast signature scores (Visium GSE220978)
# -----------------------------------------------------------------------------
sig <- readRDS(file.path(RDS, "04b_macrophage_subset_signatures.rds"))[[key]]
fib <- readRDS(file.path(RDS, "03_fibroblast.rds")); Idents(fib) <- "subtype"
fm <- FindMarkers(fib, ident.1 = target, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.25, max.cells.per.ident = 2000, verbose = FALSE)
fsig <- head(rownames(fm)[fm$pct.2 < 0.5 & !grepl("^(RP[SL]|MT-)", rownames(fm))][order(-fm$avg_log2FC[fm$pct.2 < 0.5 & !grepl("^(RP[SL]|MT-)", rownames(fm))])], 30)
rm(fib); invisible(gc())

read_visium <- function(pt, gsm) {
  pdir <- file.path(ST, "GSE220978", pt); pre <- file.path(pdir, paste0(gsm, "_", pt, "_"))
  m <- ReadMtx(paste0(pre, "matrix.mtx.gz"), paste0(pre, "barcodes.tsv.gz"), paste0(pre, "features.tsv.gz"))
  idir <- file.path(tempdir(), pt); untar(paste0(pre, "image.tar.gz"), exdir = idir); idir <- file.path(idir, paste0(pt, "_image"))
  pos <- fread(file.path(idir, "tissue_positions_list.csv"), header = FALSE, col.names = c("barcode", "in_tissue", "row", "col", "pxl_row", "pxl_col"))
  pos <- pos[in_tissue == 1 & barcode %in% colnames(m)]
  m <- m[, pos$barcode]; lc <- log1p(t(t(m) / pmax(Matrix::colSums(m), 1)) * 1e4)
  list(lc = lc, pos = pos, img = readPNG(file.path(idir, "tissue_hires_image.png")), sf = jsonlite::fromJSON(file.path(idir, "scalefactors_json.json"))$tissue_hires_scalef)
}
score <- function(lc, g) { g <- intersect(g, rownames(lc)); z <- t(scale(t(as.matrix(lc[g, , drop = FALSE])))); z[is.na(z)] <- 0; colMeans(z) }
smap <- function(v, val, title) {
  d <- data.frame(x = v$pos$pxl_col * v$sf, y = v$pos$pxl_row * v$sf, value = val)
  ggplot(d, aes(x, -y, colour = value)) + annotation_raster(v$img, 0, dim(v$img)[2], -dim(v$img)[1], 0) +
    geom_point(size = 0.8) + scale_colour_gradientn(colours = c("#3C5488", "white", "#E64B35"), name = NULL) +
    coord_fixed() + theme_void() + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5, size = 9))
}
pts <- data.frame(pt = paste0("Patient", 1:4), gsm = paste0("GSM68334", 84:87))
cors <- list(); maps <- list(); lr_maps <- list(); scat <- list()
for (i in 1:4) {
  v <- read_visium(pts$pt[i], pts$gsm[i])
  s1 <- score(v$lc, sig); s2 <- score(v$lc, fsig)
  ct <- cor.test(s1, s2, method = "spearman")
  cors[[i]] <- data.frame(sample = pts$pt[i], rho = ct$estimate, p = ct$p.value, n_spots = length(s1))
  maps[[i]] <- smap(v, s1, paste(pts$pt[i], "-", key)) + smap(v, s2, paste(pts$pt[i], "-", target))
  scat[[i]] <- ggplot(data.frame(a = s1, b = s2), aes(a, b)) + geom_point(size = 0.4, alpha = 0.4, colour = "grey40") +
    geom_smooth(method = "lm", colour = "#E64B35") + labs(x = paste(key, "score"), y = paste(target, "score"),
    title = sprintf("%s: rho = %.2f, p = %.1e", pts$pt[i], ct$estimate, ct$p.value)) + theme_classic() + umap_theme
  le <- if (lig %in% rownames(v$lc)) v$lc[lig, ] else rep(0, ncol(v$lc)); re <- if (rec %in% rownames(v$lc)) v$lc[rec, ] else rep(0, ncol(v$lc))
  lr_maps[[i]] <- smap(v, le, paste(pts$pt[i], lig)) + smap(v, re, paste(pts$pt[i], rec)) + smap(v, sqrt(le * re), paste(lig, "x", rec, "co-expression"))
}
cors <- bind_rows(cors); fwrite(cors, file.path(TAB, "Step6_spatial_keymacro_targetfib_correlation.csv"))
best_i <- which.max(cors$rho)
save_plot(maps[[best_i]], "Fig7D_Spatial_keymacro_targetfib", 10, 5)
save_plot(scat[[best_i]], "Fig7E_Scatter_spatial_correlation", 5, 4.5)
save_plot(wrap_plots(maps[-best_i], ncol = 1) | wrap_plots(scat[-best_i], ncol = 1), "FigS9_Spatial_colocalisation_other_samples", 15, 13)
save_plot(wrap_plots(lr_maps, ncol = 1), "FigS10_Spatial_ligand_receptor_coexpression", 14, 17)
print(cors)
cat("Step 6 done\n")
