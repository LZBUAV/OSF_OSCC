# Redraw figures with publication names instead of code names (deNovo_OSCC, OSF_OSCC, OSCC_NOS, T_NK, ...)
# Figures: 2C-supp, 2D, 2G, 4A-supp, 4D, 4E, 6B, 6C, 6D, S2 (dominant cell type maps + co-localisation matrix)
source("utils.R")
suppressPackageStartupMessages({ library(png); library(CellChat) })
key <- readLines(file.path(TAB, "Step3_key_subset.txt"))[1]

# ---------------- atlas UMAPs and dot plot ----------------
obj <- readRDS(file.path(RDS, "01_sc_integrated.rds"))
save_plot(umap_plot(obj, "umap", "group", pal_group, "Sample group"), "Fig2C_supp_UMAP_sample_group", 7.5, 5)
save_plot(umap_plot(obj, "umap", "celltype", setNames(unname(pal_ct_lab), names(ct_label)), "Major cell types", label = TRUE),
          "Fig2D_UMAP_major_celltypes", 7.5, 5.5)
markers <- list(Epithelial = c("KRT5", "KRT14", "KRT17"), Fibroblast = c("COL1A1", "COL1A2", "DCN"),
                Endothelial = c("PECAM1", "VWF", "CDH5"), Pericyte_SMC = c("RGS5", "ACTA2", "MYH11"),
                Myeloid = c("CD68", "CD14", "LYZ"), T_NK = c("CD3D", "CD3E", "CD2"), B_cell = c("MS4A1", "CD79A", "CD19"),
                Plasma = c("MZB1", "JCHAIN", "IGHG1"), Mast = c("TPSAB1", "CPA3", "KIT"))
obj$celltype_lab <- relab(obj$celltype)
p2g <- DotPlot(obj, features = unique(unlist(markers)), group.by = "celltype_lab", cols = c("lightgrey", "#E64B35")) +
  RotatedAxis() + theme(axis.title = element_blank()) + ggtitle("Canonical marker genes") + umap_theme
save_plot(p2g, "Fig2G_DotPlot_canonical_markers", 11, 4.5)

# data for Fig 6D (IL1B / IL1R1 per cell in GSE215403)
o215 <- subset(obj, dataset == "GSE215403"); rm(obj); invisible(gc())
fib <- readRDS(file.path(RDS, "03_fibroblast.rds")); mye <- readRDS(file.path(RDS, "04_myeloid.rds"))
lab <- setNames(as.character(o215$celltype), colnames(o215))
i <- intersect(names(lab), colnames(fib)); lab[i] <- as.character(fib$subtype[i])
i <- intersect(names(lab), colnames(mye)); lab[i] <- as.character(mye$subtype[i])
vd <- data.frame(cell_group = lab, group = relab(as.character(o215$group)),
                 IL1B = LayerData(o215, layer = "data")["IL1B", ], IL1R1 = LayerData(o215, layer = "data")["IL1R1", ])
fib_types <- levels(fib$subtype); rm(o215, fib); invisible(gc())

# ---------------- myeloid figures ----------------
save_plot(umap_plot(mye, "umap", "group", pal_group, "Sample group", size = 0.3), "Fig4A_supp_UMAP_myeloid_group", 7.5, 5.5)
sigs <- c("M1_polarization", "M2_polarization", "Hypoxia", "Glycolysis", "TCA_cycle", "OXPHOS", "Lipid_mediator",
          "Pro_inflammatory", "Interferon", "Angiogenesis", "ECM_remodeling")
sig_lab <- c(M1_polarization = "M1 polarization", M2_polarization = "M2 polarization", Hypoxia = "Hypoxia", Glycolysis = "Glycolysis",
             TCA_cycle = "TCA cycle", OXPHOS = "Oxidative phosphorylation", Lipid_mediator = "Lipid mediators",
             Pro_inflammatory = "Pro-inflammatory", Interferon = "Interferon response", Angiogenesis = "Angiogenesis",
             ECM_remodeling = "ECM remodelling")
ks <- subset(mye, subtype == key & dataset == "GSE215403"); rm(mye); invisible(gc())
avg <- aggregate(ks@meta.data[, sigs], list(sample = ks$sample), mean)
grp <- relab(as.character(ks$group[match(avg$sample, ks$sample)]))
mat <- t(scale(as.matrix(avg[, -1]))); colnames(mat) <- sub("GB_P", "Patient ", avg$sample); rownames(mat) <- sig_lab[rownames(mat)]
ann <- data.frame(Group = grp, row.names = colnames(mat)); o <- order(ann$Group != "OSF-associated OSCC")
save_pheatmap("Fig4D_Heatmap_key_macrophage_signatures", 8, 5, mat[, o], cluster_cols = FALSE, annotation_col = ann[o, , drop = FALSE],
              annotation_colors = list(Group = pal_group_lab[c("OSF-associated OSCC", "De novo OSCC")]),
              color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100), breaks = seq(-2, 2, length.out = 101),
              main = paste0(key, "s: pathway scores per patient (z-score)"), border_color = "grey70", angle_col = 90)
d4e <- reshape2::melt(ks@meta.data[, c("group", sigs)], id.vars = "group")
d4e$group <- factor(relab(as.character(d4e$group)), c("OSF-associated OSCC", "De novo OSCC"))
d4e$variable <- factor(sig_lab[as.character(d4e$variable)], sig_lab)
p4e <- ggplot(d4e, aes(group, value, fill = group)) + geom_boxplot(outlier.size = 0.2, width = 0.6) +
  facet_wrap(~variable, scales = "free_y", nrow = 2) + scale_fill_manual(values = pal_group_lab, name = NULL) +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) +
  labs(x = NULL, y = "UCell score", title = paste0(key, "s: OSF-associated vs de novo OSCC")) +
  theme_classic() + umap_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "bottom")
save_plot(p4e, "Fig4E_Boxplot_key_macrophage_signatures", 13, 6.5)

# ---------------- Fig 6B: differential interaction heat maps (own drawing, readable names) ----------------
cc <- readRDS(file.path(RDS, "06_cellchat_list.rds"))
g <- Reduce(intersect, lapply(cc, function(x) rownames(x@net$count)))
dl <- lapply(c(count = "count", weight = "weight"), function(m) cc$OSF_OSCC@net[[m]][g, g] - cc$deNovo_OSCC@net[[m]][g, g])
hm6 <- function(m, title) {
  d <- as.data.frame(as.table(relab_dimnames(m))); colnames(d) <- c("Sender", "Receiver", "diff")
  lim <- max(abs(d$diff))
  ggplot(d, aes(Receiver, Sender, fill = diff)) + geom_tile(colour = "white", linewidth = 0.2) +
    scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", limits = c(-lim, lim), name = "OSF-assoc.\nminus\nde novo") +
    labs(title = title, x = "Receiver", y = "Sender") + theme_minimal() + umap_theme +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8), axis.text.y = element_text(size = 8), panel.grid = element_blank())
}
save_plot(hm6(dl$count, "Differential number of interactions") + hm6(dl$weight, "Differential interaction strength"),
          "Fig6B_Heatmap_differential_interactions", 16, 7.5)

# ---------------- Fig 6C: L-R pairs from the key subset in both groups ----------------
lr_all <- bind_rows(lapply(names(cc), function(n) {
  d <- subsetCommunication(cc[[n]], sources.use = key, targets.use = intersect(fib_types, levels(cc[[n]]@idents)), thresh = 1)
  d$group <- n; d
}))
lrw <- fread(file.path(TAB, "Step5_keysubset_to_fibroblast_LR_OSF_vs_deNovo.csv"))
top <- head(unique(lrw[delta > 0][order(-delta)]$interaction_name), 25)
d6 <- subset(lr_all, interaction_name %in% top)
d6$pair <- paste(gsub("_", "/", d6$ligand), "->", gsub("_", "/", d6$receptor))
d6$pair <- factor(d6$pair, rev(unique(d6$pair[match(top, d6$interaction_name)])))
d6$group <- factor(relab(d6$group), c("De novo OSCC", "OSF-associated OSCC"))
d6$target <- factor(d6$target, intersect(fib_types, unique(d6$target)))
d6$sig <- ifelse(d6$pval < 0.01, "P < 0.01", ifelse(d6$pval < 0.05, "P < 0.05", "n.s."))
p6c <- ggplot(subset(d6, prob > 0), aes(target, pair)) +
  geom_point(aes(size = sig, colour = prob)) + facet_wrap(~group) +
  scale_colour_gradientn(colours = c("#4DBBD5", "#FDE725", "#E64B35"), name = "Communication\nprobability") +
  scale_size_manual(values = c(`P < 0.01` = 4, `P < 0.05` = 3, n.s. = 1), name = NULL) +
  labs(x = "Receiver fibroblast subtype", y = NULL, title = paste("L-R pairs from", key, "to fibroblast subtypes")) +
  theme_bw() + umap_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"))
save_plot(p6c, "Fig6C_Bubble_LR_keysubset_to_fibroblasts", 11, 8)

# ---------------- Fig 6D: violin plots ----------------
vd$group <- factor(vd$group, c("OSF-associated OSCC", "De novo OSCC"))
mk_v <- function(gene, types, title) {
  d <- subset(vd, cell_group %in% types); d$cell_group <- factor(d$cell_group, types)
  ggplot(d, aes(cell_group, .data[[gene]], fill = group)) + geom_violin(scale = "width", position = position_dodge(0.8), linewidth = 0.2) +
    ggpubr::stat_compare_means(aes(group = group), label = "p.signif", method = "wilcox.test") +
    scale_fill_manual(values = pal_group_lab, name = NULL) + labs(x = NULL, y = paste(gene, "expression"), title = title) +
    theme_classic() + umap_theme + theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")
}
mye_types <- intersect(c("FCN1+ Monocyte", "FOLR2+ Macrophage", "MMP12+ Macrophage", "SPP1+ Macrophage", "CXCL5+ Macrophage",
                         "IL10+ Macrophage", "KYNU+ Macrophage", "Prolif. Macrophage", "cDC2", "mregDC", "pDC", "Neutrophil"), unique(vd$cell_group))
save_plot(mk_v("IL1B", mye_types, "IL1B in myeloid populations") / mk_v("IL1R1", intersect(fib_types, unique(vd$cell_group)), "IL1R1 in fibroblast subtypes"),
          "Fig6D_Violin_ligand_receptor", 11, 9.5)

# ---------------- Supplementary Fig 2: dominant cell type maps + co-localisation matrix ----------------
rw <- fread(file.path(TAB, "Step1b_RCTD_weights_GSE220978.csv"))
cts <- setdiff(names(rw), c("sample", "barcode"))
for (i in 1:4) {
  pt <- paste0("Patient", i); gsm <- paste0("GSM68334", 83 + i)
  pre <- file.path(ST, "GSE220978", pt, paste0(gsm, "_", pt, "_")); idir <- file.path(tempdir(), pt)
  untar(paste0(pre, "image.tar.gz"), exdir = idir); idir <- file.path(idir, paste0(pt, "_image"))
  pos <- fread(file.path(idir, "tissue_positions_list.csv"), header = FALSE, col.names = c("barcode", "in_tissue", "row", "col", "pxl_row", "pxl_col"))
  img <- readPNG(file.path(idir, "tissue_hires_image.png")); sf <- jsonlite::fromJSON(file.path(idir, "scalefactors_json.json"))$tissue_hires_scalef
  w <- rw[sample == pt]; d <- merge(w, pos, by = "barcode")
  d$celltype <- factor(relab(cts[max.col(as.matrix(d[, ..cts]))]), unname(ct_label))
  p <- ggplot(d, aes(pxl_col * sf, -pxl_row * sf, colour = celltype)) + annotation_raster(img, 0, dim(img)[2], -dim(img)[1], 0) +
    geom_point(size = 0.9) + scale_colour_manual(values = pal_ct_lab, name = "Dominant\ncell type", drop = TRUE) + coord_fixed() + theme_void() +
    ggtitle(paste0("Patient ", i, ": dominant cell type per spot (RCTD)")) +
    guides(colour = guide_legend(override.aes = list(size = 3))) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(p, paste0("FigS2_RCTD_dominant_celltype_Patient", i), 7.5, 5.5)
}
cm <- relab_dimnames(cor(rw[, ..cts], method = "spearman"))
save_pheatmap("FigS2_Colocalisation_matrix_GSE220978", 7, 6.5, cm, display_numbers = TRUE, number_format = "%.2f",
              color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100), breaks = seq(-1, 1, length.out = 101),
              main = "Spatial co-localisation of RCTD proportions\n(Spearman, 17,298 spots)")
cat("relabelled figures done\n")
