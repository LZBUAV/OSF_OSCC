# =============================================================================
# Step 3: myeloid subpopulations and the key OSF-associated macrophage subset
# Reproduces: Fig 4A-G, Supplementary Table 6 (myeloid markers)
# Needs: analysis/rds/01_sc_integrated.rds
# Main contrast (oral counterpart of MASLD vs MASH): OSF-associated OSCC vs de novo OSCC (GSE215403)
# =============================================================================
source("utils.R")
suppressPackageStartupMessages({ library(harmony); library(UCell) })

f_mye <- file.path(RDS, "04_myeloid.rds")
if (!file.exists(f_mye)) {
  obj <- load_atlas()
  mye <- subset(obj, celltype == "Myeloid"); rm(obj); invisible(gc())
  mye <- NormalizeData(mye) |> FindVariableFeatures(nfeatures = 2000) |> ScaleData() |> RunPCA(npcs = 30, verbose = FALSE)
  mye <- RunHarmony(mye, group.by.vars = "sample", dims.use = 1:20, verbose = FALSE)
  mye <- FindNeighbors(mye, reduction = "harmony", dims = 1:20) |> FindClusters(resolution = 0.6)
  mye <- RunUMAP(mye, reduction = "harmony", dims = 1:20, verbose = FALSE)
  saveRDS(mye, f_mye)
} else mye <- readRDS(f_mye)

# -----------------------------------------------------------------------------
# Annotation: lineage class by module score, cluster named "<top specific marker>+ <class>"
# (same naming logic as the original paper, e.g. "DTNA+ macrophages", "VCAN+ monocytes")
# -----------------------------------------------------------------------------
mye_markers <- list(
  Monocyte    = c("VCAN", "FCN1", "S100A8", "S100A9", "CD14"),
  Macrophage  = c("C1QA", "C1QB", "C1QC", "APOE", "CD163", "MRC1"),
  cDC1        = c("CLEC9A", "XCR1", "IRF8", "CADM1"),
  cDC2        = c("CD1C", "FCER1A", "CLEC10A", "CD1E"),
  mregDC      = c("LAMP3", "CCR7", "FSCN1", "CCL19"),
  pDC         = c("LILRA4", "CLEC4C", "IL3RA", "GZMB"),
  Neutrophil  = c("CSF3R", "FCGR3B", "CXCL8", "G0S2", "AQP9"),
  Proliferating = c("MKI67", "TOP2A", "CENPF", "NUSAP1")
)
an <- annotate_clusters(mye, mye_markers)
Idents(mye) <- "seurat_clusters"
mk_cl <- FindAllMarkers(mye, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.25, max.cells.per.ident = 1000, verbose = FALSE)
spec <- mk_cl |> filter(pct.1 >= 0.3, pct.2 <= 0.25, !grepl("^(RP[SL]|MT-|HLA-|IG[HKL])", gene)) |>
  group_by(cluster) |> slice_max(avg_log2FC, n = 1, with_ties = FALSE)
an$top_gene <- spec$gene[match(an$cluster, spec$cluster)]
an$name <- ifelse(an$label %in% c("cDC1", "cDC2", "mregDC", "pDC", "Neutrophil"), an$label,
                  paste0(an$top_gene, "+ ", ifelse(an$label == "Proliferating", "Prolif. Macro", an$label)))
man <- file.path(TAB, "Step3_myeloid_subtype_manual.csv")         # optional manual override (cluster,name)
if (file.exists(man)) { m <- fread(man); an$name <- m$name[match(an$cluster, m$cluster)] }
fwrite(an, file.path(TAB, "Step3_myeloid_cluster_annotation.csv"))
mye$subtype <- an$name[match(as.character(mye$seurat_clusters), as.character(an$cluster))]
mye <- subset(mye, subtype != "Excluded")
an <- an[an$name != "Excluded", ]
mye$subtype <- factor(mye$subtype, levels = unique(an$name[order(an$label, an$name)]))
pal_mye <- setNames(rep(pal_many, 3)[seq_along(levels(mye$subtype))], levels(mye$subtype))

# Fig 4A / 4B + Supplementary Table 6
save_plot(umap_plot(mye, "umap", "subtype", pal_mye, "Myeloid subpopulations", label = TRUE, size = 0.3),
          "Fig4A_UMAP_myeloid_subtypes", 8, 6)
save_plot(umap_plot(mye, "umap", "group", pal_group, "Sample group", size = 0.3), "Fig4A_supp_UMAP_myeloid_group", 7, 5.5)
Idents(mye) <- "subtype"
mk <- FindAllMarkers(mye, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.2, max.cells.per.ident = 1000, verbose = FALSE)
fwrite(mk, file.path(TAB, "SuppTable6_myeloid_subtype_markers.csv"))
dot_genes <- unique(c("CD163", "MARCO", "C1QC", "APOE", "SPP1", "VCAN", "FCN1", "CLEC9A", "CD1C", "LAMP3",
                      "LILRA4", "CSF3R", "MKI67", na.omit(an$top_gene)))
p4b <- DotPlot(mye, features = intersect(dot_genes, rownames(mye)), cols = c("lightgrey", "#E64B35")) + RotatedAxis() +
  theme(axis.title = element_blank()) + ggtitle("Myeloid marker genes") + umap_theme
save_plot(p4b, "Fig4B_DotPlot_myeloid_markers", 13, 5)

# Fig 4C: Ro/e (main = OSF_OSCC vs deNovo in GSE215403; supp = by stage)
m215 <- subset(mye@meta.data, dataset == "GSE215403"); m215$group <- droplevels(m215$group)
r4c <- calc_roe(m215, "subtype", "group")
plot_roe(r4c, "Fig4C_RoE_myeloid_OSF_OSCC_vs_deNovo", w = 5, h = 6)
plot_roe(calc_roe(mye@meta.data, "subtype", "stage"), "Fig4C_supp_RoE_myeloid_by_stage", w = 5, h = 6)

# Per-patient proportions (patient-level test, guards against Ro/e being driven by one sample)
pp <- as.data.frame(table(m215$sample, m215$subtype)); colnames(pp) <- c("sample", "subtype", "n")
pp <- pp |> group_by(sample) |> mutate(pct = 100 * n / sum(n)) |> ungroup()
pp$group <- m215$group[match(pp$sample, m215$sample)]
ptest <- pp |> group_by(subtype) |> summarise(
  mean_OSF_OSCC = mean(pct[group == "OSF_OSCC"]), mean_deNovo = mean(pct[group == "deNovo_OSCC"]),
  p_wilcox = suppressWarnings(wilcox.test(pct[group == "OSF_OSCC"], pct[group == "deNovo_OSCC"])$p.value))
ptest$RoE_OSF_OSCC <- r4c$roe[as.character(ptest$subtype), "OSF_OSCC"]
fwrite(ptest, file.path(TAB, "Step3_myeloid_subtype_enrichment_OSF_OSCC.csv"))

# Key subset = macrophage subset most enriched in OSF-associated OSCC (can be overridden)
key_file <- file.path(TAB, "Step3_key_subset.txt")
macs <- setdiff(grep("Macrophage", as.character(ptest$subtype), value = TRUE),
                grep("Prolif", as.character(ptest$subtype), value = TRUE))   # proliferation = cell state, not identity
key <- if (file.exists(key_file)) readLines(key_file)[1] else
  as.character(ptest$subtype[ptest$subtype %in% macs][which.max(ptest$RoE_OSF_OSCC[ptest$subtype %in% macs])])
message("Key subset: ", key)
mye$is_key <- mye$subtype == key

# -----------------------------------------------------------------------------
# Signature scores (UCell)
# -----------------------------------------------------------------------------
H <- hallmark_sets(); K <- split(kegg_sets()$gene, kegg_sets()$term)
sigs <- list(
  M1_polarization  = c("IL1B", "IL6", "TNF", "CXCL9", "CXCL10", "CD80", "CD86", "IL12B", "CCL5", "IRF5", "STAT1", "SOCS3", "CCL3", "CCL4"),
  M2_polarization  = c("CD163", "MRC1", "MSR1", "CD209", "IL10", "TGFB1", "CCL18", "CCL22", "VEGFA", "MAF", "FOLR2", "LYVE1", "STAB1", "F13A1"),
  Hypoxia          = H[["HYPOXIA"]],
  Glycolysis       = H[["GLYCOLYSIS"]],
  TCA_cycle        = K[["CITRATE CYCLE TCA CYCLE"]],
  OXPHOS           = H[["OXIDATIVE PHOSPHORYLATION"]],
  Lipid_mediator   = c("PTGS2", "PTGES", "PTGES2", "ALOX5", "ALOX5AP", "LTA4H", "LTC4S", "PLA2G4A", "ALOX15", "HPGDS", "PTGDS", "TBXAS1"),
  Pro_inflammatory = H[["INFLAMMATORY RESPONSE"]],
  Interferon       = unique(c(H[["INTERFERON ALPHA RESPONSE"]], H[["INTERFERON GAMMA RESPONSE"]])),
  Angiogenesis     = H[["ANGIOGENESIS"]],
  ECM_remodeling   = c("MMP9", "MMP12", "MMP14", "SPP1", "FN1", "TGFBI", "LGMN", "CTSK", "CTSB", "TIMP1", "SPARC", "COL1A1")
)
mye <- AddModuleScore_UCell(mye, features = sigs, ncores = 4, name = "")

# Fig 4D: z-scored signature heatmap of the key subset per patient (GSE215403)
ks <- subset(mye, is_key & dataset == "GSE215403")
avg <- aggregate(ks@meta.data[, names(sigs)], list(sample = ks$sample), mean)
ann <- data.frame(group = ks$group[match(avg$sample, ks$sample)], row.names = avg$sample)
mat <- t(scale(as.matrix(avg[, -1]))); colnames(mat) <- avg$sample
mat <- mat[, order(ann$group)]
save_pheatmap("Fig4D_Heatmap_key_macrophage_signatures", 7, 5, mat, cluster_cols = FALSE,
              annotation_col = ann, annotation_colors = list(group = pal_group[c("OSF_OSCC", "deNovo_OSCC")]),
              color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100), breaks = seq(-2, 2, length.out = 101),
              main = paste(key, ": pathway z-scores per patient"), border_color = "grey70")

# Fig 4E: box plots of signature scores in the key subset, OSF_OSCC vs deNovo
d4e <- reshape2::melt(ks@meta.data[, c("group", names(sigs))], id.vars = "group")
d4e$group <- droplevels(d4e$group)
p4e <- ggplot(d4e, aes(group, value, fill = group)) + geom_boxplot(outlier.size = 0.2, width = 0.6) +
  facet_wrap(~variable, scales = "free_y", nrow = 2) + scale_fill_manual(values = pal_group) +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) +
  labs(x = NULL, y = "UCell score", title = paste(key, ": OSF-associated vs de novo OSCC")) +
  theme_classic() + umap_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p4e, "Fig4E_Boxplot_key_macrophage_signatures", 13, 6)

# Fig 4F: Hallmark activity across myeloid subsets (UCell, mean per subset, row z-score)
set.seed(1234)
cells <- unlist(lapply(split(colnames(mye), mye$subtype), function(x) sample(x, min(800, length(x)))))
ms <- subset(mye, cells = cells)
names(H) <- make.names(names(H))
ms <- AddModuleScore_UCell(ms, features = H, ncores = 4, name = "_H")
hcols <- paste0(names(H), "_H")
hm <- aggregate(ms@meta.data[, hcols], list(subtype = ms$subtype), mean)
hmat <- t(scale(as.matrix(hm[, -1]))); colnames(hmat) <- hm$subtype; rownames(hmat) <- gsub("\\.", " ", names(H))
save_pheatmap("Fig4F_Heatmap_Hallmark_myeloid_subsets", 9, 11, hmat,
              color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100), breaks = seq(-2, 2, length.out = 101),
              main = "Hallmark gene sets across myeloid subsets", angle_col = 45, border_color = "grey70")

# Fig 4G: irGSEA-style differential Hallmark activity, OSF_OSCC vs deNovo within each subset
#   (UCell scores + Wilcoxon per subset; irGSEA itself could not be installed on R 4.5)
ms215 <- subset(ms, dataset == "GSE215403")
g4 <- lapply(levels(droplevels(ms215$subtype)), function(s) {
  md <- ms215@meta.data[ms215$subtype == s, ]
  if (min(table(droplevels(md$group))) < 10 || length(unique(md$group)) < 2) return(NULL)
  do.call(rbind, lapply(hcols, function(h) {
    a <- md[md$group == "OSF_OSCC", h]; b <- md[md$group == "deNovo_OSCC", h]
    data.frame(subtype = s, pathway = gsub("\\.", " ", sub("_H$", "", h)), diff = mean(a) - mean(b),
               p = suppressWarnings(wilcox.test(a, b)$p.value))
  }))
})
g4 <- bind_rows(g4); g4$padj <- p.adjust(g4$p, "BH")
fwrite(g4, file.path(TAB, "Step3_irGSEA_style_Hallmark_OSF_OSCC_vs_deNovo.csv"))
top_pw <- g4 |> group_by(pathway) |> summarise(m = min(padj)) |> arrange(m) |> slice_head(n = 25) |> pull(pathway)
d4g <- subset(g4, pathway %in% top_pw)
d4g$sig <- cut(d4g$padj, c(-Inf, 1e-4, 1e-3, 1e-2, 0.05, Inf), labels = c("****", "***", "**", "*", ""))
p4g <- ggplot(d4g, aes(subtype, pathway, fill = diff)) + geom_tile(colour = "grey80") + geom_text(aes(label = sig), size = 3) +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", name = "Score diff\n(OSF-OSCC -\n de novo)") +
  labs(x = NULL, y = NULL, title = "Differential Hallmark activity: OSF-associated vs de novo OSCC") +
  theme_minimal() + umap_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(p4g, "Fig4G_irGSEA_Hallmark_OSF_OSCC_vs_deNovo", 10, 8)

saveRDS(mye, f_mye)
cat("Step 3 done. Key subset:", key, "\n"); print(ptest)
