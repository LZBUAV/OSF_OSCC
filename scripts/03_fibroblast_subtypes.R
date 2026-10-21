# =============================================================================
# Step 2: fibroblast / mural-cell subtypes (oral counterpart of the five HSC subtypes)
# Reproduces: Fig 3A-J, Supplementary Table 4 (subtype markers), Supplementary Table 5 (DEGs)
# Needs: analysis/rds/01_sc_integrated.rds
# =============================================================================
source("utils.R")
suppressPackageStartupMessages({ library(harmony); library(clusterProfiler) })

obj <- load_atlas()
fib <- subset(obj, celltype %in% c("Fibroblast", "Pericyte_SMC")); rm(obj); invisible(gc())

f_fib <- file.path(RDS, "03_fibroblast.rds")
if (!file.exists(f_fib)) {
  fib <- NormalizeData(fib) |> FindVariableFeatures(nfeatures = 2000) |> ScaleData() |> RunPCA(npcs = 30, verbose = FALSE)
  fib <- RunHarmony(fib, group.by.vars = "sample", dims.use = 1:20, verbose = FALSE)
  fib <- FindNeighbors(fib, reduction = "harmony", dims = 1:20) |> FindClusters(resolution = 0.4)
  fib <- RunUMAP(fib, reduction = "harmony", dims = 1:20, verbose = FALSE)
  saveRDS(fib, f_fib)
} else fib <- readRDS(f_fib)

# -----------------------------------------------------------------------------
# Subtype annotation (canonical oral fibroblast / CAF states)
# -----------------------------------------------------------------------------
fib_markers <- list(
  qFib     = c("PI16", "DPT", "MFAP5", "CD34", "SFRP2"),         # quiescent / universal fibroblast
  iCAF     = c("CXCL12", "IL6", "CXCL1", "C3", "CFD", "PLA2G2A"), # inflammatory
  myCAF    = c("POSTN", "COL11A1", "CTHRC1", "MMP11", "ACTA2", "COL10A1"),
  apCAF    = c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1"),       # antigen-presenting
  Pericyte = c("RGS5", "KCNJ8", "PDGFRB", "MCAM", "HIGD1B"),
  SMC      = c("MYH11", "DES", "CNN1", "ACTG2", "PLN")
)
an <- annotate_clusters(fib, fib_markers)
fwrite(an, file.path(TAB, "Step2_fibroblast_cluster_annotation.csv"))
man <- file.path(TAB, "Step2_fibroblast_subtype_manual.csv")   # optional manual override (cluster,label)
sub_levels <- names(fib_markers)
if (file.exists(man)) {
  m <- fread(man); an$label <- m$label[match(an$cluster, m$cluster)]
  sub_levels <- c("qFib", "iCAF", "POSTN+ myCAF", "MMP1+ myCAF", "Prolif. CAF", "Pericyte", "SMC")
}
fib$subtype <- an$label[match(as.character(fib$seurat_clusters), as.character(an$cluster))]
fib <- subset(fib, subtype != "Excluded")
fib$subtype <- droplevels(factor(fib$subtype, levels = sub_levels))
pal_sub <- setNames(pal_many[seq_along(levels(fib$subtype))], levels(fib$subtype))

# Fig 3A
save_plot(umap_plot(fib, "umap", "subtype", pal_sub, "Fibroblast subtypes", label = TRUE, size = 0.3),
          "Fig3A_UMAP_fibroblast_subtypes", 6.5, 5)
save_plot(umap_plot(fib, "umap", "seurat_clusters", NULL, "Fibroblast clusters", label = TRUE, legend = FALSE, size = 0.3),
          "Step2_check_UMAP_fibroblast_clusters", 6, 5)
save_plot(umap_plot(fib, "umap", "stage", pal_stage, "Disease status", size = 0.3), "Fig3A_supp_UMAP_fibroblast_stage", 6.5, 5)

# Fig 3B + Supplementary Table 4
Idents(fib) <- "subtype"
mk <- FindAllMarkers(fib, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.2, max.cells.per.ident = 1000, verbose = FALSE)
fwrite(mk, file.path(TAB, "SuppTable4_fibroblast_subtype_markers.csv"))
top <- mk |> group_by(cluster) |> arrange(desc(avg_log2FC), .by_group = TRUE) |> slice_head(n = 3) |> pull(gene)
dot_genes <- unique(c(unlist(lapply(fib_markers[levels(fib$subtype)], head, 2)), top))
p3b <- DotPlot(fib, features = dot_genes, cols = c("lightgrey", "#E64B35")) + RotatedAxis() +
  theme(axis.title = element_blank()) + ggtitle("Fibroblast subtype markers") + umap_theme
save_plot(p3b, "Fig3B_DotPlot_fibroblast_markers", 12, 4)

# Fig 3C: Ro/e across stages; Fig 3D: OSF-associated vs de novo OSCC (GSE215403)
plot_roe(calc_roe(fib@meta.data, "subtype", "stage"), "Fig3C_RoE_fibroblast_by_stage")
m215 <- subset(fib@meta.data, dataset == "GSE215403"); m215$group <- droplevels(m215$group)
plot_roe(calc_roe(m215, "subtype", "group"), "Fig3D_RoE_fibroblast_OSF_OSCC_vs_deNovo")

# Fig 3E: DEGs of fibroblasts. Main contrast = OSF-associated vs de novo OSCC (same batch, GSE215403,
# 3 vs 9 patients). OSCC vs OSF is supplementary: the OSF scRNA sample contains only ~46 fibroblasts.
Idents(fib) <- "group"
de <- FindMarkers(fib, ident.1 = "OSF_OSCC", ident.2 = "deNovo_OSCC", logfc.threshold = 0, min.pct = 0.1,
                  max.cells.per.ident = 3000, verbose = FALSE)
de$gene <- rownames(de); fwrite(de, file.path(TAB, "SuppTable5_fibroblast_DEG_OSF_OSCC_vs_deNovo.csv"))
save_plot(volcano_plot(de, "Fibroblasts: OSF-associated vs de novo OSCC"), "Fig3E_Volcano_fibroblast_OSF_OSCC_vs_deNovo", 6, 5.5)

Idents(fib) <- "stage"
de2 <- FindMarkers(fib, ident.1 = "OSCC", ident.2 = "OSF", logfc.threshold = 0, min.pct = 0.1,
                   max.cells.per.ident = 3000, verbose = FALSE)
de2$gene <- rownames(de2); fwrite(de2, file.path(TAB, "SuppTable5b_fibroblast_DEG_OSCC_vs_OSF.csv"))
save_plot(volcano_plot(de2, "Fibroblasts: OSCC vs OSF"), "Fig3E_supp_Volcano_OSCC_vs_OSF", 6, 5.5)

# Fig 3F: KEGG enrichment (offline MSigDB KEGG) of up / down DEGs
k2g <- kegg_sets()
enr_bar <- function(genes, title, col) {
  e <- enricher(genes, TERM2GENE = k2g, pvalueCutoff = 0.05, qvalueCutoff = 0.2, universe = rownames(fib))
  if (is.null(e) || nrow(as.data.frame(e)) == 0) return(NULL)
  d <- head(as.data.frame(e), 10); d$Description <- factor(d$Description, rev(d$Description))
  list(tab = as.data.frame(e), plot = ggplot(d, aes(-log10(p.adjust), Description)) + geom_col(fill = col, width = 0.7) +
         labs(x = "-log10 adjusted P", y = NULL, title = title) + theme_classic() + umap_theme)
}
up <- subset(de, p_val_adj < 0.05 & avg_log2FC > 0.5)$gene
dn <- subset(de, p_val_adj < 0.05 & avg_log2FC < -0.5)$gene
eu <- enr_bar(up, "KEGG: up in OSF-associated OSCC fibroblasts", "#E64B35")
ed <- enr_bar(dn, "KEGG: up in de novo OSCC fibroblasts", "#3C5488")
pl <- Filter(Negate(is.null), list(if (!is.null(ed)) ed$plot, if (!is.null(eu)) eu$plot))
if (length(pl)) save_plot(wrap_plots(pl, ncol = 1), "Fig3F_KEGG_fibroblast_DEGs", 8, 8)
if (!is.null(eu)) fwrite(eu$tab, file.path(TAB, "Step2_KEGG_up_in_OSF_OSCC_fibroblasts.csv"))
if (!is.null(ed)) fwrite(ed$tab, file.path(TAB, "Step2_KEGG_up_in_deNovo_OSCC_fibroblasts.csv"))

# -----------------------------------------------------------------------------
# Fig 3G: CytoTRACE2 differentiation potential
# -----------------------------------------------------------------------------
set.seed(1234)
keep <- unlist(lapply(split(colnames(fib), fib$subtype), function(x) sample(x, min(1000, length(x)))))
fs <- subset(fib, cells = keep)
if (requireNamespace("CytoTRACE2", quietly = TRUE)) {
  ct2 <- CytoTRACE2::cytotrace2(fs, is_seurat = TRUE, slot_type = "counts", species = "human", ncores = 4)
  fs$CytoTRACE2 <- if (inherits(ct2, "Seurat")) ct2@meta.data[colnames(fs), "CytoTRACE2_Score"] else
    ct2[colnames(fs), "CytoTRACE2_Score"]
  p3g <- ggplot(fs@meta.data, aes(reorder(subtype, -CytoTRACE2, median), CytoTRACE2, fill = subtype)) +
    geom_boxplot(outlier.size = 0.3) + scale_fill_manual(values = pal_sub) +
    labs(x = NULL, y = "CytoTRACE2 score (potency)", title = "Differentiation potential") +
    theme_classic() + umap_theme + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
  save_plot(p3g, "Fig3G_CytoTRACE2_fibroblast", 5.5, 4.5)
} else message("CytoTRACE2 not installed -> Fig 3G skipped")

# -----------------------------------------------------------------------------
# Fig 3H-J: Monocle2 pseudotime (DDRTree)
# -----------------------------------------------------------------------------
# Monocle2 runs in a separate R session (run_monocle2.R, igraph 2.0.3). Root = subtype with the
# highest median CytoTRACE2 potency (CytoTRACE "assists Monocle2 in determining the starting point")
Idents(fs) <- "subtype"
mk_o <- FindAllMarkers(fs, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.25, verbose = FALSE)
og <- unique(mk_o |> group_by(cluster) |> slice_max(avg_log2FC, n = 100) |> pull(gene))
root <- if (!is.null(fs$CytoTRACE2)) names(which.max(tapply(fs$CytoTRACE2, fs$subtype, median))) else "qFib"
message("Monocle2 root subtype: ", root)
mono_in <- file.path(RDS, "03_fibroblast_monocle_input.rds")
saveRDS(list(counts = LayerData(fs, assay = "RNA", layer = "counts"), ordering_genes = og, pal = pal_sub,
             meta = fs@meta.data[, intersect(c("sample", "group", "stage", "subtype", "CytoTRACE2"), colnames(fs@meta.data))]),
        mono_in)
system2(file.path(R.home("bin"), "Rscript.exe"),
        c(shQuote(file.path(ROOT, "scripts", "run_monocle2.R")), shQuote(mono_in), "subtype", shQuote(root), "Fig3HIJ", "Fibroblast"))

saveRDS(fib, f_fib)
cat("Step 2 done\n"); print(table(fib$subtype, fib$stage))
