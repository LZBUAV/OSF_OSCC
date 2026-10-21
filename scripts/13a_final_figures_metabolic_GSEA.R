# Regenerate figures whose titles were clipped (found by 17_check_truncated_titles.py):
#   Fig2H (both contrasts), Fig2E_supp (non-epithelial Ro/e), Fig7C_supp (complement GSEA)
source("utils.R")
suppressPackageStartupMessages(library(fgsea))

# ---- Fig 2H: metabolic GSEA heat maps from the saved result tables ----
draw_2h <- function(tag, title) {
  met <- fread(file.path(TAB, paste0("Step1b_metabolic_GSEA_", tag, ".csv")))
  top_path <- met |> group_by(pathway) |> summarise(m = min(padj, na.rm = TRUE)) |> arrange(m) |> slice_head(n = 30) |> pull(pathway)
  nes  <- reshape2::acast(subset(met, pathway %in% top_path), pathway ~ celltype, value.var = "NES")
  padj <- reshape2::acast(subset(met, pathway %in% top_path), pathway ~ celltype, value.var = "padj")
  nes[is.na(nes)] <- 0; padj[is.na(padj)] <- 1
  star <- ifelse(padj < 0.001, "***", ifelse(padj < 0.01, "**", ifelse(padj < 0.05, "*", "")))
  colnames(nes) <- relab(colnames(nes))
  save_pheatmap(paste0("Fig2H_Heatmap_metabolic_pathways_", tag), 9.5, 8, nes, display_numbers = star, number_color = "black",
                cluster_cols = FALSE, color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
                breaks = seq(-3, 3, length.out = 101), border_color = "grey70", fontsize_row = 8, main = title, angle_col = 45)
}
draw_2h("OSF_OSCC_vs_deNovo", "Metabolic pathways: OSF-associated vs de novo OSCC (GSEA NES)")
draw_2h("OSCC_vs_OSF", "Metabolic pathways: OSCC vs OSF (GSEA NES)")

# ---- Fig 2E supp: Ro/e of non-epithelial cell types by stage ----
obj <- readRDS(file.path(RDS, "01_sc_integrated.rds"))
mne <- subset(obj@meta.data, celltype != "Epithelial"); mne$celltype <- droplevels(mne$celltype); rm(obj); invisible(gc())
plot_roe(calc_roe(mne, "celltype", "stage"), "Fig2E_supp_RoE_nonepithelial_by_stage", "Ro/e, non-epithelial cells")

# ---- Fig 7C supp: GSEA plot of complement and coagulation cascades (IL1B-high vs low) ----
de <- fread(file.path(TAB, "SuppTable12_DEG_keymacro_ligand_high_vs_low.csv"))
rk <- sort(setNames(de$avg_log2FC, de$gene), decreasing = TRUE)
k2g <- kegg_sets(); ks <- split(k2g$gene, k2g$term)
gs <- fread(file.path(TAB, "Step6_GSEA_KEGG_ligand_high_vs_low.csv"))
p <- "COMPLEMENT AND COAGULATION CASCADES"
g <- plotEnrichment(ks[[p]], rk) +
  labs(title = "Complement and coagulation cascades",
       subtitle = sprintf("IL1B-high vs IL1B-low CXCL5+ macrophages; NES = %.2f, adjusted P = %.2f", gs$NES[gs$pathway == p], gs$padj[gs$pathway == p]),
       x = "Rank", y = "Enrichment score") +
  theme_classic() + umap_theme + theme(plot.subtitle = element_text(hjust = 0.5, size = 9))
save_plot(g, "Fig7C_supp_GSEA_complement_coagulation", 7, 4.8)
cat("regenerated\n")
