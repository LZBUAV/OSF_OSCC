# =============================================================================
# Step 5a: cell-cell communication (CellChat) and the TF -> ligand -> receptor axis
# Reproduces: Fig 6A-D
# Comparison: OSF-associated OSCC vs de novo OSCC (GSE215403, same batch)
# Needs: 01_sc_integrated.rds, 03_fibroblast.rds, 04_myeloid.rds, 05_scenic_results.rds
# =============================================================================
source("utils.R")
suppressPackageStartupMessages(library(CellChat))

key <- readLines(file.path(TAB, "Step3_key_subset.txt"))[1]
obj <- load_atlas(); obj <- subset(obj, dataset == "GSE215403"); invisible(gc())
fib <- readRDS(file.path(RDS, "03_fibroblast.rds")); mye <- readRDS(file.path(RDS, "04_myeloid.rds"))
lab <- as.character(obj$celltype); names(lab) <- colnames(obj)
lab[intersect(names(lab), colnames(fib))] <- as.character(fib$subtype[intersect(names(lab), colnames(fib))])
lab[intersect(names(lab), colnames(mye))] <- as.character(mye$subtype[intersect(names(lab), colnames(mye))])
# drop cells removed as doublets during sub-clustering
drop <- (obj$celltype %in% c("Fibroblast", "Pericyte_SMC") & !colnames(obj) %in% colnames(fib)) |
        (obj$celltype == "Myeloid" & !colnames(obj) %in% colnames(mye))
obj <- subset(obj, cells = colnames(obj)[!drop]); obj$cc_label <- lab[colnames(obj)]
fib_types <- levels(fib$subtype); fib_types <- fib_types[fib_types %in% obj$cc_label]
rm(fib, mye); invisible(gc())

f_cc <- file.path(RDS, "06_cellchat_list.rds")
if (!file.exists(f_cc)) {
  CellChatDB <- CellChatDB.human
  make_cc <- function(g) {
    s <- subset(obj, group == g)
    tb <- table(s$cc_label); s <- subset(s, cc_label %in% names(tb)[tb >= 10])
    cc <- createCellChat(object = LayerData(s, assay = "RNA", layer = "data"), meta = s@meta.data, group.by = "cc_label")
    cc@DB <- CellChatDB
    cc <- subsetData(cc) |> identifyOverExpressedGenes() |> identifyOverExpressedInteractions()
    cc <- computeCommunProb(cc, type = "triMean") |> filterCommunication(min.cells = 10)
    cc <- computeCommunProbPathway(cc) |> aggregateNet() |> netAnalysis_computeCentrality()
    cc
  }
  cc_list <- list(deNovo_OSCC = make_cc("deNovo_OSCC"), OSF_OSCC = make_cc("OSF_OSCC"))
  saveRDS(cc_list, f_cc)
} else cc_list <- readRDS(f_cc)
# common cell groups so both objects can be merged
grp <- Reduce(intersect, lapply(cc_list, function(x) levels(x@idents)))
cc_list <- lapply(cc_list, function(x) subsetCellChat(x, idents.use = grp))
cc <- mergeCellChat(cc_list, add.names = names(cc_list))

# Fig 6A: interaction-number networks, key subset highlighted
dev_both <- function(name, w, h, expr) {
  pdf(file.path(FIG, paste0(name, ".pdf")), w, h); expr(); dev.off()
  png(file.path(FIG, paste0(name, ".png")), w, h, units = "in", res = 300); expr(); dev.off()
}
wmax <- max(sapply(cc_list, function(x) max(x@net$count)))
dev_both("Fig6A_Network_interaction_number", 14, 7, function() {
  par(mfrow = c(1, 2), xpd = TRUE)
  for (n in names(cc_list)) {
    m <- cc_list[[n]]@net$count; m[rownames(m) != key, ] <- 0
    netVisual_circle(m, weight.scale = TRUE, edge.weight.max = wmax, label.edge = FALSE,
                     title.name = paste0(n, ": interactions sent by ", key), vertex.label.cex = 0.7)
  }
})
dev_both("Fig6A_supp_Network_all_interactions", 14, 7, function() {
  par(mfrow = c(1, 2), xpd = TRUE)
  for (n in names(cc_list)) netVisual_circle(cc_list[[n]]@net$count, weight.scale = TRUE, edge.weight.max = wmax,
                                             label.edge = FALSE, title.name = paste(n, "- number of interactions"), vertex.label.cex = 0.7)
})

# Fig 6B: differential interaction heat maps (OSF_OSCC vs deNovo)
dev_both("Fig6B_Heatmap_differential_interactions", 13, 6.5, function() {
  h1 <- netVisual_heatmap(cc, measure = "count", comparison = c(1, 2), title.name = "Differential number of interactions")
  h2 <- netVisual_heatmap(cc, measure = "weight", comparison = c(1, 2), title.name = "Differential interaction strength")
  ComplexHeatmap::draw(h1 + h2, ht_gap = grid::unit(0.8, "cm"))
})
diff_cnt <- cc_list$OSF_OSCC@net$count - cc_list$deNovo_OSCC@net$count
fwrite(data.frame(sender = rownames(diff_cnt), diff_cnt, check.names = FALSE), file.path(TAB, "Step5_diff_interaction_count_OSF_vs_deNovo.csv"))

# -----------------------------------------------------------------------------
# L-R pairs: key subset (sender) -> fibroblast subtypes (receivers), gained in OSF_OSCC
# -----------------------------------------------------------------------------
lr <- function(x, n) {
  d <- subsetCommunication(x, sources.use = key, targets.use = intersect(fib_types, levels(x@idents)), thresh = 1)
  d$group <- n; d
}
lr_all <- bind_rows(lr(cc_list$OSF_OSCC, "OSF_OSCC"), lr(cc_list$deNovo_OSCC, "deNovo_OSCC"))
lr_w <- lr_all |> group_by(interaction_name, ligand, receptor, pathway_name, target) |>
  summarise(prob_OSF = sum(prob[group == "OSF_OSCC" & pval < 0.05]), prob_deNovo = sum(prob[group == "deNovo_OSCC" & pval < 0.05]),
            .groups = "drop") |> mutate(delta = prob_OSF - prob_deNovo)

# link to SCENIC: ligand must be a target gene of a key-subset regulon TF
sc <- readRDS(file.path(RDS, "05_scenic_results.rds"))
reg <- as.data.frame(sc$regulons); reg$TF <- sub("\\(\\+\\)|\\(-\\)", "", reg$regulon)
cand_tf <- if (length(sc$shared_TFs)) sc$shared_TFs else unique(sub("\\(\\+\\)|\\(-\\)", "", head(names(sort(sc$rss[, key], TRUE)), 10)))
lig_genes <- unique(unlist(strsplit(lr_w$ligand, "_")))
tf_lig <- subset(reg, TF %in% cand_tf & gene %in% lig_genes)
lr_w$TF <- sapply(strsplit(lr_w$ligand, "_"), function(g) paste(unique(tf_lig$TF[tf_lig$gene %in% g]), collapse = ";"))
lr_w <- lr_w[order(-(lr_w$TF != ""), -lr_w$delta), ]
fwrite(lr_w, file.path(TAB, "Step5_keysubset_to_fibroblast_LR_OSF_vs_deNovo.csv"))
fwrite(tf_lig, file.path(TAB, "Step5_TF_ligand_links.csv"))

axis_file <- file.path(TAB, "Step5_axis.csv")   # optional manual choice: TF,ligand,receptor,target
best <- if (file.exists(axis_file)) fread(axis_file) else {
  b <- subset(lr_w, TF != "" & delta > 0); if (!nrow(b)) b <- subset(lr_w, delta > 0); b <- b[1, ]
  data.frame(TF = sub(";.*", "", b$TF), ligand = b$ligand, receptor = b$receptor, target = b$target)
}
fwrite(best, file.path(TAB, "Step5_axis_selected.csv"))
message("Selected axis: ", best$TF, " -> ", best$ligand, " -> ", best$receptor, " (", best$target, ")")

# Fig 6C: communication probability of L-R pairs sent by the key subset, both groups
top_lr <- head(unique(lr_w$interaction_name[lr_w$delta > 0]), 25)
pairLR <- data.frame(interaction_name = top_lr)
p6c <- netVisual_bubble(cc, sources.use = key, targets.use = intersect(fib_types, grp), pairLR.use = pairLR,
                        comparison = c(1, 2), angle.x = 45, remove.isolate = TRUE,
                        title.name = paste("L-R pairs from", key, "to fibroblasts"))
save_plot(p6c, "Fig6C_Bubble_LR_keysubset_to_fibroblasts", 10, 8)
dev_both("Fig6C_Chord_LR_keysubset_OSF_OSCC", 9, 9, function() {
  netVisual_chord_gene(cc_list$OSF_OSCC, sources.use = key, targets.use = intersect(fib_types, grp),
                       pairLR.use = pairLR, lab.cex = 0.6, legend.pos.y = 30, title.name = paste("OSF_OSCC:", key, "-> fibroblasts"))
})

# Fig 6D: violin plots of ligand (myeloid) and receptor (fibroblasts) by group
lig <- strsplit(best$ligand, "_")[[1]][1]; rec <- strsplit(best$receptor, "_")[[1]][1]
v <- subset(obj, cc_label %in% c(grp))
vd <- data.frame(cell_group = v$cc_label, group = v$group,
                 ligand = LayerData(v, layer = "data")[lig, ], receptor = LayerData(v, layer = "data")[rec, ])
vd$group <- droplevels(vd$group)
mk_v <- function(gene_col, gene, types, title) {
  d <- subset(vd, cell_group %in% types)
  ggplot(d, aes(cell_group, .data[[gene_col]], fill = group)) + geom_violin(scale = "width", position = position_dodge(0.8), linewidth = 0.2) +
    ggpubr::stat_compare_means(aes(group = group), label = "p.signif", method = "wilcox.test") +
    scale_fill_manual(values = c(OSF_OSCC = "#F39B7F", deNovo_OSCC = "#8491B4")) +
    labs(x = NULL, y = paste(gene, "expression"), title = title) + theme_classic() + umap_theme +
    theme(axis.text.x = element_text(angle = 40, hjust = 1))
}
mye_types <- grep("Macrophage|Monocyte|DC|Neutrophil", grp, value = TRUE)
save_plot(mk_v("ligand", lig, mye_types, paste(lig, "in myeloid subsets")) /
            mk_v("receptor", rec, intersect(fib_types, grp), paste(rec, "in fibroblast subtypes")), "Fig6D_Violin_ligand_receptor", 11, 9)

saveRDS(list(cc = cc, lr = lr_w, axis = best), file.path(RDS, "06_cellchat_merged.rds"))
cat("Step 5a done. Axis:", best$TF, best$ligand, best$receptor, best$target, "\n")
