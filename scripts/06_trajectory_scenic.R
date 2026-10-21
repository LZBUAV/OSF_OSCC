# =============================================================================
# Step 4: CytoTRACE2, Monocle2 pseudotime and pySCENIC regulons of myeloid cells
# Reproduces: Fig 5A-G, Fig S3A-B, Supplementary Table 7 (regulon modules / TFs)
# Needs: analysis/rds/04_myeloid.rds, tables/Step3_key_subset.txt, scenic_db/, conda env "pyscenic"
# =============================================================================
source("utils.R")

mye <- readRDS(file.path(RDS, "04_myeloid.rds"))
key <- readLines(file.path(TAB, "Step3_key_subset.txt"))[1]
pal_mye <- setNames(rep(pal_many, 3)[seq_along(levels(mye$subtype))], levels(mye$subtype))
message("Key subset: ", key)

set.seed(1234)
cells <- unlist(lapply(split(colnames(mye), mye$subtype), function(x) sample(x, min(600, length(x)))))
ms <- subset(mye, cells = cells)

# -----------------------------------------------------------------------------
# Fig 5A: CytoTRACE2 stemness across myeloid cells
# -----------------------------------------------------------------------------
f_ct2 <- file.path(TAB, "Step4_myeloid_CytoTRACE2.csv")
if (file.exists(f_ct2)) {
  ctd <- fread(f_ct2); ms <- subset(ms, cells = intersect(colnames(ms), ctd$cell))
  ms$CytoTRACE2 <- ctd$CytoTRACE2[match(colnames(ms), ctd$cell)]
} else {
  ct2 <- CytoTRACE2::cytotrace2(ms, is_seurat = TRUE, slot_type = "counts", species = "human", ncores = 1)
  ms$CytoTRACE2 <- if (inherits(ct2, "Seurat")) ct2@meta.data[colnames(ms), "CytoTRACE2_Score"] else ct2[colnames(ms), "CytoTRACE2_Score"]
}
d <- data.frame(Embeddings(ms, "umap")[, 1:2], score = ms$CytoTRACE2); colnames(d)[1:2] <- c("UMAP_1", "UMAP_2")
p5a1 <- ggplot(d[order(d$score), ], aes(UMAP_1, UMAP_2, colour = score)) + geom_point(size = 0.5) +
  scale_colour_gradientn(colours = c("#3C5488", "#4DBBD5", "#FDE725", "#E64B35"), limits = c(0, 1), name = "CytoTRACE2") +
  theme_classic() + ggtitle("CytoTRACE2 stemness (0 differentiated - 1 stem-like)") + umap_theme
p5a2 <- ggplot(ms@meta.data, aes(reorder(subtype, -CytoTRACE2, median), CytoTRACE2, fill = subtype)) +
  geom_boxplot(outlier.size = 0.2) + scale_fill_manual(values = pal_mye) + labs(x = NULL, y = "CytoTRACE2 score") +
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 40, hjust = 1))
save_plot(p5a1 + p5a2 + plot_layout(widths = c(1.1, 1)), "Fig5A_CytoTRACE2_myeloid", 13, 5.5)
fwrite(data.frame(cell = colnames(ms), subtype = ms$subtype, CytoTRACE2 = ms$CytoTRACE2), file.path(TAB, "Step4_myeloid_CytoTRACE2.csv"))

# -----------------------------------------------------------------------------
# Fig 5B-D: Monocle2 on the monocyte / macrophage lineage (separate R session, see run_monocle2.R)
# -----------------------------------------------------------------------------
lin <- grep("Monocyte|Macrophage", levels(ms$subtype), value = TRUE)
ml <- subset(ms, subtype %in% lin); ml$subtype <- droplevels(ml$subtype)
Idents(ml) <- "subtype"
mk <- FindAllMarkers(ml, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.25, verbose = FALSE)
og <- unique(mk |> group_by(cluster) |> slice_max(avg_log2FC, n = 80) |> pull(gene))
# proliferating cells score high on CytoTRACE2 by construction -> excluded from root choice
med <- tapply(ml$CytoTRACE2, ml$subtype, median); med <- med[!grepl("Prolif", names(med))]
root <- names(which.max(med))
message("Monocle2 root: ", root)
mono_in <- file.path(RDS, "05_myeloid_monocle_input.rds")
saveRDS(list(counts = LayerData(ml, assay = "RNA", layer = "counts"), ordering_genes = og, pal = pal_mye[lin],
             meta = ml@meta.data[, c("sample", "group", "stage", "subtype", "CytoTRACE2")]), mono_in)
if (!file.exists(file.path(FIG, "Fig5BCD_pseudotime_heatmap.png")))
  system2(file.path(R.home("bin"), "Rscript.exe"),
          c(shQuote(file.path(ROOT, "scripts", "run_monocle2.R")), shQuote(mono_in), "subtype", shQuote(root), "Fig5BCD", "Myeloid"))

# -----------------------------------------------------------------------------
# pySCENIC (Python) on the down-sampled myeloid cells
# -----------------------------------------------------------------------------
sc_dir <- file.path(OUT, "scenic"); dir.create(sc_dir, showWarnings = FALSE)
pre <- file.path(sc_dir, "myeloid")
if (!file.exists(paste0(pre, "_auc.csv"))) {
  cnt <- LayerData(ms, assay = "RNA", layer = "counts")
  cnt <- cnt[Matrix::rowSums(cnt > 0) >= 0.02 * ncol(cnt) & !grepl("^(MT-|RP[SL])", rownames(cnt)), ]
  ex <- as.data.frame(t(as.matrix(cnt))); ex <- cbind(cell = rownames(ex), ex)
  fwrite(ex, paste0(pre, "_expr.csv.gz")); rm(ex)
  system2(Sys.getenv("PYSCENIC_PYTHON", "python"),
          c(shQuote(file.path(ROOT, "scripts", "run_pyscenic.py")), shQuote(paste0(pre, "_expr.csv.gz")), shQuote(pre), "6"))
}
auc <- as.matrix(read.csv(paste0(pre, "_auc.csv"), row.names = 1, check.names = FALSE))
auc <- auc[colnames(ms), ]
reg <- fread(paste0(pre, "_regulons.csv"))
message("Regulons: ", ncol(auc))
lab <- ms$subtype

# -----------------------------------------------------------------------------
# Regulon modules by Connection Specificity Index (CSI), Fig 5E / S3A
# -----------------------------------------------------------------------------
pcc <- cor(auc)
n <- ncol(pcc); csi <- matrix(0, n, n, dimnames = dimnames(pcc))
for (i in seq_len(n)) for (j in seq_len(n)) {
  thr <- pcc[i, j] - 0.05
  csi[i, j] <- mean(pcc[i, -c(i, j)] < thr & pcc[j, -c(i, j)] < thr)
}
hc <- hclust(dist(csi), method = "ward.D2")
k <- max(6, min(17, round(n / 20)))
mod <- cutree(hc, k = k); mod <- setNames(paste0("M", mod), names(mod))
mods <- split(names(mod), mod)

mod_act <- sapply(mods, function(r) tapply(rowMeans(auc[, r, drop = FALSE]), lab, mean))
mod_z <- t(scale(mod_act))
ann_r <- data.frame(module = mod, row.names = names(mod))
save_pheatmap("FigS3A_CSI_regulon_modules", 10, 9, csi[hc$order, hc$order], cluster_rows = FALSE, cluster_cols = FALSE,
              show_rownames = FALSE, show_colnames = FALSE, annotation_row = ann_r,
              color = colorRampPalette(c("white", "#F39B7F", "#7E0000"))(100), main = "Regulon CSI matrix and modules")
save_pheatmap("Fig5E_Heatmap_regulon_module_activity", 8, 6, mod_z,
              color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100), breaks = seq(-2, 2, length.out = 101),
              main = "Regulon module activity (z-score)", angle_col = 45, border_color = "grey70")

# modules enriched in the key subset = module activity z-score >= 0.5 in the key subset (Fig 5E)
key_mods <- rownames(mod_z)[mod_z[, key] >= 0.5]
tf_mod <- unique(sub("\\(\\+\\)|\\(-\\)", "", unlist(mods[key_mods])))
fwrite(data.frame(regulon = names(mod), module = mod, in_key_module = mod %in% key_mods),
       file.path(TAB, "SuppTable7_regulon_modules.csv"))

# -----------------------------------------------------------------------------
# Regulon specificity score (Suo et al. 2018), Fig 5F / S3B
# -----------------------------------------------------------------------------
jsd <- function(p, q) { m <- (p + q) / 2; kl <- function(a, b) sum(ifelse(a > 0, a * log2(a / b), 0)); 0.5 * kl(p, m) + 0.5 * kl(q, m) }
rss <- sapply(levels(lab), function(s) {
  q <- as.numeric(lab == s); q <- q / sum(q)
  apply(auc, 2, function(a) { p <- a / sum(a); 1 - sqrt(jsd(p, q)) })
})
fwrite(data.frame(regulon = rownames(rss), rss, check.names = FALSE), file.path(TAB, "Step4_regulon_specificity_scores.csv"))
top10 <- lapply(colnames(rss), function(s) head(names(sort(rss[, s], decreasing = TRUE)), 10)); names(top10) <- colnames(rss)
top_all <- unique(unlist(top10))
rz <- t(scale(t(sapply(levels(lab), function(s) colMeans(auc[lab == s, top_all, drop = FALSE])))))
save_pheatmap("Fig5F_Heatmap_top10_RSS_regulons", 8, 14, t(rz), cluster_cols = FALSE,
              color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100), breaks = seq(-2.5, 2.5, length.out = 101),
              main = "Top 10 regulons (RSS) per myeloid subset", angle_col = 45, fontsize_row = 7)
rank_plot <- function(s) {
  d <- data.frame(regulon = rownames(rss), rss = rss[, s]); d <- d[order(-d$rss), ]; d$rank <- seq_len(nrow(d))
  ggplot(d, aes(rank, rss)) + geom_point(size = 0.8, colour = ifelse(d$rank <= 10, "#E64B35", "grey60")) +
    geom_text_repel(data = head(d, 10), aes(label = regulon), size = 3, max.overlaps = 20) +
    labs(title = s, x = "Rank", y = "RSS") + theme_classic() + umap_theme
}
save_plot(rank_plot(key), "Fig5F_RSS_rank_key_subset", 5, 4.5)
save_plot(wrap_plots(lapply(colnames(rss), rank_plot), ncol = 4), "FigS3B_RSS_rank_all_subsets", 18, 4 * ceiling(ncol(rss) / 4))

# -----------------------------------------------------------------------------
# Fig 5G: Venn of key-module TFs and top-10 RSS TFs of the key subset
# -----------------------------------------------------------------------------
tf_rss <- sub("\\(\\+\\)|\\(-\\)", "", top10[[key]])
shared <- intersect(tf_mod, tf_rss)
p5g <- ggVennDiagram::ggVennDiagram(list(`Key-subset regulon modules` = tf_mod, `Top10 RSS TFs` = tf_rss), label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "#E64B35") + ggtitle(paste("Shared TFs:", paste(shared, collapse = ", "))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 10))
save_plot(p5g, "Fig5G_Venn_key_TFs", 6.5, 5)
fwrite(data.frame(set = c(rep("module_TFs", length(tf_mod)), rep("RSS_top10", length(tf_rss)), rep("shared", length(shared))),
                  TF = c(tf_mod, tf_rss, shared)), file.path(TAB, "SuppTable7b_key_TFs.csv"))
saveRDS(list(auc = auc, modules = mods, key_modules = key_mods, rss = rss, shared_TFs = shared, regulons = reg),
        file.path(RDS, "05_scenic_results.rds"))
cat("Step 4 done. Key modules:", key_mods, " Shared TFs:", shared, "\n")
