# Fig 5f: top-10 RSS regulons per myeloid population, regulons as rows (previous version was transposed)
source("utils.R")
sc <- readRDS(file.path(RDS, "05_scenic_results.rds"))
c2 <- fread(file.path(TAB, "Step4_myeloid_CytoTRACE2.csv"))          # cells used for SCENIC, with their subtype
lab <- factor(c2$subtype[match(rownames(sc$auc), c2$cell)])
top10 <- lapply(colnames(sc$rss), function(s) head(names(sort(sc$rss[, s], decreasing = TRUE)), 10))
regs <- unique(unlist(top10))
m <- sapply(levels(lab), function(s) colMeans(sc$auc[lab == s, regs, drop = FALSE]))   # regulons x populations
z <- t(scale(t(m)))
ann <- data.frame(top_in = sapply(regs, function(r) colnames(sc$rss)[which.max(sc$rss[r, ])]), row.names = regs)
o <- order(match(ann$top_in, levels(lab)), -apply(z, 1, max)); z <- z[o, ]
save_pheatmap("Fig5F_Heatmap_top10_RSS_regulons", 6.5, 15, z, cluster_rows = FALSE, cluster_cols = FALSE,
              color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100), breaks = seq(-2.5, 2.5, length.out = 101),
              main = "Top-10 RSS regulons per population\n(mean AUCell, row z-score)", angle_col = 90,
              fontsize_row = 6, border_color = NA)
cat("Fig5F rows:", nrow(z), "\n")
