# =============================================================================
# Step 5b: DNA methylation (replaces the ATAC-seq analysis of the original paper)
# Reproduces: Fig S5A-C (OSF vs normal DMPs; axis-gene CpGs by genomic region), Fig S11 (promoter
#             methylation of model genes), Supplementary Tables 8-11
# Data: GSE156669 (450K, normal buccal mucosa n=5 vs OSF n=7), TCGA-HNSC 450K (oral OSCC vs normal)
# Gene lists: tables/Step5_axis_selected.csv and (optional) tables/Step7_model_genes.txt
# =============================================================================
source("utils.R")
suppressPackageStartupMessages(library(limma))

# ---- annotation (Illumina 450K manifest, hg19) ----
man_f <- file.path(METH, "GSE156669", "Platform_annotation", "GPL13534_HumanMethylation450_15017482_v.1.1.csv.gz")
f_ann <- file.path(RDS, "450k_annotation.rds")
if (!file.exists(f_ann)) {
  l <- readLines(gzfile(man_f)); s <- grep("^IlmnID", l); e <- grep("^\\[Controls\\]", l) - 1
  ann <- fread(text = l[s:e], select = c("IlmnID", "CHR", "MAPINFO", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "Relation_to_UCSC_CpG_Island"))
  saveRDS(ann, f_ann)
} else ann <- readRDS(f_ann)
ann <- ann[!CHR %in% c("X", "Y")]
# one row per probe-gene-region
ann_long <- ann[UCSC_RefGene_Name != "", {
  g <- strsplit(UCSC_RefGene_Name, ";")[[1]]; r <- strsplit(UCSC_RefGene_Group, ";")[[1]]
  .(gene = g, region = r[seq_along(g)])
}, by = .(IlmnID, CHR, MAPINFO, Relation_to_UCSC_CpG_Island)] |> unique()
region_levels <- c("TSS1500", "TSS200", "5'UTR", "1stExon", "Body", "3'UTR")
ann_long$region <- factor(ann_long$region, region_levels)

# ---- GSE156669 ----
l <- readLines(gzfile(file.path(METH, "GSE156669", "GSE156669_series_matrix.txt.gz")))
b <- grep("^!series_matrix_table_begin", l); e <- grep("^!series_matrix_table_end", l)
bm <- fread(text = l[(b + 1):(e - 1)]); ids <- bm[[1]]; beta <- as.matrix(bm[, -1]); rownames(beta) <- ids
title <- gsub('"', "", strsplit(l[grepl("^!Sample_title", l)], "\t")[[1]][-1])
grp <- factor(ifelse(grepl("OSF", title), "OSF", "Normal"), c("Normal", "OSF"))
beta <- beta[rownames(beta) %in% ann$IlmnID & rowSums(is.na(beta)) == 0, ]
M <- log2(pmax(beta, 1e-3) / pmax(1 - beta, 1e-3))
fit <- eBayes(lmFit(M, model.matrix(~grp)))
dmp <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
dmp$delta_beta <- rowMeans(beta[, grp == "OSF"]) - rowMeans(beta[, grp == "Normal"])
dmp$probe <- rownames(dmp)
dmp$status <- ifelse(dmp$P.Value < 0.05 & dmp$delta_beta > 0.1, "Hyper", ifelse(dmp$P.Value < 0.05 & dmp$delta_beta < -0.1, "Hypo", "NS"))
fwrite(dmp[dmp$status != "NS", c("probe", "delta_beta", "logFC", "P.Value", "adj.P.Val", "status")],
       file.path(TAB, "SuppTable8_DMP_OSF_vs_Normal_GSE156669.csv"))
message("DMPs: ", paste(names(table(dmp$status)), table(dmp$status), collapse = " "))

pS5a <- ggplot(dmp, aes(delta_beta, -log10(P.Value), colour = status)) + geom_point(size = 0.3, alpha = 0.5) +
  scale_colour_manual(values = c(Hyper = "#E64B35", Hypo = "#3C5488", NS = "grey80")) +
  geom_vline(xintercept = c(-0.1, 0.1), linetype = 2) + geom_hline(yintercept = -log10(0.05), linetype = 2) +
  labs(x = "Delta beta (OSF - Normal)", y = "-log10 P", colour = NULL,
       title = sprintf("Differentially methylated CpGs: OSF vs normal (hyper %d, hypo %d)", sum(dmp$status == "Hyper"), sum(dmp$status == "Hypo"))) +
  theme_classic() + umap_theme
save_plot(pS5a, "FigS5A_Volcano_DMP_OSF_vs_Normal", 7, 5.5)

# ---- axis genes / model genes ----
axis <- fread(file.path(TAB, "Step5_axis_selected.csv"))
axis_genes <- unique(c(axis$TF, strsplit(axis$ligand, "_")[[1]], strsplit(axis$receptor, "_")[[1]]))
mg_f <- file.path(TAB, "Step7_model_genes.txt")
model_genes <- if (file.exists(mg_f)) readLines(mg_f) else character(0)
genes <- unique(c(axis_genes, model_genes))
message("Genes: ", paste(genes, collapse = ", "))

ga <- merge(ann_long[gene %in% genes], data.table(IlmnID = dmp$probe, delta_beta = dmp$delta_beta, P = dmp$P.Value, status = dmp$status), by = "IlmnID")
fwrite(ga, file.path(TAB, "SuppTable9_axis_gene_CpGs_GSE156669.csv"))

# Fig S5B: genomic region distribution of CpGs (and DMPs) on axis genes
d5b <- ga[gene %in% axis_genes, .N, by = .(gene, region, dm = ifelse(status == "NS", "Not DM", status))]
pS5b <- ggplot(d5b, aes(region, N, fill = dm)) + geom_col() + facet_wrap(~gene, scales = "free_y") +
  scale_fill_manual(values = c(Hyper = "#E64B35", Hypo = "#3C5488", `Not DM` = "grey75")) +
  labs(x = NULL, y = "Number of CpGs", fill = NULL, title = "CpGs of axis genes by genomic region (OSF vs normal)") +
  theme_classic() + umap_theme + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_plot(pS5b, "FigS5B_axis_gene_CpG_regions", 10, 4.5)

# Fig S5C: delta beta of each CpG along the gene (TSS1500 -> 3'UTR) for axis genes
pS5c <- ggplot(ga[gene %in% axis_genes], aes(region, delta_beta, colour = status)) +
  geom_hline(yintercept = 0, colour = "grey60") + geom_jitter(width = 0.2, size = 2) +
  scale_colour_manual(values = c(Hyper = "#E64B35", Hypo = "#3C5488", NS = "grey60")) + facet_wrap(~gene) +
  labs(x = NULL, y = "Delta beta (OSF - Normal)", colour = NULL, title = "Methylation change relative to TSS") +
  theme_classic() + umap_theme + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_plot(pS5c, "FigS5C_axis_gene_deltabeta_by_region", 10, 4.5)

# ---- TCGA oral 450K: promoter methylation, OSCC vs normal, and correlation with expression ----
f_tm <- file.path(RDS, "TCGA_oral_450k_genes.rds")
prom_probes <- unique(ann_long[gene %in% genes & region %in% c("TSS1500", "TSS200", "5'UTR", "1stExon")]$IlmnID)
if (!file.exists(f_tm)) {
  tc <- get_TCGA_oral()
  hdr <- names(fread(file.path(ROOT, "data/TCGA", "TCGA-HNSC.methylation450.tsv.gz"), nrows = 0))
  keep_s <- hdr[substr(hdr, 1, 12) %in% tc$meta$patient & substr(hdr, 14, 15) %in% c("01", "11")]
  tm <- fread(file.path(ROOT, "data/TCGA", "TCGA-HNSC.methylation450.tsv.gz"), select = c(hdr[1], keep_s))
  tm <- tm[get(hdr[1]) %in% prom_probes]
  saveRDS(tm, f_tm)
} else tm <- readRDS(f_tm)
tb <- as.matrix(tm[, -1]); rownames(tb) <- tm[[1]]
tgrp <- ifelse(substr(colnames(tb), 14, 15) == "01", "OSCC", "Normal")

prom_mean <- function(bm, g) {
  p <- intersect(unique(ann_long[gene == g & region %in% c("TSS1500", "TSS200", "5'UTR", "1stExon")]$IlmnID), rownames(bm))
  if (!length(p)) return(rep(NA, ncol(bm))); colMeans(bm[p, , drop = FALSE], na.rm = TRUE)
}
s11 <- bind_rows(lapply(genes, function(g) bind_rows(
  data.frame(gene = g, cohort = "GSE156669 (OSF vs Normal)", group = as.character(grp), beta = prom_mean(beta, g)),
  data.frame(gene = g, cohort = "TCGA oral (OSCC vs Normal)", group = tgrp, beta = prom_mean(tb, g)))))
s11 <- s11[!is.na(s11$beta), ]
s11$group <- factor(s11$group, c("Normal", "OSF", "OSCC"))
fwrite(s11, file.path(TAB, "SuppTable10_promoter_methylation_genes.csv"))
pS11 <- ggplot(s11, aes(group, beta, fill = group)) + geom_boxplot(outlier.size = 0.3, width = 0.6) +
  ggpubr::stat_compare_means(method = "wilcox.test", label = "p.signif", label.x = 1.5) +
  facet_grid(cohort ~ gene, scales = "free_x") + scale_fill_manual(values = pal_stage) +
  labs(x = NULL, y = "Promoter beta value (TSS1500/TSS200/5'UTR/1stExon)", title = "Promoter methylation of axis / model genes") +
  theme_bw() + umap_theme + theme(legend.position = "none")
save_plot(pS11, "FigS11_promoter_methylation_genes", max(6, 1.6 * length(genes)), 6.5)

# methylation-expression correlation in TCGA oral tumours (Supplementary Table 11)
tc <- get_TCGA_oral()
cs <- intersect(colnames(tb), colnames(tc$expr))
mc <- bind_rows(lapply(intersect(genes, rownames(tc$expr)), function(g) {
  b <- prom_mean(tb[, cs, drop = FALSE], g); if (all(is.na(b))) return(NULL)
  ct <- suppressWarnings(cor.test(b, tc$expr[g, cs], method = "spearman"))
  data.frame(gene = g, rho_meth_expr = ct$estimate, p = ct$p.value, n = sum(!is.na(b)))
}))
fwrite(mc, file.path(TAB, "SuppTable11_methylation_expression_correlation_TCGA.csv"))
print(mc)
cat("Step 5b done\n")
