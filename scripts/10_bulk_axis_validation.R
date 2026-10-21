# =============================================================================
# Step 5c: bulk validation of the TF -> ligand -> receptor axis, lactate genes, drugs
# Reproduces: Fig 6E (fibrosis score / TF vs ligand), Fig 6F (key macrophage vs target fibroblast
#             correlation across bulk cohorts), Fig S4 (other cohorts), Fig S6 (lactate genes),
#             Fig S7 (validation scRNA GSE103322), Fig S8 (CTD drug candidates)
# =============================================================================
source("utils.R")
source(file.path(ROOT, "scripts", "bulk_loaders.R"))
suppressPackageStartupMessages(library(GSVA))

key  <- readLines(file.path(TAB, "Step3_key_subset.txt"))[1]
axis <- fread(file.path(TAB, "Step5_axis_selected.csv"))
TF <- axis$TF; lig <- strsplit(axis$ligand, "_")[[1]][1]; rec <- strsplit(axis$receptor, "_")[[1]][1]; target <- axis$target
message("Axis: ", TF, " -> ", lig, " -> ", rec, " (", key, " -> ", target, ")")

bulk <- readRDS(file.path(RDS, "10_bulk_cohorts.rds"))
extra <- load_all_bulk(c("GSE319719", "GSE20170")); bulk <- c(bulk, extra)

# signatures: myeloid subsets (Step 3b) + fibroblast subtypes
msig <- readRDS(file.path(RDS, "04b_macrophage_subset_signatures.rds"))
f_fs <- file.path(RDS, "07_fibroblast_subtype_signatures.rds")
if (!file.exists(f_fs)) {
  fib <- readRDS(file.path(RDS, "03_fibroblast.rds")); Idents(fib) <- "subtype"
  fsig <- lapply(setNames(levels(fib$subtype), levels(fib$subtype)), function(s) {
    m <- FindMarkers(fib, ident.1 = s, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.25, max.cells.per.ident = 2000, verbose = FALSE)
    m <- m[m$pct.2 < 0.5 & !grepl("^(RP[SL]|MT-)", rownames(m)), ]; head(rownames(m)[order(-m$avg_log2FC)], 30)
  })
  saveRDS(fsig, f_fs); rm(fib)
} else fsig <- readRDS(f_fs)
H <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2")
fibrosis <- unique(H$gene_symbol[H$gs_name == "NABA_CORE_MATRISOME"])
stopifnot(length(fibrosis) > 100)
L <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
lactate <- unique(c(L$gene_symbol[L$gs_name == "GOBP_LACTATE_METABOLIC_PROCESS"],
                    "LDHA", "LDHB", "SLC16A1", "SLC16A3", "SLC16A4", "SLC16A7", "HIF1A", "PDK1", "PKM", "HK2", "BSG"))
sets <- c(setNames(msig, paste0("Mac:", names(msig))), setNames(fsig, paste0("Fib:", names(fsig))), list(Fibrosis = fibrosis))
scores <- lapply(bulk, function(b) gsva(ssgseaParam(b$expr, sets, minSize = 5)))
saveRDS(scores, file.path(RDS, "07_bulk_ssgsea.rds"))

cor_scatter <- function(x, y, xl, yl, title) {
  ct <- cor.test(x, y, method = "pearson")
  ggplot(data.frame(x, y), aes(x, y)) + geom_point(colour = "#3C5488", size = 2) + geom_smooth(method = "lm", colour = "#E64B35", fill = "grey85") +
    labs(x = xl, y = yl, title = sprintf("%s\nPearson R = %.2f, %s", title, ct$estimate, fmt_p(ct$p.value))) + theme_classic() + umap_theme
}
# Fig 6E (GSE125866) and Fig S4 (all other cohorts)
pan <- list(); stat <- list()
for (cn in names(bulk)) {
  e <- bulk[[cn]]$expr; s <- scores[[cn]]
  if (!all(c(lig, TF) %in% rownames(e))) next
  ttl <- sprintf("%s (n = %d)", relab(cn), ncol(e))
  pan[[cn]] <- cor_scatter(s["Fibrosis", ], e[lig, ], "Fibrosis score (core matrisome ssGSEA)", paste(lig, "expression"), ttl) +
               cor_scatter(e[TF, ], e[lig, ], paste(TF, "expression"), paste(lig, "expression"), ttl)
  stat[[cn]] <- data.frame(cohort = cn, n = ncol(e),
                           r_fibrosis_ligand = cor(s["Fibrosis", ], e[lig, ]), p1 = cor.test(s["Fibrosis", ], e[lig, ])$p.value,
                           r_TF_ligand = cor(e[TF, ], e[lig, ]), p2 = cor.test(e[TF, ], e[lig, ])$p.value)
}
fwrite(bind_rows(stat), file.path(TAB, "Step5c_axis_bulk_correlations.csv"))
# main panel = largest cohort (TCGA oral, n = 262); all other cohorts in Fig S4
save_plot(pan[["TCGA_oral"]], "Fig6E_Scatter_fibrosis_TF_ligand_TCGA_oral", 10, 4.5)
save_plot(wrap_plots(pan[setdiff(names(pan), "TCGA_oral")], ncol = 1), "FigS4_Scatter_axis_other_cohorts", 10, 4.5 * (length(pan) - 1))

# Fig 6F: correlation between key-macrophage and fibroblast-subtype scores in every cohort
cm <- bind_rows(lapply(names(scores), function(cn) {
  s <- scores[[cn]]; bind_rows(lapply(paste0("Fib:", names(fsig)), function(f) {
    ct <- cor.test(s[paste0("Mac:", key), ], s[f, ]); data.frame(cohort = cn, fib = sub("Fib:", "", f), r = ct$estimate, p = ct$p.value)
  }))
}))
fwrite(cm, file.path(TAB, "Step5c_keymac_fib_correlation_bulk.csv"))
cm$lab <- paste0(sprintf("%.2f", cm$r), ifelse(cm$p < 0.001, "\n***", ifelse(cm$p < 0.01, "\n**", ifelse(cm$p < 0.05, "\n*", ""))))
nmap <- sapply(bulk, function(b) ncol(b$expr))
cm$cohort_lab <- factor(sprintf("%s (n = %d)", relab(cm$cohort), nmap[cm$cohort]),
                        rev(sprintf("%s (n = %d)", relab(names(scores)), nmap[names(scores)])))
cm$fib <- factor(cm$fib, names(fsig))
p6f <- ggplot(cm, aes(fib, cohort_lab, fill = r)) + geom_tile(colour = "white") + geom_text(aes(label = lab), size = 3) +
  scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, fill = "Pearson R", title = paste(key, "score vs fibroblast-subtype scores (bulk cohorts)")) +
  theme_minimal() + umap_theme + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_plot(p6f, "Fig6F_Correlation_matrix_keymac_fibroblast_bulk", 8, 5.5)

# Fig S6: lactate metabolism genes vs TF / ligand / receptor
s6 <- bind_rows(lapply(c("GSE125866", "TCGA_oral"), function(cn) {
  e <- bulk[[cn]]$expr; lg <- intersect(lactate, rownames(e)); ax <- intersect(c(TF, lig, rec), rownames(e))
  bind_rows(lapply(ax, function(a) bind_rows(lapply(lg, function(g) {
    ct <- suppressWarnings(cor.test(e[a, ], e[g, ], method = "spearman")); data.frame(cohort = cn, axis_gene = a, lactate_gene = g, rho = ct$estimate, p = ct$p.value)
  }))))
}))
fwrite(s6, file.path(TAB, "Step5c_lactate_correlation.csv"))
s6$star <- ifelse(s6$p < 0.001, "***", ifelse(s6$p < 0.01, "**", ifelse(s6$p < 0.05, "*", "")))
s6$cohort <- relab(s6$cohort)
pS6 <- ggplot(s6, aes(lactate_gene, axis_gene, fill = rho)) + geom_tile(colour = "white") + geom_text(aes(label = star), size = 3) +
  facet_wrap(~cohort, ncol = 1) + scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35", limits = c(-1, 1)) +
  labs(x = NULL, y = NULL, fill = "Spearman rho", title = "Axis genes vs lactate-metabolism genes") +
  theme_minimal() + umap_theme + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
save_plot(pS6, "FigS6_Heatmap_axis_lactate_genes", 13, 5)

# Fig S7: validation scRNA cohort GSE103322 (HNSCC, Smart-seq2, log2(TPM/10+1))
f103 <- file.path(RDS, "GSE103322_axis.rds")
if (!file.exists(f103)) {
  d <- fread(file.path(SC, "GSE103322", "GSE103322_HNSCC_all_data.txt.gz"))
  meta <- as.data.frame(t(d[1:5, -1])); colnames(meta) <- d[[1]][1:5]
  g <- gsub("'", "", d[[1]][-(1:5)]); keep <- g %in% c(TF, lig, rec, msig[[key]])
  m <- as.matrix(d[-(1:5)][keep, -1]); rownames(m) <- g[keep]
  meta$celltype <- ifelse(meta[["classified  as cancer cell"]] == "1", "Malignant", meta[["non-cancer cell type"]])
  saveRDS(list(m = m, meta = meta), f103)
} else { x <- readRDS(f103); m <- x$m; meta <- x$meta }
storage.mode(m) <- "numeric"   # columns were read as text because of the 5 annotation rows
v <- data.frame(celltype = meta$celltype, t(m[intersect(c(TF, lig, rec), rownames(m)), , drop = FALSE]), check.names = FALSE) |>
  tidyr::pivot_longer(-celltype, names_to = "gene", values_to = "expr")
v <- subset(v, celltype != "0" & celltype != "-Fibroblast")
pS7a <- ggplot(v, aes(celltype, expr, fill = celltype)) + geom_violin(scale = "width") + facet_wrap(~gene, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = "log2(TPM/10 + 1)", title = "Validation cohort GSE103322: axis genes by cell type") + theme_classic() + umap_theme +
  theme(legend.position = "none", axis.text.x = element_text(angle = 40, hjust = 1))
mac <- meta$celltype == "Macrophage"
ks <- colMeans(m[intersect(msig[[key]], rownames(m)), mac, drop = FALSE])
pS7b <- cor_scatter(ks, m[TF, mac], paste(key, "signature (macrophages)"), TF, "GSE103322 macrophages")
save_plot(pS7a | pS7b, "FigS7_Validation_GSE103322_axis", 13, 8)

# Fig S8: CTD chemicals interacting with all three axis genes (+ any oral/fibrosis/cancer disease link)
ctd_f <- file.path(ROOT, "data", "CTD", "CTD_chem_gene_ixns.tsv.gz")
read_ctd <- function(f, select = NULL) {   # CTD files: comment block, "# Fields:", "# <header>", "#", data
  h <- readLines(gzfile(f), n = 60); i <- grep("^# Fields:", h)
  cols <- strsplit(sub("^# ", "", h[i + 1]), "\t")[[1]]
  d <- fread(f, skip = i + 2, header = FALSE, sep = "\t", quote = "", fill = TRUE, select = select)
  setnames(d, cols[if (is.null(select)) seq_len(ncol(d)) else select]); d
}
if (file.exists(ctd_f)) {
  ctd <- read_ctd(ctd_f)
  ctd <- ctd[Organism == "Homo sapiens" | is.na(Organism) | Organism == ""]
  hits <- ctd[GeneSymbol %in% c(TF, lig, rec), .(n_genes = uniqueN(GeneSymbol), genes = paste(sort(unique(GeneSymbol)), collapse = ";"),
                                                n_refs = .N), by = .(ChemicalName, ChemicalID)]
  all3 <- hits[n_genes == length(unique(c(TF, lig, rec)))][order(-n_refs)]
  cd_f <- file.path(ROOT, "data", "CTD", "CTD_chemicals_diseases.tsv.gz")
  if (file.exists(cd_f) && nrow(all3)) {
    cd <- read_ctd(cd_f, select = 1:6)
    rel <- "^(Fibrosis|Oral Submucous Fibrosis|Mouth Neoplasms|Tongue Neoplasms|Carcinoma, Squamous Cell|Squamous Cell Carcinoma of Head and Neck|Head and Neck Neoplasms|Liver Cirrhosis|Pulmonary Fibrosis|Idiopathic Pulmonary Fibrosis)$"
    cd <- cd[DirectEvidence == "therapeutic" & grepl(rel, DiseaseName)]
    all3$diseases <- sapply(all3$ChemicalName, function(c) paste(unique(cd$DiseaseName[cd$ChemicalName == c]), collapse = "; "))
  }
  fwrite(all3, file.path(TAB, "Step5c_CTD_chemicals_axis_all.csv"))
  drugs <- all3[nchar(diseases) > 0][order(-n_refs)]
  fwrite(drugs, file.path(TAB, "Step5c_CTD_drugs_axis.csv"))
  top <- head(drugs, 20); top$ChemicalName <- factor(top$ChemicalName, rev(top$ChemicalName))
  top$label <- ifelse(nchar(top$diseases) > 60, paste0(substr(top$diseases, 1, 57), "..."), top$diseases)
  pS8 <- ggplot(top, aes(n_refs, ChemicalName)) + geom_col(fill = "#E64B35") +
    geom_text(aes(label = label), hjust = -0.02, size = 2.6) + scale_x_continuous(expand = expansion(mult = c(0, 1.2))) +
    labs(x = "CTD chemical-gene interactions with all three axis genes", y = NULL,
         title = paste("Drugs modulating", TF, "/", lig, "/", rec, "with therapeutic evidence for fibrosis / OSCC")) +
    theme_classic() + umap_theme
  save_plot(pS8, "FigS8_CTD_candidate_drugs", 9, 6.5)
}
cat("Step 5c done\n")
