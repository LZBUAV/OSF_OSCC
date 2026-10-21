# =============================================================================
# Step 3b: multi-evidence selection of the key macrophage subset
#   evidence 1: Ro/e in OSF-associated vs de novo OSCC (scRNA, GSE215403)
#   evidence 2: spatial co-localisation of the subset signature with fibroblasts (Visium GSE220978)
#   evidence 3: signature change from OSF to paired OSCC (bulk GSE125866, 8 pairs)
# Output: tables/Step3b_key_subset_evidence.csv, tables/Step3_key_subset.txt
# =============================================================================
source("utils.R")
suppressPackageStartupMessages(library(GSVA))

mye <- readRDS(file.path(RDS, "04_myeloid.rds"))
enr <- fread(file.path(TAB, "Step3_myeloid_subtype_enrichment_OSF_OSCC.csv"))
macs <- setdiff(grep("Macrophage", levels(mye$subtype), value = TRUE), "Prolif. Macrophage")

# subset signatures: top 30 specific markers vs all other myeloid cells
Idents(mye) <- "subtype"
sig <- lapply(setNames(macs, macs), function(s) {
  m <- FindMarkers(mye, ident.1 = s, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.25, max.cells.per.ident = 2000, verbose = FALSE)
  m <- m[m$pct.2 < 0.4 & !grepl("^(RP[SL]|MT-|HB[AB]|IG[HKL])", rownames(m)), ]
  head(rownames(m)[order(-m$avg_log2FC)], 30)
})
saveRDS(sig, file.path(RDS, "04b_macrophage_subset_signatures.rds"))
fwrite(data.table(subset = rep(names(sig), lengths(sig)), gene = unlist(sig)), file.path(TAB, "Step3b_macrophage_subset_signatures.csv"))

# evidence 2: Visium
rw <- fread(file.path(TAB, "Step1b_RCTD_weights_GSE220978.csv"))
sp <- list()
for (i in 1:4) {
  pt <- paste0("Patient", i); gsm <- paste0("GSM68334", 83 + i); pre <- file.path(ST, "GSE220978", pt, paste0(gsm, "_", pt, "_"))
  m <- ReadMtx(paste0(pre, "matrix.mtx.gz"), paste0(pre, "barcodes.tsv.gz"), paste0(pre, "features.tsv.gz"))
  w <- rw[sample == pt]; m <- m[, w$barcode]
  lc <- log1p(t(t(m) / Matrix::colSums(m)) * 1e4)
  for (s in names(sig)) {
    g <- intersect(sig[[s]], rownames(lc)); sc <- Matrix::colMeans(lc[g, , drop = FALSE])
    ct <- cor.test(sc, w$Fibroblast, method = "spearman")
    sp[[paste(pt, s)]] <- data.frame(sample = pt, subset = s, rho_fibroblast = ct$estimate, p = ct$p.value)
  }
}
sp <- bind_rows(sp)
fwrite(sp, file.path(TAB, "Step3b_spatial_subset_fibroblast_coloc.csv"))

# evidence 3: GSE125866 bulk (counts), OSF vs paired OSCC
read_bulk125866 <- function() {
  fs <- list.files(file.path(BULK, "GSE125866"), "FPKM.xls.gz$", recursive = TRUE, full.names = TRUE)
  l <- lapply(fs, function(f) { d <- fread(f, header = FALSE); setNames(d$V2, d$V1) })
  g <- Reduce(intersect, lapply(l, names))
  cnt <- sapply(l, function(x) x[g]); colnames(cnt) <- sub("\\.gene.*", "", basename(fs))
  meta <- data.frame(sample = colnames(cnt), group = basename(dirname(fs)),
                     id = sub("^GSM\\d+_", "", colnames(cnt)))
  meta$patient <- ifelse(meta$group == "Normal", meta$id, sub("[NT]$", "", meta$id))
  list(counts = cnt, meta = meta)
}
b <- read_bulk125866()
saveRDS(b, file.path(RDS, "bulk_GSE125866.rds"))
lcpm <- log2(t(t(b$counts) / colSums(b$counts)) * 1e6 + 1)
ss <- gsva(ssgseaParam(lcpm, sig))
bk <- lapply(names(sig), function(s) {
  d <- data.frame(score = ss[s, ], b$meta)
  pr <- intersect(d$patient[d$group == "OSF"], d$patient[d$group == "OSCC"])
  a <- d$score[d$group == "OSCC"][match(pr, d$patient[d$group == "OSCC"])]
  o <- d$score[d$group == "OSF"][match(pr, d$patient[d$group == "OSF"])]
  data.frame(subset = s, n_pairs = length(pr), mean_diff_OSCC_minus_OSF = mean(a - o),
             p_paired_wilcox = wilcox.test(a, o, paired = TRUE)$p.value)
})
bk <- bind_rows(bk)

ev <- data.frame(subset = macs, RoE_OSF_OSCC = enr$RoE_OSF_OSCC[match(macs, enr$subtype)]) |>
  left_join(sp |> group_by(subset) |> summarise(mean_rho_fibroblast = mean(rho_fibroblast)), by = "subset") |>
  left_join(bk, by = "subset")
ev$rank_sum <- rank(-ev$RoE_OSF_OSCC) + rank(-ev$mean_rho_fibroblast) + rank(-ev$mean_diff_OSCC_minus_OSF)
ev <- ev[order(ev$rank_sum), ]
fwrite(ev, file.path(TAB, "Step3b_key_subset_evidence.csv"))
print(ev)
writeLines(ev$subset[1], file.path(TAB, "Step3_key_subset.txt"))
cat("Key subset (multi-evidence):", ev$subset[1], "\n")
