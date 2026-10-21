# Bulk transcriptome loaders. Each returns list(expr = log-scale gene x sample matrix (symbols),
# meta = data.frame(sample, group, ...), counts = raw counts if available, name)
suppressPackageStartupMessages({ library(data.table); library(Matrix) })

.ens2sym <- function() {
  pm <- fread(file.path(ROOT, "data/TCGA", "gencode.v36.annotation.gtf.gene.probemap"))
  setNames(pm$gene, sub("\\..*", "", pm$id))
}
.collapse <- function(m, sym) {               # rows -> gene symbols, keep the row with highest mean
  ok <- !is.na(sym) & sym != ""; m <- m[ok, , drop = FALSE]; sym <- sym[ok]
  o <- order(rowMeans(m), decreasing = TRUE); m <- m[o, , drop = FALSE]; sym <- sym[o]
  m <- m[!duplicated(sym), , drop = FALSE]; rownames(m) <- sym[!duplicated(sym)]; m
}
.logcpm <- function(cnt) log2(t(t(cnt) / colSums(cnt)) * 1e6 + 1)
.series_matrix <- function(f) {
  l <- readLines(gzfile(f)); b <- grep("^!series_matrix_table_begin", l); e <- grep("^!series_matrix_table_end", l)
  m <- fread(text = l[(b + 1):(e - 1)]); ids <- m[[1]]; m <- as.matrix(m[, -1]); rownames(m) <- ids
  hdr <- l[grepl("^!Sample_(title|geo_accession|characteristics_ch1|source_name_ch1)", l)]
  list(mat = m, hdr = hdr)
}
.hdr_field <- function(hdr, key) gsub('"', "", strsplit(hdr[grepl(key, hdr)][1], "\t")[[1]][-1])
.illumina_ann <- function(f) {
  l <- readLines(gzfile(f)); s <- grep("^\\[Probes\\]", l) + 1; e <- grep("^\\[Controls\\]", l) - 1
  a <- fread(text = l[s:e], fill = TRUE, quote = ""); setNames(a$Symbol, a$Probe_Id)
}

load_GSE125866 <- function() {
  fs <- list.files(file.path(BULK, "GSE125866"), "FPKM.xls.gz$", recursive = TRUE, full.names = TRUE)
  l <- lapply(fs, function(f) { d <- fread(f, header = FALSE); setNames(d$V2, d$V1) })
  g <- Reduce(intersect, lapply(l, names)); cnt <- sapply(l, function(x) x[g]); rownames(cnt) <- g
  colnames(cnt) <- sub("\\.gene.*", "", basename(fs))
  meta <- data.frame(sample = colnames(cnt), group = basename(dirname(fs)), id = sub("^GSM\\d+_", "", colnames(cnt)))
  meta$patient <- ifelse(meta$group == "Normal", meta$id, sub("[NT]$", "", meta$id))
  cnt <- cnt[rowSums(cnt >= 10) >= 3, ]
  list(name = "GSE125866", counts = cnt, expr = .logcpm(cnt), meta = meta)
}
load_GSE274203 <- function() {
  d <- fread(file.path(BULK, "GSE274203", "GSE274203_RNA_seq_Read_counts.txt.gz"))
  cnt <- as.matrix(d[, -1]); cnt <- .collapse(cnt, .ens2sym()[d[[1]]]); cnt <- cnt[rowSums(cnt >= 10) >= 2, ]
  grp <- c(NOM = "Normal", OSF = "OSF", OSCC = "OSCC")[sub("_\\d+$", "", colnames(cnt))]
  list(name = "GSE274203", counts = cnt, expr = .logcpm(cnt), meta = data.frame(sample = colnames(cnt), group = grp))
}
load_GSE319719 <- function() {
  d <- fread(file.path(BULK, "GSE319719", "GSE319719_Raw_counts.txt.gz"))
  cnt <- as.matrix(d[, -1]); cnt <- .collapse(cnt, .ens2sym()[d[[1]]]); cnt <- cnt[rowSums(cnt >= 10) >= 3, ]
  grp <- ifelse(grepl("OSMF", colnames(cnt)), "OSF_OSCC", "deNovo_OSCC")
  list(name = "GSE319719", counts = cnt, expr = .logcpm(cnt), meta = data.frame(sample = colnames(cnt), group = grp))
}
load_GSE64216 <- function() {
  s <- .series_matrix(file.path(BULK, "GSE64216", "GSE64216_series_matrix.txt.gz"))
  ann <- .illumina_ann(file.path(BULK, "GSE64216", "Platform_annotation", "GPL10558_HumanHT-12_V4_0_R2_15002873_B.txt.gz"))
  m <- log2(pmax(s$mat, 1)); m <- limma::normalizeBetweenArrays(m, "quantile"); m <- .collapse(m, ann[rownames(m)])
  title <- .hdr_field(s$hdr, "Sample_title")
  grp <- c(NOM = "Normal", OSFWT = "OSF", OSFWD = "OSF_dysplasia", OSCC = "OSF_OSCC")[sub("2$", "", title)]
  list(name = "GSE64216", expr = m, meta = data.frame(sample = colnames(m), group = grp, title = title))
}
load_GSE37991 <- function() {
  s <- .series_matrix(file.path(BULK, "GSE37991", "GSE37991_series_matrix.txt.gz"))
  ann <- .illumina_ann(file.path(BULK, "GSE37991", "Platform_annotation", "GPL6883_HumanRef-8_V3_0_R0_11282963_A.bgx.gz"))
  m <- s$mat; if (max(m, na.rm = TRUE) > 50) m <- log2(pmax(m, 1)); m <- .collapse(m, ann[rownames(m)])
  title <- .hdr_field(s$hdr, "Sample_title")
  list(name = "GSE37991", expr = m, meta = data.frame(sample = colnames(m), group = ifelse(grepl("^T", title), "OSCC", "Normal"),
                                                      patient = sub("^[TN](\\d+).*", "\\1", title)))
}
load_GSE20170 <- function() {   # Agilent two-colour: OSF vs pooled normal log10 ratios -> log2
  fs <- list.files(file.path(BULK, "GSE20170", "OSF_vs_PooledNormal"), "txt.gz$", full.names = TRUE)
  l <- lapply(fs, function(f) {
    x <- readLines(gzfile(f)); h <- grep("^FEATURES", x)
    d <- fread(text = x[h:length(x)], select = c("ControlType", "GeneName", "LogRatio"))
    d <- d[ControlType == 0]; v <- tapply(d$LogRatio * log2(10), d$GeneName, mean); v
  })
  g <- Reduce(intersect, lapply(l, names)); m <- sapply(l, function(v) v[g]); rownames(m) <- g
  colnames(m) <- sub("\\.txt\\.gz$", "", basename(fs))
  list(name = "GSE20170", expr = m, meta = data.frame(sample = colnames(m), group = "OSF_vs_Normal_ratio"))
}
load_TCGA_oral <- function() {
  td <- file.path(ROOT, "data/TCGA")
  cl <- fread(file.path(td, "TCGA-HNSC.clinical.tsv.gz"))
  site_col <- intersect(c("tissue_or_organ_of_origin.diagnoses", "tissue_or_organ_of_origin"), names(cl))[1]
  oral_sites <- c("Tongue, NOS", "Border of tongue", "Ventral surface of tongue, NOS", "Anterior 2/3 of tongue, NOS",
                  "Floor of mouth, NOS", "Anterior floor of mouth", "Overlapping lesion of floor of mouth",
                  "Cheek mucosa", "Mouth, NOS", "Overlapping lesion of other and unspecified parts of mouth",
                  "Gum, NOS", "Lower gum", "Upper gum", "Upper Gum", "Hard palate", "Retromolar area", "Lip, NOS",
                  "Overlapping lesion of lip", "Lower lip, inner aspect", "Upper lip, NOS", "Mucosa of lower lip",
                  "Palate, NOS", "Vestibule of mouth")
  # "Overlapping lesion of lip, oral cavity and pharynx" (C14.8) and base of tongue (oropharynx) are excluded
  oral <- unique(substr(cl$sample[cl[[site_col]] %in% oral_sites], 1, 12))
  tp <- fread(file.path(td, "TCGA-HNSC.star_tpm.tsv.gz"))
  m <- as.matrix(tp[, -1]); m <- .collapse(m, .ens2sym()[sub("\\..*", "", tp[[1]])])
  keep <- substr(colnames(m), 1, 12) %in% oral & substr(colnames(m), 14, 15) %in% c("01", "11")
  m <- m[, keep]
  list(name = "TCGA_oral", expr = m,
       meta = data.frame(sample = colnames(m), group = ifelse(substr(colnames(m), 14, 15) == "01", "OSCC", "Normal"),
                         patient = substr(colnames(m), 1, 12)))
}
load_all_bulk <- function(which = c("GSE125866", "GSE64216", "GSE274203", "GSE319719", "GSE37991", "GSE20170", "TCGA_oral")) {
  fun <- list(GSE125866 = load_GSE125866, GSE64216 = load_GSE64216, GSE274203 = load_GSE274203, GSE319719 = load_GSE319719,
              GSE37991 = load_GSE37991, GSE20170 = load_GSE20170, TCGA_oral = load_TCGA_oral)
  out <- list()
  for (w in which) { r <- try(fun[[w]](), silent = TRUE); if (inherits(r, "try-error")) message("skip ", w, ": ", r) else out[[w]] <- r }
  out
}
