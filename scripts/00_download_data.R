# =============================================================================
# Download every public input used by the analysis and arrange it in the folder layout the scripts expect.
# Run from the scripts/ folder:  Rscript 00_download_data.R
# Total download ~6 GB (largest: TCGA-HNSC 450K methylation, 1.9 GB). Re-running skips completed files.
# =============================================================================
options(timeout = 7200)
ROOT <- Sys.getenv("OSF_OSCC_ROOT", normalizePath(".."))
D <- file.path(ROOT, "data")

geo <- function(gse, sub, file) sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/%s/%s",
                                        substr(gse, 1, nchar(gse) - 3), gse, sub, file)
fetch <- function(url, dest, tries = 5) {
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(dest) && file.size(dest) > 0) return(invisible(dest))
  for (i in seq_len(tries)) {
    ok <- tryCatch({ download.file(url, paste0(dest, ".part"), mode = "wb", quiet = FALSE); TRUE }, error = function(e) FALSE)
    if (ok) { file.rename(paste0(dest, ".part"), dest); return(invisible(dest)) }
    message("retry ", i, ": ", url); Sys.sleep(5)
  }
  stop("download failed: ", url)
}
# extract a GEO RAW.tar and move each member into a sub-folder chosen by `rule(filename)` (NA = skip)
sort_tar <- function(tar, outdir, rule) {
  tmp <- file.path(tempdir(), basename(tar)); dir.create(tmp, showWarnings = FALSE)
  untar(tar, exdir = tmp)
  for (f in list.files(tmp, full.names = TRUE)) {
    sub <- rule(basename(f)); if (is.na(sub)) next
    dir.create(file.path(outdir, sub), recursive = TRUE, showWarnings = FALSE)
    file.copy(f, file.path(outdir, sub, basename(f)), overwrite = TRUE)
  }
  unlink(tmp, recursive = TRUE)
}
raw <- function(gse, area) { f <- file.path(D, area, gse, paste0(gse, "_RAW.tar")); fetch(geo(gse, "suppl", basename(f)), f); f }

# ---------------- scRNA-seq ----------------
sort_tar(raw("GSE262446", "scRNA"), file.path(D, "scRNA", "GSE262446"),
         function(f) if (grepl("^GSM8169899", f)) "Control" else if (grepl("^GSM8169898", f)) "OSF" else NA)
sort_tar(raw("GSE215403", "scRNA"), file.path(D, "scRNA", "GSE215403"), function(f) "OSCC")
pre <- c("GSM5851565", "GSM5851567", "GSM5851569", "GSM5851571")          # pre-treatment scRNA-seq only
sort_tar(raw("GSE195832", "scRNA"), file.path(D, "scRNA", "GSE195832"),
         function(f) if (sub("_.*", "", f) %in% pre) "scRNA_Pre" else NA)
fetch(geo("GSE103322", "suppl", "GSE103322_HNSCC_all_data.txt.gz"), file.path(D, "scRNA", "GSE103322", "GSE103322_HNSCC_all_data.txt.gz"))

# ---------------- spatial (GSE310797 also holds the scRNA-seq part of that study) ----------------
sort_tar(raw("GSE220978", "spatial"), file.path(D, "spatial", "GSE220978"), function(f) sub(".*_(Patient\\d)_.*", "\\1", f))
sort_tar(raw("GSE310797", "spatial"), file.path(D, "spatial", "GSE310797"), function(f) if (grepl("_SP_", f)) "Spatial" else "scRNA")

# ---------------- bulk transcriptomes ----------------
sort_tar(raw("GSE125866", "bulk"), file.path(D, "bulk", "GSE125866"), function(f) {
  id <- sub("^GSM\\d+_([^.]+)\\..*", "\\1", f)
  if (id %in% c("0N1", "ON2")) "Normal" else if (grepl("N$", id)) "OSF" else "OSCC" })
for (g in c("GSE64216", "GSE37991")) {
  fetch(geo(g, "matrix", paste0(g, "_series_matrix.txt.gz")), file.path(D, "bulk", g, paste0(g, "_series_matrix.txt.gz")))
  sort_tar(raw(g, "bulk"), file.path(D, "bulk", g), function(f) if (grepl("^GPL", f)) "Platform_annotation" else NA)
}
fetch(geo("GSE274203", "suppl", "GSE274203_RNA_seq_Read_counts.txt.gz"), file.path(D, "bulk", "GSE274203", "GSE274203_RNA_seq_Read_counts.txt.gz"))
fetch(geo("GSE319719", "suppl", "GSE319719_Raw_counts.txt.gz"), file.path(D, "bulk", "GSE319719", "GSE319719_Raw_counts.txt.gz"))
sort_tar(raw("GSE20170", "bulk"), file.path(D, "bulk", "GSE20170"), function(f) "OSF_vs_PooledNormal")

# ---------------- DNA methylation ----------------
fetch(geo("GSE156669", "matrix", "GSE156669_series_matrix.txt.gz"), file.path(D, "methylation", "GSE156669", "GSE156669_series_matrix.txt.gz"))
sort_tar(raw("GSE156669", "methylation"), file.path(D, "methylation", "GSE156669"),
         function(f) if (grepl("^GPL13534", f)) "Platform_annotation" else NA)

# ---------------- TCGA-HNSC (UCSC Xena GDC hub) ----------------
x <- "https://gdc-hub.s3.us-east-1.amazonaws.com/download/"
for (f in c("TCGA-HNSC.clinical.tsv.gz", "TCGA-HNSC.star_tpm.tsv.gz", "TCGA-HNSC.methylation450.tsv.gz", "gencode.v36.annotation.gtf.gene.probemap"))
  fetch(paste0(x, f), file.path(D, "TCGA", f))

# ---------------- Comparative Toxicogenomics Database ----------------
for (f in c("CTD_chem_gene_ixns.tsv.gz", "CTD_chemicals_diseases.tsv.gz"))
  fetch(paste0("https://ctdbase.org/reports/", f), file.path(D, "CTD", f))

# ---------------- cisTarget resources for pySCENIC ----------------
cis <- "https://resources.aertslab.org/cistarget"
fetch(paste0(cis, "/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"),
      file.path(D, "cisTarget", "hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"))
fetch(paste0(cis, "/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"),
      file.path(D, "cisTarget", "hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather"))
fetch(paste0(cis, "/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"), file.path(D, "cisTarget", "motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"))
fetch("https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_tfs.txt", file.path(D, "cisTarget", "hs_hgnc_tfs.txt"))

# ---------------- manual annotation / decision files used by the scripts ----------------
tab <- file.path(ROOT, "output", "tables"); dir.create(tab, recursive = TRUE, showWarnings = FALSE)
file.copy(list.files(file.path(ROOT, "config"), full.names = TRUE), tab, overwrite = TRUE)
cat("All inputs downloaded to", D, "\n")
