# Figure fixes after visual QC (titles, legends, labels); regenerates figures from saved objects/tables only
source("utils.R")
suppressPackageStartupMessages({ library(png); library(UCell) })

# ---- Ro/e heat maps (two-line titles, P < 2.2e-16 instead of "p = 0") ----
obj <- readRDS(file.path(RDS, "01_sc_integrated.rds")); md <- obj@meta.data; rm(obj); invisible(gc())
plot_roe(calc_roe(md, "celltype", "stage"), "Fig2E_RoE_major_celltypes_by_stage", "Ro/e of major cell types")
m215 <- subset(md, dataset == "GSE215403"); m215$group <- droplevels(m215$group)
plot_roe(calc_roe(m215, "celltype", "group"), "Fig2F_RoE_OSF_OSCC_vs_deNovo_OSCC", "Ro/e, GSE215403")
mne <- subset(md, celltype != "Epithelial"); mne$celltype <- droplevels(mne$celltype)
plot_roe(calc_roe(mne, "celltype", "stage"), "Fig2E_supp_RoE_nonepithelial_by_stage", "Ro/e, non-epithelial cells")

fib <- readRDS(file.path(RDS, "03_fibroblast.rds"))
plot_roe(calc_roe(fib@meta.data, "subtype", "stage"), "Fig3C_RoE_fibroblast_by_stage", "Ro/e of fibroblast subtypes")
f215 <- subset(fib@meta.data, dataset == "GSE215403"); f215$group <- droplevels(f215$group)
plot_roe(calc_roe(f215, "subtype", "group"), "Fig3D_RoE_fibroblast_OSF_OSCC_vs_deNovo", "Ro/e, GSE215403")
rm(fib); invisible(gc())

mye <- readRDS(file.path(RDS, "04_myeloid.rds"))
y215 <- subset(mye@meta.data, dataset == "GSE215403"); y215$group <- droplevels(y215$group)
plot_roe(calc_roe(y215, "subtype", "group"), "Fig4C_RoE_myeloid_OSF_OSCC_vs_deNovo", "Ro/e, GSE215403", w = 5, h = 6)
plot_roe(calc_roe(mye@meta.data, "subtype", "stage"), "Fig4C_supp_RoE_myeloid_by_stage", "Ro/e of myeloid populations", w = 5, h = 6)

# ---- Fig 4F: Hallmark heat map, vertical column labels so no label is clipped ----
set.seed(1234)
cells <- unlist(lapply(split(colnames(mye), mye$subtype), function(x) sample(x, min(800, length(x)))))
ms <- subset(mye, cells = cells); rm(mye); invisible(gc())
H <- hallmark_sets(); names(H) <- make.names(names(H))
ms <- AddModuleScore_UCell(ms, features = H, ncores = 4, name = "_H")
hcols <- paste0(names(H), "_H")
hm <- aggregate(ms@meta.data[, hcols], list(subtype = ms$subtype), mean)
hmat <- t(scale(as.matrix(hm[, -1]))); colnames(hmat) <- as.character(hm$subtype); rownames(hmat) <- gsub("\\.", " ", names(H))
save_pheatmap("Fig4F_Heatmap_Hallmark_myeloid_subsets", 9, 12, hmat,
              color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100), breaks = seq(-2, 2, length.out = 101),
              main = "Hallmark gene sets across myeloid populations", angle_col = 90, border_color = "grey70")

# ---- Fig 5G: Venn with complete set names ----
tf <- fread(file.path(TAB, "SuppTable7b_key_TFs.csv"))
sets <- list(`TFs of enriched modules` = tf$TF[tf$set == "module_TFs"], `Top-10 RSS TFs` = tf$TF[tf$set == "RSS_top10"])
shared <- tf$TF[tf$set == "shared"]
p5g <- ggVennDiagram::ggVennDiagram(sets, label_alpha = 0, set_size = 4.2) +
  scale_fill_gradient(low = "white", high = "#E64B35", name = "TFs") +
  scale_x_continuous(expand = expansion(mult = 0.25)) +
  labs(title = "Key transcription factors of CXCL5+ macrophages", subtitle = paste("Shared:", paste(shared, collapse = ", "))) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), plot.subtitle = element_text(hjust = 0.5, size = 10))
save_plot(p5g, "Fig5G_Venn_key_TFs", 8, 6)

# ---- Fig 2I and supplementary maps: RCTD proportions, vertical layout, named colour bar ----
rw <- fread(file.path(TAB, "Step1b_RCTD_weights_GSE220978.csv"))
read_img <- function(pt, gsm) {
  pre <- file.path(ST, "GSE220978", pt, paste0(gsm, "_", pt, "_")); idir <- file.path(tempdir(), pt)
  untar(paste0(pre, "image.tar.gz"), exdir = idir); idir <- file.path(idir, paste0(pt, "_image"))
  pos <- fread(file.path(idir, "tissue_positions_list.csv"), header = FALSE, col.names = c("barcode", "in_tissue", "row", "col", "pxl_row", "pxl_col"))
  list(pos = pos, img = readPNG(file.path(idir, "tissue_hires_image.png")),
       sf = jsonlite::fromJSON(file.path(idir, "scalefactors_json.json"))$tissue_hires_scalef)
}
smap <- function(v, w, k, title) {
  d <- merge(w[, c("barcode", k), with = FALSE], v$pos, by = "barcode"); setnames(d, k, "value")
  d$x <- d$pxl_col * v$sf; d$y <- d$pxl_row * v$sf
  ggplot(d, aes(x, -y, colour = value)) + annotation_raster(v$img, 0, dim(v$img)[2], -dim(v$img)[1], 0) +
    geom_point(size = 0.8) + scale_colour_gradientn(colours = c("grey95", "#FDE725", "#E64B35", "#7E0000"),
                                                    limits = c(0, max(d$value)), name = "RCTD\nproportion") +
    coord_fixed(xlim = range(d$x) + c(-20, 20), ylim = -rev(range(d$y)) + c(-20, 20)) + theme_void() + ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, size = 11, face = "bold"))
}
pts <- data.frame(pt = paste0("Patient", 1:4), gsm = paste0("GSM68334", 84:87))
key_ct <- c("Myeloid", "Fibroblast", "Epithelial")
allmaps <- list()
for (i in 1:4) {
  v <- read_img(pts$pt[i], pts$gsm[i]); w <- rw[sample == pts$pt[i]]
  ps <- lapply(key_ct, function(k) smap(v, w, k, paste0("Patient ", i, ": ", relab(k))))
  if (i == 1) save_plot(wrap_plots(ps, ncol = 1), "Fig2I_Spatial_RCTD_Patient1", 5, 13)
  allmaps[[i]] <- wrap_plots(ps, nrow = 1)
}
save_plot(wrap_plots(allmaps, ncol = 1), "Fig2I_Spatial_RCTD_Myeloid_Fibroblast_Epithelial", 13, 15)
cat("figures fixed\n")
