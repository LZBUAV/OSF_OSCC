# Redraw Monocle2 trajectory plots coloured by subtype with the legend on the right (legend was clipped)
suppressPackageStartupMessages({ library(monocle); library(ggplot2); library(patchwork) })
OUT <- file.path(Sys.getenv("OSF_OSCC_ROOT", normalizePath("..")), "output"); FIG <- file.path(OUT, "figures"); RDS <- file.path(OUT, "rds")
save_plot <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}
th <- theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "right",
            legend.text = element_text(size = 9), legend.title = element_blank())
for (x in list(c("03_fibroblast_monocle_input", "Fig3HIJ", "Fibroblast subtypes"),
               c("05_myeloid_monocle_input", "Fig5BCD", "Monocyte/macrophage populations"))) {
  cds <- readRDS(file.path(RDS, paste0(x[1], "_cds.rds"))); pal <- readRDS(file.path(RDS, paste0(x[1], ".rds")))$pal
  p <- plot_cell_trajectory(cds, color_by = "subtype", cell_size = 0.6) + scale_colour_manual(values = pal) +
    guides(colour = guide_legend(ncol = 1, override.aes = list(size = 3))) + ggtitle(x[3]) + th
  save_plot(p, paste0(x[2], "_trajectory_subtype"), 8.5, 5.5)
}
cat("trajectory plots fixed\n")
