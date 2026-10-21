# =============================================================================
# Monocle2 (DDRTree) runner, used for Fig 3H-J (fibroblasts) and Fig 5B-D (myeloid)
# Monocle2 calls igraph::nei() inside vertex indexing, which is defunct in igraph >= 2.1; the calls are
# rewritten to the supported .nei() at load time.
# Usage: Rscript run_monocle2.R <input.rds> <label_col> <root_label> <fig_prefix> <title>
#   input.rds = list(counts = dgCMatrix, meta = data.frame, ordering_genes = character, pal = named colours)
# =============================================================================
suppressPackageStartupMessages({ library(igraph); library(monocle); library(ggplot2); library(patchwork); library(data.table) })
ns <- asNamespace("monocle")
for (fn in ls(ns)) {
  f <- get(fn, envir = ns)
  if (is.function(f) && any(grepl("suppressWarnings(nei(", deparse(f), fixed = TRUE))) {
    body(f) <- parse(text = gsub("suppressWarnings(nei(", "(.nei(", deparse(body(f)), fixed = TRUE))[[1]]
    environment(f) <- ns
    assignInNamespace(fn, f, "monocle"); message("patched monocle::", fn)
  }
}

a <- commandArgs(trailingOnly = TRUE)
inp <- a[1]; label_col <- a[2]; root_label <- a[3]; prefix <- a[4]; title <- a[5]
OUT <- file.path(Sys.getenv("OSF_OSCC_ROOT", normalizePath("..")), "output"); FIG <- file.path(OUT, "figures"); TAB <- file.path(OUT, "tables")
save_plot <- function(p, name, w, h) {
  ggsave(file.path(FIG, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(FIG, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}
th <- theme(plot.title = element_text(hjust = 0.5, face = "bold"))

x <- readRDS(inp)
pd <- new("AnnotatedDataFrame", data = x$meta)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(x$counts), row.names = rownames(x$counts)))
cds <- newCellDataSet(x$counts, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds); cds <- estimateDispersions(cds)
cds <- setOrderingFilter(cds, x$ordering_genes)
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- orderCells(cds)
# root = state holding most cells of the progenitor-like population (highest CytoTRACE2 / quiescent)
st <- names(which.max(table(pData(cds)$State[pData(cds)[[label_col]] == root_label])))
cds <- orderCells(cds, root_state = st)

pS <- plot_cell_trajectory(cds, color_by = "State", cell_size = 0.6) + ggtitle(paste(title, "- cell states")) + th
pP <- plot_cell_trajectory(cds, color_by = "Pseudotime", cell_size = 0.6) + scale_colour_viridis_c() + ggtitle("Pseudotime") + th
pL <- plot_cell_trajectory(cds, color_by = label_col, cell_size = 0.6) + scale_colour_manual(values = x$pal) +
  ggtitle(paste(title, "- subtypes")) + th
save_plot(pS + pP, paste0(prefix, "_trajectory_state_pseudotime"), 11, 5)
save_plot(pL, paste0(prefix, "_trajectory_subtype"), 6.5, 5.5)
save_plot(pL + facet_wrap(as.formula(paste("~ `", label_col, "`", sep = "")), nrow = 2), paste0(prefix, "_trajectory_subtype_facet"), 12, 7)
if (!is.null(pData(cds)$CytoTRACE2)) {
  pC <- plot_cell_trajectory(cds, color_by = "CytoTRACE2", cell_size = 0.6) + scale_colour_gradientn(colours = c("#3C5488", "#FDE725", "#E64B35")) +
    ggtitle("CytoTRACE2 potency") + th
  save_plot(pC, paste0(prefix, "_trajectory_CytoTRACE2"), 6.5, 5.5)
}

dg <- differentialGeneTest(cds[x$ordering_genes, ], fullModelFormulaStr = "~sm.ns(Pseudotime)", cores = 1)
dg <- dg[order(dg$qval), ]
fwrite(dg, file.path(TAB, paste0(prefix, "_pseudotime_DEG.csv")))
sig <- head(dg$gene_short_name, 100)
hm <- plot_pseudotime_heatmap(cds[sig, ], num_clusters = 4, cores = 1, show_rownames = TRUE, return_heatmap = TRUE)
pdf(file.path(FIG, paste0(prefix, "_pseudotime_heatmap.pdf")), 7, 13); grid::grid.newpage(); grid::grid.draw(hm$gtable); dev.off()
png(file.path(FIG, paste0(prefix, "_pseudotime_heatmap.png")), 7, 13, units = "in", res = 300); grid::grid.newpage(); grid::grid.draw(hm$gtable); dev.off()
gene_cl <- data.frame(gene = names(cutree(hm$tree_row, 4)), module = cutree(hm$tree_row, 4))
fwrite(gene_cl, file.path(TAB, paste0(prefix, "_pseudotime_gene_modules.csv")))

fwrite(cbind(cell = colnames(cds), pData(cds)[, c("State", "Pseudotime", label_col)]), file.path(TAB, paste0(prefix, "_pseudotime.csv")))
saveRDS(cds, sub("\\.rds$", "_cds.rds", inp))
cat("Monocle2 done:", prefix, "\n")
