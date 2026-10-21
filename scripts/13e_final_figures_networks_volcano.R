# Redraw circle networks with radial (non-overlapping) labels, a readable chord diagram, and volcano plots
# with non-overlapping gene labels. Replaces: Fig6A_Network_interaction_number, Fig6A_supp_Network_all_interactions,
# Fig6C_Chord_LR_keysubset_OSF_OSCC, Fig3E_*, Fig7A_*
source("utils.R")
suppressPackageStartupMessages(library(CellChat))

key <- readLines(file.path(TAB, "Step3_key_subset.txt"))[1]
cc_list <- readRDS(file.path(RDS, "06_cellchat_list.rds"))
groups <- Reduce(union, lapply(cc_list, function(x) rownames(x@net$count)))
pal <- setNames(rep(c(pal_many, "#A20056", "#008B45", "#631879", "#5F559B", "#BB0021", "#008280",
                      "#925E9F", "#FDAF91", "#AD002A", "#ADB6B6", "#1B1919"), 3)[seq_along(groups)], groups)

# ---- circle network: nodes on a circle, labels rotated along the radius ----
circle_net <- function(m, title, only_sender = NULL, wmax = NULL) {
  n <- length(groups); th <- pi / 2 - 2 * pi * (seq_len(n) - 1) / n
  nd <- data.frame(g = groups, x = cos(th), y = sin(th), th = th)
  nd$size <- sapply(groups, function(g) if (g %in% rownames(m)) sum(m[g, ]) + sum(m[, g]) else 0)
  nd$ang <- th * 180 / pi; nd$left <- cos(th) < 0
  nd$ang <- ifelse(nd$left, nd$ang + 180, nd$ang)
  ed <- as.data.frame(as.table(m)); colnames(ed) <- c("from", "to", "w")
  ed <- subset(ed, w > 0 & from != to)
  if (!is.null(only_sender)) ed <- subset(ed, from == only_sender)
  ed <- cbind(ed, setNames(nd[match(ed$from, nd$g), c("x", "y")], c("x", "y")),
              setNames(nd[match(ed$to, nd$g), c("x", "y")], c("xend", "yend")))
  wmax <- if (is.null(wmax)) max(ed$w) else wmax
  ggplot() +
    geom_curve(data = ed, aes(x, y, xend = xend, yend = yend, colour = from, linewidth = w),
               curvature = 0.25, alpha = if (is.null(only_sender)) 0.35 else 0.8,
               arrow = arrow(length = unit(0.12, "cm"), type = "closed"), show.legend = FALSE) +
    geom_point(data = nd, aes(x, y, fill = g, size = size), shape = 21, colour = "grey20", show.legend = FALSE) +
    geom_text(data = nd, aes(x * 1.12, y * 1.12, label = relab(g), angle = ang, hjust = ifelse(left, 1, 0)), size = 3.2) +
    scale_colour_manual(values = pal) + scale_fill_manual(values = pal) +
    scale_linewidth(range = c(0.2, 2.5), limits = c(0, wmax)) + scale_size(range = c(2, 7)) +
    coord_fixed(xlim = c(-1.9, 1.9), ylim = c(-1.9, 1.9)) + theme_void() + ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
}
grp_name <- c(deNovo_OSCC = "De novo OSCC", OSF_OSCC = "OSF-associated OSCC")
wk <- max(sapply(cc_list, function(x) max(x@net$count[key, ])))
pa <- lapply(names(cc_list), function(n) circle_net(cc_list[[n]]@net$count,
        sprintf("%s: interactions sent by %s\n(n = %d)", grp_name[n], key, sum(cc_list[[n]]@net$count[key, ])), key, wk))
save_plot(wrap_plots(pa, nrow = 1), "Fig6A_Network_interaction_number", 15, 7.5)
wa <- max(sapply(cc_list, function(x) max(x@net$count)))
pb <- lapply(names(cc_list), function(n) circle_net(cc_list[[n]]@net$count,
        sprintf("%s: all inferred interactions\n(n = %d)", grp_name[n], sum(cc_list[[n]]@net$count)), NULL, wa))
save_plot(wrap_plots(pb, nrow = 1), "Fig6A_supp_Network_all_interactions", 15, 7.5)

# ---- dot plot of the 15 most strongly gained L-R pairs (replaces the crowded chord diagram) ----
lr <- fread(file.path(TAB, "Step5_keysubset_to_fibroblast_LR_OSF_vs_deNovo.csv"))
top <- head(unique(lr[delta > 0][order(-delta)]$interaction_name), 15)
fib_order <- c("qFib", "iCAF", "POSTN+ myCAF", "MMP1+ myCAF", "Prolif. CAF", "Pericyte", "SMC")
dd <- lr[interaction_name %in% top]
dd$pair_lab <- paste(gsub("_", "/", dd$ligand), "->", gsub("_", "/", dd$receptor))
lab_order <- dd$pair_lab[match(top, dd$interaction_name)]
dd$pair <- factor(dd$pair_lab, levels = rev(unique(lab_order)))
dd$target <- factor(dd$target, intersect(fib_order, unique(dd$target)))
dd$prob_OSF[dd$prob_OSF == 0] <- NA
pdot <- ggplot(dd, aes(target, pair)) +
  geom_point(aes(size = prob_OSF, colour = prob_OSF)) +
  geom_point(data = subset(dd, prob_deNovo > 0), aes(target, pair), shape = 4, size = 2.5, colour = "grey30") +
  scale_colour_gradient(low = "#FDE0D9", high = "#B2182B", name = "Probability\n(OSF-assoc.)", na.value = NA) +
  scale_size(range = c(1.5, 7), name = "Probability\n(OSF-assoc.)") + guides(size = "none") +
  labs(x = "Receiver fibroblast subtype", y = NULL,
       title = paste("Top 15 L-R pairs gained in OSF-associated OSCC\n", key, "-> fibroblasts"),
       caption = "x = also detected in de novo OSCC; empty = not significant in OSF-associated OSCC") +
  theme_bw() + umap_theme + theme(axis.text.x = element_text(angle = 40, hjust = 1), plot.caption = element_text(size = 8))
save_plot(pdot, "Fig6C_Chord_LR_keysubset_OSF_OSCC", 8, 7)

# ---- volcano plots with non-overlapping labels ----
volcano2 <- function(de, title, n_lab = 8, fc = 0.5, p = 0.05) {
  de$sig <- ifelse(de$p_val_adj < p & de$avg_log2FC > fc, "Up", ifelse(de$p_val_adj < p & de$avg_log2FC < -fc, "Down", "NS"))
  de$y <- -log10(pmax(de$p_val_adj, 1e-300))
  amb <- "^(HB[ABDGMQZ]|IG[HKL]|JCHAIN|MT-|RP[SL]|GNLY|NKG7|A[CLP]\\d{6}|LINC\\d)"
  lab <- bind_rows(de |> filter(sig == "Up", !grepl(amb, gene)) |> arrange(p_val_adj, desc(avg_log2FC)) |> head(n_lab),
                   de |> filter(sig == "Down", !grepl(amb, gene)) |> arrange(p_val_adj, avg_log2FC) |> head(n_lab))
  ggplot(de, aes(avg_log2FC, y, colour = sig)) + geom_point(size = 0.8, alpha = 0.7) +
    scale_colour_manual(values = c(Up = "#E64B35", Down = "#3C5488", NS = "grey75"), name = NULL) +
    geom_vline(xintercept = c(-fc, fc), linetype = 2, colour = "grey50") + geom_hline(yintercept = -log10(p), linetype = 2, colour = "grey50") +
    geom_text_repel(data = lab, aes(label = gene), size = 3, colour = "black", box.padding = 0.6, point.padding = 0.3,
                    min.segment.length = 0, max.overlaps = Inf, force = 3, seed = 1) +
    labs(x = "log2 fold change", y = "-log10 adjusted P", title = title) + theme_classic() + umap_theme
}
d <- fread(file.path(TAB, "SuppTable5_fibroblast_DEG_OSF_OSCC_vs_deNovo.csv"))
save_plot(volcano2(d, "Fibroblasts: OSF-associated vs de novo OSCC"), "Fig3E_Volcano_fibroblast_OSF_OSCC_vs_deNovo", 6.5, 5.5)
d <- fread(file.path(TAB, "SuppTable5b_fibroblast_DEG_OSCC_vs_OSF.csv"))
save_plot(volcano2(d, "Fibroblasts: OSCC vs OSF"), "Fig3E_supp_Volcano_OSCC_vs_OSF", 6.5, 5.5)
d <- fread(file.path(TAB, "SuppTable12_DEG_keymacro_ligand_high_vs_low.csv"))
save_plot(volcano2(d, paste0(key, ": IL1B-high vs IL1B-low samples")), "Fig7A_Volcano_ligand_high_vs_low", 6.5, 5.5)
cat("networks, chord and volcano plots redrawn\n")
