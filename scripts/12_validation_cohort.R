# Fig S1: validation scRNA cohort GSE103322 (HNSCC, Smart-seq2) - same cell types and marker patterns
# + Supplementary Table 2 (all datasets used)
source("utils.R")

d <- fread(file.path(SC, "GSE103322", "GSE103322_HNSCC_all_data.txt.gz"))
meta <- as.data.frame(t(d[1:5, -1])); colnames(meta) <- c("maxima", "lymph_node", "cancer", "noncancer", "celltype_raw")
m <- as.matrix(d[-(1:5), -1]); storage.mode(m) <- "numeric"; rownames(m) <- make.unique(gsub("'", "", d[[1]][-(1:5)])); rm(d)
meta$celltype <- ifelse(meta$cancer == "1", "Malignant epithelial", meta$celltype_raw)
meta$celltype <- recode(meta$celltype, `-Fibroblast` = "Fibroblast", Macrophage = "Myeloid", Dendritic = "Myeloid",
                        `B cell` = "B_cell", `T cell` = "T_NK", myocyte = "Myocyte", Mast = "Mast")
meta$tumor <- sub("_.*", "", rownames(meta))
keep <- meta$celltype != "0" & meta$lymph_node == "0"            # primary tumours, annotated cells
so <- CreateSeuratObject(counts = as(m[, keep], "dgCMatrix"), meta.data = meta[keep, ])
so <- SetAssayData(so, layer = "data", new.data = as(m[, keep], "dgCMatrix"))
so <- FindVariableFeatures(so, nfeatures = 3000) |> ScaleData() |> RunPCA(npcs = 30, verbose = FALSE)
so <- RunUMAP(so, dims = 1:20, verbose = FALSE)

pal <- setNames(pal_many[seq_along(unique(so$celltype))], sort(unique(so$celltype)))
pA <- umap_plot(so, "umap", "celltype", pal, "Validation cohort GSE103322: cell types", label = TRUE, size = 0.6)
mk <- c("KRT5", "KRT14", "KRT17", "COL1A1", "DCN", "LUM", "PECAM1", "VWF", "RGS5", "ACTA2", "CD68", "CD14", "LYZ", "CD163",
        "CD3D", "CD2", "NKG7", "MS4A1", "CD79A", "TPSAB1", "CPA3", "CXCL5", "IL1B", "ETS2", "MMP1", "IL1R1")
pB <- DotPlot(so, features = intersect(mk, rownames(so)), group.by = "celltype", cols = c("lightgrey", "#E64B35")) + RotatedAxis() +
  theme(axis.title = element_blank()) + ggtitle("Canonical markers and axis genes") + umap_theme
pr <- as.data.frame(table(so$tumor, so$celltype)); colnames(pr) <- c("tumor", "celltype", "n")
pC <- ggplot(pr, aes(tumor, n, fill = celltype)) + geom_col(position = "fill") + scale_fill_manual(values = pal) +
  labs(x = NULL, y = "Proportion", title = "Cell-type composition per tumour") + theme_classic() + umap_theme +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
save_plot(pA, "FigS1A_UMAP_validation_GSE103322", 7, 5.5)
save_plot(pB, "FigS1B_DotPlot_validation_GSE103322", 12, 4.5)
save_plot(pC, "FigS1C_Composition_validation_GSE103322", 9, 4.5)

st2 <- data.frame(
  dataset = c("GSE262446", "GSE215403", "GSE195832", "GSE310797", "GSE103322", "GSE220978", "GSE310797 (spatial)", "GSE295691",
              "GSE125866", "GSE64216", "GSE274203", "GSE319719", "GSE37991", "GSE20170", "TCGA-HNSC (oral)", "GSE156669", "TCGA-HNSC 450K (oral)"),
  data_type = c(rep("scRNA-seq (10x)", 4), "scRNA-seq (Smart-seq2)", "Visium", "Visium", "Visium", "RNA-seq", "Microarray (Illumina HT-12)",
                "RNA-seq", "RNA-seq", "Microarray (Illumina Ref-8)", "Microarray (Agilent 2-colour)", "RNA-seq (STAR TPM)",
                "Methylation 450K", "Methylation 450K"),
  samples = c("Normal 1, OSF 1", "OSF-associated OSCC 3, de novo OSCC 9", "OSCC pre-treatment 4 (post-treatment excluded)",
              "OSCC 15 (PT13 excluded: truncated file on GEO)", "HNSCC primary tumours", "OSF-associated OSCC 4", "OSCC 6", "OSCC 10 (BQ+ 4, BQ- 6)",
              "Normal 2, OSF 8, OSCC 8 (6 OSF/OSCC pairs)", "Normal 2, OSF 2, OSF+dysplasia 2, OSF-OSCC 2", "Normal 2, OSF 2, OSCC 2",
              "OSF-derived OSCC 5, de novo OSCC 8", "OSCC 40, adjacent normal 40", "OSF 10 vs pooled normal", "OSCC 245, normal 17",
              "Normal 5, OSF 7", "Oral OSCC / normal"),
  role = c("Discovery (Normal / OSF)", "Discovery (main contrast)", "Discovery (atlas)", "Discovery (atlas)", "Validation (Fig S1, S7)",
           "Spatial (Fig 2I, 7D-E, S2, S9-S10)", "Spatial (supplementary)", "Spatial (supplementary)", "Training set (ML), Fig 6E",
           "Validation (ML, Fig 6F)", "Validation (ML, Fig 6F)", "Validation (Fig 6F)", "Validation (ML, Fig 6F)", "Validation (Fig 6F)",
           "Validation (ML, Fig 6F, S11)", "Methylation (Fig S5, S11)", "Methylation (Fig S11)"))
fwrite(st2, file.path(TAB, "SuppTable2_all_datasets.csv"))
cat("Fig S1 done\n")
