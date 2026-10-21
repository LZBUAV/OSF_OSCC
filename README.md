# Analysis code: CXCL5+ macrophages and an ETS2–IL1B–IL1R1 axis in the malignant transformation of oral submucous fibrosis

This repository contains the code used for the analyses in
*"Integrated single-cell, spatial and bulk multi-omics with machine learning identifies CXCL5+ macrophages and an
ETS2–IL1B–IL1R1 axis associated with the malignant transformation of oral submucous fibrosis"* [[authors, journal, year]].

All input data are public (GEO, TCGA via the UCSC Xena GDC hub, CTD, cisTarget). No new data were generated.
`scripts/00_download_data.R` downloads every input and arranges it in the folder layout expected by the scripts.

## Repository structure

```
code_submission/
├── README.md
├── scripts/            analysis scripts (run in numerical order, from this folder)
├── config/             manual curation files used by the scripts (copied to output/tables by 00_download_data.R)
├── environment/        R / Python versions and package versions used in the study
├── data/               created by 00_download_data.R (public inputs, ~6 GB)
└── output/             created by the scripts: figures/, tables/, rds/ (checkpoints), scenic/, logs/
```

## Requirements

- R 4.5.0 (package versions in `environment/R_packages.csv`; installer: `scripts/00_install_packages.R`).
  On Windows, Rtools is needed for packages installed from GitHub.
- Python 3.10 environment for pySCENIC 0.12.1 (`environment/python_pyscenic_requirements.txt`; key pins:
  numpy 1.23.5, pandas 1.5.3, numba 0.56.4, dask/distributed 2023.3.2, and `setuptools<70`):
  ```
  conda create -n pyscenic python=3.10
  conda activate pyscenic
  pip install -r environment/python_pyscenic_requirements.txt "setuptools<70"
  ```
- The study was run on Windows 11 with 64 GB RAM and 24 logical cores. The single-cell integration step needs about
  15–20 GB of RAM; the complete pipeline takes several hours.

## Environment variables (optional)

| Variable | Meaning | Default |
|---|---|---|
| `OSF_OSCC_ROOT` | repository root | parent folder of `scripts/` |
| `PYSCENIC_PYTHON` | Python executable of the pySCENIC environment | `python` |
| `SCENIC_DB` | folder with cisTarget databases | `data/cisTarget` |

## How to run

Run every script from the `scripts/` folder, e.g. `Rscript 01_scRNA_integration.R`.

| Script | Content | Main outputs |
|---|---|---|
| 00_install_packages.R | install R packages | |
| 00_download_data.R | download and arrange all public inputs; copy `config/` files | `data/` |
| 01_scRNA_integration.R | per-sample MAD quality control, Harmony integration, clustering, major cell-type annotation, Ro/e | Fig. 2a–g; Supplementary Tables 1, 3 |
| 02_metabolism_spatial_RCTD.R | KEGG metabolic GSEA per cell type; RCTD deconvolution of Visium sections | Fig. 2h,i; Supplementary Fig. 2 |
| 03_fibroblast_subtypes.R | fibroblast/mural sub-clustering, DEGs, KEGG, CytoTRACE2, Monocle2 | Fig. 3; Supplementary Tables 4, 5 |
| 04_myeloid_subtypes.R | myeloid sub-clustering, Ro/e, UCell signatures, Hallmark analysis | Fig. 4; Supplementary Table 6 |
| 05_key_subset_evidence.R | multi-evidence ranking of macrophage subsets (Ro/e, paired bulk RNA-seq, spatial) | Supplementary Table 15 |
| 06_trajectory_scenic.R | CytoTRACE2, Monocle2 and pySCENIC (calls `run_pyscenic.py`), CSI modules, RSS | Fig. 5; Supplementary Fig. 3; Supplementary Table 7 |
| 07_cellchat_axis.R | CellChat in OSF-associated vs de novo OSCC; TF–ligand–receptor axis | Fig. 6a–d; Supplementary Table 14 |
| 08_ligand_spatial.R | IL1B-high vs IL1B-low macrophages (DEG, GSEA); spatial co-localisation | Fig. 7; Supplementary Figs. 9, 10; Supplementary Table 12 |
| 09_machine_learning.R | 113 machine-learning models, ROC, permutation SHAP | Fig. 8, 9; Supplementary Table 13 |
| 10_bulk_axis_validation.R | bulk correlations, lactate genes, validation cohort GSE103322, CTD drugs | Fig. 6e,f; Supplementary Figs. 4, 6–8; Supplementary Table 16 |
| 11_methylation.R | 450K differential methylation (GSE156669) and TCGA promoter methylation | Supplementary Figs. 5, 11; Supplementary Tables 8–11 |
| 12_validation_cohort.R | validation scRNA-seq cohort GSE103322; dataset table | Supplementary Fig. 1; Supplementary Table 2 |
| 13a–13f_final_figures_*.R | redraw the published versions of figures from saved results (publication labels, legends, layout); no new analysis | final panels in `output/figures/` |

Helper files: `utils.R` (paths, plotting helpers, Ro/e, display labels), `bulk_loaders.R` (one loader per bulk dataset),
`run_monocle2.R` (Monocle2 runner, called by scripts 03 and 06), `run_pyscenic.py` (GRNBoost2–cisTarget–AUCell).

## Notes on reproducibility

- Random seeds are fixed (`set.seed(1234)`; pySCENIC `seed=1234`).
- Manual curation is recorded in `config/`: cluster-to-cell-type assignments after marker review
  (`Step1_celltype_manual.csv`, `Step2_fibroblast_subtype_manual.csv`, `Step3_myeloid_subtype_manual.csv`, with marker
  evidence for each decision; clusters labelled `Excluded` were removed as doublets or low-quality cells) and the key
  macrophage subset (`Step3_key_subset.txt`, selected by `05_key_subset_evidence.R`).
- Sample GSM9309226 (GSE310797, PT13) is excluded because its deposited count matrix is truncated.
- OSF status of GSE215403 tumours follows the original publication (concomitant OSF in patients P8, P10, P12).
- Monocle2 calls `igraph::nei()`, which is defunct in igraph >= 2.1; `run_monocle2.R` rewrites these calls to the
  equivalent `.nei()` at run time.
- Results are checkpointed in `output/rds/`; delete a checkpoint to recompute that step.

## License

[[Choose a license, e.g. MIT, before making the repository public.]]

## Contact

[[Corresponding author, e-mail]]
