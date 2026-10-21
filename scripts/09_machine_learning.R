# =============================================================================
# Step 7-8: 113 machine-learning models on the key-macrophage signature + SHAP interpretation
# Reproduces: Fig 8A-H, Fig 9A-D
# Task: OSCC (1) vs non-malignant mucosa (OSF + normal, 0) - malignant transformation of OSF
# Training: GSE125866 (OSF 8 / OSCC 8 / normal 2); validation: GSE64216, GSE274203, GSE37991, TCGA oral
# Features: CXCL5+ (key subset) signature genes, z-scored within each cohort
# 113 models = 22 single-algorithm variants + 91 combinations (7 feature selectors x 13 classifiers)
# =============================================================================
source("utils.R")
source(file.path(ROOT, "scripts", "bulk_loaders.R"))
suppressPackageStartupMessages({ library(glmnet); library(pROC); library(randomForest); library(xgboost); library(gbm)
  library(e1071); library(MASS); library(naivebayes); library(mboost); library(plsRglm) })
set.seed(1234)

key <- readLines(file.path(TAB, "Step3_key_subset.txt"))[1]
sig <- readRDS(file.path(RDS, "04b_macrophage_subset_signatures.rds"))[[key]]
mye <- readRDS(file.path(RDS, "04_myeloid.rds")); Idents(mye) <- "subtype"
mk <- FindMarkers(mye, ident.1 = key, only.pos = TRUE, logfc.threshold = 0.5, min.pct = 0.25, max.cells.per.ident = 2000, verbose = FALSE)
mk <- mk[!grepl("^(RP[SL]|MT-|HB[AB]|IG[HKL])", rownames(mk)), ]
cand <- unique(c(sig, head(rownames(mk)[order(-mk$avg_log2FC)], 60))); rm(mye); invisible(gc())

f_b <- file.path(RDS, "10_bulk_cohorts.rds")
if (!file.exists(f_b)) { bulk <- load_all_bulk(c("GSE125866", "GSE64216", "GSE274203", "GSE37991")); bulk$TCGA_oral <- get_TCGA_oral(); saveRDS(bulk, f_b) } else bulk <- readRDS(f_b)
lab_fun <- function(g) ifelse(g %in% c("OSCC", "OSF_OSCC"), 1, ifelse(g %in% c("Normal", "OSF", "OSF_dysplasia"), 0, NA))
genes <- Reduce(intersect, c(list(cand), lapply(bulk, function(b) rownames(b$expr))))
message("Candidate genes present in all cohorts: ", length(genes))
prep <- function(b) {
  y <- lab_fun(b$meta$group); keep <- !is.na(y)
  x <- t(b$expr[genes, keep, drop = FALSE]); x <- scale(x); x[is.na(x)] <- 0
  list(x = x, y = y[keep], group = b$meta$group[keep])
}
D <- lapply(bulk, prep)
cohorts <- c("GSE125866", "GSE64216", "GSE274203", "GSE37991", "TCGA_oral")
sapply(D, function(d) table(d$y))
tr <- D$GSE125866

# ---- RFE on the training set (random forest importance, repeated) ----
imp <- rowMeans(sapply(1:50, function(i) importance(randomForest(tr$x, factor(tr$y), ntree = 500), type = 2)[, 1]))
imp <- sort(imp, decreasing = TRUE)
rfe_auc <- sapply(seq(3, min(30, length(imp))), function(k) {
  g <- names(imp)[1:k]; f <- cv.glmnet(tr$x[, g], tr$y, family = "binomial", alpha = 0.5, nfolds = 5, type.measure = "deviance")
  as.numeric(auc(tr$y, as.numeric(predict(f, tr$x[, g], s = "lambda.min")), quiet = TRUE))
})
k_best <- seq(3, min(30, length(imp)))[which.max(rfe_auc - 0.002 * seq_along(rfe_auc))]
feat <- names(imp)[1:max(k_best, 8)]
message("RFE retained ", length(feat), " genes: ", paste(feat, collapse = ", "))
fwrite(data.frame(gene = names(imp), RF_importance = imp, selected = names(imp) %in% feat), file.path(TAB, "Step7_RFE_importance.csv"))

# xgboost >= 3.0 no longer accepts xgboost(data=, label=) for binary:logistic -> use xgb.train
fit_xgb <- function(x, y) xgb.train(params = list(objective = "binary:logistic", max_depth = 2, eta = 0.1, nthread = 1),
                                    data = xgb.DMatrix(x, label = y), nrounds = 100, verbose = 0)

# ---- learners: fit(x, y) -> predict(newx) returning probability-like scores ----
fit_learner <- function(name, x, y) {
  yf <- factor(y, c(0, 1)); df <- data.frame(y = yf, x, check.names = TRUE)
  if (grepl("^Enet", name) || name %in% c("Lasso", "Ridge")) {
    a <- switch(name, Lasso = 1, Ridge = 0, as.numeric(sub("Enet\\[a=(.*)\\]", "\\1", name)))
    f <- cv.glmnet(x, y, family = "binomial", alpha = a, nfolds = min(5, min(table(y))), type.measure = "deviance")
    return(function(nx) as.numeric(predict(f, nx, s = "lambda.min", type = "response")))
  }
  if (grepl("^Stepglm", name)) {
    dir <- sub("Stepglm\\[(.*)\\]", "\\1", name)
    f0 <- suppressWarnings(glm(y ~ ., data = df, family = binomial))
    f <- suppressWarnings(stepAIC(if (dir == "forward") glm(y ~ 1, data = df, family = binomial) else f0,
                                  scope = list(upper = formula(f0), lower = ~1), direction = dir, trace = 0, steps = 50))
    return(function(nx) suppressWarnings(as.numeric(predict(f, data.frame(nx, check.names = TRUE), type = "response"))))
  }
  switch(name,
    glmBoost = { f <- glmboost(y ~ ., data = df, family = Binomial(), control = boost_control(mstop = 200))
                 function(nx) as.numeric(predict(f, data.frame(nx, check.names = TRUE), type = "response")) },
    plsRglm = { f <- plsRglm(y, x, nt = min(3, ncol(x)), modele = "pls-glm-logistic", verbose = FALSE)
                function(nx) as.numeric(predict(f, newdata = nx, type = "response")) },
    RF = { f <- randomForest(x, yf, ntree = 500); function(nx) predict(f, nx, type = "prob")[, "1"] },
    GBM = { f <- gbm.fit(x, y, distribution = "bernoulli", n.trees = 300, interaction.depth = 2, shrinkage = 0.05,
                         n.minobsinnode = 2, bag.fraction = 0.8, verbose = FALSE)
            function(nx) predict(f, nx, n.trees = 300, type = "response") },
    XGBoost = { f <- fit_xgb(x, y); function(nx) predict(f, xgb.DMatrix(nx[, colnames(x), drop = FALSE])) },
    SVM = { f <- svm(x, yf, probability = TRUE, kernel = "radial")
            function(nx) attr(predict(f, nx, probability = TRUE), "probabilities")[, "1"] },
    LDA = { f <- lda(x, yf); function(nx) predict(f, nx)$posterior[, "1"] },
    NaiveBayes = { f <- naive_bayes(x, yf, usekernel = FALSE); function(nx) predict(f, nx, type = "prob")[, "1"] })
}
selected_genes <- function(sel, x, y) {        # feature-selecting algorithms
  g <- switch(sel,
    Lasso = { f <- cv.glmnet(x, y, family = "binomial", alpha = 1, nfolds = 5); b <- coef(f, s = "lambda.min")[-1, 1]; names(b)[b != 0] },
    `Enet[a=0.5]` = { f <- cv.glmnet(x, y, family = "binomial", alpha = 0.5, nfolds = 5); b <- coef(f, s = "lambda.min")[-1, 1]; names(b)[b != 0] },
    `Stepglm[both]` = { df <- data.frame(y = factor(y), x, check.names = TRUE); f0 <- suppressWarnings(glm(y ~ ., data = df, family = binomial))
                        f <- suppressWarnings(stepAIC(f0, direction = "both", trace = 0, steps = 50)); colnames(x)[make.names(colnames(x)) %in% attr(terms(f), "term.labels")] },
    glmBoost = { f <- glmboost(factor(y) ~ ., data = data.frame(y = y, x, check.names = TRUE), family = Binomial(), control = boost_control(mstop = 200))
                 v <- names(coef(f))[-1]; colnames(x)[make.names(colnames(x)) %in% v] },
    RF = { im <- importance(randomForest(x, factor(y), ntree = 500), type = 2)[, 1]; names(sort(im, TRUE))[1:max(3, ceiling(length(im) / 2))] },
    XGBoost = { xgb.importance(model = fit_xgb(x, y))$Feature },
    plsRglm = { f <- plsRglm(y, x, nt = min(3, ncol(x)), modele = "pls-glm-logistic", verbose = FALSE)
                w <- abs(f$wwetoile[, 1]); names(sort(w, TRUE))[1:max(3, ceiling(length(w) / 2))] })
  g <- intersect(g, colnames(x)); if (length(g) < 2) g <- colnames(x)[1:min(5, ncol(x))]; g
}

singles <- c("Lasso", "Ridge", paste0("Enet[a=", seq(0.1, 0.9, 0.1), "]"), "Stepglm[both]", "Stepglm[backward]", "Stepglm[forward]",
             "glmBoost", "plsRglm", "RF", "GBM", "XGBoost", "SVM", "LDA", "NaiveBayes")                       # 22
selectors <- c("Lasso", "Enet[a=0.5]", "Stepglm[both]", "glmBoost", "RF", "XGBoost", "plsRglm")                # 7
classifiers <- c("Lasso", "Ridge", "Enet[a=0.5]", "Enet[a=0.9]", "Stepglm[both]", "glmBoost", "plsRglm",
                 "RF", "GBM", "XGBoost", "SVM", "LDA", "NaiveBayes")                                            # 13
models <- c(lapply(singles, function(s) list(name = s, sel = NA, clf = s)),
            unlist(lapply(selectors, function(s) lapply(classifiers, function(c) list(name = paste(s, "+", c), sel = s, clf = c))), recursive = FALSE))
stopifnot(length(models) == 113)

auc_of <- function(y, p) if (length(unique(y)) < 2 || any(!is.finite(p))) NA else as.numeric(auc(y, p, quiet = TRUE, direction = "<"))
cv_auc <- function(m, x, y, k = 5, reps = 5) {          # repeated stratified CV on the training set
  mean(sapply(1:reps, function(r) {
    set.seed(r); fold <- integer(length(y))
    for (cl in 0:1) { i <- which(y == cl); fold[i] <- sample(rep(1:k, length.out = length(i))) }
    p <- numeric(length(y))
    for (f in 1:k) {
      te <- fold == f; g <- if (is.na(m$sel)) colnames(x) else selected_genes(m$sel, x[!te, ], y[!te])
      pr <- try(fit_learner(m$clf, x[!te, g, drop = FALSE], y[!te]), silent = TRUE)
      p[te] <- if (inherits(pr, "try-error")) NA else pr(x[te, g, drop = FALSE])
    }
    auc_of(y, p)
  }), na.rm = TRUE)
}
f_res <- file.path(RDS, "10_ml_results.rds")
if (!file.exists(f_res)) {
  res <- lapply(models, function(m) {
    message(m$name)
    g <- if (is.na(m$sel)) feat else tryCatch(selected_genes(m$sel, tr$x[, feat], tr$y), error = function(e) feat)
    pr <- try(fit_learner(m$clf, tr$x[, g, drop = FALSE], tr$y), silent = TRUE)
    a <- sapply(cohorts[-1], function(cn) if (inherits(pr, "try-error")) NA else auc_of(D[[cn]]$y, tryCatch(pr(D[[cn]]$x[, g, drop = FALSE]), error = function(e) NA)))
    cv <- tryCatch(cv_auc(m, tr$x[, feat], tr$y), error = function(e) NA)
    list(name = m$name, genes = g, auc = c(GSE125866_CV = cv, a))
  })
  saveRDS(res, f_res)
} else res <- readRDS(f_res)
A <- do.call(rbind, lapply(res, function(r) r$auc)); rownames(A) <- sapply(res, `[[`, "name")
A <- cbind(A, Mean = rowMeans(A, na.rm = TRUE)); A <- A[order(-A[, "Mean"]), ]
fwrite(data.frame(model = rownames(A), round(A, 3), n_genes = sapply(res, function(r) length(r$genes))[match(rownames(A), sapply(res, `[[`, "name"))], check.names = FALSE),
       file.path(TAB, "Step7_113_models_AUC.csv"))

# Fig 8A: AUC heat map of the 113 models
lab_mat <- matrix(sprintf("%.2f", A), nrow(A)); lab_mat[is.na(A)] <- ""
A_plot <- A; colnames(A_plot) <- relab(colnames(A_plot))
save_pheatmap("Fig8A_Heatmap_113_models_AUC", 9, 26, A_plot, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = lab_mat,
              fontsize_number = 6, fontsize_row = 6.5, number_color = "black", gaps_col = ncol(A) - 1,
              color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100), breaks = seq(0.4, 1, length.out = 101),
              main = sprintf("AUC of 113 models\nbest: %s, mean AUC = %.3f", rownames(A)[1], A[1, "Mean"]), angle_col = 90)

best <- res[[which(sapply(res, `[[`, "name") == rownames(A)[1])]]
writeLines(best$genes, file.path(TAB, "Step7_model_genes.txt"))
message("Best model: ", best$name, " genes: ", paste(best$genes, collapse = ", "))
mb <- models[[which(sapply(models, `[[`, "name") == best$name)]]
pred <- fit_learner(mb$clf, tr$x[, best$genes, drop = FALSE], tr$y)

# Fig 8B-G: ROC of model score and of each model gene in every cohort (+ OSF vs OSCC only)
roc_panel <- function(cn, osf_only = FALSE) {
  d <- D[[cn]]; keep <- if (osf_only) d$group %in% c("OSF", "OSF_dysplasia", "OSCC", "OSF_OSCC") else rep(TRUE, length(d$y))
  if (length(unique(d$y[keep])) < 2) return(NULL)
  sc <- c(list(Model = pred(d$x[keep, best$genes, drop = FALSE])), lapply(setNames(best$genes, best$genes), function(g) d$x[keep, g]))
  rl <- lapply(sc, function(s) roc(d$y[keep], s, quiet = TRUE, direction = "auto"))
  au <- sapply(rl, function(r) as.numeric(r$auc))
  dd <- bind_rows(lapply(names(rl), function(n) data.frame(var = sprintf("%s (AUC %.2f)", n, au[n]),
                                                           fpr = 1 - rl[[n]]$specificities, tpr = rl[[n]]$sensitivities)))
  ggplot(dd, aes(fpr, tpr, colour = var)) + geom_path(linewidth = 0.7) + geom_abline(linetype = 2, colour = "grey60") +
    scale_colour_manual(values = rep(pal_many, 3)) + coord_equal() +
    labs(x = "1 - Specificity", y = "Sensitivity", colour = NULL, title = paste0(relab(cn), if (osf_only) ", OSF vs OSCC" else "", " (n = ", sum(keep), ")")) +
    theme_classic() + umap_theme + theme(legend.text = element_text(size = 7))
}
rp <- Filter(Negate(is.null), c(lapply(cohorts, roc_panel), list(roc_panel("GSE125866", TRUE), roc_panel("GSE64216", TRUE), roc_panel("GSE274203", TRUE))))
save_plot(wrap_plots(rp, ncol = 3), "Fig8BG_ROC_model_genes_all_cohorts", 18, 6 * ceiling(length(rp) / 3))

# Fig 8H: expression of model genes, OSF vs OSCC (GSE125866) and normal/OSCC (TCGA)
bx <- bind_rows(lapply(c("GSE125866", "TCGA_oral", "GSE37991"), function(cn) {
  b <- bulk[[cn]]; data.frame(cohort = cn, group = b$meta$group, t(b$expr[best$genes, , drop = FALSE]), check.names = FALSE)
})) |> tidyr::pivot_longer(-c(cohort, group), names_to = "gene", values_to = "expr")
bx$group <- factor(bx$group, c("Normal", "OSF", "OSCC"))
bx$cohort <- factor(relab(bx$cohort), relab(c("GSE125866", "GSE37991", "TCGA_oral")))
p8h <- ggplot(bx, aes(group, expr, fill = group)) + geom_boxplot(outlier.size = 0.3, width = 0.6) +
  ggpubr::stat_compare_means(comparisons = list(c("OSF", "OSCC"), c("Normal", "OSCC")), method = "wilcox.test", label = "p.signif", size = 3) +
  facet_grid(cohort ~ gene, scales = "free") + scale_fill_manual(values = pal_stage) +
  labs(x = NULL, y = "log2 expression") + theme_bw() + theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1))
save_plot(p8h, "Fig8H_Boxplot_model_genes", max(8, 2 * length(best$genes)), 9)

# -----------------------------------------------------------------------------
# Fig 9: SHAP (Kernel SHAP via fastshap, plots via shapviz)
# -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(shapviz))
X <- as.data.frame(rbind(tr$x[, best$genes, drop = FALSE], D$GSE64216$x[, best$genes, drop = FALSE], D$GSE274203$x[, best$genes, drop = FALSE]))
pf <- function(object, newdata) pred(as.matrix(newdata))
# Monte-Carlo permutation SHAP (Strumbelj & Kononenko 2014), background = all rows of X; additive by construction
perm_shap <- function(X, f, nsim = 200) {
  X <- as.matrix(X); p <- ncol(X); S <- matrix(0, nrow(X), p, dimnames = dimnames(X))
  for (i in seq_len(nrow(X))) for (s in seq_len(nsim)) {
    o <- sample(p); z <- X[sample(nrow(X), 1), ]; cur <- z; prev <- f(NULL, t(cur))
    for (j in o) { cur[j] <- X[i, j]; new <- f(NULL, t(cur)); S[i, j] <- S[i, j] + (new - prev); prev <- new }
  }
  S / nsim
}
set.seed(1234)
sh <- perm_shap(X, pf, nsim = 200)
sv <- shapviz(sh, X = X, baseline = mean(pf(NULL, X)))
fwrite(data.frame(gene = colnames(sh), mean_abs_SHAP = colMeans(abs(sh))), file.path(TAB, "Step8_SHAP_importance.csv"))
save_plot(sv_importance(sv, kind = "bar", fill = "#E64B35") + ggtitle("Mean |SHAP|") + theme_classic() + umap_theme, "Fig9A_SHAP_bar", 6, 5)
save_plot(sv_importance(sv, kind = "beeswarm") + ggtitle("SHAP value distribution") + theme_classic() + umap_theme, "Fig9B_SHAP_beeswarm", 7, 5)
i_case <- which.max(pf(NULL, X))
save_plot(sv_waterfall(sv, row_id = i_case) + ggtitle(paste("Waterfall plot:", rownames(X)[i_case])) + umap_theme, "Fig9C_SHAP_waterfall", 7, 5)
save_plot(sv_force(sv, row_id = i_case) + ggtitle(paste("Force plot:", rownames(X)[i_case])) + umap_theme, "Fig9D_SHAP_force", 9, 3.5)
cat("Step 7-8 done. Best:", best$name, "\n")
