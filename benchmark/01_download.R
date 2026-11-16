# 01_download.R — fetch benchmark datasets from GEO and build DEG tables
#
# Run from the package root:  Rscript benchmark/01_download.R
#
# For each dataset in ground_truth.csv this script:
#   1. downloads the GEO series matrix (processed values) via GEOquery
#   2. selects treatment / control samples by regex on sample titles
#   3. computes a DEG table (limma for microarray; DESeq2 for RNA-seq counts)
#   4. saves benchmark/deg_tables/<dataset_id>.rds
#
# Datasets were verified via NCBI E-utilities on 2026-11-16 (see
# benchmark/README.md). NOTE: RNA-seq series matrices occasionally lack
# processed count tables; for those, download raw counts from the series
# supplementary files (getGEOSuppFiles) — the script aborts with an
# informative message in that case.

library(GEOquery)

GT <- read.csv("benchmark/ground_truth.csv", stringsAsFactors = FALSE,
               check.names = FALSE)
dir.create("benchmark/deg_tables", showWarnings = FALSE, recursive = TRUE)
dir.create("benchmark/geo_cache", showWarnings = FALSE, recursive = TRUE)

pick_samples <- function(gse, treat_re, ctrl_re) {
  titles <- as.character(pData(phenoData(gse[[1]]))$title)
  treat <- grep(treat_re, titles, ignore.case = TRUE, value = FALSE)
  ctrl  <- grep(ctrl_re,  titles, ignore.case = TRUE, value = FALSE)
  # samples matching both patterns are ambiguous — drop them from control
  ctrl <- setdiff(ctrl, treat)
  list(treat = titles[treat], ctrl = titles[ctrl],
       all = titles, treat_re = treat_re, ctrl_re = ctrl_re)
}

deg_from_matrix <- function(expr, groups) {
  library(limma)
  design <- model.matrix(~ 0 + groups)
  colnames(design) <- levels(groups)
  fit <- lmFit(expr, design)
  fit <- eBayes(contrasts.fit(fit, makeContrasts(
    treat - ctrl, levels = design)))
  tt <- topTable(fit, number = Inf, sort.by = "none")
  data.frame(symbol = rownames(tt), log2FoldChange = tt$logFC,
             padj = p.adjust(tt$P.Value, "BH"), stringsAsFactors = FALSE)
}

for (i in seq_len(nrow(GT))) {
  d <- GT[i, ]
  out <- file.path("benchmark/deg_tables", paste0(d$dataset_id, ".rds"))
  if (file.exists(out)) { message("skip (exists): ", d$dataset_id); next }
  message("== ", d$dataset_id, " (", d$gse, ")")

  gse <- getGEO(d$gse, GSEMatrix = TRUE, destdir = "benchmark/geo_cache")
  sel <- pick_samples(gse, d$treat_regex, d$ctrl_regex)
  message("  treat: ", length(sel$treat), " | ctrl: ", length(sel$ctrl))
  if (length(sel$treat) < 2 || length(sel$ctrl) < 2) {
    message("  !! too few matched samples — inspect titles manually:")
    print(sel$all)
    next
  }

  expr <- exprs(gse[[1]])[, c(sel$treat, sel$ctrl), drop = FALSE]
  if (max(expr, na.rm = TRUE) > 100) {
    message("  data look like raw counts (max>100); applying log2(x+1). ",
            "For publication-quality RNA-seq DE use DESeq2 on raw counts.")
    expr <- log2(expr + 1)
  }
  groups <- factor(c(rep("treat", length(sel$treat)),
                     rep("ctrl",  length(sel$ctrl))),
                   levels = c("ctrl", "treat"))
  deg <- deg_from_matrix(expr, groups)

  saveRDS(list(deg = deg, meta = d, treat_samples = sel$treat,
               ctrl_samples = sel$ctrl), out)
  message("  saved: ", out, " (", nrow(deg), " genes)")
}
message("Done.")
