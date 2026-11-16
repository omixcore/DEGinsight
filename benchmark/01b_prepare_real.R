# 01b_prepare_real.R — dataset-specific DEG preparation for the live run
#
# Why this exists: GEO series matrices for RNA-seq contain no expression
# table (verified 2026-11-16), so each series needs its own supplementary
# file handling. Run after data-raw/install_deps.R and the curl downloads
# in benchmark/README.md.
#
# Inputs (benchmark/geo_cache/):
#   GSE318529_Raw_Counts.csv.gz                  (rapamycin, raw counts)
#   GSE304202_normalized_counts.csv.gz           (dexamethasone, normalized)
#   GSE304202_gene_annotation.csv.gz             (ENSEMBL -> Symbol)
#   GSE321679_Normalised_counts_all_samples.xlsx (imatinib; via python csv)
# Outputs: benchmark/deg_tables/<dataset_id>.rds (DEG tables via limma)

lib <- "G:/WorkBuddy/2026-11-16-20-16-33/Rlib"
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
if (requireNamespace("limma", quietly = TRUE)) library(limma)

CACHE <- "benchmark/geo_cache"
dir.create("benchmark/deg_tables", showWarnings = FALSE, recursive = TRUE)

ann <- read.csv(file.path(CACHE, "GSE304202_gene_annotation.csv.gz"))
ens2sym <- setNames(ann$Symbol, ann$ENSEMBL)

limma_deg <- function(mat, groups) {
  # mat: gene x sample numeric matrix (counts or normalized values)
  voom_input <- log2(mat + 1)
  treat <- voom_input[, groups == "treat", drop = FALSE]
  ctrl  <- voom_input[, groups == "ctrl",  drop = FALSE]
  if (requireNamespace("limma", quietly = TRUE)) {
    design <- model.matrix(~ 0 + groups)
    colnames(design) <- levels(groups)
    fit <- eBayes(contrasts.fit(lmFit(voom_input, design),
                                makeContrasts(treat - ctrl, levels = design)))
    tt <- topTable(fit, number = Inf, sort.by = "none")
    data.frame(symbol = rownames(tt), log2FoldChange = tt$logFC,
               padj = p.adjust(tt$P.Value, "BH"), stringsAsFactors = FALSE)
  } else {
    message("limma unavailable — falling back to Welch t-test + BH ",
            "(document this in the paper; use limma for final runs).")
    # drop near-constant genes (t.test would return NaN)
    vars <- apply(voom_input, 1, var)
    keep <- is.finite(vars) & vars > 1e-8
    voom_input <- voom_input[keep, ]
    treat <- treat[keep, ]; ctrl <- ctrl[keep, ]
    p <- apply(voom_input, 1, function(g) {
      tryCatch(t.test(g[groups == "treat"], g[groups == "ctrl"])$p.value,
               error = function(e) NA_real_)
    })
    p[is.na(p)] <- 1
    data.frame(symbol = rownames(voom_input),
               log2FoldChange = rowMeans(treat) - rowMeans(ctrl),
               padj = p.adjust(p, "BH"), stringsAsFactors = FALSE)
  }
}

collapse_symbols <- function(counts, ids) {
  sym <- ens2sym[sub("\\..*$", "", ids)]
  keep <- !is.na(sym) & nzchar(sym)
  counts <- counts[keep, ]; sym <- sym[keep]
  counts <- rowsum(counts, sym)  # sum duplicate symbols
  counts
}

# ---- rapamycin: GSE318529, Stage 4, 15-C, mTORi vs DMSO -------------------
raw <- read.csv(file.path(CACHE, "GSE318529_Raw_Counts.csv.gz"),
                row.names = 1, check.names = FALSE)
treat <- grep("^S4_15C_mTORi", colnames(raw), value = TRUE)
ctrl  <- grep("^S4_15C_DMSO",  colnames(raw), value = TRUE)
cat("rapamycin treat:", length(treat), "ctrl:", length(ctrl), "\n")
mat <- collapse_symbols(as.matrix(raw[, c(treat, ctrl)]), rownames(raw))
mat <- mat[rowSums(mat >= 10) >= 3, ]  # filter low counts
groups <- factor(c(rep("treat", length(treat)), rep("ctrl", length(ctrl))),
                 levels = c("ctrl", "treat"))
deg <- limma_deg(mat, groups)
saveRDS(list(deg = deg, treat_samples = treat, ctrl_samples = ctrl),
        "benchmark/deg_tables/rapamycin_nk.rds")
cat("rapamycin DEGs (padj<0.05, |logFC|>1):",
    sum(deg$padj < 0.05 & abs(deg$log2FoldChange) > 1, na.rm = TRUE), "\n")

# ---- dexamethasone: GSE304202, Dexa vs Control ----------------------------
nc <- read.csv(file.path(CACHE, "GSE304202_normalized_counts.csv.gz"),
               row.names = 1, check.names = FALSE)
treat <- grep("Dexa", colnames(nc), value = TRUE)
ctrl  <- grep("Control", colnames(nc), value = TRUE)
cat("dexamethasone treat:", length(treat), "ctrl:", length(ctrl), "\n")
mat <- collapse_symbols(as.matrix(nc[, c(treat, ctrl)]), rownames(nc))
mat <- mat[rowMeans(mat) > 1, ]
groups <- factor(c(rep("treat", length(treat)), rep("ctrl", length(ctrl))),
                 levels = c("ctrl", "treat"))
deg <- limma_deg(mat, groups)
saveRDS(list(deg = deg, treat_samples = treat, ctrl_samples = ctrl),
        "benchmark/deg_tables/dexamethasone_organoid.rds")
cat("dexamethasone DEGs:", sum(deg$padj < 0.05 & abs(deg$log2FoldChange) > 1, na.rm = TRUE), "\n")

# ---- imatinib: GSE321679, K562 WT imatinib vs DMSO ------------------------
# xlsx converted to csv beforehand (python openpyxl, see README)
im_csv <- file.path(CACHE, "GSE321679_counts.csv")
if (file.exists(im_csv)) {
  # sheet 'K562-WT- DMSO VS Imatinib': ENSEMBL, SYMBOL, 3 imatinib, 3 DMSO
  im <- read.csv(im_csv, check.names = FALSE)
  counts <- as.matrix(im[, grep("WT-", colnames(im))])
  rownames(counts) <- make.unique(im$SYMBOL)
  counts <- counts[!is.na(im$SYMBOL) & nzchar(im$SYMBOL), ]
  treat <- grep("imatinib", colnames(counts), ignore.case = TRUE, value = TRUE)
  ctrl  <- grep("DMSO", colnames(counts), ignore.case = TRUE, value = TRUE)
  cat("imatinib treat:", length(treat), "ctrl:", length(ctrl), "\n")
  mat <- counts[, c(treat, ctrl)]
  storage.mode(mat) <- "numeric"
  mat <- mat[rowMeans(mat) > 1, ]
  groups <- factor(c(rep("treat", length(treat)), rep("ctrl", length(ctrl))),
                   levels = c("ctrl", "treat"))
  deg <- limma_deg(mat, groups)
  saveRDS(list(deg = deg, treat_samples = treat, ctrl_samples = ctrl),
          "benchmark/deg_tables/imatinib_k562.rds")
  cat("imatinib DEGs:", sum(deg$padj < 0.05 & abs(deg$log2FoldChange) > 1, na.rm = TRUE), "\n")
} else cat("imatinib csv not ready, skipping\n")

# ---- vorinostat: GSE319257, mock arm only, SAHA vs mock ------------------
# Combined by 01c_vorinostat_xlsx.py from the per-sample GEO xlsx files.
vo_csv <- file.path(CACHE, "GSE319257_counts.csv")
if (file.exists(vo_csv)) {
  vo <- read.csv(vo_csv, check.names = FALSE, stringsAsFactors = FALSE)
  counts <- as.matrix(vo[, setdiff(colnames(vo), "gene")])
  rownames(counts) <- make.unique(vo$gene)
  storage.mode(counts) <- "numeric"
  treat <- grep("mock_vorinostat", colnames(counts), value = TRUE)
  ctrl  <- grep("mock_control",    colnames(counts), value = TRUE)
  cat("vorinostat treat:", length(treat), "ctrl:", length(ctrl), "\n")
  mat <- counts[, c(treat, ctrl), drop = FALSE]
  keep <- rowSums(mat >= 10) >= 3
  mat <- rowsum(mat[keep, , drop = FALSE], rownames(mat)[keep])  # dedupe symbols
  groups <- factor(c(rep("treat", length(treat)), rep("ctrl", length(ctrl))),
                   levels = c("ctrl", "treat"))
  deg <- limma_deg(mat, groups)
  saveRDS(list(deg = deg, treat_samples = treat, ctrl_samples = ctrl),
          "benchmark/deg_tables/vorinostat_a498.rds")
  cat("vorinostat DEGs:", sum(deg$padj < 0.05 & abs(deg$log2FoldChange) > 1, na.rm = TRUE), "\n")
} else cat("vorinostat csv not ready, skipping\n")

# ---- LPS: GSE297509, siRNA-ctrl arm, LPS vs mock, 2 donors ----------------
# Per-sample files carry: Gene ID, Name (symbol), ..., unique reads (count).
lps_dir <- file.path(CACHE, "lps")
lps_map <- c(
  donor1_mock = "GSM8993640_A_C_minus.txt.gz",
  donor1_LPS  = "GSM8993641_A_C_LPS.txt.gz",
  donor2_mock = "GSM8993646_B_C_minus.txt.gz",
  donor2_LPS  = "GSM8993647_B_C_LPS.txt.gz")
if (all(file.exists(file.path(lps_dir, lps_map)))) {
  read_lps <- function(f) {
    x <- read.delim(gzfile(file.path(lps_dir, f)), check.names = FALSE,
                    stringsAsFactors = FALSE)
    setNames(x[["unique reads"]], x[["Name"]])
  }
  lst <- lapply(lps_map, read_lps)
  mat <- do.call(cbind, lst)
  rownames(mat) <- sub("\\..*$", "", rownames(mat))
  mat <- rowsum(mat, rownames(mat))
  treat <- grep("_LPS$", colnames(mat), value = TRUE)
  ctrl  <- grep("_mock$", colnames(mat), value = TRUE)
  cat("LPS treat:", length(treat), "ctrl:", length(ctrl), "\n")
  mat <- mat[rowSums(mat >= 10) >= 2, , drop = FALSE]
  groups <- factor(c(rep("treat", length(treat)), rep("ctrl", length(ctrl))),
                   levels = c("ctrl", "treat"))
  deg <- limma_deg(mat[, c(treat, ctrl), drop = FALSE], groups)
  saveRDS(list(deg = deg, treat_samples = treat, ctrl_samples = ctrl),
          "benchmark/deg_tables/lps_macrophage.rds")
  cat("LPS DEGs:", sum(deg$padj < 0.05 & abs(deg$log2FoldChange) > 1, na.rm = TRUE), "\n")
} else cat("LPS per-sample files not ready, skipping\n")

# ---- TNF-alpha: GSE305796, iPSC astrocytes, TNF vs untreated, 3 donors ----
# Salmon transcript counts mapped to symbols by 01d_tnf_quant.py (mygene.info).
tnf_csv <- file.path(CACHE, "GSE305796_counts.csv")
if (file.exists(tnf_csv)) {
  tf <- read.csv(tnf_csv, check.names = FALSE, stringsAsFactors = FALSE)
  counts <- as.matrix(tf[, setdiff(colnames(tf), "gene")])
  rownames(counts) <- tf$gene
  storage.mode(counts) <- "numeric"
  treat <- grep("^T", colnames(counts), value = TRUE)
  ctrl  <- grep("^C", colnames(counts), value = TRUE)
  cat("TNF treat:", length(treat), "ctrl:", length(ctrl), "\n")
  mat <- counts[, c(treat, ctrl), drop = FALSE]
  mat <- mat[rowSums(mat >= 10) >= 3, , drop = FALSE]
  groups <- factor(c(rep("treat", length(treat)), rep("ctrl", length(ctrl))),
                   levels = c("ctrl", "treat"))
  deg <- limma_deg(mat, groups)
  saveRDS(list(deg = deg, treat_samples = treat, ctrl_samples = ctrl),
          "benchmark/deg_tables/tnf_astrocyte.rds")
  cat("TNF DEGs:", sum(deg$padj < 0.05 & abs(deg$log2FoldChange) > 1, na.rm = TRUE), "\n")
} else cat("TNF counts not ready, skipping\n")

cat("01b done\n")
