# Run from the repository root after 00_fetch_geo.sh and 01_convert_inputs.py.
stopifnot(requireNamespace("limma", quietly = TRUE),
          requireNamespace("edgeR", quietly = TRUE))
cache <- "benchmark/geo_cache"
output <- "benchmark/deg_tables"
dir.create(output, showWarnings = FALSE, recursive = TRUE)

read_input <- function(name, ...) {
  read.csv(file.path(cache, name), check.names = FALSE, ...)
}

fit_deg <- function(mat, samples, symbols, type, paired = FALSE) {
  if (!identical(colnames(mat), samples$sample) || anyDuplicated(samples$sample) ||
      is.null(rownames(mat)) || anyDuplicated(rownames(mat)) || !is.numeric(mat) ||
      any(!is.finite(mat)) || (type != "logcpm" && any(mat < 0))) {
    stop("Invalid expression matrix or sample alignment.")
  }
  samples$condition <- factor(samples$condition, levels = c("control", "treated"))
  if (anyNA(samples$condition) || min(table(samples$condition)) < 2L) {
    stop("Expected at least two samples per condition.")
  }
  if ("block" %in% names(samples)) {
    samples$block <- factor(samples$block)
    if (paired && !all(table(samples$block, samples$condition) == 1L)) {
      stop("Expected exactly one sample per block and condition.")
    }
    design <- model.matrix(~ block + condition, samples)
  } else {
    design <- model.matrix(~ condition, samples)
  }
  if (qr(design)$rank != ncol(design) || nrow(design) <= ncol(design)) {
    stop("The design has confounding or insufficient residual degrees of freedom.")
  }
  if (type %in% c("counts", "scaled_tpm_counts")) {
    dge <- edgeR::DGEList(mat)
    keep <- edgeR::filterByExpr(dge, group = samples$condition)
    normalize <- if (utils::packageVersion("edgeR") >= "4.4.0" &&
                     "normLibSizes" %in% getNamespaceExports("edgeR")) edgeR::normLibSizes else edgeR::calcNormFactors
    dge <- normalize(dge[keep, , keep.lib.sizes = FALSE], method = "TMM")
    expression <- limma::voom(dge, design, plot = FALSE)
    method <- "filterByExpr + TMM + limma-voom"
  } else if (type == "normalized_counts") {
    keep <- rowSums(mat >= 1) >= min(table(samples$condition))
    expression <- log2(mat[keep, , drop = FALSE] + 1)
    method <- "log2(normalized counts + 1) + limma-trend"
  } else {
    expression <- mat
    method <- "submitted logCPM + limma-trend (original voom weights unavailable)"
  }
  fit <- limma::eBayes(limma::lmFit(expression, design),
                       trend = !type %in% c("counts", "scaled_tpm_counts"))
  tt <- limma::topTable(fit, coef = "conditiontreated", number = Inf, sort.by = "none")
  deg <- data.frame(gene_id = rownames(tt), symbol = unname(symbols[rownames(tt)]),
                    log2FoldChange = tt$logFC, p_value = tt$P.Value,
                    padj = tt$adj.P.Val, statistic = tt$t, mean_expression = tt$AveExpr)
  list(deg = deg, samples = samples, input_type = type, method = method,
       design = design, tested_genes = nrow(deg), residual_df = nrow(design) - ncol(design),
       pipeline = "limma-v2")
}

save_dataset <- function(id, expression, sample_data, symbols, type, paired = FALSE) {
  obj <- fit_deg(expression, sample_data, symbols, type, paired)
  obj$preparation_sha256 <- digest::digest(file = "benchmark/01b_prepare_real.R", algo = "sha256")
  obj$sources_sha256 <- digest::digest(file = "benchmark/sources.tsv", algo = "sha256")
  saveRDS(obj, file.path(output, paste0(id, ".rds")))
  write.csv(obj$samples, file.path(output, paste0(id, "_samples.csv")), row.names = FALSE)
  message(id, ": ", obj$tested_genes, " tested genes, residual df = ", obj$residual_df)
}

select_samples <- function(samples, treated, control, n_per_arm) {
  treat <- grepl(treated, samples)
  ctrl <- grepl(control, samples)
  stopifnot(!any(treat & ctrl), sum(treat) == n_per_arm, sum(ctrl) == n_per_arm)
  keep <- treat | ctrl
  data.frame(sample = samples[keep], condition = ifelse(treat[keep], "treated", "control"))
}

converted <- function(name) {
  tab <- read_input(name)
  mat <- as.matrix(tab[setdiff(names(tab), c("gene_id", "symbol"))])
  rownames(mat) <- tab$gene_id
  list(mat = mat, symbols = setNames(tab$symbol, tab$gene_id))
}

raw <- as.matrix(read_input("GSE318529_Raw_Counts.csv.gz", row.names = 1))
ann <- read_input("ensembl97_symbols.csv")
samples <- select_samples(colnames(raw), "^S4_15C_mTORi_", "^S4_15C_DMSO_", 3)
samples$block <- sub("^.*_", "", samples$sample)
save_dataset("rapamycin_nk", raw[, samples$sample], samples,
             setNames(ann$symbol, ann$gene_id), "counts", paired = TRUE)

norm <- as.matrix(read_input("GSE304202_normalized_counts.csv.gz", row.names = 1))
ann <- read_input("GSE304202_gene_annotation.csv.gz")
meta <- read_input("GSE304202_sample_metadata.csv.gz", row.names = 1)
samples <- select_samples(colnames(norm), "_Dexa_", "_Control_", 6)
stopifnot(!anyDuplicated(meta$Sample), setequal(meta$Sample, samples$sample))
index <- match(samples$sample, meta$Sample)
stopifnot(all(meta$Condition[index] == ifelse(samples$condition == "treated", "Dexa", "Control")))
samples$submitted_id <- rownames(meta)[index]
samples$sex <- meta$Sex[index]
# The study used each patient's sample as its own control (doi:10.1016/j.redox.2026.104057).
# The deposited sample IDs share a suffix between matched vehicle and treated samples.
stopifnot(all(grepl("^[CD][0-9]+$", samples$submitted_id)),
          all(substr(samples$submitted_id, 1, 1) == ifelse(samples$condition == "treated", "D", "C")))
samples$block <- sub("^[CD]", "", samples$submitted_id)
stopifnot(all(vapply(split(samples$sex, samples$block), function(x) length(unique(x)) == 1L, logical(1))))
save_dataset("dexamethasone_organoid", norm[, samples$sample], samples,
             setNames(ann$Symbol, ann$ENSEMBL), "normalized_counts", paired = TRUE)

im <- converted("GSE321679_logcpm.csv")
samples <- select_samples(colnames(im$mat), "^WT-imatinib-", "^WT-DMSO-", 3)
samples$block <- sub("^.*-", "", samples$sample)
save_dataset("imatinib_k562", im$mat[, samples$sample], samples,
             im$symbols, "logcpm", paired = TRUE)

vo <- converted("GSE319257_counts.csv")
samples <- select_samples(colnames(vo$mat), "_mock_vorinostat_", "_mock_control_", 6)
save_dataset("vorinostat_a498", vo$mat[, samples$sample], samples, vo$symbols, "counts")

lps <- converted("GSE297509_counts.csv")
samples <- select_samples(colnames(lps$mat), "_treat$", "_ctrl$", 2)
samples$block <- sub("^GSM[0-9]+_(donor[12])_.*$", "\\1", samples$sample)
save_dataset("lps_macrophage", lps$mat[, samples$sample], samples,
             lps$symbols, "counts", paired = TRUE)

tnf <- converted("GSE305796_scaled_tpm_counts.csv")
samples <- select_samples(colnames(tnf$mat), "_T[0-9]+$", "_C[0-9]+$", 3)
samples$block <- sub("^.*_[TC]", "", samples$sample)
save_dataset("tnf_astrocyte", tnf$mat[, samples$sample], samples, tnf$symbols,
             "scaled_tpm_counts", paired = TRUE)

writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
