# 02_run_benchmark.R — run DEGinsight over datasets x models
#
# Run from the package root:  Rscript benchmark/02_run_benchmark.R
# Requires: credentials for at least one model (see benchmark/models.csv)
#   LLM_API_KEY / LLM_BASE_URL          -> DeepSeek models
#   MOONSHOT_API_KEY                    -> kimi-k2       (optional)
#   OPENAI_API_KEY                      -> gpt-4o-mini   (optional)
# Results: benchmark/results/<dataset>__<model>.rds  (safe to re-run;
#          completed combos are skipped; LLM responses are disk-cached)

lib <- "G:/WorkBuddy/2026-11-16-20-16-33/Rlib"
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
if (requireNamespace("DEGinsight", quietly = TRUE)) {
  library(DEGinsight)
} else {
  for (f in list.files("R", full.names = TRUE)) source(f)
}

# ---- configuration --------------------------------------------------------
GENE_SETS_GMT <- Sys.getenv("BENCH_GMT", "benchmark/gene_sets.gmt")
PADJ <- 0.05; LOGFC <- 1; TOP_N <- 20

# Models to evaluate. models.csv keeps the per-vendor endpoint and the
# environment variable that holds each key, so adding a provider is a
# one-line edit rather than a code change.
mods <- read.csv("benchmark/models.csv", stringsAsFactors = FALSE,
               fileEncoding = "UTF-8-BOM")
pick <- Sys.getenv("BENCH_MODELS", "")
if (nzchar(pick)) {
  mods <- mods[mods$model_id %in% strsplit(pick, ",")[[1]], ]
}
mods$base_url <- vapply(seq_len(nrow(mods)), function(i)
  Sys.getenv(mods$base_url_env[i], mods$default_base_url[i]), "")
mods$api_key <- vapply(mods$api_key_env, function(e) Sys.getenv(e, ""), "")
mods <- mods[nzchar(mods$api_key), ]
if (!nrow(mods)) stop("No usable model credentials found (see benchmark/models.csv)")
message("models: ", paste(mods$model_id, collapse = ", "))

GT <- read.csv("benchmark/ground_truth.csv", stringsAsFactors = FALSE,
               check.names = FALSE, fileEncoding = "UTF-8-BOM")
dir.create("benchmark/results", showWarnings = FALSE, recursive = TRUE)

gene_sets <- if (file.exists(GENE_SETS_GMT)) {
  read_gmt(GENE_SETS_GMT)
} else {
  data(example_data); example_gene_sets
}

# Adaptive DEG selection: strict DE thresholds first; when the design is
# underpowered (small n), fall back to top-N genes by p-value so ORA still
# has input. The fallback is logged and must be reported in the paper.
select_degs <- function(deg, padj = 0.05, logfc = 1, top_n = 150) {
  ok <- !is.na(deg$padj) & !is.na(deg$log2FoldChange)
  sig <- ok & deg$padj < padj & abs(deg$log2FoldChange) > logfc
  if (sum(sig) >= 20) return(list(genes = deg$symbol[sig], mode = "strict"))
  deg <- deg[ok, ]
  deg <- deg[order(deg$padj), ]
  message("  strict DEGs = ", sum(sig),
          " -> falling back to top ", top_n, " by p-value")
  up <- deg$symbol[deg$log2FoldChange > 0][seq_len(min(top_n / 2, sum(deg$log2FoldChange > 0)))]
  dn <- deg$symbol[deg$log2FoldChange < 0][seq_len(min(top_n / 2, sum(deg$log2FoldChange < 0)))]
  list(genes = intersect(unique(c(up, dn)), deg$symbol), mode = "topN")
}

for (i in seq_len(nrow(GT))) {
  d <- GT[i, ]
  deg_file <- file.path("benchmark/deg_tables", paste0(d$dataset_id, ".rds"))
  if (!file.exists(deg_file)) { message("no DEG table, skip: ", d$dataset_id); next }
  obj <- readRDS(deg_file)
  deg <- obj$deg

  contrast <- paste0(d$perturbagen, "-treated vs control in ", d$cell_type)
  sel <- select_degs(deg)
  ora <- tryCatch(
    run_ora(sel$genes, deg$symbol, gene_sets),
    error = function(e) { message("  ORA failed: ", e$message); NULL }
  )
  if (is.null(ora)) next
  message("== ", d$dataset_id, " (", sel$mode, " DEGs, n=", length(sel$genes),
          ", ", nrow(ora), " pathways)")

  for (j in seq_len(nrow(mods))) {
    m <- mods$model_id[j]
    out <- file.path("benchmark/results",
                     paste0(d$dataset_id, "__", gsub("[^A-Za-z0-9]", "-", m), ".rds"))
    if (file.exists(out)) { message("  skip (exists): ", basename(out)); next }

    t0 <- Sys.time()
    interp <- tryCatch(
      interpret_enrichment(ora, contrast = contrast,
                           species = "Homo sapiens", tissue = d$cell_type,
                           model = m, base_url = mods$base_url[j],
                           api_key = mods$api_key[j]),
      error = function(e) { message("  LLM failed (", m, "): ", e$message); NULL }
    )
    if (is.null(interp)) next
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

    # deterministic keyword scoring (no extra API spend)
    kw <- score_recovery_keywords(interp, strsplit(d$gt_components, "\\|")[[1]])

    saveRDS(list(dataset = d$dataset_id, model = m,
                 model_label = mods$label[j],
                 deg_mode = sel$mode, n_degs = length(sel$genes),
                 ora = ora, interpretation = interp, keyword_score = kw,
                 elapsed_sec = elapsed), out)
    message("  ", m, " keyword recovery: ", kw$recovery_rate,
            " | hallucination: ", interp$validation$hallucination_rate,
            " | ", round(elapsed, 1), "s")
  }
}
message("Benchmark run complete.")
