# 02_run_benchmark.R — run DEGinsight over datasets x models
#
# Run from the package root:  Rscript benchmark/02_run_benchmark.R
# Requires: Sys.setenv(LLM_API_KEY=..., LLM_BASE_URL=...)
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
MODELS <- strsplit(Sys.getenv("BENCH_MODELS", "deepseek-chat"), ",")[[1]]
GENE_SETS_GMT <- Sys.getenv("BENCH_GMT", "benchmark/gene_sets.gmt")
# Provide a GMT (e.g. MSigDB Hallmark h.all.*.symbols.gmt) at the path above,
# or leave empty to use the small bundled example_gene_sets for a dry run.

PADJ <- 0.05; LOGFC <- 1; TOP_N <- 20

GT <- read.csv("benchmark/ground_truth.csv", stringsAsFactors = FALSE,
               check.names = FALSE)
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
  sig <- deg$padj < padj & abs(deg$log2FoldChange) > logfc
  if (sum(sig) >= 20) return(list(genes = deg$symbol[sig], mode = "strict"))
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

  contrast <- paste0(d$perturbagen, "-treated vs control ", d$cell_type)
  sel <- select_degs(deg)
  ora <- tryCatch(
    run_ora(sel$genes, deg$symbol, gene_sets),
    error = function(e) { message("  ORA failed: ", e$message); NULL }
  )
  if (is.null(ora)) next

  for (m in MODELS) {
    out <- file.path("benchmark/results",
                     paste0(d$dataset_id, "__", gsub("[^A-Za-z0-9]", "-", m), ".rds"))
    if (file.exists(out)) { message("skip (exists): ", basename(out)); next }
    message("== ", d$dataset_id, " x ", m)

    interp <- tryCatch(
      interpret_enrichment(ora, contrast = contrast,
                           species = "Homo sapiens", tissue = d$cell_type,
                           model = m),
      error = function(e) { message("  LLM failed: ", e$message); NULL }
    )
    if (is.null(interp)) next

    # deterministic keyword scoring (no extra API spend)
    kw <- score_recovery_keywords(interp, strsplit(d$gt_components, "\\|")[[1]])

    saveRDS(list(dataset = d$dataset_id, model = m, ora = ora,
                 interpretation = interp, keyword_score = kw), out)
    message("  keyword recovery: ", kw$recovery_rate,
            " | hallucination: ", interp$validation$hallucination_rate)
  }
}
message("Benchmark run complete.")
