# 03_run_ablation.R — ablation grid: few-shot x top_n x tissue context
#
# Run from the package root:  Rscript benchmark/03_run_ablation.R
# Design: full factorial on a subset of datasets and one reference model.
# Factors:
#   few_shot : FALSE / TRUE        (prompting strategy)
#   top_n    : 5 / 10 / 20 / 30    (number of pathways shown to the model)
#   tissue   : TRUE / FALSE        (cell-type context injected or withheld)
# Results: benchmark/ablation/<dataset>__fs<0|1>__n<k>__ctx<0|1>.rds

lib <- "G:/WorkBuddy/2026-11-16-20-16-33/Rlib"
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
if (requireNamespace("DEGinsight", quietly = TRUE)) {
  library(DEGinsight)
} else {
  for (f in list.files("R", full.names = TRUE)) source(f)
}

MODEL   <- Sys.getenv("ABLATION_MODEL", "deepseek-chat")
DATASETS <- strsplit(Sys.getenv("ABLATION_DATASETS",
                                "rapamycin_nk,imatinib_k562,dexamethasone_organoid"), ",")[[1]]
GRID <- expand.grid(few_shot = c(FALSE, TRUE),
                    top_n = c(5, 10, 20, 30),
                    tissue = c(TRUE, FALSE),
                    stringsAsFactors = FALSE)

GT <- read.csv("benchmark/ground_truth.csv", stringsAsFactors = FALSE,
               check.names = FALSE)
GENE_SETS_GMT <- Sys.getenv("BENCH_GMT", "benchmark/gene_sets.gmt")
gene_sets <- if (file.exists(GENE_SETS_GMT)) read_gmt(GENE_SETS_GMT) else {
  data(example_data); example_gene_sets
}
dir.create("benchmark/ablation", showWarnings = FALSE, recursive = TRUE)

select_degs <- function(deg, padj = 0.05, logfc = 1, top_n = 150) {
  sig <- deg$padj < padj & abs(deg$log2FoldChange) > logfc
  if (sum(sig) >= 20) return(deg$symbol[sig])
  deg <- deg[order(deg$padj), ]
  up <- deg$symbol[deg$log2FoldChange > 0][seq_len(min(top_n / 2, sum(deg$log2FoldChange > 0)))]
  dn <- deg$symbol[deg$log2FoldChange < 0][seq_len(min(top_n / 2, sum(deg$log2FoldChange < 0)))]
  unique(c(up, dn))
}

for (ds in DATASETS) {
  d <- GT[GT$dataset_id == ds, ]
  obj <- readRDS(file.path("benchmark/deg_tables", paste0(ds, ".rds")))
  deg <- obj$deg
  ora <- run_ora(select_degs(deg), deg$symbol, gene_sets)
  contrast <- paste0(d$perturbagen, "-treated vs control ", d$cell_type)
  gt_comp <- strsplit(d$gt_components, "\\|")[[1]]

  for (g in seq_len(nrow(GRID))) {
    cfg <- GRID[g, ]
    out <- file.path("benchmark/ablation",
      sprintf("%s__fs%d__n%d__ctx%d.rds", ds, cfg$few_shot, cfg$top_n, cfg$tissue))
    if (file.exists(out)) next
    message("== ", basename(out))

    interp <- tryCatch(
      interpret_enrichment(
        ora, contrast = contrast, species = "Homo sapiens",
        tissue = if (cfg$tissue) d$cell_type else NULL,
        top_n = cfg$top_n, few_shot = cfg$few_shot, model = MODEL),
      error = function(e) { message("  failed: ", e$message); NULL })
    if (is.null(interp)) next

    saveRDS(list(dataset = ds, config = cfg, interpretation = interp,
                 keyword_score = score_recovery_keywords(interp, gt_comp)),
            out)
  }
}
message("Ablation run complete.")
