# 03_run_ablation.R — ablation grid: few-shot x top_n x tissue context
#
# Run from the package root:  Rscript benchmark/03_run_ablation.R
# Design: full factorial on a subset of datasets and one reference model.
# Factors:
#   few_shot : FALSE / TRUE        (prompting strategy)
#   top_n    : 5 / 10 / 20 / 30    (number of pathways shown to the model)
#   tissue   : TRUE / FALSE        (cell-type context injected or withheld)
# Results: benchmark/ablation/<dataset>__fs<0|1>__n<k>__ctx<0|1>.rds

library(DEGinsight)

MODEL   <- Sys.getenv("ABLATION_MODEL", "deepseek-chat")
DATASETS <- strsplit(Sys.getenv("ABLATION_DATASETS",
                                "rapamycin_nk,imatinib_k562,tnf_hela"), ",")[[1]]
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

for (ds in DATASETS) {
  d <- GT[GT$dataset_id == ds, ]
  obj <- readRDS(file.path("benchmark/deg_tables", paste0(ds, ".rds")))
  deg <- obj$deg
  sig <- deg$padj < 0.05 & abs(deg$log2FoldChange) > 1
  ora <- run_ora(deg$symbol[sig], deg$symbol, gene_sets)
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
