source("benchmark/common.R")
set.seed(101)
out <- file.path(BENCH, "scores")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
results <- lapply(GT$dataset_id, function(id) {
  message("Preranked GSEA: ", id)
  deg <- get_dataset(id)$deg
  result <- run_gsea(setNames(deg$statistic, deg$symbol), gene_sets)
  cbind(dataset = rep(id, nrow(result)), result)
})
write.csv(do.call(rbind, results), file.path(out, "baseline_gsea.csv"), row.names = FALSE)
