source("benchmark/common.R")
out <- file.path(BENCH, "scores")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

expected_terms <- function(d) {
  terms <- unique(strsplit(d$expected_gmt, "|", fixed = TRUE)[[1]])
  stopifnot(length(terms) > 0L, all(terms %in% names(gene_sets)))
  terms
}

rows <- list()
for (i in seq_len(nrow(GT))) {
  d <- GT[i, ]
  obj <- get_dataset(d$dataset_id)
  for (direction in c("mixed", "up", "down")) {
    b <- benchmark_ora(obj$deg, direction)
    exp <- expected_terms(d)
    hits <- b$ora$pathway %in% exp
    rows[[length(rows) + 1L]] <- data.frame(dataset = d$dataset_id, direction = direction,
      method = obj$method, residual_df = obj$residual_df, n_tested = obj$tested_genes,
      n_symbols = nrow(obj$deg), n_degs = b$n_degs, pathways_with_overlap = nrow(b$ora),
      pathways_fdr_0_05 = sum(b$ora$p_adjust < 0.05),
      expected_at_5 = sum(head(hits, 5L)), expected_at_10 = sum(head(hits, 10L)),
      expected_at_20 = sum(head(hits, 20L)), n_expected = length(exp))
  }
}
write.csv(do.call(rbind, rows), file.path(out, "baseline_ora.csv"), row.names = FALSE)
message("Wrote ORA baselines for ", nrow(GT), " datasets.")
