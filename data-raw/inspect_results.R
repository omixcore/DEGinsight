lib <- "G:/WorkBuddy/2026-11-16-20-16-33/Rlib"
.libPaths(c(lib, .libPaths()))
GT <- read.csv("benchmark/ground_truth.csv", stringsAsFactors = FALSE, check.names = FALSE)
for (f in list.files("benchmark/results", pattern = "rds$", full.names = TRUE)) {
  r <- readRDS(f)
  d <- GT[GT$dataset_id == r$dataset, ]
  cat("\n=====", r$dataset, "x", r$model, "=====\n")
  cat("-- ORA top terms:\n")
  print(head(r$ora[, c("pathway", "overlap", "p_adjust")], 8))
  cat("-- hypotheses:\n")
  for (h in r$interpretation$hypotheses) {
    cat(sprintf("[%s] %s\n   %s\n", h$evidence_grade, h$name,
                substr(h$mechanism, 1, 220)))
  }
  cat("-- keyword components:\n")
  print(r$keyword_score$per_component[, c("component", "recovered")])
  cat("-- summary:", substr(r$interpretation$overall_summary, 1, 260), "\n")
}
