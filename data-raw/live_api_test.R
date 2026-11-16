lib <- "G:/WorkBuddy/2026-11-16-20-16-33/Rlib"
.libPaths(c(lib, .libPaths()))
for (f in list.files("R", full.names = TRUE)) source(f)
load("data/example_data.rda")

res <- interpret_enrichment(
  example_ora,
  contrast = "TNF-stimulated vs control HeLa cells",
  species = "Homo sapiens",
  tissue = "HeLa",
  direction = "up",
  model = "deepseek-chat"
)
print(res)
cat("\n--- hypotheses ---\n")
for (h in res$hypotheses) {
  cat("\n## ", h$name, " [", h$evidence_grade, "]\n", sep = "")
  cat(h$mechanism, "\n")
  cat("pathways:", paste(unlist(h$supporting_pathways), collapse = "; "), "\n")
}
cat("\nsummary:", res$overall_summary, "\n")
cat("tokens:", unlist(res$tokens), "\n")
cat("hallucination rate:", res$validation$hallucination_rate, "\n")

# judge round-trip on the same output
gt <- c("NF-kB activation", "inflammatory cytokine induction", "IL-17 signaling")
kw <- score_recovery_keywords(res, gt)
print(kw$per_component)
saveRDS(res, "benchmark/results/example__deepseek-chat.rds")
cat("LIVE API TEST OK\n")
