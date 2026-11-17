# 04b_baseline_ora.R — no-LLM control: how much does plain ORA recover?
#
# Reviewers always ask "does the LLM add anything over the enrichment table it
# was given?". This script answers that with a free, deterministic baseline:
# treat the top-k ORA terms themselves as the answer, and score them against
# the same ground truth used for the LLM outputs.
#
# Both arms are scored on the SAME projection of the ground truth onto the 49
# terms of benchmark/gene_sets.gmt (column `expected_gmt` in
# ground_truth.csv). Using one frozen projection for both arms is what makes
# the comparison fair: any bias in the projection cancels out.
#
# Run from the package root:  Rscript benchmark/04b_baseline_ora.R
# Produces:
#   benchmark/scores/baseline_ora.csv    (per dataset: does top-k ORA hit?)
#   benchmark/scores/pathway_recall.csv  (ORA vs each LLM, per dataset)

lib <- "G:/WorkBuddy/2026-11-16-20-16-33/Rlib"
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
if (requireNamespace("DEGinsight", quietly = TRUE)) {
  library(DEGinsight)
} else {
  for (f in list.files("R", full.names = TRUE)) source(f)
}

norm_path <- function(x) {
  x <- gsub("\u03ba", "k", x)            # NF-κB -> NF-kB
  x <- tolower(x)
  x <- gsub("[^a-z0-9]", "", x)
  x
}

GT <- read.csv("benchmark/ground_truth.csv", stringsAsFactors = FALSE,
               check.names = FALSE, fileEncoding = "UTF-8-BOM")
gene_sets <- read_gmt("benchmark/gene_sets.gmt")
gmt_names <- names(gene_sets)

# sanity: every expected term must exist in the GMT, otherwise the projection
# silently scores 0 for both arms
for (i in seq_len(nrow(GT))) {
  exp <- strsplit(GT$expected_gmt[i], "\\|")[[1]]
  miss <- setdiff(norm_path(exp), norm_path(gmt_names))
  if (length(miss)) warning(GT$dataset_id[i], " unknown GMT terms: ",
                            paste(miss, collapse = ", "))
}

# same adaptive DEG rule as 02_run_benchmark.R, so the ORA input is identical
select_degs <- function(deg, padj = 0.05, logfc = 1, top_n = 150) {
  ok <- !is.na(deg$padj) & !is.na(deg$log2FoldChange)
  sig <- ok & deg$padj < padj & abs(deg$log2FoldChange) > logfc
  if (sum(sig) >= 20) return(list(genes = deg$symbol[sig], mode = "strict"))
  deg <- deg[ok, ]
  deg <- deg[order(deg$padj), ]
  up <- deg$symbol[deg$log2FoldChange > 0][seq_len(min(top_n / 2, sum(deg$log2FoldChange > 0)))]
  dn <- deg$symbol[deg$log2FoldChange < 0][seq_len(min(top_n / 2, sum(deg$log2FoldChange < 0)))]
  list(genes = intersect(unique(c(up, dn)), deg$symbol), mode = "topN")
}

recall_at <- function(hit) {
  c(r5 = any(hit[seq_len(min(5, length(hit)))]),
    r10 = any(hit[seq_len(min(10, length(hit)))]),
    n_hit = sum(hit))
}

rows <- list()
for (i in seq_len(nrow(GT))) {
  d <- GT[i, ]
  f <- file.path("benchmark/deg_tables", paste0(d$dataset_id, ".rds"))
  if (!file.exists(f)) { message("no DEG table: ", d$dataset_id); next }
  obj <- readRDS(f)
  sel <- select_degs(obj$deg)
  ora <- run_ora(sel$genes, obj$deg$symbol, gene_sets)
  if (is.null(ora) || !nrow(ora)) next
  exp <- norm_path(strsplit(d$expected_gmt, "\\|")[[1]])
  hit <- norm_path(ora$pathway) %in% exp
  rr <- recall_at(hit)
  first <- if (any(hit)) which(hit)[1] else NA_integer_
  rows[[length(rows) + 1]] <- data.frame(
    dataset = d$dataset_id, arm = "ORA_only", model = "none",
    deg_mode = sel$mode, n_degs = length(sel$genes), n_ora = nrow(ora),
    hit_at_5 = as.integer(rr[["r5"]]), hit_at_10 = as.integer(rr[["r10"]]),
    n_expected_hit = as.integer(rr[["n_hit"]]), n_expected = length(exp),
    rank_first_hit = first, stringsAsFactors = FALSE)
  message(sprintf("%-24s ORA: hit@5=%d hit@10=%d (%d/%d expected) first rank=%s",
                  d$dataset_id, rr[["r5"]], rr[["r10"]], rr[["n_hit"]],
                  length(exp), first))
}

# ---- LLM arm, from whatever has been run so far ---------------------------
files <- list.files("benchmark/results", pattern = "\\.rds$", full.names = TRUE)
for (f in files) {
  r <- readRDS(f)
  d <- GT[GT$dataset_id == r$dataset, ]
  if (!nrow(d)) next
  exp <- norm_path(strsplit(d$expected_gmt, "\\|")[[1]])
  hyps <- r$interpretation$hypotheses
  pw <- unlist(lapply(hyps, function(h) unlist(h$supporting_pathways)))
  pw <- unique(pw[nzchar(pw)])
  hit <- norm_path(pw) %in% exp
  rr <- recall_at(hit)
  rows[[length(rows) + 1]] <- data.frame(
    dataset = r$dataset, arm = "LLM", model = r$model,
    deg_mode = r$deg_mode, n_degs = r$n_degs, n_ora = nrow(r$ora),
    hit_at_5 = as.integer(rr[["r5"]]), hit_at_10 = as.integer(rr[["r10"]]),
    n_expected_hit = as.integer(rr[["n_hit"]]), n_expected = length(exp),
    rank_first_hit = NA_integer_, stringsAsFactors = FALSE)
}

out <- do.call(rbind, rows)
dir.create("benchmark/scores", showWarnings = FALSE, recursive = TRUE)
write.csv(out, "benchmark/scores/pathway_recall.csv", row.names = FALSE)

# per-dataset side-by-side, for the report table
ora <- out[out$arm == "ORA_only", ]
cat("\n== ORA-only baseline (all datasets) ==\n")
print(ora[, c("dataset", "deg_mode", "n_degs", "n_ora", "hit_at_5",
              "hit_at_10", "n_expected_hit", "n_expected")], row.names = FALSE)
cat("\nmean hit@5  =", round(mean(ora$hit_at_5), 3),
    "\nmean hit@10 =", round(mean(ora$hit_at_10), 3), "\n")

if (any(out$arm == "LLM")) {
  llm <- out[out$arm == "LLM", ]
  cat("\n== LLM pathway recovery (completed runs) ==\n")
  print(llm[, c("dataset", "model", "hit_at_5", "n_expected_hit", "n_expected")],
        row.names = FALSE)
}
