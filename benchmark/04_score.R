# 04_score.R — aggregate benchmark / ablation results into tables
#
# Run from the package root:  Rscript benchmark/04_score.R
# Produces:
#   benchmark/scores/benchmark_summary.csv   (per dataset x model)
#   benchmark/scores/ablation_summary.csv
#   benchmark/scores/model_ranking.csv       (per-model means)
# Optional LLM-judge scoring: set RUN_JUDGE=1 (uses API; keyword scores
# are always computed for free).

lib <- "G:/WorkBuddy/2026-11-16-20-16-33/Rlib"
if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
if (requireNamespace("DEGinsight", quietly = TRUE)) {
  library(DEGinsight)
} else {
  for (f in list.files("R", full.names = TRUE)) source(f)
}

GT <- read.csv("benchmark/ground_truth.csv", stringsAsFactors = FALSE,
               check.names = FALSE)
dir.create("benchmark/scores", showWarnings = FALSE, recursive = TRUE)
RUN_JUDGE <- identical(Sys.getenv("RUN_JUDGE"), "1")
JUDGE_MODEL <- Sys.getenv("JUDGE_MODEL", "deepseek-chat")

cost_of <- function(interp, model) {
  usage <- interp$tokens
  if (is.null(usage)) return(NA_real_)
  price_table <- if (exists("default_prices", inherits = TRUE)) {
    default_prices()
  } else {
    DEGinsight:::default_prices()
  }
  prices <- price_table[[model]]
  if (is.null(prices)) return(NA_real_)
  (usage$prompt_tokens %||% 0) / 1e6 * prices$input +
    (usage$completion_tokens %||% 0) / 1e6 * prices$output
}

score_file <- function(f, run_judge = FALSE) {
  r <- readRDS(f)
  d <- GT[GT$dataset_id == r$dataset, ]
  gt_comp <- strsplit(d$gt_components, "\\|")[[1]]
  # recompute keyword score against the CURRENT ground truth so that
  # vocabulary refinements apply without re-calling the LLM (responses
  # are disk-cached anyway)
  kw <- score_recovery_keywords(r$interpretation, gt_comp)
  gt_context <- paste0(
    "Perturbation: ", d$perturbagen, " in ", d$cell_type, " (", d$gse,
    "). Established mechanism components: ",
    paste(gt_comp, collapse = "; "),
    ". Canonical pathways: ", gsub("\\|", ", ", d$gt_pathways), ".")

  judge <- NULL
  if (run_judge) {
    judge <- tryCatch(
      judge_mechanism(r$interpretation, gt_comp, gt_context,
                      model = JUDGE_MODEL),
      error = function(e) { message("judge failed: ", e$message); NULL })
  }

  data.frame(
    dataset = r$dataset,
    model = r$model %||% NA_character_,
    recovery_keyword = kw$recovery_rate,
    recovery_judge = if (!is.null(judge)) judge$recovery_rate else NA_real_,
    hallucination_rate = r$interpretation$validation$hallucination_rate,
    n_hypotheses = length(r$interpretation$hypotheses),
    prompt_tokens = r$interpretation$tokens$prompt_tokens %||% NA_integer_,
    completion_tokens = r$interpretation$tokens$completion_tokens %||% NA_integer_,
    cost_usd = cost_of(r$interpretation, r$model %||% NA_character_),
    stringsAsFactors = FALSE
  )
}

# ---- main benchmark -------------------------------------------------------
files <- list.files("benchmark/results", pattern = "\\.rds$", full.names = TRUE)
if (length(files)) {
  bench <- do.call(rbind, lapply(files, score_file, run_judge = RUN_JUDGE))
  write.csv(bench, "benchmark/scores/benchmark_summary.csv", row.names = FALSE)

  ranking <- aggregate(
    cbind(recovery_keyword, recovery_judge, hallucination_rate, cost_usd) ~ model,
    bench, mean, na.rm = TRUE, na.action = na.pass)
  ranking <- ranking[order(-ranking$recovery_keyword), ]
  write.csv(ranking, "benchmark/scores/model_ranking.csv", row.names = FALSE)
  print(ranking)
}

# ---- ablation -------------------------------------------------------------
afiles <- list.files("benchmark/ablation", pattern = "\\.rds$", full.names = TRUE)
if (length(afiles)) {
  abl <- do.call(rbind, lapply(afiles, function(f) {
    r <- readRDS(f)
    data.frame(dataset = r$dataset,
               few_shot = r$config$few_shot, top_n = r$config$top_n,
               tissue = r$config$tissue,
               recovery_keyword = r$keyword_score$recovery_rate,
               hallucination_rate = r$interpretation$validation$hallucination_rate,
               prompt_tokens = r$interpretation$tokens$prompt_tokens %||% NA_integer_,
               stringsAsFactors = FALSE)
  }))
  write.csv(abl, "benchmark/scores/ablation_summary.csv", row.names = FALSE)
  cat("\nAblation means by factor level:\n")
  print(aggregate(recovery_keyword ~ few_shot, abl, mean))
  print(aggregate(recovery_keyword ~ top_n, abl, mean))
  print(aggregate(recovery_keyword ~ tissue, abl, mean))
  print(aggregate(hallucination_rate ~ tissue, abl, mean))
}
