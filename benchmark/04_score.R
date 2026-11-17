source("benchmark/common.R")
out <- file.path(BENCH, "scores")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
use_judge <- identical(Sys.getenv("RUN_JUDGE"), "1")
if (use_judge) {
  judge <- model_config(Sys.getenv("JUDGE_MODEL"))
  if (!nzchar(judge$api_key)) stop("Set credentials for JUDGE_MODEL.")
}
score <- function(r) {
  d <- GT[match(r$dataset, GT$dataset_id), ]
  components <- strsplit(d$gt_components, "|", fixed = TRUE)[[1]]
  exp <- unique(strsplit(d$expected_gmt, "|", fixed = TRUE)[[1]])
  stopifnot(all(exp %in% names(gene_sets)))
  i <- r$interpretation
  grounded <- i$validation$grounded
  n <- length(grounded)
  kw <- score_recovery_keywords(i, components)$recovery_rate
  verdict <- NA_real_
  if (use_judge && r$status == "completed") {
    if (identical(judge$model, r$model)) stop("Use a fixed judge distinct from all evaluated models.")
    verdict <- do.call(judge_mechanism, c(list(interpretation = i,
      gt_components = components, gt_context = d$gt_pathways), judge))$recovery_rate
  }
  usage <- i$tokens
  prices <- MODELS[match(r$model, MODELS$model_id), ]
  prompt_tokens <- usage$prompt_tokens %||% NA_integer_
  completion_tokens <- usage$completion_tokens %||% NA_integer_
  cost <- (prompt_tokens * prices$price_input_per_1m +
           completion_tokens * prices$price_output_per_1m) / 1e6
  data.frame(dataset = r$dataset, model = r$model, endpoint = r$config$base_url,
    status = r$status, few_shot = r$config$few_shot, top_n = r$config$top_n,
    tissue = r$config$tissue, n_degs = r$n_degs, recovery_keyword = kw,
    recovery_judge = verdict, judge_model = if (use_judge) judge$model else NA_character_,
    prompt_tokens = prompt_tokens, completion_tokens = completion_tokens,
    estimated_cost_usd = cost, elapsed_sec = r$elapsed_sec,
    n_cited = length(i$validation$cited), n_grounded = n,
    unrecognized_pathway_rate = i$validation$hallucination_rate,
    n_unsupported_drivers = length(i$validation$unsupported_driver_genes),
    expected_cited = sum(exp %in% grounded),
    expected_in_ora_at_same_budget = sum(exp %in% head(i$input_pathways, n)),
    expected_in_prompt = sum(exp %in% i$input_pathways), n_expected = length(exp))
}
for (kind in c("results", "ablation")) {
  rows <- read_runs(kind)
  summary_file <- file.path(out, paste0(kind, "_summary.csv"))
  if (length(rows)) {
    scores <- do.call(rbind, lapply(rows, score))
    write.csv(scores, summary_file, row.names = FALSE)
  } else if (file.exists(summary_file)) {
    unlink(summary_file)
  }
  message(kind, ": ", sum(vapply(rows, function(r) r$status == "completed", logical(1))),
          " model calls completed")
}
