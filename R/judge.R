#' Judge mechanism recovery against ground truth
#'
#' Two complementary scorers used in the DEGinsight benchmark:
#'
#' \itemize{
#'   \item \code{judge_mechanism()}: LLM-as-judge. An independent judge model
#'     reads the interpretation plus the ground-truth mechanism description
#'     and scores each ground-truth component as recovered (1) or not (0),
#'     with a short justification. Uses a different model than the one being
#'     evaluated wherever possible to reduce self-preference bias.
#'   \item \code{score_recovery_keywords()}: deterministic keyword fallback
#'     (case-insensitive matching of ground-truth keywords against
#'     hypothesis text). Reported alongside the judge score as a robustness
#'     check.
#' }
#'
#' @param interpretation A `"deginsight_interpretation"` object (or a raw
#'   list with a `hypotheses` element).
#' @param gt_components Character vector of ground-truth mechanism
#'   components to check, e.g. `c("mTORC1 inhibition", "translational
#'   repression")`.
#' @param gt_context Character scalar describing the perturbation and its
#'   established mechanism (shown to the judge).
#' @param model,base_url,api_key Judge model settings (see
#'   \code{\link{call_llm}}).
#' @param ... Passed to \code{\link{call_llm}}.
#'
#' @return `judge_mechanism()` returns a list with `per_component`
#'   (data.frame: component, recovered, justification), `recovery_rate`
#'   (fraction of components recovered), `judge_model`, and `tokens`.
#'   `score_recovery_keywords()` returns the same structure with
#'   `judge_model = "keyword"`.
#'
#' @export
judge_mechanism <- function(interpretation,
                            gt_components,
                            gt_context,
                            model = "deepseek-chat",
                            base_url = Sys.getenv("LLM_BASE_URL", "https://api.deepseek.com/v1"),
                            api_key = Sys.getenv("LLM_API_KEY"),
                            ...) {
  hyp_text <- hypotheses_to_text(interpretation)

  prompt <- paste0(
    "## Ground truth\n", gt_context, "\n\n",
    "## Mechanism components to verify\n",
    paste(sprintf("%d. %s", seq_along(gt_components), gt_components),
          collapse = "\n"), "\n\n",
    "## Interpretation under evaluation\n", hyp_text, "\n\n",
    "## Task\n",
    "For each ground-truth component, decide whether the interpretation ",
    "recovered it (recovered = true) or not (recovered = false). ",
    "A component counts as recovered if the interpretation states the same ",
    "mechanism, even with different wording. Partial or vague mentions ",
    "count as NOT recovered. Be strict."
  )
  schema <- list(
    verdicts = "array of objects: {component: string, recovered: boolean, justification: string}"
  )

  raw <- call_llm(prompt,
                  system = "You are a strict, unbiased evaluator of mechanistic hypotheses.",
                  model = model, base_url = base_url, api_key = api_key,
                  response_schema = schema, ...)

  verdicts <- raw$verdicts %||% list()
  per_component <- do.call(rbind, lapply(verdicts, function(v) {
    data.frame(component = v$component %||% NA_character_,
               recovered = isTRUE(v$recovered),
               justification = v$justification %||% NA_character_,
               stringsAsFactors = FALSE)
  }))

  list(per_component = per_component,
       recovery_rate = if (nrow(per_component)) mean(per_component$recovered) else NA_real_,
       judge_model = model,
       tokens = attr(raw, "usage"))
}

#' @rdname judge_mechanism
#' @export
score_recovery_keywords <- function(interpretation, gt_components) {
  text <- tolower(hypotheses_to_text(interpretation))
  per_component <- do.call(rbind, lapply(gt_components, function(comp) {
    # synonyms separated by "|" or "/" (e.g. "mTOR|mTORC1",
    # "translational/ribosomal suppression")
    kws <- trimws(strsplit(comp, "[|/]")[[1]])
    hit <- any(vapply(tolower(kws), function(k) grepl(k, text, fixed = TRUE),
                      logical(1)))
    data.frame(component = comp, recovered = hit,
               justification = if (hit) "keyword match" else "no keyword match",
               stringsAsFactors = FALSE)
  }))
  list(per_component = per_component,
       recovery_rate = mean(per_component$recovered),
       judge_model = "keyword",
       tokens = NULL)
}

#' Flatten hypotheses into plain text for scoring
#' @noRd
hypotheses_to_text <- function(interpretation) {
  hyps <- interpretation$hypotheses %||% list()
  parts <- vapply(hyps, function(h) {
    paste(c(h$name, h$mechanism,
            unlist(h$supporting_pathways), unlist(h$driver_genes),
            h$follow_up_experiment), collapse = " | ")
  }, character(1))
  paste(c(parts, interpretation$overall_summary), collapse = "\n")
}
