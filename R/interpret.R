#' Interpret enrichment results with a single LLM
#'
#' Runs the full single-model pipeline: prompt construction -> LLM call ->
#' ontology validation -> hallucination scoring.
#'
#' @param enrichment Enrichment results data.frame (see
#'   \code{\link{build_interpretation_prompt}}).
#' @param contrast,species,tissue,direction,top_n Passed to
#'   \code{\link{build_interpretation_prompt}}.
#' @param few_shot Logical; prepend a worked example to the prompt
#'   (ablation factor, see \code{\link{build_interpretation_prompt}}).
#' @param model,base_url,api_key Passed to \code{\link{call_llm}}.
#' @param vocabulary Optional character vector of valid term names used to
#'   detect hallucinated pathways (e.g. all GO BP term names). If NULL,
#'   validation is performed against the input pathway names only.
#' @param ... Further arguments passed to \code{\link{call_llm}}.
#'
#' @return An object of class `"deginsight_interpretation"`: a list with
#'   `hypotheses`, `flagged_pathways`, `overall_summary`,
#'   `validation` (see \code{\link{validate_terms}}), `model`,
#'   `tokens`, and the raw `prompt`.
#'
#' @examples
#' \dontrun{
#' res <- interpret_enrichment(example_ora,
#'                             contrast = "TNF vs control",
#'                             species = "Homo sapiens",
#'                             model = "deepseek-chat")
#' res$hypotheses
#' }
#'
#' @export
interpret_enrichment <- function(enrichment,
                                 contrast,
                                 species,
                                 tissue = NULL,
                                 direction = c("mixed", "up", "down"),
                                 top_n = 20,
                                 model = "deepseek-chat",
                                 base_url = Sys.getenv("LLM_BASE_URL", "https://api.deepseek.com/v1"),
                                 api_key = Sys.getenv("LLM_API_KEY"),
                                 vocabulary = NULL,
                                 few_shot = FALSE,
                                 ...) {
  p <- build_interpretation_prompt(enrichment, contrast, species, tissue,
                                   match.arg(direction), top_n,
                                   few_shot = few_shot)
  raw <- call_llm(p$prompt, system = p$system, model = model,
                  base_url = base_url, api_key = api_key,
                  response_schema = p$schema, ...)

  valid_terms <- unique(c(vocabulary, enrichment$pathway))
  validation <- validate_terms(raw, reference_terms = valid_terms,
                               input_pathways = enrichment$pathway)

  out <- list(
    hypotheses = raw$hypotheses %||% list(),
    flagged_pathways = raw$flagged_pathways %||% list(),
    overall_summary = raw$overall_summary %||% NA_character_,
    validation = validation,
    model = attr(raw, "model") %||% model,
    tokens = attr(raw, "usage"),
    prompt = p$prompt
  )
  class(out) <- "deginsight_interpretation"
  out
}

#' End-to-end interpretation of a DEG table
#'
#' Convenience wrapper: ORA (via \code{\link{run_ora}}) followed by
#' LLM interpretation (single model or ensemble vote).
#'
#' @param deg A data.frame with at least a gene column and (optionally)
#'   logFC / adjusted p-value columns for pre-filtering.
#' @param gene_col,padj_col,logfc_col Column names in `deg`.
#' @param padj_cutoff,logfc_cutoff DEG filtering thresholds.
#' @param universe Background genes; defaults to all rows of `deg`.
#' @param gene_sets Named list of gene sets (see \code{\link{run_ora}}).
#' @param models Character vector of one or more model IDs. If length > 1,
#'   an ensemble interpretation is produced via
#'   \code{\link{interpret_ensemble}}.
#' @param ... Passed to \code{\link{interpret_enrichment}} /
#'   \code{\link{interpret_ensemble}}.
#'
#' @return A list with elements `ora` (data.frame) and `interpretation`.
#'
#' @examples
#' \dontrun{
#' gs <- read_gmt("hallmark.gmt")
#' out <- interpret_deg(my_deg, gene_col = "symbol", padj_col = "padj",
#'                      logfc_col = "log2FoldChange", gene_sets = gs,
#'                      contrast = "treated vs control", species = "Homo sapiens",
#'                      models = c("deepseek-chat", "kimi-k2"))
#' }
#'
#' @export
interpret_deg <- function(deg,
                          gene_col,
                          padj_col = NULL,
                          logfc_col = NULL,
                          padj_cutoff = 0.05,
                          logfc_cutoff = 1,
                          universe = NULL,
                          gene_sets,
                          models = "deepseek-chat",
                          contrast,
                          species,
                          ...) {
  sig <- rep(TRUE, nrow(deg))
  if (!is.null(padj_col)) sig <- sig & deg[[padj_col]] < padj_cutoff
  if (!is.null(logfc_col)) sig <- sig & abs(deg[[logfc_col]]) > logfc_cutoff
  genes <- deg[[gene_col]][sig]
  universe <- universe %||% deg[[gene_col]]

  ora <- run_ora(genes, universe, gene_sets)

  interp <- if (length(models) == 1) {
    interpret_enrichment(ora, contrast = contrast, species = species,
                         model = models, ...)
  } else {
    interpret_ensemble(ora, contrast = contrast, species = species,
                       models = models, ...)
  }

  list(ora = ora, interpretation = interp)
}

#' @export
print.deginsight_interpretation <- function(x, ...) {
  cat("DEGinsight interpretation (model:", x$model, ")\n")
  cat("  hypotheses:", length(x$hypotheses),
      " | hallucination rate:",
      sprintf("%.1f%%", 100 * x$validation$hallucination_rate), "\n")
  for (h in x$hypotheses) {
    cat(sprintf("- [%s] %s\n", h$evidence_grade %||% "?", h$name %||% "?"))
  }
  invisible(x)
}
