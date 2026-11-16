#' Build an ontology-constrained interpretation prompt
#'
#' Converts enrichment results into a structured prompt for mechanism
#' interpretation. The prompt design is a core methodological component of
#' DEGinsight:
#'
#' \itemize{
#'   \item \strong{Grounding}: the model is given the exact enriched term
#'     names, statistics, and overlapping genes — it is explicitly forbidden
#'     from inventing pathways.
#'   \item \strong{Context injection}: experimental contrast, species and
#'     tissue are provided to anchor biological reasoning.
#'   \item \strong{Evidence grading}: the model must grade each hypothesis
#'     (strong / moderate / speculative) based only on the supplied data.
#'   \item \strong{Structured output}: a JSON schema constrains the answer
#'     for downstream validation (\code{\link{validate_terms}}).
#' }
#'
#' @param enrichment A data.frame of enrichment results. Required columns:
#'   `pathway`, `p_adjust`; optional: `overlap`, `size`, `genes`.
#'   Output of \code{\link{run_ora}} or a `clusterProfiler` result coerced
#'   with \code{as.data.frame()}.
#' @param contrast Character scalar describing the comparison, e.g.
#'   `"drug-treated vs control HepG2 cells"`.
#' @param species Character scalar, e.g. `"Homo sapiens"`.
#' @param tissue Optional tissue / cell-line context.
#' @param direction Character: are these up-regulated, down-regulated or
#'   mixed DEGs? One of `"up"`, `"down"`, `"mixed"`.
#' @param top_n Number of top pathways to include.
#' @param max_genes Maximum overlapping genes listed per pathway.
#'
#' @return A list with elements `system`, `prompt` and `schema`
#'   (the JSON response schema).
#'
#' @examples
#' data(example_ora, package = "DEGinsight")
#' p <- build_interpretation_prompt(example_ora,
#'                                  contrast = "TNF-stimulated vs control",
#'                                  species = "Homo sapiens")
#' cat(substr(p$prompt, 1, 200))
#'
#' @export
build_interpretation_prompt <- function(enrichment,
                                        contrast,
                                        species,
                                        tissue = NULL,
                                        direction = c("mixed", "up", "down"),
                                        top_n = 20,
                                        max_genes = 30) {
  direction <- match.arg(direction)
  df <- utils::head(enrichment[order(enrichment$p_adjust), ], top_n)

  pathway_lines <- vapply(seq_len(nrow(df)), function(i) {
    genes <- if ("genes" %in% names(df) && !is.na(df$genes[i])) {
      g <- strsplit(df$genes[i], ",")[[1]]
      paste0(utils::head(g, max_genes), collapse = ", ")
    } else ""
    sprintf("%d. %s (adj.P = %.2e, overlap = %s/%s)%s",
            i, df$pathway[i], df$p_adjust[i],
            if ("overlap" %in% names(df)) df$overlap[i] else "?",
            if ("size" %in% names(df)) df$size[i] else "?",
            if (nzchar(genes)) paste0("\n   Genes: ", genes) else "")
  }, character(1))

  system <- paste(
    "You are a careful molecular biologist specializing in transcriptomics.",
    "You interpret pathway enrichment results strictly based on the data",
    "provided. You NEVER invent pathway names, genes, or citations.",
    "When evidence is weak you say so explicitly."
  )

  prompt <- paste0(
    "## Experiment\n",
    "- Contrast: ", contrast, "\n",
    "- Species: ", species, "\n",
    if (!is.null(tissue)) paste0("- Tissue/cells: ", tissue, "\n") else "",
    "- Direction of DEGs analysed: ", direction, "-regulated\n\n",
    "## Enriched pathways (top ", nrow(df), " by adjusted p-value)\n",
    paste(pathway_lines, collapse = "\n"), "\n\n",
    "## Task\n",
    "1. Synthesize these pathways into 1-5 coherent mechanistic hypotheses ",
    "explaining the biological response in this contrast.\n",
    "2. For each hypothesis: name it, describe the mechanism, list the ",
    "supporting pathway names EXACTLY as written above, list key driver ",
    "genes (only from the genes shown), and grade the evidence as ",
    "'strong', 'moderate' or 'speculative'.\n",
    "3. Flag any pathway that looks like a false positive or artifact.\n",
    "4. Suggest one concrete follow-up experiment per hypothesis.\n",
    "Use only pathway names and genes from the lists above."
  )

  schema <- list(
    hypotheses = paste(
      "array of objects: {name: string, mechanism: string,",
      "supporting_pathways: array of string, driver_genes: array of string,",
      "evidence_grade: 'strong'|'moderate'|'speculative',",
      "follow_up_experiment: string}"
    ),
    flagged_pathways = "array of objects: {pathway: string, reason: string}",
    overall_summary = "string"
  )

  list(system = system, prompt = prompt, schema = schema)
}
