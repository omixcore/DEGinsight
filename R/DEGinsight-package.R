#' DEGinsight: LLM-Assisted Mechanism Interpretation for Bulk RNA-seq DEGs
#'
#' DEGinsight combines classical enrichment analysis with large language model
#' (LLM) reasoning to turn differential expression results into
#' evidence-graded, ontology-grounded mechanistic hypotheses.
#'
#' @section Main workflow:
#' \enumerate{
#'   \item \code{\link{run_ora}} — over-representation analysis of DEGs
#'   \item \code{\link{build_interpretation_prompt}} — ontology-constrained prompt
#'   \item \code{\link{call_llm}} — unified OpenAI-compatible LLM backend
#'   \item \code{\link{interpret_enrichment}} / \code{\link{interpret_deg}} —
#'         single-model and ensemble interpretation
#'   \item \code{\link{validate_terms}} — hallucination detection against
#'         controlled vocabularies
#'   \item \code{\link{render_report}} — reproducible HTML report
#' }
#'
#' @keywords internal
"_PACKAGE"
NULL
