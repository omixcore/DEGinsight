#' Call an LLM through a unified OpenAI-compatible backend
#'
#' Sends a prompt to any OpenAI-compatible chat-completions endpoint
#' (OpenAI, DeepSeek, Moonshot/Kimi, local Ollama, etc.) and returns the
#' parsed JSON object requested via `response_schema`.
#'
#' Reliability design (see the accompanying paper):
#' \itemize{
#'   \item temperature defaults to 0 for reproducibility;
#'   \item JSON-mode / structured output is requested when supported;
#'   \item transient failures are retried with exponential backoff;
#'   \item responses are cached on disk keyed by a hash of the full request,
#'         so repeated runs are free and byte-identical.
#' }
#'
#' @param prompt Character scalar. The user prompt.
#' @param system Character scalar or NULL. Optional system prompt.
#' @param model Model identifier, e.g. `"deepseek-flash"`, `"gpt-4o-mini"`,
#'   `"qwen2.5:14b"` (Ollama).
#' @param base_url API endpoint. Defaults to the environment variable
#'   `LLM_BASE_URL`, falling back to DeepSeek's public endpoint.
#' @param api_key API key. Defaults to `Sys.getenv("LLM_API_KEY")`. For local
#'   backends (e.g. Ollama) any non-empty string works.
#' @param response_schema Named list describing the expected JSON object
#'   (used both for prompt augmentation and response validation). If NULL the
#'   raw text is returned.
#' @param temperature Sampling temperature. Default 0.
#' @param max_retry Number of retries on transient HTTP errors.
#' @param use_cache Logical; read/write the on-disk response cache.
#'
#' @return A named list parsed from the model's JSON output (or a character
#'   scalar if `response_schema` is NULL), with attributes `model`,
#'   `usage` (token counts) and `cached`.
#'
#' @examples
#' \dontrun{
#' Sys.setenv(LLM_API_KEY = "sk-...", LLM_BASE_URL = "https://api.deepseek.com")
#' call_llm("Reply with {\"ok\": true}",
#'          response_schema = list(ok = "logical"))
#' }
#'
#' @importFrom httr2 request req_headers req_body_json req_retry req_perform resp_body_json
#' @importFrom jsonlite fromJSON toJSON
#' @export
call_llm <- function(prompt,
                     system = NULL,
                     model = "deepseek-flash",
                     base_url = Sys.getenv("LLM_BASE_URL", "https://api.deepseek.com"),
                     api_key = Sys.getenv("LLM_API_KEY"),
                     response_schema = NULL,
                     temperature = 0,
                     max_retry = 3,
                     use_cache = TRUE) {
  if (!nzchar(api_key)) {
    stop("No API key found. Set Sys.setenv(LLM_API_KEY = '...') first.",
         call. = FALSE)
  }

  if (!is.null(response_schema)) {
    prompt <- paste0(
      prompt,
      "\n\nRespond with a single JSON object (json only, no prose) ",
      "matching this schema (field: type):\n",
      jsonlite::toJSON(response_schema, auto_unbox = TRUE, pretty = TRUE)
    )
  }

  messages <- c(
    if (!is.null(system)) list(list(role = "system", content = system)),
    list(list(role = "user", content = prompt))
  )
  body <- list(model = model, messages = messages, temperature = temperature)

  if (!is.null(response_schema)) body$response_format <- list(type = "json_object")
  endpoint <- paste0(sub("/+$", "", base_url), "/chat/completions")
  cache_key <- digest::digest(list(endpoint = endpoint, body = body,
                                   schema = response_schema), algo = "sha256")
  if (use_cache) {
    hit <- cache_get(cache_key)
    if (!is.null(hit)) {
      attr(hit, "cached") <- TRUE
      return(hit)
    }
  }

  req <- httr2::request(endpoint) |>
    httr2::req_headers(
      Authorization = paste("Bearer", api_key),
      `Content-Type` = "application/json"
    ) |>
    httr2::req_body_json(body) |>
    httr2::req_timeout(180) |>
    httr2::req_retry(max_tries = max_retry, retry_on_failure = TRUE,
                     is_transient = function(resp) httr2::resp_status(resp) %in%
                       c(408L, 429L, 500L, 502L, 503L, 504L))

  resp <- httr2::req_perform(req)
  payload <- httr2::resp_body_json(resp)

  choice <- if (length(payload$choices)) payload$choices[[1]] else NULL
  text <- choice$message$content
  if (!is.character(text) || length(text) != 1L || !nzchar(trimws(text)) ||
      identical(choice$finish_reason, "length")) stop("Empty or truncated LLM response.")
  out <- if (is.null(response_schema)) text else parse_json_safely(text)
  if (!is.null(response_schema) &&
      (is.null(names(out)) || !all(names(response_schema) %in% names(out)))) {
    stop("LLM response is missing required JSON fields.")
  }
  if (identical(names(response_schema), c("hypotheses", "flagged_pathways", "overall_summary"))) {
    validate_interpretation_response(out)
  }

  attr(out, "model") <- model
  attr(out, "usage") <- payload$usage
  attr(out, "cached") <- FALSE

  if (use_cache) cache_set(cache_key, out)
  out
}

#' Parse JSON from model output, tolerating code fences
#' @noRd
parse_json_safely <- function(text) {
  text <- gsub("^```(json)?\\s*|\\s*```$", "", trimws(text))
  tryCatch(
    jsonlite::fromJSON(text, simplifyVector = FALSE),
    error = function(e) {
      # last resort: extract the outermost {...} block
      m <- regmatches(text, regexpr("\\{.*\\}", text, perl = TRUE))
      if (length(m) == 0 || !nzchar(m)) {
        stop("LLM output is not parseable JSON.", call. = FALSE)
      }
      jsonlite::fromJSON(m, simplifyVector = FALSE)
    }
  )
}

#' @noRd
validate_interpretation_response <- function(out) {
  scalar_string <- function(x) is.character(x) && length(x) == 1L && !is.na(x)
  array <- function(x) is.list(x) && is.null(names(x))
  strings <- function(x) array(x) && all(vapply(x, scalar_string, logical(1)))
  hypothesis <- function(h) {
    is.list(h) && all(vapply(h[c("name", "mechanism", "follow_up_experiment")],
                             scalar_string, logical(1))) &&
      strings(h$supporting_pathways) && strings(h$driver_genes) &&
      scalar_string(h$evidence_grade) &&
      h$evidence_grade %in% c("strong", "moderate", "speculative")
  }
  flag <- function(x) is.list(x) && scalar_string(x$pathway) && scalar_string(x$reason)
  if (!scalar_string(out$overall_summary) || !array(out$hypotheses) ||
      length(out$hypotheses) > 5L || !all(vapply(out$hypotheses, hypothesis, logical(1))) ||
      !array(out$flagged_pathways) || !all(vapply(out$flagged_pathways, flag, logical(1)))) {
    stop("Invalid LLM interpretation response structure.", call. = FALSE)
  }
  invisible(out)
}
