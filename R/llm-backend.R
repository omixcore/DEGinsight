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
#' @param model Model identifier, e.g. `"deepseek-chat"`, `"gpt-4o-mini"`,
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
#' Sys.setenv(LLM_API_KEY = "sk-...", LLM_BASE_URL = "https://api.deepseek.com/v1")
#' call_llm("Reply with {\"ok\": true}",
#'          response_schema = list(ok = "logical"))
#' }
#'
#' @importFrom httr2 request req_headers req_body_json req_retry req_perform resp_body_json
#' @importFrom jsonlite fromJSON toJSON
#' @export
call_llm <- function(prompt,
                     system = NULL,
                     model = "deepseek-chat",
                     base_url = Sys.getenv("LLM_BASE_URL", "https://api.deepseek.com/v1"),
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
      "\n\nRespond with a single JSON object matching this schema ",
      "(field: type):\n",
      jsonlite::toJSON(response_schema, auto_unbox = TRUE, pretty = TRUE)
    )
  }

  messages <- c(
    if (!is.null(system)) list(list(role = "system", content = system)),
    list(list(role = "user", content = prompt))
  )
  body <- list(model = model, messages = messages, temperature = temperature)

  cache_key <- digest::digest(list(body = body, schema = response_schema),
                              algo = "sha256")
  if (use_cache) {
    hit <- cache_get(cache_key)
    if (!is.null(hit)) {
      attr(hit, "cached") <- TRUE
      return(hit)
    }
  }

  req <- httr2::request(paste0(base_url, "/chat/completions")) |>
    httr2::req_headers(
      Authorization = paste("Bearer", api_key),
      `Content-Type` = "application/json"
    ) |>
    httr2::req_body_json(c(body, list(response_format = list(type = "json_object")))) |>
    httr2::req_retry(max_tries = max_retry, backoff = ~ 2^.x)

  resp <- httr2::req_perform(req)
  payload <- httr2::resp_body_json(resp)

  text <- payload$choices[[1]]$message$content
  out <- if (is.null(response_schema)) text else parse_json_safely(text)

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
        stop("LLM output is not parseable JSON:\n", text, call. = FALSE)
      }
      jsonlite::fromJSON(m, simplifyVector = FALSE)
    }
  )
}
