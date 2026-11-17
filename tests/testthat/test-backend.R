test_that("LLM cache is separated by endpoint and obeys the requested directory", {
  dir <- tempfile()
  withr::local_options(list(DEGinsight.cache_dir = NULL))
  DEGinsight_cache(dir)
  count <- 0L
  response <- list(choices = list(list(message = list(content = "ok"))))
  testthat::local_mocked_bindings(req_perform = function(req) {
    count <<- count + 1L
    req
  }, resp_body_json = function(resp) response, .package = "httr2")
  first <- call_llm("prompt", base_url = "https://one.example/v1", api_key = "mock")
  repeat_call <- call_llm("prompt", base_url = "https://one.example/v1/", api_key = "mock")
  second <- call_llm("prompt", base_url = "https://two.example/v1", api_key = "mock")
  expect_false(attr(first, "cached"))
  expect_true(attr(repeat_call, "cached"))
  expect_false(attr(second, "cached"))
  expect_equal(count, 2L)
  expect_length(list.files(dir), 2L)
})

test_that("interpretation structure requires appropriate types and enumerations", {
  expect_error(validate_interpretation_response(list(
    overall_summary = "ok", hypotheses = list(list(name = "missing fields")),
    flagged_pathways = list())), "Invalid")
  expect_silent(validate_interpretation_response(list(
    overall_summary = "ok", hypotheses = list(), flagged_pathways = list())))
})
