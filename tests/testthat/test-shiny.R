test_that("Shiny prepares evidence locally and audits an explicitly requested result", {
  skip_if_not_installed("shiny")
  requests <- 0L
  local_mocked_bindings(call_llm = function(...) {
    requests <<- requests + 1L
    list(hypotheses = list(list(
      name = "Example", mechanism = "A test mechanism",
      supporting_pathways = c("Pilot pathway", "Unsupported pathway"),
      driver_genes = c("A", "B"), evidence_grade = "hypothesis",
      follow_up_experiment = "Repeat the experiment"
    )), overall_summary = "Example result.")
  })
  withr::local_envvar(LLM_API_KEY = NA_character_)
  path <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(data.frame(pathway = "Pilot pathway", p_adjust = 0.01, genes = "A"),
                   path, row.names = FALSE)

  shiny::testServer(shiny_server, {
    session$setInputs(mode = "enrichment", table = list(name = "input.csv", datapath = path),
                      model = "test-model", contrast = "case versus control",
                      species = "Homo sapiens", tissue = "", top_n = 10L, prepare = 1L)
    expect_match(prepared()$prompt$prompt, "Pilot pathway")
    expect_identical(requests, 0L)
    expect_null(result())
    session$setInputs(interpret = 1L)
    expect_identical(requests, 1L)
    expect_identical(result()$validation$grounded, "Pilot pathway")
    expect_identical(result()$validation$unsupported_driver_genes, "B")
    session$setInputs(table = list(name = "input.txt", datapath = path), prepare = 2L)
    expect_null(prepared())
    expect_null(result())
  })
})
