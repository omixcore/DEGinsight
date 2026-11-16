test_that("prompt contains grounding constraints and pathway names", {
  data(example_data)
  p <- build_interpretation_prompt(example_ora, contrast = "TNF vs control",
                                   species = "Homo sapiens", top_n = 3)
  expect_match(p$prompt, "TNF signaling via NF-kB", fixed = TRUE)
  expect_match(p$prompt, "NEVER invent")
  expect_match(p$system, "NEVER invent")
  expect_true(all(c("system", "prompt", "schema") %in% names(p)))
})

test_that("top_n limits included pathways", {
  data(example_data)
  p <- build_interpretation_prompt(example_ora, contrast = "x", species = "y",
                                   top_n = 1)
  expect_match(p$prompt, "TNF signaling via NF-kB", fixed = TRUE)
  expect_false(grepl("IL-17 signaling", p$prompt, fixed = TRUE))
})
