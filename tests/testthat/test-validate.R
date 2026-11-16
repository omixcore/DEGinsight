fake_interpretation <- function() {
  list(hypotheses = list(
    list(name = "NF-kB activation",
         supporting_pathways = list("TNF signaling via NF-kB", "Made-up pathway XYZ"),
         evidence_grade = "strong"),
    list(name = "Immune recruitment",
         supporting_pathways = list("IL-17 signaling"),
         evidence_grade = "moderate")
  ))
}

test_that("hallucinated terms are detected and scored", {
  data(example_data)
  v <- validate_terms(fake_interpretation(),
                      reference_terms = character(0),
                      input_pathways = example_ora$pathway)
  expect_equal(v$hallucinated, "Made-up pathway XYZ")
  expect_equal(v$hallucination_rate, 1 / 3)
  expect_equal(nrow(v$per_hypothesis), 2)
})

test_that("clean output yields zero hallucination rate", {
  data(example_data)
  good <- list(hypotheses = list(
    list(name = "x", supporting_pathways = list(example_ora$pathway[1]))
  ))
  v <- validate_terms(good, character(0), example_ora$pathway)
  expect_equal(v$hallucination_rate, 0)
})
