test_that("majority_vote keeps only widely supported elements", {
  sets <- list(m1 = c("A", "B"), m2 = c("A", "C"), m3 = c("A", "B"))
  expect_setequal(majority_vote(sets, min_agreement = 0.5), c("A", "B"))
  expect_equal(majority_vote(sets, min_agreement = 0.7), "A")
  expect_setequal(majority_vote(sets, min_agreement = 1 / 3), c("A", "B", "C"))
})

test_that("estimate_cost prices known models", {
  data(example_data)
  cost <- estimate_cost(example_ora, contrast = "x", species = "y",
                        models = c("deepseek-chat", "gpt-4o-mini"))
  expect_equal(nrow(cost), 2)
  expect_true(all(cost$cost_usd > 0))
  expect_true(cost$cost_usd[1] < 0.01)  # a single run should cost fractions of a cent
})
