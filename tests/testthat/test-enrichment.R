test_that("run_ora recovers the planted pathway", {
  data(example_data)
  sig <- example_deg$padj < 0.05 & abs(example_deg$log2FoldChange) > 1
  res <- run_ora(example_deg$symbol[sig], example_deg$symbol,
                 example_gene_sets, min_size = 3, p_cutoff = 0.05)
  expect_s3_class(res, "deginsight_ora")
  expect_equal(res$pathway[1], "TNF signaling via NF-kB")
  expect_true(all(res$p_adjust <= 0.05))
})

test_that("run_ora errors cleanly on empty input", {
  data(example_data)
  expect_error(run_ora("NOT_A_GENE", example_deg$symbol, example_gene_sets))
})

test_that("read_gmt round-trips", {
  f <- tempfile(fileext = ".gmt")
  writeLines(c("setA\tdesc\tG1\tG2\tG3", "setB\tdesc\tG2\tG4"), f)
  gs <- read_gmt(f)
  expect_equal(names(gs), c("setA", "setB"))
  expect_equal(gs$setA, c("G1", "G2", "G3"))
})

test_that("ORA adjusts over zero-hit pathways and deduplicates the universe", {
  universe <- letters[1:10]
  sets <- list(hit = c("a", "b", "b"), miss = c("c", "d", "e"))
  res <- run_ora(c("a", "b", "b"), c(universe, "a", NA, ""), sets, min_size = 2)
  expected <- stats::phyper(1, 2, 8, 2, lower.tail = FALSE)
  expect_equal(res$p_value, c(expected, 1))
  expect_equal(res$p_adjust, c(2 * expected, 1))
  expect_equal(run_ora(character(), universe, sets, min_size = 2)$p_value, c(1, 1))
  expect_equal(nrow(run_ora("a", universe, sets, min_size = 5)), 0L)
})

test_that("GSEA preserves both enrichment directions", {
  skip_if_not_installed("fgsea")
  set.seed(420)
  ranks <- setNames(seq(3, -2, length.out = 150), paste0("g", 1:150))
  out <- run_gsea(ranks, list(positive = names(ranks)[1:20],
                             negative = names(ranks)[131:150]))
  expect_gt(out$nes[out$pathway == "positive"], 0)
  expect_lt(out$nes[out$pathway == "negative"], 0)
  expect_match(build_interpretation_prompt(out, "t vs c", "human")$prompt, "NES")
})
