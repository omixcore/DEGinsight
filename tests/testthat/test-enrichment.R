test_that("run_ora recovers the planted pathway", {
  data(example_data)
  sig <- example_deg$padj < 0.05 & abs(example_deg$log2FoldChange) > 1
  res <- run_ora(example_deg$symbol[sig], example_deg$symbol,
                 example_gene_sets, min_size = 3)
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
