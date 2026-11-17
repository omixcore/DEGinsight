test_that("directional ORA separates up and down lists and handles empty arms", {
  local_mocked_bindings(call_llm = function(...) list(hypotheses = list()))
  deg <- data.frame(gene = letters[1:10], adj = c(0.01, 0.01, rep(0.8, 7), NA),
                    lfc = c(2, -2, rep(0, 8)))
  sets <- list(up = letters[c(1, 3:6)], down = letters[c(2, 3:6)])
  out <- interpret_deg(deg, "gene", "adj", "lfc", gene_sets = sets,
                       contrast = "t vs c", species = "human", direction = "split")
  expect_equal(out$up$ora$genes[out$up$ora$pathway == "up"], "a")
  expect_equal(out$down$ora$genes[out$down$ora$pathway == "down"], "b")
  out <- interpret_deg(deg, "gene", "adj", "lfc", gene_sets = sets,
                       contrast = "t vs c", species = "human", logfc_cutoff = 4)
  expect_length(out$interpretation$hypotheses, 0L)
})

test_that("prompt grounding audits only pathways and genes actually shown", {
  local_mocked_bindings(call_llm = function(...) list(hypotheses = list(
    list(name = "Test", supporting_pathways = list("A", "B", "c"),
         driver_genes = list("g1", "g2")))))
  x <- data.frame(pathway = c("A", "B", "C"), p_adjust = c(0.01, 0.02, 0.03),
                  overlap = c(1, 1, 0), genes = c("g1", "g2", ""))
  out <- interpret_enrichment(x, "t vs c", "human", top_n = 1L, vocabulary = "B")
  expect_equal(out$validation$grounded, "A")
  expect_equal(out$validation$unsupported_pathways, c("B", "c"))
  expect_equal(out$validation$hallucinated, "c")
  expect_equal(out$validation$unsupported_driver_genes, "g2")
})

test_that("consensus requires independent backing for the same grounded claim", {
  local_mocked_bindings(call_llm = function(..., model) {
    list(hypotheses = list(list(name = "Claim", mechanism = "Shared response.",
                                supporting_pathways = if (model == "m1") "A" else c("A", "B"),
                                driver_genes = "g1", evidence_grade = "speculative")))
  })
  x <- data.frame(pathway = c("A", "B"), p_adjust = c(0.01, 0.02), genes = "g1")
  out <- interpret_ensemble(x, "t vs c", "human", models = c("m1", "m2"))
  expect_equal(out$consensus$pathways, "A")
  expect_equal(out$consensus$mechanisms[[1]]$pathways, "A")
  expect_equal(out$consensus$supporting_models$A, c("m1", "m2"))
})

test_that("opposing mechanism notation never forms a shared claim", {
  local_mocked_bindings(call_llm = function(..., model) {
    mechanism <- if (model == "m1") "A -> B" else "A -| B"
    list(hypotheses = list(list(name = "Regulation", mechanism = mechanism,
                               supporting_pathways = "Pathway",
                               driver_genes = "g1", evidence_grade = "speculative")))
  })
  x <- data.frame(pathway = "Pathway", p_adjust = 0.01, genes = "g1")
  out <- interpret_ensemble(x, "t vs c", "human", models = c("m1", "m2"))
  expect_equal(out$consensus$pathways, "Pathway")
  expect_length(out$consensus$mechanisms, 0L)
})
