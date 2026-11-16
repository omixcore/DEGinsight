# Smoke test: parse + source all R files, exercise non-LLM functions.
for (f in list.files("R", full.names = TRUE)) { parse(f); cat("parsed:", f, "\n") }
for (f in list.files("R", full.names = TRUE)) source(f)

load("data/example_data.rda")
sig <- example_deg$padj < 0.05 & abs(example_deg$log2FoldChange) > 1
ora <- run_ora(example_deg$symbol[sig], example_deg$symbol, example_gene_sets, min_size = 3)
print(ora[, c("pathway", "p_adjust")])

p <- build_interpretation_prompt(ora, contrast = "TNF vs control",
                                 species = "Homo sapiens", top_n = 3)
cat(nchar(p$prompt), "chars in prompt\n")

v <- validate_terms(
  list(hypotheses = list(list(name = "x",
                              supporting_pathways = list("Apoptosis", "Fake pathway")))),
  character(0), ora$pathway)
cat("hallucination rate:", v$hallucination_rate, "\n")

print(majority_vote(list(m1 = c("A", "B"), m2 = c("A", "C"), m3 = c("A", "B"))))
print(estimate_cost(ora, contrast = "x", species = "y",
                    models = c("deepseek-chat", "gpt-4o-mini")))
cat("SMOKE TEST OK\n")
