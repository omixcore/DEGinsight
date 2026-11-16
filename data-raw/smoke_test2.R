# Smoke test 2: benchmark scoring logic (no API calls)
for (f in list.files("R", full.names = TRUE)) { parse(f); source(f) }
load("data/example_data.rda")

# --- few-shot prompt variant ----------------------------------------------
p0 <- build_interpretation_prompt(example_ora, "x", "y", few_shot = FALSE)
p1 <- build_interpretation_prompt(example_ora, "x", "y", few_shot = TRUE)
stopifnot(!grepl("## Example", p0$prompt, fixed = TRUE))
stopifnot(grepl("## Example", p1$prompt, fixed = TRUE))
stopifnot(nchar(p1$prompt) > nchar(p0$prompt))
cat("few-shot toggle OK (", nchar(p0$prompt), "->", nchar(p1$prompt), "chars)\n")

# --- mock interpretation resembling an LLM output --------------------------
mock <- list(
  hypotheses = list(
    list(name = "mTORC1 inhibition blocks translation",
         mechanism = "mTOR signaling suppression with ribosomal program shutdown indicates mTORC1-driven translational repression.",
         supporting_pathways = list("TNF signaling via NF-kB", "Ribosome biogenesis"),
         driver_genes = list("RPTOR", "RPS6"),
         evidence_grade = "strong",
         follow_up_experiment = "phospho-S6K western"),
    list(name = "Immune note", mechanism = "minor", 
         supporting_pathways = list("IL-17 signaling"), driver_genes = list(),
         evidence_grade = "speculative", follow_up_experiment = "qPCR")
  ),
  flagged_pathways = list(),
  overall_summary = "Dominated by mTORC1 translational shutdown."
)
class(mock) <- "deginsight_interpretation"
attr(mock, "model") <- "mock-model"
attr(mock, "usage") <- list(prompt_tokens = 900, completion_tokens = 600)

# keyword scoring: synonym groups
gt <- c("mTOR|mTORC1 inhibition", "translational/ribosomal suppression",
        "autophagy induction")
kw <- score_recovery_keywords(mock, gt)
print(kw$per_component)
stopifnot(kw$recovery_rate == 2/3)
cat("keyword recovery OK:", kw$recovery_rate, "\n")

# hallucination validation on mock
v <- validate_terms(mock, reference_terms = character(0),
                    input_pathways = example_ora$pathway)
cat("hallucination rate:", v$hallucination_rate, "\n")  # Ribosome biogenesis not in input -> flagged
stopifnot(v$hallucination_rate > 0)

# cost accounting mirroring 04_score.R
usage <- attr(mock, "usage")
prices <- default_prices()[["deepseek-chat"]]
cost <- usage$prompt_tokens / 1e6 * prices$input + usage$completion_tokens / 1e6 * prices$output
cat("mock cost USD:", cost, "\n")
stopifnot(cost > 0 & cost < 0.01)

# ablation grid sanity (as in 03_run_ablation.R)
grid <- expand.grid(few_shot = c(FALSE, TRUE), top_n = c(5, 10, 20, 30),
                    tissue = c(TRUE, FALSE))
stopifnot(nrow(grid) == 16)
cat("ablation grid OK: 16 configs\n")

cat("SMOKE TEST 2 OK\n")
