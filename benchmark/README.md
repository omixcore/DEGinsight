# DEGinsight benchmark

Goal: can LLM-assisted interpretation **recover a known mechanism** from a
bulk RNA-seq DE table — beyond what the enrichment table already contains?
Plus hallucination auditing, cost accounting, an ORA-only control, and a
3-factor ablation (few-shot, number of pathways, tissue context).

The written-up results live in two places:

| output | what it is |
|---|---|
| [`05_report.Rmd`](05_report.Rmd) | benchmark report; renders itself from `scores/*.csv` |
| [`../paper/preprint.Rmd`](../paper/preprint.Rmd) | bioRxiv-style manuscript draft |

Both are generated from the CSVs on disk, so they never claim more coverage
than the runs actually completed.

## Datasets (frozen ground truth)

Six GEO series whose mechanism of action is known independently of the
transcriptome (`ground_truth.csv`):

| dataset | GSE | perturbagen | established MoA | DEGs (strict) |
|---|---|---|---|---|
| rapamycin_nk | GSE318529 | rapamycin | mTORC1 inhibition | 0 → top-150 fallback |
| dexamethasone_organoid | GSE304202 | dexamethasone | glucocorticoid receptor agonist | 49 |
| vorinostat_a498 | GSE319257 | vorinostat/SAHA | HDAC inhibition | 1150 |
| imatinib_k562 | GSE321679 | imatinib | BCR-ABL TKI | 121 |
| tnf_astrocyte | GSE305796 | TNF-alpha | TNFR1 → NF-kB | 792 |
| lps_macrophage | GSE297509 | LPS | TLR4 → NF-kB/IFN | 453 |

Columns of `ground_truth.csv`:

- `gt_components` — mechanism components the model must recover. Synonyms are
  separated by `/` or `|`; multiple acceptable components are separated by
  `|`. Frozen before the first model call (see git history).
- `expected_gmt` — the projection of the mechanism onto the 49 terms of
  `gene_sets.gmt`. Used so that the ORA-only control and the LLM arm are
  scored through the *same* mapping (§2 of the report).
- `treat_regex` / `ctrl_regex` — sample selection, matched against GEO titles.

## How to run

```bash
export LLM_API_KEY=sk-...            # any OpenAI-compatible endpoint
export LLM_BASE_URL=https://api.deepseek.com/v1
export RSTUDIO_PANDOC=/path/to/pandoc   # only for rendering

bash   benchmark/00_fetch_geo.sh        # resumable, size-pinned downloads
Rscript benchmark/01b_prepare_real.R    # GEO -> DEG tables (limma / Welch t)
Rscript benchmark/02_run_benchmark.R    # datasets x models -> interpretations
Rscript benchmark/04b_baseline_ora.R    # ORA-only control (needs NO API key)
RUN_JUDGE=1 Rscript benchmark/04_score.R  # aggregate + cross-model judging
Rscript benchmark/03_run_ablation.R     # few-shot x top_n x tissue grid
Rscript benchmark/04_score.R            # re-aggregate including ablation
```

Render the two documents:

```bash
cd benchmark && Rscript -e 'rmarkdown::render("05_report.Rmd")'
cd paper     && Rscript -e 'rmarkdown::render("preprint.Rmd")'
```

All scripts are idempotent (existing results are skipped) and every LLM
response is disk-cached by prompt hash, so a re-run costs nothing and a
re-render is byte-identical.

`04b_baseline_ora.R` deliberately needs no key: the control that decides
whether the model adds anything must be reproducible by a reviewer who has no
API budget.

## Metrics

- **Mechanism recovery** — fraction of pre-registered components recovered,
  scored two ways: a deterministic synonym-aware matcher (free) and an
  independent LLM judge with a strict rubric (`judge_mechanism()`), always a
  *different* model from the one being scored.
- **Hallucination rate** — fraction of cited pathways absent from both the
  input table and the controlled vocabulary (`validate_terms()`). This is an
  audit output, not a property of the model.
- **Pathway recall (ORA vs LLM)** — hit@5 / hit@10 against `expected_gmt`.
- **Cost** — tokens actually billed x per-model prices recorded in
  `models.csv`, so the accounting can be recomputed later.

## Current status

| item | state |
|---|---|
| DEG tables | 6/6 datasets complete |
| ORA-only control | 6/6 datasets complete |
| LLM runs | 3/12 dataset-model cells (API credit exhausted mid-run) |
| Ablation grid | 48/48 configurations (3 datasets, reference model) |
| Cross-model judge | not run (requires credit) |

Remaining cells are filled by re-running `02_run_benchmark.R`; it skips what
exists and resumes. Nothing in the report or the manuscript is imputed — the
coverage table at the top of `05_report.html` states the count.

## Notes / caveats

- The rapamycin dataset yields **zero** genes at padj < 0.05 and |log2FC| > 1;
  the adaptive rule falls back to the 150 smallest-p genes, balanced by
  direction, and records `deg_mode` in every result file. Two of six datasets
  use the fallback.
- Sample selection uses regexes on GEO sample titles — inspect the matched
  titles on first run, the annotations are inconsistent between submitters.
- The vorinostat table covers ~8,000 genes (panel-limited platform), so its
  DEG count is not comparable in absolute terms; all reported metrics are
  ratios.
- The judge model must differ from the model being evaluated
  (`JUDGE_MODEL` / `JUDGE_MODE=cross`) to limit self-preference bias.
- `02_run_benchmark.R` silently uses bundled toy gene sets if
  `benchmark/gene_sets.gmt` is missing — that mode is for pipeline testing
  only.

## What past runs showed

- 48-config ablation: few-shot **hurts** (0.80 → 0.64) while costing 31% more
  prompt tokens; recovery improves monotonically with wider pathway context
  (0.64 @ top_n=10 → 0.81 @ 30); tissue context is marginal; hallucination 0%
  everywhere.
- ORA alone recovers a canonical pathway in the top five for 4/6 datasets and
  in the top ten for 5/6 — which is why both documents report the LLM against
  that control rather than against nothing.
- Qualitative outputs were textbook-correct where they were correct
  (rapamycin → mTORC1-driven lipid/hypoxia reprogramming; dexamethasone → GR
  activation with NF-kB repression and FKBP5 induction; imatinib → release of
  the erythroid differentiation block in K562). Early keyword misses were
  vocabulary problems ("NF-kB" vs "NF-kappaB"), fixed in
  `normalize_omics_text()` — always report keyword and judge scores together.
