# DEGinsight benchmark

Goal: can LLM-assisted interpretation **recover known drug mechanisms** from
bulk RNA-seq enrichment results? Plus hallucination-rate and cost
accounting, and a 3-factor ablation (few-shot, number of pathways, tissue
context).

## Datasets (ground truth)

Seven GEO series with unambiguous mechanisms of action, verified via NCBI
E-utilities on 2026-11-16 (`ground_truth.csv`):

| dataset | GSE | perturbagen | established MoA |
|---|---|---|---|
| rapamycin_nk | GSE318529 | rapamycin | mTORC1 inhibition |
| dexamethasone_organoid | GSE304202 | dexamethasone | glucocorticoid receptor agonist |
| vorinostat_a498 | GSE319257 | vorinostat/SAHA | HDAC inhibition |
| imatinib_k562 | GSE321679 | imatinib | BCR-ABL TKI |
| tnf_hela | GSE244770 | TNF-α | TNFR → NF-κB |
| lps_monocyte | GSE289434 | LPS | TLR4 → NF-κB/IFN |
| wortmannin_hepg2 | GSE299063 | wortmannin | PI3K inhibition |

`gt_components` = mechanism components the model must recover
(synonyms separated by `|` for keyword scoring).

## How to run

```bash
# 0. prerequisites: GEOquery, limma; an OpenAI-compatible API key
export LLM_API_KEY=sk-...
export LLM_BASE_URL=https://api.deepseek.com/v1
# optional: place an MSigDB Hallmark GMT at benchmark/gene_sets.gmt

Rscript benchmark/01_download.R        # GEO -> DEG tables (limma)
Rscript benchmark/02_run_benchmark.R   # datasets x models -> interpretations
RUN_JUDGE=1 Rscript benchmark/04_score.R   # + LLM-judge scoring; writes CSVs
Rscript benchmark/03_run_ablation.R    # few-shot x top_n x tissue grid
Rscript benchmark/04_score.R           # aggregate ablation too
# render benchmark/05_report.Rmd for figures
```

All scripts are idempotent (skip existing results) and all LLM responses
are disk-cached by DEGinsight, so re-runs cost nothing.

## Metrics

- **Mechanism recovery rate** — fraction of ground-truth components
  recovered, scored two ways: deterministic keyword match (free) and an
  independent LLM judge with a strict rubric (`judge_mechanism()`).
- **Hallucination rate** — fraction of cited pathway names absent from the
  input table and reference vocabulary (`validate_terms()`).
- **Cost** — actual token usage from API responses x per-model prices.

## Notes / caveats

- Sample selection uses regexes on sample titles (`treat_regex` /
  `ctrl_regex` in `ground_truth.csv`); **inspect matched titles on first
  run** — GEO annotations are inconsistent.
- For publication-quality DE on RNA-seq, replace the log-transform + limma
  shortcut in `01_download.R` with DESeq2 on raw counts from the series
  supplementary files.
- Judge model should differ from the models being evaluated
  (`JUDGE_MODEL` env var) to limit self-preference bias.
- `02_run_benchmark.R` uses the bundled toy gene sets if no GMT is found —
  that mode is only for pipeline testing, not for real results.
