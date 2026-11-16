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
#   (e.g. from Enrichr: https://maayanlab.cloud/Enrichr/geneSetLibrary?mode=text&libraryName=MSigDB_Hallmark_2020)

Rscript benchmark/01_download.R        # GEO -> DEG tables (limma)
Rscript benchmark/02_run_benchmark.R   # datasets x models -> interpretations
RUN_JUDGE=1 Rscript benchmark/04_score.R   # + LLM-judge scoring; writes CSVs
Rscript benchmark/03_run_ablation.R    # few-shot x top_n x tissue grid
Rscript benchmark/04_score.R           # aggregate ablation too
# render benchmark/05_report.Rmd for figures
```

Network note: if NCBI TLS fails ("SSL connect error"), download series
matrices / supplementary files with `curl -k` into `benchmark/geo_cache/`
manually, then use `01b_prepare_real.R` (dataset-specific supplementary-file
handling — RNA-seq series matrices contain no expression table).

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

## Live run (2026-11-16)

First end-to-end execution on 3 datasets (rapamycin_nk, dexamethasone_organoid,
imatinib_k562) x deepseek-chat, Hallmark gene sets, DEG tables via
`01b_prepare_real.R` (limma; Welch t-test fallback where limma was
unavailable during the run — flagged in the log, rerun with limma for the
paper). Results in `scores/`:

| metric | value |
|---|---|
| mechanism recovery (keyword) | 0.722 |
| mechanism recovery (LLM judge) | 0.778 |
| hallucination rate | 0.00 |
| mean cost per dataset | $0.002 |

Ablation (48 configs, 3 datasets):

| factor | finding |
|---|---|
| few-shot | hurts recovery on these data (0.80 without vs 0.64 with) — the worked example appears to bias framing; report as a negative result |
| top_n pathways | recovery improves with more pathways (0.64 @10 -> 0.81 @30) |
| tissue context | marginal (0.71 vs 0.73) |
| hallucination | 0% in all 48 configs |

Model quality note: the qualitative interpretations were textbook-correct
(rapamycin -> lipid metabolic reprogramming via mTORC1; dexamethasone -> GR
activation with NF-kB repression and FKBP5 induction; imatinib -> release of
erythroid differentiation block in K562). The remaining keyword misses were
vocabulary mismatches (e.g. "NF-kB" vs "NF-kappaB"), fixed by
`normalize_omics_text()` — always report keyword + judge scores together.

Pending: tnf_hela (GSE244770 stores only bigWigs — needs recount),
lps_monocyte (suppl lacks the paired untreated-vs-LPS contrast tables),
vorinostat_a498 / wortmannin_hepg2 (RAW.tar extraction pending).
