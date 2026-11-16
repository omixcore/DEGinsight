# DEGinsight 0.1.0

* Initial release.
* Core pipeline: `run_ora()` -> `build_interpretation_prompt()` ->
  `interpret_enrichment()` / `interpret_ensemble()` -> `validate_terms()`.
* Unified OpenAI-compatible LLM backend with structured JSON output,
  retry, and on-disk response caching.
* Hallucination auditing against controlled vocabularies
  (`validate_terms()`, `go_vocabulary()`).
* Multi-model majority-vote ensembles with agreement matrices.
* `estimate_cost()` budgeting utility.
* Reproducible HTML reports (`render_report()`).
* Bundled simulated TNF-stimulation example dataset.
