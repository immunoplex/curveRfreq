# curveRfreq 0.4.0 (2026-07-29)

* Lockstep version bump — **no functional changes**. Released so the curveR stack
  shares a version and the worker image can pin `curveRfreq@v0.4.0` for
  reproducible builds. Rebuilt/verified against curveRcore 0.4.0.

# curveRfreq 0.1.0

* Initial release.
* `fit_calibration_freq()` — single-curve NLS calibration with
  multi-start Levenberg-Marquardt and AIC + eligibility-gate selection.
* `fit_calibration_freq_multiplate()` — multi-curve wrapper that splits
  by `curve_id` and handles per-curve errors gracefully.
* Per-model precision grids: every converged model gets a full `pcov`
  profile, not just the selected best.
* Four eligibility gates: `at_bound`, `vcov_condition`, `rel_se`,
  `dynamic_range` — intercept unidentified models before AIC ranking.
* `summary_table()` and `collect_samples()` for tidy extraction from
  multi-curve results.
* `bead_assay_example` synthetic dataset for two antigens × three plates.

# curveRfreq <next>

## Verified compatible with curveRcore 0.3.0 (mask-aware preprocessing)

* No code changes required. `fit_calibration_freq*()` receive standards and
  blanks already preprocessed by `curveRcore::preprocess_standards()`, and the
  worker passes only the *included* subset (masked rows are filtered out before
  the fit). Verified there is no internal call to `preprocess_standards()`,
  `correct_prozone()`, `perform_blank_operation()`, or `compute_log_response()`,
  no recomputation of set-level statistics, and no database reads in the fit
  path. Blanks are stored in `result$blanks` for QA/plotting only and do not
  enter the frequentist fit, so masked points cannot influence it.
