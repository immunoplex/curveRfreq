# curveRfreq 0.4.3 (2026-09-26)

* **`summary_table()$converged` now means the selected model actually fitted.**
  It was `!is.na(selection$best_model_name)`, and the eligibility fallback in
  `curveRcore::select_best_eligible()` always returns a name -- so a plate on
  which no `nls()` succeeded was still reported as converged. This was masked
  until now by a crash (below) that made the whole per-plate call fail and
  record the curve as failed for the wrong reason.
* **New additive `fallback` column** in `summary_table()`, recording whether the
  selection came from the eligibility fallback. "Fitted cleanly" and "fitted but
  passed no gate" were previously conflated; on real plates the latter is common.
* **Fixed three NULL-fit dereferences.** `extract_best_parameters()` guarded the
  model *name* but not the fit and hit `summary(NULL)$coefficients`
  (`$ operator is invalid for atomic vectors`); the grid fallback and the sample
  path in `fit_calibration()` had the same shape and reached `vcov(NULL)`. All
  three now guard the fit, matching the check step 6a already performed.
* Validated across 1,706 plates / 119 antigen units in four studies: no
  convergence count changed, and the 3-level case that per-plate NLS genuinely
  cannot fit is now reported as clean non-convergence rather than a caught crash.


# curveRfreq 0.4.2

* Created a new pcov_gate_class and changed the basis for pcov_pass classifcations.


# curveRfreq 0.4.1

* Fix: predict_samples_freq() now gates pcov_pass on pcov_threshold (not cv_x_max), matching predict_grid_freq(). Previously frequentist test samples with pcov between pcov_threshold and cv_x_max were incorrectly marked pcov_pass = TRUE.


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
