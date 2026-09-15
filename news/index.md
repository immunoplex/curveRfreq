# Changelog

## curveRfreq 0.4.2

- Created a new pcov_gate_class and changed the basis for pcov_pass
  classifcations.

## curveRfreq 0.4.1

- Fix: predict_samples_freq() now gates pcov_pass on pcov_threshold (not
  cv_x_max), matching predict_grid_freq(). Previously frequentist test
  samples with pcov between pcov_threshold and cv_x_max were incorrectly
  marked pcov_pass = TRUE.

## curveRfreq 0.4.0 (2026-07-29)

- Lockstep version bump — **no functional changes**. Released so the
  curveR stack shares a version and the worker image can pin
  `curveRfreq@v0.4.0` for reproducible builds. Rebuilt/verified against
  curveRcore 0.4.0.

## curveRfreq 0.1.0

- Initial release.
- [`fit_calibration_freq()`](https://immunoplex.github.io/curveRfreq/reference/fit_calibration_freq.md)
  — single-curve NLS calibration with multi-start Levenberg-Marquardt
  and AIC + eligibility-gate selection.
- [`fit_calibration_freq_multiplate()`](https://immunoplex.github.io/curveRfreq/reference/fit_calibration_freq_multiplate.md)
  — multi-curve wrapper that splits by `curve_id` and handles per-curve
  errors gracefully.
- Per-model precision grids: every converged model gets a full `pcov`
  profile, not just the selected best.
- Four eligibility gates: `at_bound`, `vcov_condition`, `rel_se`,
  `dynamic_range` — intercept unidentified models before AIC ranking.
- [`summary_table()`](https://immunoplex.github.io/curveRfreq/reference/summary_table.md)
  and
  [`collect_samples()`](https://immunoplex.github.io/curveRfreq/reference/collect_samples.md)
  for tidy extraction from multi-curve results.
- `bead_assay_example` synthetic dataset for two antigens × three
  plates.
