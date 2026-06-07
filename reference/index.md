# Package index

## Main Entry Points

The two functions most users call directly. Both expect preprocessed
data on the fitting scale (use curveRcore::preprocess_standards()
upstream).

- [`bead_assay_example`](https://immunoplex.github.io/curveRfreq/reference/bead_assay_example.md)
  : Synthetic Multi-Plate Bead Assay Example Dataset
- [`fit_calibration_freq()`](https://immunoplex.github.io/curveRfreq/reference/fit_calibration_freq.md)
  : Fit a Frequentist Calibration Curve (Single Curve)
- [`fit_calibration_freq_multiplate()`](https://immunoplex.github.io/curveRfreq/reference/fit_calibration_freq_multiplate.md)
  : Fit Frequentist Calibration Curves Across Multiple Curves

## Multi-Curve Helpers

Convenience extractors that operate on calibration_result_multiplate
objects.

- [`collect_samples()`](https://immunoplex.github.io/curveRfreq/reference/collect_samples.md)
  : Extract All Sample Predictions from Multi-Curve Results
- [`summary_table()`](https://immunoplex.github.io/curveRfreq/reference/summary_table.md)
  : Extract a Combined Summary Table from Multi-Curve Results

## NLS Fitting Engine

Lower-level functions that implement multi-start Levenberg-Marquardt
fitting for one plate. Called internally by fit_calibration_freq() but
exported for users who need fine-grained control.

- [`compute_model_constraints()`](https://immunoplex.github.io/curveRfreq/reference/compute_model_constraints.md)
  : Compute Parameter Bounds for All Candidate Models
- [`fit_ensemble_nls()`](https://immunoplex.github.io/curveRfreq/reference/fit_ensemble_nls.md)
  : Fit Ensemble of NLS Models for One Plate
- [`generate_start_lists()`](https://immunoplex.github.io/curveRfreq/reference/generate_start_lists.md)
  : Generate Multi-Start Lists for NLS Fitting

## Model Selection

AIC-based ranking and best-parameter extraction. select_best_aic() is a
pure ranking step; eligibility gating is handled upstream by
curveRcore::select_best_eligible().

- [`summarise_ensemble()`](https://immunoplex.github.io/curveRfreq/reference/summarise_ensemble.md)
  : Summarise Ensemble Fit Statistics
- [`extract_best_parameters()`](https://immunoplex.github.io/curveRfreq/reference/extract_best_parameters.md)
  : Extract Parameters from the Best Fit
- [`select_best_aic()`](https://immunoplex.github.io/curveRfreq/reference/select_best_aic.md)
  : Select the Best Model by AIC

## Grid Prediction & Uncertainty

Evaluates the best-fit model across the prediction grid, computes
delta-method confidence intervals on the response, back-calculates
concentration, and propagates uncertainty to se_concentration and pcov.

- [`predict_grid_freq()`](https://immunoplex.github.io/curveRfreq/reference/predict_grid_freq.md)
  : Predict Grid with Uncertainty (Frequentist)

## Sample Back-Calculation

Applies the analytical model inverse to test-sample responses,
propagates uncertainty via the delta method, and multiplies by the
dilution factor to produce final_concentration.

- [`predict_samples_freq()`](https://immunoplex.github.io/curveRfreq/reference/predict_samples_freq.md)
  : Predict Concentration for Test Samples
