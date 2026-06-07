# curveRfreq

Frequentist nonlinear calibration curves for the curveR suite.

## Installation

```r
# Install curveRcore first
devtools::install_github("immunoplex/curveRcore")
# Then curveRfreq
devtools::install_github("immunoplex/curveRfreq")
```

## Quick start

```r
library(curveRfreq)

data(bead_assay_example, package = "curveRfreq")

# Preprocess standards (upstream of fitting)
std_pre <- curveRcore::preprocess_standards(...)

# Fit all curves
mp <- fit_calibration_freq_multiplate(
  standards      = std_pre,
  samples        = bead_assay_example$samples,
  response_var   = "mfi",
  std_curve_conc = 30
)

# One-row-per-curve summary
summary_table(mp)

# All sample back-calculations
collect_samples(mp)
```

---

## Quick Reference

| Task | Command |
|---|---|
| Regenerate man pages | `devtools::document()` |
| Full site rebuild | `pkgdown::build_site(dest_dir = "docs")` |
| Reference pages only | `pkgdown::build_reference()` |
| Vignette only | `pkgdown::build_articles()` |
| Home page only | `pkgdown::build_home()` |
| Check no topics are unassigned | `setdiff(pkg$topics$name, unlist(lapply(cfg$reference, "[[", "contents")))` |
| Preview locally | Open `docs/index.html` |
| Push and deploy | `git add docs/ && git commit -m "..." && git push` |

---

## Troubleshooting

| Problem | Fix |
|---|---|
| `curveRcore` not found during CI | Add the `install_github("immunoplex/curveRcore")` step before `setup-r-dependencies` in the workflow |
| Vignette fails during `R CMD check` | Add `eval = requireNamespace("curveRcore", quietly = TRUE)` to the setup chunk |
| `bead_assay_example` not found | Run `usethis::use_data(bead_assay_example)` and add a `R/data.R` roxygen block |
| Function missing from Reference | Check `@export` is present; re-run `devtools::document()` |
| LaTeX equations not rendering | Add `mathjax: true` under `template:` in `_pkgdown.yml` |
| `docs/` not served by GitHub | Settings → Pages: branch `main`, folder `/docs` |

---

See `vignette("frequentist-quickstart", package = "curveRfreq")` for
the full worked example.
