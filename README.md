
# curveRfreq

`curveRfreq` is an R package designed for fast, transparent standard
curve fitting for immunoassay data including data from bead array and
ELISA assays. It leverages nonlinear least squares (Levenberg-Marquart)
with automatic model selection via AIC. It provides direct control over
asymotote constraints, and estimates sample concentrations by inversion
of the fitted equation. In addition quality control measures including
limits of detection and limits of quantification can be calculated as
well as sample measurement error. This package provides a pracrical
default for routine single-plate analysis where computational speed and
reliability matter most.

## Installation

The package is not yet on CRAN but can be installed directly from the
repository.

``` r
devtools::install_github("immunoplex/curveRfreq")
```
