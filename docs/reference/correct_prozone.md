# Correct the Prozone (Hook) Effect in Standard Curve Data

The prozone effect occurs when analyte concentration is so high that it
saturates antibody binding sites, causing the measured signal to
decrease at very high concentrations. This function dampens the apparent
decrease beyond the signal peak by compressing the post-peak delta
toward the peak value.

## Usage

``` r
correct_prozone(
  stdframe = NULL,
  prop_diff = NULL,
  dil_scale = NULL,
  response_variable = "mfi",
  independent_variable = "concentration",
  verbose = TRUE
)
```

## Arguments

- stdframe:

  Data frame of standard curve data for a single dilution series. Must
  contain `response_variable` and `independent_variable` columns with no
  NA values.

- prop_diff:

  Numeric. Dampening factor applied to the post-peak signal delta (e.g.
  `0.1`).

- dil_scale:

  Numeric. Dilution scale factor used in the dampening formula (e.g.
  `2`).

- response_variable:

  Character. Name of the response (y) column (default `"mfi"`).

- independent_variable:

  Character. Name of the concentration (x) column (default
  `"concentration"`).

- verbose:

  Logical. If `TRUE` (default), prints diagnostics to the console.

## Value

`stdframe` with post-peak response values adjusted.

## Details

The method follows the computational approach described in:

- ACS Meas. Sci. Au 2024, 4, 4, 452–458.

- Sensors and Actuators B: Chemical, 324, 128756 (2020).

- Sensors and Actuators B: Chemical, 304, 127408 (2020).
