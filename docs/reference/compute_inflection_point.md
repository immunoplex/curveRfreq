# Compute Inflection Point of Fitted Curve

Calculates the inflection point (x and y coordinates) for a fitted
nonlinear model. The x-coordinate is computed analytically based on the
model type, and the y-coordinate is obtained by evaluating the fitted
model at that x-value (via
[`predict()`](https://rdrr.io/r/stats/predict.html) when possible).

## Usage

``` r
compute_inflection_point(
  model_name,
  fit,
  fixed_a_result,
  independent_variable,
  verbose = TRUE
)
```

## Arguments

- model_name:

  Character string specifying the model type. Supported values include:

  - `"logistic5"`

  - `"logistic4"`

  - `"loglogistic5"`

  - `"loglogistic4"`

  - `"gompertz4"`

- fit:

  A fitted model object with coefficients accessible via
  [`coef()`](https://rdrr.io/r/stats/coef.html) and compatible with
  [`predict()`](https://rdrr.io/r/stats/predict.html).

- fixed_a_result:

  Optional numeric value for parameter `a`. If provided, this overrides
  the estimated `a` from the model coefficients.

- independent_variable:

  Character string giving the name of the predictor variable used in the
  model. This is required to construct `newdata` for prediction.

- verbose:

  Logical; if `TRUE`, prints a message including the computed inflection
  point coordinates.

## Value

A named list with:

- `inflect_x` Numeric x-coordinate of the inflection point

- `inflect_y` Numeric y-coordinate of the inflection point

## Details

Supported model types include 4- and 5-parameter logistic and
log-logistic models, as well as the 4-parameter Gompertz model.

The inflection point is defined as the point on the curve where the
second derivative equals zero. For supported models, the x-coordinate is
computed analytically:

- Logistic 5-parameter: \\x = c - b \log(g)\\

- Logistic 4-parameter: \\x = c\\

- Log-logistic 5-parameter: \\x = c + \log(g)/b\\

- Log-logistic 4-parameter: \\x = c\\

- Gompertz 4-parameter: \\x = c\\

The y-coordinate is computed by evaluating the fitted model at the
inflection x-value. If
[`predict()`](https://rdrr.io/r/stats/predict.html) fails, an analytical
expression is used as a fallback.

When `verbose = TRUE`, the function prints the inflection point as a
coordinate:

    [compute_inflection_point] Inflection point: (x = ..., y = ...)

## Examples

``` r
if (FALSE) { # \dontrun{
result <- compute_inflection_point(
  model_name = "logistic5",
  fit = fitted_model,
  fixed_a_result = NULL,
  independent_variable = "log_dilution"
)

result$inflect_x
result$inflect_y
} # }
```
