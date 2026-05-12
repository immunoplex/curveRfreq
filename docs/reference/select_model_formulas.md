# Select Model Formulas for Standard Curve Fitting

Returns a named list of
[`stats::nls()`](https://rdrr.io/r/stats/nls.html)-compatible formulas
for all candidate models (logistic5, loglogistic5, logistic4,
loglogistic4, gompertz4). When `fixed_constraint` is supplied, parameter
`a` (lower asymptote) is substituted as a numeric constant in every
formula; otherwise `a` is a free parameter.

## Usage

``` r
select_model_formulas(
  fixed_constraint,
  response_variable,
  is_log_response,
  model_names
)
```

## Arguments

- fixed_constraint:

  Numeric or `NULL`. If non-`NULL`, the lower asymptote `a` is fixed to
  this value in every formula. Typically derived from blank-well
  measurements. Values that are non-positive or non-finite are rejected
  with a message, falling back to free-`a` formulas.

- response_variable:

  Character. Column name for the assay response (e.g., `"mfi"`, `"od"`).
  Used as the left-hand side of each formula.

- is_log_response:

  Logical. If `TRUE`, `fixed_constraint` is \\\log\_{10}\\-transformed
  before substitution (with a small epsilon to avoid \\\log(0)\\).

- model_names:

  a vector of model formulas to be considered.

## Value

Named list of formulas keyed by model name: `"logistic5"`,
`"loglogistic5"`, `"logistic4"`, `"loglogistic4"`, `"gompertz4"`.

## Details

The concentration column is hard-coded as `concentration` in all
formulas. The [`I()`](https://rdrr.io/r/base/AsIs.html) wrappers protect
complex sub-expressions from formula parsing by
[`stats::nls()`](https://rdrr.io/r/stats/nls.html).

## See also

[`logistic4()`](https://immunoplex.github.io/curveRfreq/reference/logistic4.md),
[`logistic5()`](https://immunoplex.github.io/curveRfreq/reference/logistic5.md),
[`loglogistic4()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic4.md),
[`loglogistic5()`](https://immunoplex.github.io/curveRfreq/reference/loglogistic5.md),
[`gompertz4()`](https://immunoplex.github.io/curveRfreq/reference/gompertz4.md)
for the underlying model functions;
[`compute_robust_curves()`](https://immunoplex.github.io/curveRfreq/reference/compute_robust_curves.md)
which consumes these formulas.

## Examples

``` r
# Fixed lower asymptote
forms <- select_model_formulas(fixed_constraint = 50,
                               response_variable = "mfi",
                               is_log_response = FALSE,
                               model_names = c("logistic4", "logistic5",
                               "loglogistic4", "loglogistic5",
                               "gompertz4"))
#> Lower asymptote is fixed at50
names(forms)
#> [1] "logistic4"    "logistic5"    "loglogistic4" "loglogistic5" "gompertz4"   
forms$logistic4
#> mfi ~ d + (((50) - d)/I((1 + exp((concentration - c)/b))))
#> <environment: 0x12ea97278>

# Free lower asymptote
forms_free <- select_model_formulas(fixed_constraint = NULL,
                                    response_variable = "mfi",
                                    is_log_response = FALSE,
                                    model_names = c("logistic4", "logistic5",
                                     "loglogistic4", "loglogistic5",
                                      "gompertz4"))
forms_free$logistic4
#> mfi ~ d + (a - d)/I(1 + exp((concentration - c)/b))
#> <environment: 0x12ea97278>
```
