# Identify Free Parameters in a List of Model Formulae

Identify Free Parameters in a List of Model Formulae

## Usage

``` r
obtain_free_variables(formulas, dep = "mfi", indep = "concentration")
```

## Arguments

- formulas:

  Named list of `formula` objects.

- dep:

  Character. Dependent variable name (default `"mfi"`).

- indep:

  Character. Independent variable name (default `"concentration"`).

## Value

Named list of character vectors, one per formula, containing the free
parameter names in alphabetical order.
