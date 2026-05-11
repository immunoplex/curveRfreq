# Extract the Response Variable Name from a List of Formulae

Extract the Response Variable Name from a List of Formulae

## Usage

``` r
obtain_response_variable(formulas)
```

## Arguments

- formulas:

  Named list of `formula` objects.

## Value

Character scalar: the response (LHS) variable name. Throws an error if
the formulas disagree on the response variable.
