# Safely Collapse Unique Values

Returns the unique non-missing values of a vector as a single
semicolon-separated string. If no non-missing values are present,
returns `NA_character_`.

## Usage

``` r
safe_unique(x)
```

## Arguments

- x:

  A vector of values (character, numeric, or factor).

## Value

A character string containing unique non-missing values separated by
semicolons, or `NA_character_` if no such values exist.

## Details

This is useful for summarizing identifiers or metadata fields that may
contain repeated or missing values within grouped data.

The function removes `NA` values before computing uniqueness. If
multiple unique values remain, they are concatenated using a semicolon
(`";"`) separator.

## Examples

``` r
safe_unique(c("A", "A", "B", NA))
#> [1] "A;B"
# "A;B"

safe_unique(c(NA, NA))
#> [1] NA
# NA

safe_unique(1:3)
#> [1] "1;2;3"
# "1;2;3"
```
