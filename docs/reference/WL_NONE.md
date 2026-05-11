# Sentinel value for missing wavelength

Constant used to represent a missing or undefined wavelength in grouping
keys. This avoids the use of `NA` or `NULL`, ensuring compatibility with
SQL `UNIQUE` constraints and stable joins in R workflows, particularly
for bead-array data where wavelength may not be defined.

## Usage

``` r
WL_NONE
```

## Format

A character scalar.

## Details

This value must match the sentinel used for assays that do not have a
wavelength (bead arrays)

## Examples

``` r
WL_NONE
#> [1] "__none__"
```
