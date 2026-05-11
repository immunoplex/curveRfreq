# Null-coalescing operator

Returns the left-hand side if it is not `NULL`, otherwise returns the
right-hand side.

## Usage

``` r
a %||% b
```

## Arguments

- a:

  First value.

- b:

  Fallback value if `a` is `NULL`.

## Value

The first non-NULL value.

## Examples

``` r
NULL %||% 5
#> [1] 5
10 %||% 5
#> [1] 10
```
