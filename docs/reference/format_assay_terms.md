# Format common assay variable names for display

Normalises and maps assay-related variable name strings (e.g. column
names such as `"mfi"` or `"concentration"`) to their conventional
display forms. Matching is case-insensitive. Unrecognized values are
returned unchanged.

## Usage

``` r
format_assay_terms(x)
```

## Arguments

- x:

  Character vector of variable names to format.

## Value

Character vector the same length as `x` with display-ready names. Known
mappings are: `"mfi"` → `"MFI"`, `"absorbance"` → `"Absorbance"`,
`"fluorescence"` → `"Fluorescence"`, `"od"` → `"OD"`, `"concentration"`
→ `"Concentration"`.
