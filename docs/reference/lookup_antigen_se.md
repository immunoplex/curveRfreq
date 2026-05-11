# Look Up SE for a Specific Antigen from the SE Table

Look Up SE for a Specific Antigen from the SE Table

## Usage

``` r
lookup_antigen_se(
  se_table,
  study_accession,
  experiment_accession,
  source,
  antigen,
  feature
)
```

## Arguments

- se_table:

  data.frame from compute_antigen_se_table()

- study_accession:

  study identifier

- experiment_accession:

  experiment identifier

- source:

  source identifier

- antigen:

  antigen identifier

- feature:

  feature identifier

## Value

numeric SE value, or NA if not found
