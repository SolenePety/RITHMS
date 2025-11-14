# Rarefaction step within read_input_data

Rarefaction step within read_input_data

## Usage

``` r
rarefied_microbiome(microbiome_matrix, id_column = "ind_id", threshold = 0.05)
```

## Arguments

- microbiome_matrix:

  Microbiome counting table, OTUs in columns and a column specifying the
  IDs

- id_column:

  String specifying the name of the column containing the individual
  IDs. Default is "ind_id"

- threshold:

  Minimum prevalence threshold

## Value

The filtered microbiome matrix, without OTUs whose prevalence is lower
than the set threshold

## See also

[`read_input_data()`](read_input_data.md),
[`generate_founder()`](generate_founder.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  microbiome <- data.table::fread("/path/to/microbiome.txt")
  microbiome_filtered <- rarefied_microbiome(microbiome, id_columns = "sample_id", threshold = 0.05)
} # }
```
