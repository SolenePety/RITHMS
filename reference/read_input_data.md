# Formatting data from file paths to base population object

Formatting data from file paths to base population object

## Usage

``` r
read_input_data(
  path_to_microbiome,
  path_to_genotype,
  file_type = "pedmap",
  biome_id_column = "ind_id",
  threshold = 0.05,
  ind_selected = NULL
)
```

## Arguments

- path_to_microbiome:

  String giving the path to count table file.
  `"path/to/microbiome.{extension}"`

- path_to_genotype:

  String giving the path and prefix to the genotype file. This should be
  one of `"pedmap"` or `"vcf"`. (default: `"pedmap"`), see `file_type`
  param. For a `"example_pedmap.ped"` file enter only
  `"path/to/example_pedmap"`.

- file_type:

  String specifying the file type used to load genotype data. This
  should be one of `"pedmap"` or `"vcf"`. (default: `"pedmap"`).

- biome_id_column:

  String specifying the name of the column containing the individual IDs
  for the microbiome matrix. Default is "ind_id".

- threshold:

  Threshold for rarefaction, DEFAULT = 0.05

- ind_selected:

  Vector of string values with individuals to keep, have to match
  rownames of count table file, DEFAULT = NULL

## Value

A `list` containing two object. A `matrix` containing the
microbiome(`microbiome`). A `list` from
[MoBPS](https://pubmed.ncbi.nlm.nih.gov/32229505/) package, containing
the genotypes (`population`). Individuals are in columns and taxa or
SNPs are in rows.

## Note

The number of individuals in the microbiome an genotype data must be
equivalent.

## See also

[`generate_founder()`](generate_founder.md),
[`rarefied_microbiome()`](rarefied_microbiome.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create founder object from PED/MAP set
founder_object <- read_input_data(path_to_microbiome = "/path/to/microbiome.txt",
                                  path_to_pedmap = "/path/to/pedmap/'prefix'",
                                  biome_id_column = "ind_id",
                                  threshold = 0.05,
                                  ind_selected = NULL,
                                  file_type = "pedmap")

# Create founder object from VCF
founder_object <- read_input_data(path_to_microbiome = "/path/to/microbiome.txt",
                                  path_to_pedmap = "/path/to/vcf/'prefix'",
                                  biome_id_column = "ind_id",
                                  threshold = 0.05,
                                  ind_selected = NULL,
                                  file_type = "vcf")
} # }
```
