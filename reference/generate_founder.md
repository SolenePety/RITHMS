# Formatting of ped/map files or vcf into haplotypes table and add it inside the `founder_object` list.

Formatting of ped/map files or vcf into haplotypes table and add it
inside the `founder_object` list.

## Usage

``` r
generate_founder(path = NULL, microbiome_matrix, file_type = "pedmap")
```

## Arguments

- path:

  String giving the path and prefix to ped and map file

- microbiome_matrix:

  Filtered microbiome matrix, with OTUs in columns and individuals in
  rows

- file_type:

  String specifying the file type used to load genotype data. This
  should be one of `"pedmap"` or `"vcf"`. (default: `"pedmap"`).

## Value

A `list` of two `matrix`. The `microbiome` matrix contains taxa (in
columns) accross individuals (in rows). The `population` object contains
genotypes, encoded as 0,1,2 and generated using the
[MoBPS](https://pubmed.ncbi.nlm.nih.gov/32229505/) package.

## See also

[`read_input_data()`](read_input_data.md),
[`rarefied_microbiome()`](rarefied_microbiome.md)

## Examples

``` r
if (FALSE) { # \dontrun{
  # Generate founder object from PED/MAP set
  microbiome <- data.table::fread("/path/to/microbiome.txt")
  founder_object <- generate_founder(path = "/path/to/pedmap", microbiome, file_type = "pedmap")
  
  # Generate founder object from VCF set
  microbiome <- data.table::fread("/path/to/microbiome.txt")
  founder_object <- generate_founder(path = "/path/to/vcf", microbiome, file_type = "vcf") 
} # }
```
