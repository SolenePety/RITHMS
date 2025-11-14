# Convert a 0/1/2 genotype matrix into a VCF-like format

This function converts a genotype matrix encoded as 0,1,2 into a
VCF-like format. It can either return the VCF content as a `data.frame`,
write it to a .vcf file or do both.

## Usage

``` r
transform_geno_into_vcf(
  geno_matrix,
  output_type = c("file", "dataframe", "both"),
  output_path = NULL
)
```

## Arguments

- geno_matrix:

  A matrix of genotypes with values 0, 1, 2. SNPs are in rows and
  individuals in columns.

- output_type:

  A character, that specifies the output type. Choose between `"file"`
  (write `.vcf`), `"dataframe"` (return `data.frame`), or `"both"` (do
  both).

- output_path:

  A character string that specifies the output path. Requiered if
  `output_type` is "file" or "both".

## Value

A `data.frame` if `output_type = "dataframe" or "both"`, or just write
the `.vcf` file if `output_type = "file" or "both"`.

## Examples

``` r
n_ind <- 50
n_snp <- 50

# Create a small genotype matrix
set.seed(123)
geno_matrix <- matrix(sample(0:2, size = n_ind*n_snp, replace = TRUE),
                      nrow = n_snp, ncol = n_ind)

# Transform genotype matrix into VCF
vcf <- transform_geno_into_vcf(geno_matrix, output_type = "dataframe")
head(vcf[1:3,1:7])
#>   CHROM POS  ID REF ALT QUAL FILTER
#> 1     1   1 rs1   A   T    .      .
#> 2     1   2 rs2   A   T    .      .
#> 3     1   3 rs3   A   T    .      .
```
