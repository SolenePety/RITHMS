# Package index

## Formatting data functions

Set of useful functions for preparing input data

- [`read_input_data()`](read_input_data.md) : Formatting data from file
  paths to base population object

- [`rarefied_microbiome()`](rarefied_microbiome.md) : Rarefaction step
  within read_input_data

- [`generate_founder()`](generate_founder.md) :

  Formatting of ped/map files or vcf into haplotypes table and add it
  inside the `founder_object` list.

## Taxa assignation function

Assign all taxa to a cluster

- [`assign_taxa()`](formatting_data.md) : Assign all taxa to a cluster,
  eventually under genetic control using hclust

## Calibration of genetic effect

- [`gen_effect_calibration()`](gen_effect_calibration.md) : Calibration
  of genetic effect from founder population data

## HoloSimu

The main function

- [`holo_simu()`](holo_simu.md) : Main function to generate
  transgenerational hologenomic data

## Useful function relative to HoloSimu

Set of functions used within HoloSimu to generate the components of each
generation, such as phenotypes and microbiomes.

- [`compute_beta_g()`](compute_beta_g.md) : Compute product of matrixes
  based on few parameters
- [`compute_current_microbiome()`](compute_current_microbiome.md) :
  Compute microbiomes for all individuals of current generation
  gathering all relevant object already computed.
- [`compute_mean_microbiome()`](compute_mean_microbiome.md) : Compute
  the mean for each taxa across population
- [`compute_phenotypes()`](compute_phenotypes.md) : Compute phenotype
  values based on generated objects of the current generation
- [`calibrate_params_phenotypes()`](calibrate_params_phenotypes.md) :
  Calibration coefficients to compute phenotypes
- [`compute_beta_matrix_cluster()`](compute_beta_matrix_cluster.md) :
  Generate beta matrix giving genetic effect per SNP on taxa abundances.

## Get functions

Useful get functions for accessing the components of each generation
from HoloSimu output.

- [`get_microbiomes()`](get_microbiomes.md) : Useful function to extract
  microbiomes from RITHMS output
- [`get_mean_phenotypes()`](get_mean_phenotypes.md) : Useful function to
  extract mean phenotype values from RITHMS output
- [`get_phenotypes_value()`](get_phenotypes_value.md) : Useful function
  to extract phenotype values from RITHMS output ("gq" by default)
- [`get_phenotypes()`](get_phenotypes.md) : Useful function to extract
  all phenotype values as data frame from RITHMS output
- [`get_om_beta_g()`](get_om_beta_g.md) : Useful function to extract
  omega beta G values from RITHMS output
- [`get_selected_ind()`](get_selected_ind.md) : Useful function to
  extract ID of selected individuals from RITHMS output
- [`get_mean_diversity()`](get_mean_diversity.md) : Useful function to
  extract Shannon diversity from RITHMS output

## Estimate diversity metrics

- [`richness_from_abundances_gen()`](richness_from_abundances_gen.md) :
  Estimate diversity metrics from relative abundances

## Another little useful functions across the HoloSimu process

- [`replace_zero()`](replace_zero.md) : Replace zero in the matrix by
  1e-15 to avoid infinite values when CLR transformation is computed
- [`root()`](root.md) : Compute root value for beta matrix construction
- [`select_individual()`](select_individual.md) : Obtain id of the
  parents selected for the next generation based on the criteria chosen
  by the user
- [`transform_geno_into_vcf()`](transform_geno_into_vcf.md) : Convert a
  0/1/2 genotype matrix into a VCF-like format

## Dataset

- [`Deru`](Deru.md) : Deru Dataset
