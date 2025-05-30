# WARNING - Generated by {fusen} from dev/flat_first.Rmd: do not edit by hand # nolint: line_length_linter.

#' Useful function to extract mean phenotype values from RITHMS output
#'
#' The gets functions use the output of [holo_simu()] to extract the information of interest from a given generation.
#' `get_mean_phenotypes`extract average phenotype from a generation object.
#' 
#' @importFrom purrr pluck
#' @importFrom magrittr %>%
#' 
#' @param data List corresponding to one generation, as returned by [holo_simu()]. Containing simulation output.
#' @return The mean phenotype value for a given generation.
#' @examples
#' \dontrun{
#' library(purrr)
#' library(magrittr)
#' data("Deru")
#' ToyData <- Deru
#' taxa_assign_g <- assign_taxa(founder_object = ToyData)
#' generations_simu <- holo_simu(h2 = 0.25, b2 = 0.25, founder_object = ToyData,
#'                               n_clust = taxa_assign_g, n_ind = 500,
#'                               verbose = FALSE, seed = 1234)
#'                               
#' # Extract mean phenotype value for G1 generation
#' G1_mean_phenotype <- get_mean_phenotypes(generations_simu$G1)
#'
#' # Extract mean phenotype values of all generations
#' ## Don't forget to substract the metadata
#' mean_phenotypes <- generations_simu[-1] %>% map(get_mean_phenotypes)
#' }
#'
#' @seealso [get_microbiomes()], [get_phenotypes_value()], [get_om_beta_g()], [get_selected_ind()], [get_phenotypes()]
#' @rdname get_mean_phenotypes
#' @export
get_mean_phenotypes <- function(data) {
  data |> pluck("phenotypes","y") |> mean()
}

