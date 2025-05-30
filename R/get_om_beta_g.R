# WARNING - Generated by {fusen} from dev/flat_first.Rmd: do not edit by hand # nolint: line_length_linter.

#' Useful function to extract omega beta G values from RITHMS output
#' 
#' The gets functions use the output of [holo_simu()] to extract the information of interest from a given generation.
#' `get_om_beta_g`extract omega beta G values from a generation object.
#' 
#' @importFrom purrr pluck
#' @importFrom magrittr %>%
#' 
#' @param data List corresponding to one generation, as returned by [holo_simu()]. Containing simulation output.
#' @return A `data.frame` of omega beta G values for each individuals of a given generation.
#' 
#' @examples
#' \dontrun{
#' library(magrittr)
#' library(purrr)
#' data("Deru")
#' ToyData <- Deru
#' taxa_assign_g <- assign_taxa(founder_object = ToyData)
#' generations_simu <- holo_simu(h2 = 0.25, b2 = 0.25, founder_object = ToyData,
#'                               n_clust = taxa_assign_g, n_ind = 500,
#'                               verbose = FALSE, seed = 1234)
#'                               
#' #Extract omega beta G values for each individuals of G1 generation
#' G1_om_beta_g <- get_om_beta_g(generations_simu$G1)
#'  
#' #Extract omega beta G values for each individuals of all generations
#' ## Don't forget to substract the metadata
#' om_beta_g <- generations_simu[-1] %>% map(get_om_beta_g)
#'}
#'  
#' @seealso [get_mean_phenotypes()], [get_phenotypes_value()], [get_microbiomes()], [get_selected_ind()], [get_phenotypes()]
#' @rdname get_om_beta_g
#' @export
get_om_beta_g <- function(data) {
  data |> pluck("omega_beta_g")
}
