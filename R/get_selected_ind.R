# WARNING - Generated by {fusen} from dev/flat_first.Rmd: do not edit by hand # nolint: line_length_linter.

#' Useful function to extract ID of selected individuals from RITHMS output
#'
#'The gets functions use the output of [holo_simu()] to extract the information of interest from a given generation.
#' `get_selected_ind`extract selectd individuals IDs from a generation object.
#' 
#' @importFrom purrr pluck
#' @importFrom magrittr %>%
#' 
#' @param data List corresponding to one generation, as returned by [holo_simu()]. Containing simulation output.
#' 
#' @return A`list` of the selected individuals of a given generation.
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
#' #Extract selected individuals IDs for the G1 generation
#' G1_selected_ind <- get_selected_ind(generations_simu$G1)
#'  
#' #Extract selected individuals IDs for all generations
#' ## Don't forget to substract the metadata
#' selected_ind <- generations_simu[-1] %>% map(get_selected_ind)
#' }
#'
#' @seealso [get_mean_phenotypes()], [get_phenotypes_value()], [get_om_beta_g()], [get_microbiomes()], [get_phenotypes()]
#' @rdname get_selected_ind
#' @export
get_selected_ind <- function(data) {
  data |> pluck("ID_selected")
}
