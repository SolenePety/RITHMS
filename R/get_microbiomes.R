# WARNING - Generated by {fusen} from dev/flat_first.Rmd: do not edit by hand # nolint: line_length_linter.

#' Useful function to extract microbiomes from RITHMS output
#' 
#' The gets functions use the output of [holo_simu()] to extract the information of interest from a given generation.
#' `get_microbiomes` extract the microbiome abundance matrix from a generation object, with or without CLR transformation or transposition.
#' 
#' @importFrom compositions clr
#' @importFrom purrr pluck
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' 
#' @param data List corresponding to one generation, as returned by `holo_simu()`. Containing simulation output.
#' @param transpose Logical; if `TRUE`, transpose the microbiome matrix (OTUs in rows, individuals in columns).
#' @param CLR Logical; if `TRUE`, applies a CLR transformation to the abundance data. This transformation requires `transpose = TRUE`.
#' 
#' @return A `data.frame`containing the microbiome abundances of individuals. Default, individuals are in rows and OTUs in columns. Change `transpose`parameter if needed.
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
#' # Extract microbiome matrix for G1 generation
#' G1_microbiome <- get_microbiomes(generations_simu$G1)
#'
#' # Extract with transposition
#' G1_t_microbiome <- get_microbiomes(generations_simu$G1, transpose = TRUE)
#'
#' # Extract with CLR transformation
#' G1_CLR_microbiome <- get_microbiomes(generations_simu$G1, transpose = TRUE, CLR = TRUE)
#'
#' # Extract all microbiome matrices of all generations
#' # substract metadata
#' microbiomes <- generations_simu[-1] %>% map(get_microbiomes)
#' }
#'
#' @seealso [get_mean_phenotypes()], [get_phenotypes_value()], [get_om_beta_g()], [get_selected_ind()], [get_phenotypes()]
#' @rdname get_microbiomes
#' @export
get_microbiomes <- function(data, transpose = F, CLR = F) {
  if(transpose){
    if(CLR){
      return(data |> pluck("microbiome") |> t() |> replace_zero() |> clr() |> as.data.frame())
    }else{
      return(data |> pluck("microbiome") |> t() |> as.data.frame())
    }
  }
  data |> pluck("microbiome") |> as.data.frame()
}
