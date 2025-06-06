# WARNING - Generated by {fusen} from dev/flat_first.Rmd: do not edit by hand # nolint: line_length_linter.

#' Main function to generate transgenerational hologenomic data
#' 
#' @description Simulates hologenomic data over multiple generations with genetic, microbial, and environmental effects.
#'
#' @importFrom MoBPS get.geno breeding.diploid get.pedigree
#' @importFrom compositions clr clrInv
#' @importFrom glue glue
#' @importFrom phyloseq phyloseq estimate_richness otu_table
#' @importFrom purrr list_c
#' @importFrom magrittr %>%
#' 
#' @param h2 direct heritability value, between 0 and 1. 
#' @param b2 microbiability value, between 0 and 1.
#' @param founder_object output of `generate_founder()` function.
#' @param n_ind number of individual per generation.
#' @param n_clust vector with taxa assignment, typical output of `assign_taxa()`
#' @param n_gen number of generation, DEFAULT = 5
#' @param qtn_y number of causal SNPs for the phenotypes.
#' @param correlation Correlation between taxa within the same cluster, value between 0 and 1, DEFAULT = 0.5
#' @param otu_g percentage of taxa under genetic control, between 0 and 1, DEFAULT = 0.05
#' @param lambda proportion of microbiome of offspring coming from vertical transmission, DEFAULT = 0.5
#' @param effect.size Vector giving the size of genetic effect to try.
#' @param mix.params Vector of two numeric values giving the weights for the regularisation of the base population microbiome. `mix.params[1]` = weight for raw microbiome and `mix.params[2]` = weight for mean microbiome. DEFAULT = c(0.75,0.25).
#' @param mix.params.M A vector of two numeric values specifying weights between Dirichlet samples and the original mean. DEFAULT = c(0.75,0.25). With `mix.params.M[1]` the dirichlet microbiome coefficient and `mix.params.M[2]` the mean microbiome coefficient.
#' @param noise.microbiome sd of microbiome noise, DEFAULT = 0.1
#' @param dir Logical; Mentions if the ambient microbiome is generated via a Dirichlet law or a only the `mean_microbiome`. DEFAULT = T
#' @param ao A numeric scalar used as the concentration parameter for the Dirichlet distribution.
#' @param size_rmultinom Integer; specifying the total number of object for the multinomial sampling(default: 10000, according to DeruPop.rds dataset).
#' @param selection bool, if selection process needed, DEFAULT = FALSE
#' @param size_selection_F percentage of female to select.
#' @param size_selection_M percentage of male to select.
#' @param selection_type mode of selection to be used, value in ("GB", "B", "G", "diversity", "div.GB"), DEFAULT = "GB"
#' @param w.param in case div.GB selection mode is chosen.
#' @param thetaX Optionnal matrix specifying environmental effects applied to the microbiome, such as antibiotic treatment. This matrix should be of dimension `n_taxa x n_individuals`, and typically constructed as the product of:
#' - a vector of taxa-specific effects `theta` (negative or positive values for all taxa).
#' - a binary vector `X` encoding individual exposure (1 for treated, 0 ofr untreated).
#' If `NULL` (default value), no environmental effect is applied.
#' Example use case: refer to the vignette on [Generate figures](https://solenepety.github.io/RITHMS/articles/generate-figures.html#introduction-of-transient-perturbations-of-the-microbiota)
#' @param env_gen vector of booleans.
#' @param seed seed value for samplings in the function.
#' @param verbose bool, DEFAULT = T
#'
#' @section Required parameters:
#' - `h2`
#' - `b2`
#' - `founder_object`
#' - `n_ind`
#' - `n_clust`
#'
#' @section Generation parameters:
#' - `n_gen`
#'
#' @section Genotype-related parameters:
#' - `qtn_y`
#'
#' @section Microbiome-related parameters:
#' - `correlation`
#' - `otu_g`
#' - `lambda`
#' - `effect.size`
#' - `mix.params`
#' - `mix.params.M`
#' - `noise.microbiome`
#' - `dir`
#' - `ao`
#'
#' @section Selection-related parameters:
#' - `selection` 
#' - `size_selection_F`, `size_selection_M` 
#' - `selection_type` 
#'
#' @section Environmental effects:
#' - `thetaX`
#' - `env_gen`
#'
#' @section Misc:
#' - `seed`
#' - `verbose` 
#' 
#' @return
#' An important list with several different objects :
#' - `parameters`: the set of parameters used for the simulation when calling `holo_simu()`
#' - `metadata`: metadata info such as beta matrix details
#' - `G0` to `G5` (5 generations by default): computational characteristics of generations. For each generation, genotypes, microbiomes, phenotypes, pedigree and individuals selected can be reachable.
#' 
#' @seealso [compute_beta_matrix_cluster()], [compute_mean_microbiome()], [compute_current_microbiome()], [compute_phenotypes()]
#' @rdname holo_simu
#' @export
#' @examples
#' \dontrun{
#'   data("Deru")
#'   ToyData <- Deru
#'   taxa_assign_g <- assign_taxa(founder_object = ToyData)
#'   generations_simu <- holo_simu(h2 = 0.25,
#'                                  b2 = 0.25,
#'                                  founder_object = ToyData,
#'                                  n_clust = taxa_assign_g,
#'                                  n_ind = 500,
#'                                  verbose = FALSE,
#'                                  noise.microbiome = 0.5,
#'                                  effect.size = 0.3,
#'                                  lambda = 0.5,
#'                                  dir = TRUE,
#'                                  selection = FALSE,
#'                                  seed = 1234)
#' }  
holo_simu <- function(### Required parameters ###
                      h2,
                      b2,
                      founder_object,
                      n_ind = NULL,
                      n_clust = NULL,
                      ###
                      
                      ### Generations parameters ###
                      n_gen = 5,
                      ###
                      
                      ### Genotypes parameters ###
                      qtn_y = NULL,        
                      ###
                      
                      ### Microbiome-related parameters ###
                      correlation = 0.5,
                      otu_g = 0.05,
                      lambda = 0.5,
                      effect.size = 0.1,
                      mix.params = c(0.75,0.25),
                      mix.params.M = c(0.75,0.25),
                      noise.microbiome = 0.1,
                      dir = T,
                      ao = 25,
                      ###
                      
                      ### Selection-related parameters ###
                      size_rmultinom = 10000,
                      selection = F,       
                      size_selection_F = NULL,
                      size_selection_M = NULL,
                      selection_type = "GB",
                      w.param = c(0.5,0.5),
                      ###
                      
                      ### Environmental-effects-related parameters ###
                      thetaX = NULL,
                      env_gen = NULL,
                      ###
                      
                      ### Additional parameters ###
                      seed=1234,
                      verbose = T
                      ###
                      ){
  
  parameters <- as.list(environment()) %>% within(., rm(founder_object))
  
  set.seed(seed)
  
  population <- founder_object$population
  #obtain genotypes for all generation using MoBPS breeding.diploid
  if(is.null(n_ind)){
    n_ind <- ncol(founder_object$microbiome)
  }
  
  ######
  #call to compute beta matrix function, common for all generations
  ######
  microbiome_filtered <- founder_object$microbiome |> t()
  n_g <- nrow(get.geno(population, gen=1))
  
  #####################
  ### beta matrix generation
  #####################
    n_clust_g <- n_clust[n_clust != 0] %>% unique() %>% length()
    n_qtl_o <- ((n_g * 0.2)/n_clust_g)  %>% round()
    beta <- compute_beta_matrix_cluster(n_b = nrow(microbiome_filtered), 
                                        n_g = n_g, 
                                        n_clust = n_clust, 
                                        n_qtl_o = n_qtl_o, 
                                        n_otus = sum(n_clust != 0), 
                                        correlation = correlation, 
                                        effect_size = effect.size)
  beta_info <- beta$sim_params
  otu_list <- list_c(beta_info$id_otu)
  
  ######
  #loop on generation to obtain microbiome
  ######
  ## ML estimate of founder pop compositions
  microbiome_gen <- apply(microbiome_filtered, 2, \(x) x/sum(x))
  ## shrinkage estimator: weighted average of ML composition and pop-level mean composition followed by clr transform
  microbiome_gen <- (mix.params[1] * microbiome_gen + mix.params[2] * matrix(rowMeans(microbiome_gen),nrow=nrow(microbiome_gen),ncol=ncol(microbiome_gen),byrow=F)) |> t() |> clr() |> t()
  # microbiome_gen <- microbiome_gen |> replace_zero() |> t() |> compositions::clr() |> t()
  #abundances
  taxa_scale <- microbiome_gen |> apply(1, sd)
  colnames(microbiome_gen) <- get.geno(population, gen = 1) %>% colnames()
  
  if(is.null(qtn_y)){
    Nqtl_y = 100
  }else{
    Nqtl_y = qtn_y
  }
  
  #####################
  ### founder phenotypes and parameters calibration
  #####################
  params_phenotypes <- calibrate_params_phenotypes(X0 = get.geno(population,gen = 1),
                                                   B0 = microbiome_gen, #CLR abundances ,
                                                   h2 = h2,
                                                   b2 = b2,
                                                   otu_list = otu_list,
                                                   Nqtl_y = Nqtl_y,
                                                   Notu_y = length(otu_list))
  
  phenotypes_founder <- compute_phenotypes(X = get.geno(population,gen = 1), 
                                           B = microbiome_gen, #CLR abundances 
                                           Nqtl_y = Nqtl_y, 
                                           Notu_y = length(otu_list), 
                                           otu_list = otu_list,
                                           qtl_list = params_phenotypes$qtl_list,
                                           beta_qtl = params_phenotypes$alpha,
                                           beta_otu = params_phenotypes$omega,
                                           se = params_phenotypes$se)
  
  
  #####################
  ### first elements output
  #####################
  metadata <- list(Beta_matrix = beta,
                   QTN_y = params_phenotypes$qtl_list, 
                   G_coeff = params_phenotypes$alpha, 
                   B_coeff = params_phenotypes$omega)
  G0 <- list(microbiome = microbiome_gen |>  t() |> clrInv() |> t(), 
             genotypes = get.geno(population,gen = 1), 
             phenotypes = phenotypes_founder)
  list_output <- list(parameters = parameters, metadata = metadata, G0 = G0)
  
  
  #####################
  #Generation for loop
  ####################
  for(i in 1:n_gen){
    start_g <- proc.time()
    if(verbose){
      print(glue("Generation : {i}"))
    }
    #Warning : outside MoBPS functions, index of generation is i+1
    ######
    # Step 1 : get_current_genotype
    ######
    if(i==1){
      #renames microbiome matrix with individuals names
      population <- breeding.diploid(population,
                                     breeding.size = c(n_ind/2,n_ind/2),
                                     verbose = F)
    }else{
      population <- breeding.diploid(population,
                                     breeding.size = c(n_ind/2,n_ind/2), 
                                     selection.m.database = cbind(i,1,as.numeric(ID_selected$M_id)), 
                                     selection.f.database = cbind(i,2,as.numeric(ID_selected$F_id)),
                                     verbose = F)
    }
    
    ######
    # Step 2 : get_current_microbiome
    ######
    #Extract microbiomes of the mothers
    microbiome_gen <- microbiome_gen |>  t() |> clrInv() |> t()
    current_mother_microbiomes <- sapply(get.pedigree(population, gen=i+1)[,"mother"],function(x){
      microbiome_gen[,x]
    })
    
    mean_microbiome <- compute_mean_microbiome(microbiome_gen, 
                                               dir = dir,
                                               n_ind = n_ind,
                                               ao = ao,
                                               mix.params = mix.params.M)
    
    if(!is.null(env_gen)){
      if(env_gen[i]){
        microbiome_gen <- compute_current_microbiome(beta = beta$matrix,
                                                     current_genotypes = get.geno(population,gen = i+1),
                                                     mother_microbiomes = current_mother_microbiomes,
                                                     mean_microbiome = mean_microbiome,
                                                     noise = noise.microbiome, taxa_scale = taxa_scale, 
                                                     lambda = lambda,
                                                     dir = dir,
                                                     thetaX = thetaX) #not mandatory if we check dimension of mean_microbiome
      }else{
        microbiome_gen <- compute_current_microbiome(beta = beta$matrix,
                                                     current_genotypes = get.geno(population,gen = i+1),
                                                     mother_microbiomes = current_mother_microbiomes,
                                                     mean_microbiome = mean_microbiome,
                                                     noise = noise.microbiome, taxa_scale = taxa_scale, 
                                                     lambda = lambda,
                                                     dir = dir,
                                                     thetaX = NULL)
      }
    }else{
      microbiome_gen <- compute_current_microbiome(beta = beta$matrix,
                                                   current_genotypes = get.geno(population,gen = i+1),
                                                   mother_microbiomes = current_mother_microbiomes,
                                                   mean_microbiome = mean_microbiome,
                                                   noise = noise.microbiome, taxa_scale = taxa_scale, 
                                                   lambda = lambda,
                                                   dir = dir,
                                                   thetaX = NULL)
    }
    ######
    # Step 3 : get_current_phenotype
    ######
    phenotypes <- compute_phenotypes(X = get.geno(population,gen = i+1), 
                                     B = microbiome_gen, #CLR abundances 
                                     Nqtl_y = Nqtl_y, 
                                     Notu_y=length(otu_list), 
                                     otu_list = otu_list,
                                     qtl_list = params_phenotypes$qtl_list,
                                     beta_qtl = params_phenotypes$alpha,
                                     beta_otu = params_phenotypes$omega,
                                     se = params_phenotypes$se)
    
    
    ######
    # Step 4 : compute TBV to select individuals
    ######
    
    ID_selected <- select_individual(phenotypes = phenotypes,
                                     microbiomes = microbiome_gen |>  t() |> clrInv() |> t(),
                                     genotypes = get.geno(population,gen = 1+i),
                                     beta = beta$matrix,
                                     beta_otu = params_phenotypes$omega,
                                     selection = selection,
                                     size_selection_F = size_selection_F,
                                     size_selection_M = size_selection_M,
                                     selection_type = selection_type,
                                     size_rmultinom = size_rmultinom,
                                     w.param = w.param) 

    ###
    #update output
    ###
    list_output[[glue("G{i}")]] <- list(microbiome = microbiome_gen |>  t() |> clrInv() |> t(),
                                        genotypes = get.geno(population,gen = 1+i),
                                        phenotypes = phenotypes,
                                        pedigree = get.pedigree(population,gen = i+1),
                                        mean_microbiome = mean_microbiome,
                                        noise_microbiome = attr(microbiome_gen,"noise_microbiome"),
                                        ID_selected = c(glue("F_{ID_selected$F_id}"),glue("M_{ID_selected$M_id}")),
                                        omega_beta_g = as.vector(params_phenotypes$omega %*% (beta$matrix[rowSums(beta$matrix) != 0, ] %*% get.geno(population,gen = 1+i))))
    if(verbose){
      end_g <- proc.time() - start_g
      print(glue("Executed in {round(end_g[1] |> as.numeric(),3)} seconds"))
    }
  }
  return(list_output)
}



#' Estimate diversity metrics from relative abundances
#' 
#' This function estimate diversity metrics (Observed, Shannon, Inverse Simpson) from the matrix of relative abundances (see [get_microbiomes()]).
#' It uses multinomial sampling to simulate read counts from abundances, and computes diversity metrics across `n_loop` in order to obtain robust estimation. This function is particularly useful when selection is based on diversity.
#' 
#' @importFrom phyloseq phyloseq estimate_richness otu_table
#' 
#' @param microbiome_matrix A matrix of relative abundances (individuals in rows and OTUs in columns, see [get_microbiomes()] output).
#' @param size_rmultinom Integer; specifying the total number of object for the multinomial sampling(default: 10000, according to DeruPop.rds dataset).
#' @param n_loop Integer; number of multinomial resampling iterations to perform (default: 10).
#' @param plot Logical; not currently implemented
#' 
#' @return A `data.frame`of average diversity metrics (Observed, Shannon, Inverse Simpson) for each sample.
#'
#' @seealso [get_microbiomes()], [phyloseq::estimate_richness()]
#' @export
#' @examples
#' \dontrun{
#' library(magrittr)
#' library(purrr)
#' data("Deru")
#' ToyData <- Deru
#' taxa_assign_g <- assign_taxa(founder_object = ToyData)
#' generations_simu <- holo_simu(h2 = 0.25, b2 = 0.25, founder_object = ToyData,
#'                                n_clust = taxa_assign_g, n_ind = 500,
#'                                verbose = FALSE, seed = 1234)
#'                                
#' # Extract microbiomes matrix for each generations
#' microbiomes <- generations_simu[-1] %>% map(get_microbiomes)
#'  
#' # Estimate diversity metrics
#' richness_from_abundances <- microbiomes %>% map(richness_from_abundances_gen, size_rmultinom = 10000) 
#' ## size_rmultinom = 10000 according to DeruPops dataset
#' }
richness_from_abundances_gen <- function(microbiome_matrix, size_rmultinom = 10000, n_loop = 10, plot=T){
  microbiome_matrix[microbiome_matrix<0] <- 0
  for(i in 1:n_loop){
    microbiome_matrix <- apply(microbiome_matrix, 2, function(x){
      rmultinom(1,size_rmultinom,as.vector(x))
    })
    tmp_physeq <- phyloseq(otu_table(microbiome_matrix,taxa_are_rows = T))
    if(i==1){
      richness <- estimate_richness(tmp_physeq,measures = c("Observed","Shannon","InvSimpson"))
    }else{
      richness <- richness + estimate_richness(tmp_physeq,measures = c("Observed","Shannon","InvSimpson"))
    }
  }
  return(richness/n_loop)
}

