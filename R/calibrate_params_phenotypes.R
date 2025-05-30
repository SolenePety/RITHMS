# WARNING - Generated by {fusen} from dev/flat_first.Rmd: do not edit by hand # nolint: line_length_linter.

#' Calibration coefficients to compute phenotypes
#' 
#' This function is part of the first step of [holo_simu()], focused on calibrating phynotypes parameters before looping on the generations. All coefficents are rescale during the process to ensure that standard deviation = 1 and the variance of microbiote effect and genetic effect satisfy the target microbiability and direct heritability.
#' 
#' @inheritParams holo_simu
#' 
#' @param X0 Matrix of the primary generation genotypes, given from \link[MoBPS:get.geno]{get.geno()} MoBPS function. SNPs are in rows and individuals in columns.
#' @param B0 Matrix of the primary generation microbiomes, as CLR abundances. OTUs are in rows and individuals in columns.
#' @param otu_list List of causal OTUs for the phenotypes.
#' @param Nqtl_y Integer; number of causal SNPs for the phenotypes.
#' @param Notu_y Integer; number of causal OTUs for the phenotypes.
#' 
#' @examples
#' library(magrittr)
#' n_ind <- 10
#' n_snp <- 50
#' n_otu <- 50
#'
#' # Simulate a small genotype matrix 
#' set.seed(123)
#' X0 <- matrix(sample(0:2, n_snp * n_ind, replace = TRUE), nrow = n_snp, ncol = n_ind)
#'
#' # Simulate a small microbiome counts matrix
#' B0_counts_table <- matrix(abs(rnorm(n_snp * n_ind, mean = 10, sd = 3)), nrow = n_otu, ncol = n_ind)
#'
#' # Transform to relative abundances per individuals and then apply CLR transformation
#' B0_abund <- apply(B0_counts_table, 2, function(x) x/sum(x))
#' B0 <- compositions::clr(t(B0_abund)) %>% t()
#' 
#' # Randomly select causal OTUs
#' rownames(B0) <- paste0("OTU", 1:n_otu)
#' otu_list <- c(y = sample(rownames(B0), 10))
#' 
#' 
#' params <- suppressWarnings(RITHMS:::calibrate_params_phenotypes(X0 = X0,
#'                                       B0 = B0,
#'                                       h2 = 0.25,
#'                                       b2 = 0.25,
#'                                       otu_list = otu_list,
#'                                       Nqtl_y = 10,
#'                                       Notu_y = length(otu_list)))
#' str(params)
#'  
#' @return A `list` of phenotype parameters such as alpha, omega, list of causal SNPs for the phenotypes, standard deviation (=1) of phenotypes
#' 
#' @seealso [compute_phenotypes()]
#' @rdname calibrate_params_phenotypes

calibrate_params_phenotypes <- function(X0,
                                        B0,
                                        h2,
                                        b2,
                                        otu_list,
                                        Nqtl_y,
                                        Notu_y = length(otu_list)){
  if (Notu_y == 0) b2 = 0
  if (Nqtl_y == 0) h2 = 0
  
  ## heritability and microbiability should be non-negative and the sum should be positive
  stopifnot(h2 >= 0, b2 >= 0, h2 + b2 > 0)
  
  #init 
  Nind = ncol(X0)
  Nsnp = nrow(X0)
  Notu = nrow(B0)
  Bclust=seq(Notu)
  Nclust = max(Bclust)
  gb = rep(0,Nind)
  gq = gb
  
  ## list of snps directly influencing y :
  qtl_list  = sample(seq(Nsnp), size = Nqtl_y)
  beta_qtl = rgamma(Nqtl_y, shape = 0.4, scale = 5) * sample(c(1,-1), size=Nqtl_y, replace = T)
  beta_otu = rgamma(Notu_y, shape = 1.4, scale = 3.8) * sample(c(1,-1), size = Notu_y, replace = T)
  
  # OTU part  
  if(Notu_y>0){ #arrow from B to y
    B <- as.matrix(B0)
    gb = t(B[otu_list, , drop = FALSE]) %*% beta_otu
    varb <- var(gb)
  }
  
  ## If not genetic effect:
  ## - no scaling for microbiote parameters
  ## - scale noise variance using microbiote as a reference
  if (h2 == 0) {
    
    ## genetic effects
    beta_qtl[] <- 0
    gq <- rep(0, Nind)
    
    ## noise scaling
    se <- sqrt(varb * (1-b2) / b2)
  } else {
    ## Scale microbiote and noise variance using genetic effects as reference
    
    ## genetic effects
    stopifnot(Nqtl_y > 0) ## Nqtl_y should be positive if h2 > 0
    X <- as.matrix(X0)
    gq = t(X[qtl_list, , drop = FALSE]) %*% beta_qtl
    varg <- var(gq)
    
    ## microbiote and noise scaling
    microbiote_scaling <- sqrt( (b2 * varg) / (h2 * varb) )
    beta_otu[] <- microbiote_scaling * beta_otu
    se <- sqrt(varg * (1- h2 - b2) / h2)
  }
  
  ## Rescale all coefficents to ensure that se = 1 and the variances of gb and gq satisfy the 
  ## target microbiability and "heritability" 
  return(list('alpha' = beta_qtl / c(se), 
              'omega' = beta_otu / c(se),
              'qtl_list' = qtl_list,
              'se' = se / se ## = 1
  ))
}

