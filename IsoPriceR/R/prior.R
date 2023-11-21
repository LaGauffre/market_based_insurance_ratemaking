# titre: prior.r
# auteur: Pierre-O Goffard
# date: 05/08/2023
# Functions to set up the prior distribution over one parameters

#' prior_dist is an object that corresponds to the prior distribution of some parameter
#'
#' @param dist_name string the name of the distribution  (e.g. uniform)
#' @param parm_name string name of the model parameter
#' @param hyper_parms vector value of tyhe hyper parameter of the prior distribution
#'
#' @return a prior distribution for a parameter
#' @export
#'
#' @examples
#' alpha_prior <- prior_dist("uniform", "alpha", c(0,10))
prior_dist <- function(dist_name, parm_name, hyper_parms){
  characteristics <- list(dist_name = dist_name, parm_name = parm_name,
                          hyper_parms = hyper_parms)
  class(characteristics) <- "prior_dist"
  return(characteristics)
}



#' sample_prior generate sample from teh prior distribution
#'
#' @param prior_dist object of type prior distribution
#' @param R scalar size of the prior sample
#'
#' @return a sample generated fgrom the prior distribution
#' @export
#'
#' @examples
# alpha_prior <- prior_dist("uniform", "alpha", c(0,10))
# sample_prior(alpha_prior, 1)
sample_prior <- function(prior_dist, R){
  if(prior_dist$dist_name == "uniform"){
    runif(R, min = prior_dist$hyper_parm[1], max= prior_dist$hyper_parm[2])
  }else if(prior_dist$dist_name == "constant"){
    rep(prior_dist$hyper_parms[1], R)
  }
  else{
    print("That prior distribution is not considered")
  }
}


#' logd_prior compute the probability density function associated to the prior distribution
#'
#' @param prior_dist object of type prior distribution
#' @param x vector points at which the pdf is being evaluated
#'
#' @return value of the pdf of the prior distribution
#' @export
#'
#' @examples
#' alpha_prior <- prior_dist("uniform", "alpha", c(0,10))
#' alphas <- sample_prior(alpha_prior, 10)
#' logd_prior(alpha_prior, alphas)
logd_prior <- function(prior_dist, x){
  if(prior_dist$dist_name == "uniform"){
    dunif(x, min = prior_dist$hyper_parm[1], max= prior_dist$hyper_parm[2], log)
  }else if(prior_dist$dist_name == "constant"){
      rep(0, length(x))
  }else{
    print("That prior distribution is not considered")
  }
}



#' independent_priors function to set up the independent priors over all the parameters
#'
#' @param prior_distributions list list of prior distribution sto be combined
#'
#' @return an object of class independent priors
#' @export
#'
#' @examples
#' alpha_prior <- prior_dist("uniform", "alpha", c(0,10))
#' beta_prior <- prior_dist("constant", "beta", c(1))
#' independent_priors(list(alpha_prior, beta_prior))
independent_priors <- function(prior_distributions){
  characteristics <- list(
    prior_distributions = prior_distributions,
    names = sapply(prior_distributions, function(prior){prior$dist_name}),
    hyper_parms = lapply(prior_distributions, function(prior){prior$hyper_parms}),
    parms_name = sapply(prior_distributions, function(prior){prior$parm_name})
  )
  class(characteristics) <- "independent_priors"
  return(characteristics)
}

#' sample_independent_priors generate value from the prior distribution
#'
#' @param independent_priors object of the class independent priors
#' @param R scalar Sample size
#'
#' @return A sample from the prior distribution
#' @export
#'
#' @examples
#' alpha_prior <- prior_dist("uniform", "alpha", c(0,10))
#' beta_prior <- prior_dist("constant", "beta", c(1))
#' model_prior <- independent_priors(list(alpha_prior, beta_prior))
#' sample_independent_priors(model_prior, 1)
#' sample_independent_priors(model_prior, 10)
sample_independent_priors <- function(model_prior, R){
  n_priors <- length(model_prior$names)
  if(R==1){
    mat <- t(as.matrix(
      sapply(1:length(model_prior$names), function(k)
        sample_prior(model_prior$prior_distributions[[k]], R)),
      ncol = n_priors))
    colnames(mat) <- model_prior$parms_name
  }else{
    mat <- sapply(1:length(model_prior$names), function(k)
        sample_prior(model_prior$prior_distributions[[k]], R))
    colnames(mat) <- model_prior$parms_name
  }
  mat
}


#' pdf_independent_priors compute the density of the prior distribution
#'
#' @param independent_priors object of the class independent priors
#' @param x matrix, points at which the log probabilities are evaluated
#'
#' @return value of the pdf of a set of parameter from the prior distribution viewpoint
#' @export
#'
#' @examples
#' alpha_prior <- prior_dist("uniform", "alpha", c(0,10))
#' beta_prior <- prior_dist("constant", "beta", c(1))
#' model_prior <- independent_priors(list(alpha_prior, beta_prior))
#' cloud <- sample_independent_priors(model_prior, 10)
#' logd_independent_priors(model_prior, cloud)
logd_independent_priors <- function(model_prior, cloud){
  n_priors <- length(model_prior$names)
  if(nrow(cloud)==1){
    sum(sapply(1:length(model_prior$names), function(k)
      logd_prior(model_prior$prior_distributions[[k]], cloud[,k])))
  }else{
    rowSums(sapply(1:length(model_prior$names), function(k)
      logd_prior(model_prior$prior_distributions[[k]], cloud[,k])))
  }

}


