# titre: premium.r
# auteur: Pierre-O Goffard
# date: 05/08/2023
# Functions to compute insurance premium

#' pure_premium creates an instance of pure premium
#'
#' @param coverage_type string, type of coverage like excess of loss ("xs")
#' @param parm_names vector, names of the parameters of the coverage type.
#' In the case pf excess of loss, we have rate, deductible and limit
#'
#' @return an instance of the class pure premium
#' @export
#'
#' @examples
#' pure_premium("xs", c("r", "d", "l"))
pure_premium <- function(coverage_type, parm_names){
  characteristics <- list(coverage_type = coverage_type,
                          parm_names = parm_names)
  class(characteristics) <- "pure_premium"
  return(characteristics)
}

#' compute_pure_premia_mat compute all the premiums at once using MC simulations for
#' matrix of insurance coverage parameter
#'
#' @param X_mat matrix samples of aggregate claims
#' @param coverage_type string Type of insurance coverage considered
#' @param th_premium matrix parameter of the premium
#'
#' @return the value of the pure premia
#' @export
#'
#' @examples
#' freq_dist <- model("negative binomial", c(1, 0.3))
# sev_dist <- model("gamma", c("alpha", "beta"))
# l_m <- loss_model(freq_dist, sev_dist)
# l_m$sev_dist$parm_names
# freq = l_m$freq_dist$dist_name; sev = l_m$sev_dist$dist_name
#
# N_prior <- prior_dist("constant", "N", c(1))
# p_prior <- prior_dist("uniform", "p", c(0, 1))
# alpha_prior <- prior_dist("uniform", "alpha", c(0,10))
# beta_prior <- prior_dist("constant", "beta", c(1))
# model_prior <- independent_priors(list(N_prior, p_prior, alpha_prior, beta_prior))
# cloud <- sample_independent_priors(model_prior, 1000)
#
# R <- 1000
# X <- sample_X(l_m, cloud, R)
# p_p <- pure_premium("xs", c("r", "d", "l"))
# th_premium <-  as.matrix(expand.grid( c(0.5, 0.75),  rep(30, 50),  c(1e3)))
# colnames(th_premium) <- p_p$parm_names
# coverage_type <- p_p$coverage_type
# compute_pure_premia(X, coverage_type , th_premium)
compute_pure_premia <- function(X,
                                coverage_type = "xs",
                                th_premium =
                                  as.matrix(expand.grid( c(0.75),  c(30, 50),  c(1e3)))){
  if(coverage_type == "xs"){
    if(ncol(X)>1){
      t(
        sapply(1:nrow(th_premium), function(k)
          colMeans(
            pmin(pmax(
              th_premium[k,1] * X - th_premium[k,2], 0), th_premium[k, 3]
              )
          )
        )
      )
    }else{
      as.matrix(sapply(1:nrow(th_premium), function(k)
        mean(pmin(pmax(th_premium[k,1] * X - th_premium[k,2], 0), th_premium[k, 3]))), ncol = 1)
    }
  }
}

sd_pure_premia <- function(X,
                                coverage_type = "xs",
                                th_premium =
                                  as.matrix(expand.grid( c(0.75),  c(30, 50),  c(1e3)))){
  if(coverage_type == "xs"){
    if(ncol(X)>1){
      t(
        sapply(1:nrow(th_premium), function(k)
          resample::colStdevs(
            pmin(pmax(
              th_premium[k,1] * X - th_premium[k,2], 0), th_premium[k, 3]
            )
          )
        )
      )
    }else{
      as.matrix(sapply(1:nrow(th_premium), function(k)
        sd(pmin(pmax(th_premium[k,1] * X - th_premium[k,2], 0), th_premium[k, 3]))), ncol = 1)
    }
  }
}




