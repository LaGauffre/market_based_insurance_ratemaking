# titre: loss_model.r
# auteur: Pierre-O Goffard
# date: 21/09/2023
# Functions to set up the loss models

#' model Set up a distribution for the frequency or severity for instance
#'
#' @param model_name string, name of the frequency distribution
#' @param parm_names vector, name of the parameters of the frequency distribution
#'
#' @return an object of the type model
#' @export
#'
#' @examples
model <- function(dist_name, parm_names){
  characteristics <- list(dist_name = dist_name,
                          parm_names = parm_names
                          )
  class(characteristics) <- "model"
  return(characteristics)
}

#' loss_model
#'
#' @param freq_dist object of the type model
#' @param sev_dist object of the type model
#'
#' @return an object of type loss model
#' @export
#'
#' @examples
#' freq_dist <- model("none", c())
#' sev_dist <- model("gamma", c("alpha", "beta"))
#' l_m <- loss_model(freq_dist, sev_dist)
#' l_m$sev_dist$parm_names
loss_model <- function(freq_dist, sev_dist, model_name){
  characteristics <- list(freq_dist = freq_dist,
                          sev_dist = sev_dist,
                          model_name = model_name ,
                          parm_names = c(freq_dist$parm_names, sev_dist$parm_names)
  )
  class(characteristics) <- "loss_model"
  return(characteristics)
}

#' sample_X generates sample from the risk X for each set of parameters
#'
#' @param l_m object of type loss_model
#' @param cloud matrix, cloud of particles (stes of parameters for the loss
#' model)
#' @param R scalar number of MC replications
#'
#' @return the value of the pure premium
#' @export
#'
#' @examples
sample_X <- function(l_m, cloud, R = 3){
  freq = l_m$freq_dist$dist_name; sev = l_m$sev_dist$dist_name
  if(freq == "binomial"){
    d_freq <- 2
    th_freq <- cloud[, 1:d_freq, drop = FALSE]
    N <- sapply(1:nrow(cloud), function(k) rbinom(n = R, round(th_freq[k, 1]), th_freq[k, 2]))
  }else if(freq == "negative binomial"){
    d_freq <- 2
    th_freq <- cloud[, 1:d_freq,  drop = FALSE]
    N <- sapply(1:nrow(cloud), function(k) rnbinom(n = R, th_freq[k, 1], th_freq[k, 2]))
  }else if(freq == "poisson"){
    d_freq <- 1
    N <- sapply(1:nrow(cloud), function(k) rpois(n = R, cloud[k, 1:d_freq, drop = FALSE]))
  }else if(freq == "none"){
    d_freq <- 0
  }

  if(freq == "none"){
    if(sev == "lognormal"){
      d_freq <- 0
      th_sev <- cloud[, (d_freq + 1): (d_freq + 2), drop = FALSE]
      X <- sapply(1:nrow(cloud), function(k) rlnorm(R, th_sev[k, 1], th_sev[k, 2]))
    }else if(sev == "gamma"){
      th_sev <- cloud[, (d_freq + 1): (d_freq + 2), drop = FALSE]
      X <- sapply(1:nrow(cloud), function(k) rgamma(R, th_sev[k, 1], th_sev[k, 2]))
    }
  }else{
    if(R == 1){
      N <- t(N)
    }
    if(sev == "lognormal"){
      th_sev <- cloud[, (d_freq + 1): (d_freq + 2), drop = FALSE]
      X <- sapply(1:nrow(cloud), function(k)
        sapply(1:nrow(N), function(l)
          sum(rlnorm(N[l, k], th_sev[k, 1], th_sev[k, 2]))
        )
      )

    }else if(sev == "gamma"){
      th_sev <- cloud[, (d_freq + 1): (d_freq + 2),  drop = FALSE]
      X <- sapply(1:nrow(cloud), function(k)
        sapply(1:nrow(N), function(l)
          sum(rgamma(N[l, k], th_sev[k, 1], th_sev[k, 2]))
        )
      )
    }
  }
  X
}

#' sample_freq_sev generates sample of frequency and severity for each set of parameters
#'
#' @param l_m object of type loss_model
#' @param cloud matrix, cloud of particles (stes of parameters for the loss
#' model)
#' @param R scalar number of MC replications
#'
#' @return the value of the pure premium
#' @export
#'
#' @examples
sample_freq_sev <- function(l_m, cloud, R = 3){
  freq = l_m$freq_dist$dist_name; sev = l_m$sev_dist$dist_name
  if(freq == "binomial"){
    d_freq <- 2
    th_freq <- cloud[, 1:d_freq, drop = FALSE]
    N <- sapply(1:nrow(cloud), function(k) rbinom(n = R, round(th_freq[k, 1]), th_freq[k, 2]))
  }else if(freq == "negative binomial"){
    d_freq <- 2
    th_freq <- cloud[, 1:d_freq,  drop = FALSE]
    N <- sapply(1:nrow(cloud), function(k) rnbinom(n = R, th_freq[k, 1], th_freq[k, 2]))
  }else if(freq == "poisson"){
    d_freq <- 1
    N <- sapply(1:nrow(cloud), function(k) rpois(n = R, cloud[k, 1:d_freq, drop = FALSE]))
  }else if(freq == "none"){
    d_freq <- 0
  }
  if(sev == "lognormal"){
    th_sev <- cloud[, (d_freq + 1): (d_freq + 2), drop = FALSE]
    U <- sapply(1:nrow(cloud), function(k) rlnorm(R, th_sev[k, 1], th_sev[k, 2]))
    }else if(sev == "gamma"){
    th_sev <- cloud[, (d_freq + 1): (d_freq + 2), drop = FALSE]
    U <- sapply(1:nrow(cloud), function(k) rgamma(R, th_sev[k, 1], th_sev[k, 2]))
    }
  res <- list(N, U)
  names(res) <- c("freqs", "sevs")
  return(res)
}




