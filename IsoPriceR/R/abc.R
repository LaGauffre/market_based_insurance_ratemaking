# # abc.r
# auteur: Pierre-O Goffard
# date: 21/09/2023
# Abc method


#' sample_first Initialize the ABC algorithm.
#'
#' @param data dataframe, contains the premium parameter and the value of the commercial premium
#' @param model_prior object of type independent_prior
#' @param l_m object oftype loss_model
#' @param p_p objecte of type pure_premium
#' @param popSize Size of the cloud of particles
#' @param MC_R NUmber of replication foer the M%onte Carlo estimation of tghe pure premiums
#'
#' @return A list containing the initial cloud,
#' the distance wrt the real data,
#' the indices of the particles to be retained to build the next cloud of particles
#' the tolerance level, the Effective Sample Size, the KDE estimator to represent the particles distribution
#' teh weights of each particles and the empirical variance-covariance matrix.
#' @export
#'
#' @examples
sample_first <- function(data, model_prior, l_m, p_p, popSize, MC_R, sp_bounds){
  # premium parameter as a matrix
  ths_premium <- as.matrix(data[p_p$parm_names]); premium_type <- p_p$premium_type
  cps <- data$x
  # which parameters are random (to be inferred)
  random_parms <- l_m$parm_names[sapply(1:length(model_prior$prior_distributions),
                                       function(k) model_prior$prior_distributions[[k]]$dist_name != "constant")]

  cloud <- sample_independent_priors(model_prior, R = popSize)
  w <- rep(1, nrow(cloud))
  if(length(random_parms) == 1){
    w.cov <- sd(cloud[, random_parms])

  }else{
    w.cov <- cov.wt(cloud, wt = w, cor = FALSE, center = TRUE)$cov
  }
  X <- sample_X(l_m , cloud, R = MC_R)
  pp_fake <- compute_pure_premia(X, coverage_type = "xs", ths_premium)
  ds1 <- sapply(1:ncol(pp_fake), function(k) sqrt(mean((pmax(cps - pp_fake[,k] /sp_bounds[1], 0) )^2)))

  ds2 <- sapply(1:ncol(pp_fake), function(k) sqrt(mean((pmax( pp_fake[,k] / sp_bounds[2] - cps , 0) )^2)))
  # colMeans(sapply(1:popSize, function(k) (pp_fake[, k] - pps)^2))

  ress_iso <- lapply(1:ncol(pp_fake), function(k) isoreg(pp_fake[, k], cps))

  fs <- lapply(ress_iso, function(res_iso) stepfun(sort(res_iso$x), c(min(res_iso$yf) /2 ,
                                                                      res_iso$yf), f = 1))
  ds <- sapply(1:nrow(cloud), function(k) mean((fs[[k]](pp_fake[, k]) - data$x)^2)^(1/2)) + ds1 + ds2
  sort_index <- order(ds)
  ESSs <- sapply(1:length(w), function(k) 1 / sum((w[sort_index][1:k] / sum(w[sort_index][1:k]))^2) )
  indices <- order(ds)[1:which((ESSs > popSize/2))[1]]
  # kde <- ks::kde(x = cloud[indices, random_parms],
  #                H =  2 * w.cov[random_parms,random_parms],
  #                eval.points = cloud[indices, random_parms],
  #                w = w[indices])
  # kde$estimate

  return(list(cloud = cloud, ds = ds,  indices = indices, eps = max(ds[indices]), ESS = sum(w[indices])^2/sum(w[indices]^2), w = w, w.cov = w.cov))

}

#' sample_g Sample generation g of the Population Monte carlo sampler in ABC
#'
#' @param data dataframe, contains the premium parameter and the value of the commercial premium
#' @param model_prior object of type independent_prior
#' @param l_m object oftype loss_model
#' @param p_p objecte of type pure_premium
#' @param popSize Size of the cloud of particles
#' @param sample_res list, result from sample_first or sample_g
#' @param MC_R Number of replication foer the Monte Carlo estimation of tghe pure premiums
#'
#' @return A list containing the initial cloud,
#' the distance wrt the real data,
#' the indices of the particles to be retained to build the next cloud of particles
#' the tolerance level, the Effective Sample Size, the KDE estimator to represent the particles distribution
#' the weights of each particles and the empirical variance-covariance matrix.
#'
#' @return
#' @export
#'
#' @examples
sample_g <- function(data, model_prior, l_m, p_p, popSize, sample_res, MC_R, sp_bounds){
  # premium parameter as a matrix
  ths_premium <- as.matrix(data[p_p$parm_names]); premium_type <- p_p$premium_type
  cps <- data$x
  # which parameters are random (to be infered)
  random_parms <- l_m$parm_names[sapply(1:length(model_prior$prior_distributions),
                                       function(k) model_prior$prior_distributions[[k]]$dist_name != "constant")]

  old_cloud <- sample_res$cloud[, l_m$parm_names]; indices <- sample_res$indices; w <- sample_res$w

  eps <- sample_res$eps ; w.cov <- sample_res$w.cov
  cloud <- matrix(nrow = 0, ncol = ncol(old_cloud)); ds <- c()
  pb = utils::txtProgressBar(min = 0, max = popSize, initial = 0)
  while(nrow(cloud) < popSize){

    selected_indices <- sample(indices, size = popSize - nrow(cloud), replace = T, prob = w[indices])
    if(length(random_parms) == 1){
      new_cloud <- old_cloud[selected_indices, , drop = FALSE] + rnorm(n = popSize - nrow(cloud),
                                                                                  mean = 0,
                                                                                  sd = 2 * w.cov)


    }else{
      new_cloud <- old_cloud[selected_indices, , drop = FALSE] + mvtnorm::rmvnorm(n = popSize - nrow(cloud),
                                                                                  mean = rep(0,length(model_prior$parms_name)),
                                                                                  sigma = 2 * w.cov)
    }

    # We blend the new particles with the old ones
    # new_cloud <- rbind(new_cloud[s, ], cloud[selected_indices, ][ !s,])
    # new_cloud <- new_cloud[(logd_independent_priors(model_prior, new_cloud) != - Inf) &
    #                          (predict(kde, x = new_cloud[, random_parms]) >0 ), , drop = FALSE]
    new_cloud <- new_cloud[logd_independent_priors(model_prior, new_cloud) != - Inf, , drop = FALSE]
    if(nrow(new_cloud) >0){
      X <- sample_X(l_m , new_cloud, R = MC_R)
      pp_fake <- compute_pure_premia(X, coverage_type = "xs", ths_premium)
      ds1 <- sapply(1:ncol(pp_fake), function(k) sqrt(mean((pmax(cps - pp_fake[,k] /sp_bounds[1], 0) )^2)))

      ds2 <- sapply(1:ncol(pp_fake), function(k) sqrt(mean((pmax( pp_fake[,k] / sp_bounds[2] - cps , 0) )^2)))

      ress_iso <- lapply(1:ncol(pp_fake), function(k) isoreg(pp_fake[, k], cps))

      fs <- lapply(ress_iso, function(res_iso) stepfun(sort(res_iso$x), c(min(res_iso$yf) /2 ,
                                                                          res_iso$yf), f = 1))
      ds_temp <- sapply(1:nrow(new_cloud), function(k) mean((fs[[k]](pp_fake[, k]) - data$x)^2)^(1/2)) + ds1 + ds2

      cloud <- rbind(cloud, new_cloud[ds_temp <= eps , ])
      ds <- c(ds, ds_temp[ds_temp <= eps ])

    }else{
      cloud <- cloud; ds <- ds
    }
    utils::setTxtProgressBar(pb,nrow(cloud))

  }
  close(pb)
  prior_prob <- exp(logd_independent_priors(model_prior, cloud))
  # kde <- ks::kde(x = old_cloud[indices, random_parms],
  #                H =  2 * w.cov[random_parms,random_parms],
  #                eval.points = cloud[, random_parms],
  #                w = w[indices])
  # inter_prob <- kde$estimate
  if(length(random_parms) == 1){
    inter_prob <- ks::kde(cloud[, random_parms],eval.points = cloud[, random_parms])$estimate


    w.cov <- sd(cloud[, random_parms])

  }else{
  inter_prob <- ks::dmvnorm.mixt(cloud[, random_parms],
               mus = old_cloud[indices, random_parms],
               Sigmas=matrix(rep(2 * w.cov[random_parms,random_parms], length(indices)), ncol = length(random_parms), byrow = TRUE),
               props = w[indices] / sum(w[indices]))
  w.cov <- cov.wt(cloud[, l_m$parm_names], wt = w, cor = FALSE, center = TRUE)$cov
  }




  w <- (prior_prob / inter_prob)
  # w[is.infinite(w)] <- 0
  w <- w / sum(w)

  # w.cov <- tryCatch(
  #   {
  #     cov.wt(cloud[, l_m$parm_names], wt = w, cor = FALSE, center = TRUE)$cov
  #   },
  #   error = function(e) {
  #     browser()
  #   }
  # )
  # w.cov <- cov.wt(cloud[, l_m$parm_names], wt = w, cor = FALSE, center = TRUE)$cov
  sort_index <- order(ds)
  ESSs <- sapply(1:length(w), function(k) 1 / sum((w[sort_index][1:k] / sum(w[sort_index][1:k]))^2) )
  if(sum(ESSs > popSize / 2) == 0){
    indices <- sort_index
  }else{
    indices <- sort_index[1:which(ESSs > popSize/2)[1]]
  }
  eps <- max(max(ds[indices]) - 10^(-8), eps_min)


  # kde <- suppressWarnings(
  #   ks::kde(x = cloud[indices, random_parms], H =  2 * w.cov[random_parms,random_parms],
  #           # eval.points = cloud[indices, random_parms] ,
  #           w = w[indices]))
  return(list(cloud = cloud, ds = ds, indices = indices, eps = eps, ESS = sum(w[indices])^2/sum(w[indices]^2), w = w, w.cov = w.cov))
}


#' abc ABC algorithm based on a Population Monte Carlo Sampler
#'
#' @param data dataframe, contains the premium parameter and the value of the commercial premium
#' @param model_prior object of type independent_prior
#' @param l_m object oftype loss_model
#' @param p_p objecte of type pure_premium
#' @param popSize Size of the cloud of particles
#' @param MC_R Number of replication foer the Monte Carlo estimation of tghe pure premiums
#' @param acc scalar, level of accuracy that is the difference between to consecutive tolerance level (epsilon)
#'
#' @return A list containing all the clouds and the distance to the real data for each particle.
#' @export
#'
#' @examples
abc <- function(data, model_prior, l_m, p_p, popSize, MC_R, acc = 0.01, sp_bounds, eps_min = Inf){
  g<- 1; Epss <- c(Inf); ESS <- c(popSize)
  print(paste("Sampling generation of particles ", g, "; eps = ", Epss[g], "; ESS = ", ESS[g]))
  sample_res <- sample_first(data, model_prior, l_m, p_p, popSize, MC_R, sp_bounds)
  Clouds <- list(sample_res$cloud); dss <- list(sample_res$ds)
  g <- 2
  Epss[g] <- sample_res$eps; ESS[g] <- sample_res$ESS
  print(paste("Sampling generation of particles ", g, "; eps = ", Epss[g], "; ESS = ", ESS[g]))
  sample_res <- sample_g(data, model_prior, l_m, p_p, popSize, sample_res, MC_R, sp_bounds)
  Clouds[[g]] <- sample_res$cloud; dss[[g]] <- sample_res$ds
  while( (Epss[g-1] - Epss[g] > acc) | (Epss[g] > eps_min)){
    g <- g + 1
    Epss[g] <- sample_res$eps; ESS[g] <- sample_res$ESS
    print(paste("Sampling generation of particles ", g, "; eps = ", Epss[g], "; ESS = ", ESS[g]))
    sample_res <- sample_g(data, model_prior, l_m, p_p, popSize, sample_res, MC_R, sp_bounds)
    Clouds[[g]] <- sample_res$cloud; dss[[g]] <- sample_res$ds
  }
  return(list(Clouds = Clouds, ds = dss))
}

