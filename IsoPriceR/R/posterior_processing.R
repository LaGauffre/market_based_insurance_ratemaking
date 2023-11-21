# posterior_processing.r
# auteur: Pierre-O Goffard
# date: 21/09/2023
# Plots for the posterior distribution

#' posterior_plots plot the posterior distribution of each parameter
#'
#' @param cloud matrix, The cloud of particles
#' @param ref_line boolean, do we have reference value for the parameter
#' @param ref_matrix matrix, reference value for the parameter
#' @param font_size integer, size of the font
#'
#' @return a grid of plots with the posterior density of each parameter
#' @export
#'
#' @examples
posterior_plots <- function(cloud, model_prior, l_m, ref_line = TRUE, ref_matrix = true_parms, font_size = 20 ){
  popSize <- nrow(cloud)
  df_post <- data.frame(parm = as.vector(cloud),
                        parm_name = as.vector(
                          sapply(
                            colnames(cloud), function(parm_name) rep(parm_name, popSize))
                        )
  )
  random_parms <- l_m$parm_names[sapply(1:length(model_prior$prior_distributions),
                                        function(k) model_prior$prior_distributions[[k]]$dist_name != "constant")]
  df_post <- dplyr::filter(df_post, parm_name %in% random_parms)

  g <- ggplot2::ggplot(df_post, ggplot2::aes(x=parm))+
    ggplot2::geom_density(fill="gray") +
    ggplot2::facet_wrap(. ~ as.factor(parm_name), nrow = 2, scales='free') + ggplot2::labs(x= "") +
    ggplot2::theme_classic(base_size = font_size) + ggplot2::theme(legend.position="none")

  if(ref_line){
    df_ref <- data.frame(parm = true_parms[1, ],
                         parm_name = colnames(true_parms))
    df_ref <- dplyr::filter(df_ref, parm_name %in% random_parms)
    g <- g +   ggplot2::geom_vline(data=df_ref, ggplot2::aes(xintercept=parm, color="red"),
                                   linetype="dashed")
  }
  g

}

#' last_cloud_extract extract the last generation of particles
#'
#' @param res_abc list, output of the abc function
#'
#' @return matrix, last generation of particles
#' @export
#'
#' @examples
extract_last_cloud <- function(res_abc){
  res_abc$Clouds[[length(res_abc$Clouds)]]
}

extract_mode <- function(res_abc){
  cloud <- res_abc$Clouds[[length(res_abc$Clouds)]]
  ks::kde
}

#' best_cloud_extract extract the best particles over all the generation
#'
#' @param res_abc list, output of the abc function
#'
#' @return matrix, best particles over all the generations
#' @export
#'
#' @examples
extract_best_cloud <- function(res_abc){
  popSize <- nrow(res_abc$Clouds[[1]])
  all_clouds <- do.call(rbind, res_abc$Clouds); all_ds <- do.call(c, res_abc$ds)
  return(all_clouds[order(all_ds), ][1:popSize,])
}

#' best_particle_extract extract the best particle over all the generations
#'
#' @param res_abc list, output of the abc function
#'
#' @return vector, best particle over all the generations
#' @export
#'
#' @examples
extract_best_particle <- function(res_abc){
  all_clouds <- do.call(rbind, res_abc$Clouds); all_ds <- do.call(c, res_abc$ds)
  return(all_clouds[all_ds == min(all_ds), ])
}

#' map_extract extract mean a posteriori
#'
#' @param res_abc list, output of the abc function
#'
#' @return vector, mean a posteriori
#' @export
#'
#' @examples
extract_map <- function(res_abc){
  return(colMeans(res_abc$Clouds[[length(res_abc$Clouds)]]))
}

#' extract_mode extrcat the mode of the posterior distribution
#'
#' @param res_abc object returned by the abc function
#' @param l_m object of the type loss model
#' @param model_prior object of th etype prior distribution
#'
#' @return vector, mode a posteriori
#' @export
#'
#' @examples
extract_mode <- function(res_abc, l_m, model_prior){
  last_cloud <- extract_last_cloud(res_abc)
  random_parms <- l_m$parm_names[sapply(1:length(model_prior$prior_distributions),
                                        function(i) model_prior$prior_distributions[[i]]$dist_name != "constant")]
  map <- colMeans(last_cloud[,random_parms])
  f_hat <- function(x) -log(ks::dmvnorm.mixt(x,
                                             mus = last_cloud[,random_parms],
                                             Sigmas=matrix(rep(2 * cov(last_cloud[,random_parms]), nrow(last_cloud)), ncol = 2, byrow = TRUE)))
  mode <- colMeans(last_cloud)
  mode[random_parms] <- optim(map, f_hat)$par
  return(mode)
}


#' compute_metric Provide an estimation of the pure premium, the probability of no claims,
#'  the average Loss ratio, claim frequency, claim sizes, aggregated claim sizes
#'
#' @param data
#' @param p_p
#' @param l_m
#' @param particle
#' @param method
#'
#' @return
#' @export
#'
#' @examples
# particle <-true_parms
# particle <- true_parms

compute_metric_particle <- function(data, p_p, l_m, particle, method = "none"){
  ths_premium <- as.matrix(data[p_p$parm_names]); premium_type <- p_p$premium_type
  X <- sample_X(l_m, particle, R = 10000)
  freqs_sevs <- sample_freq_sev(l_m, particle, R = 10000)

  pp <- as.vector(compute_pure_premia(X, coverage_type = "xs", ths_premium))
  estimates <-data.frame(model = rep(l_m$model_name, length(particle) + 5),
                         method = rep(method, length(particle) + 5),
                         parm_names = c(colnames(particle), "Prob_0", "LR", "Claim_Frequency", "Claim_Amount", "Total_Claim_Amount"),
                         parm_values = c(particle,
                                         mean(X==0),
                                         mean(pp / data$x),
                                         mean(freqs_sevs$freqs),
                                         mean(freqs_sevs$sevs),
                                         mean(X)
                                         )
                         )
  res <- list(estimates, pp)
  names(res) <- c("estimates", "pp")
  return(res)
}



#' compute_metric_cloud Compute the pure premium and the metrics for a whole cloud
#'
#' @param data dataframe, commercial premiums and premium parameters
#' @param p_p object of type pure premium
#' @param l_m object of the type loss model
#' @param cloud matrix, cloud of particles
#'
#' @return a list that containss the pure premium for each particle of the cloud
#' and an estimation of the quantity of interest
#' @export
#'
#' @examples
# k <- 9
# l_m <- loss_models[[k]]
# l_m$model_name
# cloud <- extract_last_cloud(res_abc_list[[k]])
# freqs_sevs <- sample_freq_sev(l_m, cloud, R = 10000)
# colMeans(freqs_sevs$freqs)
compute_metric_cloud <- function(data, p_p, l_m, cloud){
  popSize <- nrow(cloud)
  ths_premium <- as.matrix(data[p_p$parm_names]); premium_type <- p_p$premium_type
  X <- sample_X(l_m, cloud, R = 10000)
  freqs_sevs <- sample_freq_sev(l_m, cloud, R = 10000)
  pp <- compute_pure_premia(X, coverage_type = "xs", ths_premium)
  df <- data.frame(
          model = rep(l_m$model_name, popSize),
          Claim_Frequency = colMeans(freqs_sevs$freqs),
          Claim_Amount = colMeans(freqs_sevs$sevs),
          Prob_0 = colMeans(X == 0),
          Total_Claim_Amount = colMeans(X),
          LR = sapply(1:popSize, function(k) mean(pp[,k] / data$x))
          )

  res <- list(pp, df)
  names(res) <- c("pps", "estimates")
  return(res)
}
