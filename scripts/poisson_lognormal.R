# poisson_lognormal.r
# Simulation study over the Poisson Lognormal model

# Libraries----
library(ggplot2)

library(dplyr)
library(rstudioapi)
library(IsoPriceR)
library(readr)
library(parallel)
setwd(dirname(getActiveDocumentContext()$path))
getwd()



# Data ----
pet_insurance_FR <- read_csv("../data/pet_insurance_FR.csv")
data_pet_insurance <- data.frame(carrier = pet_insurance_FR$Company, race = pet_insurance_FR$Race,
                                 r = as.numeric(sub("%", "", pet_insurance_FR$rate))/100,
                                 d = pet_insurance_FR$deductible,
                                 l = pet_insurance_FR$limit,
                                 x = pet_insurance_FR$Premium)

# Sample sizes----
ns <- c( 50, 100, 250, 500)


# Premium an loss model Settings ----
p_p <- pure_premium("xs", c("r", "d", "l"))
freq_dist <- model("poisson", c("lambda"))
sev_dist <- model("lognormal", c("mu", "sigma"))
l_m <- loss_model(freq_dist, sev_dist, "Lognormal(mu, sigma = 1)-Poisson(lambda)")

# True parameters of the loss model ----
lambda_true <- 0.58; mu_true <- 5.75; sigma_true <- 1
true_parms <- t(as.matrix(c(lambda_true, mu_true, sigma_true)))
colnames(true_parms) <- l_m$parm_names

# Prior assumptions----
prior_lambda <- prior_dist("uniform", "lambda",c(0, 10))
prior_mu <- prior_dist("uniform", "mu",c(-10, 10))
prior_sigma <- prior_dist("constant", "sigma", c(1))
prior_distributions <- list(prior_lambda, prior_mu, prior_sigma)
model_prior <- independent_priors(prior_distributions)

# ABC parameters----
simulation_abc <- function(k, popSize = 100,  MC_R = 1000, acc = 1, sp_bounds =c(0.4, 0.7)){
  set.seed(k)
  print(k)
  # Fake Data
  d <- sample(data_pet_insurance$d, max(ns), replace = TRUE)
  l <- sample(data_pet_insurance$l, max(ns), replace = TRUE)
  r <- sample(data_pet_insurance$r, max(ns), replace = TRUE)
  th_premium <- matrix(c(r, d, l), nrow = max(ns))
  data <- as.data.frame(th_premium); colnames(data) <- p_p$parm_names
  X_true <- sample_X(l_m, true_parms, R = 10000)
  data["pp"] <- compute_pure_premia(X_true, coverage_type = "xs", th_premium)
  lr <- runif(max(ns) ,0.4, 0.7)
  data["x"] <- data["pp"] / lr
  # ABC
  res <- data.frame()
  n <- ns[1]
  for(n in ns){

    sub_data <- data[1:n,]
    res_abc <- abc(sub_data, model_prior, l_m, p_p,
                   popSize,  MC_R, acc = 1, sp_bounds = c(0.4, 0.7))
    map <- t(as.matrix(colMeans(extract_last_cloud(res_abc))))
    mode <- t(as.matrix(extract_mode(res_abc, l_m, model_prior)))

    res <- rbind(res, compute_metric_particle(data, p_p, l_m, map, method = "map")$estimates,
                 compute_metric_particle(data, p_p, l_m, mode, method = "mode")$estimates)
  }

  res['ss'] <- as.vector(sapply(ns, function(n) rep(n, nrow(res)/length(ns))))
  return(res)

}


n.cores_max <- detectCores()
n.cores <- n.cores_max -1
Run <- 7


cl <- makeCluster(n.cores)
clusterExport(cl, c(unclass(lsf.str(envir = asNamespace("IsoPriceR"), all = T)),
                    'set.seed', 'simulation_abc','ns','data_pet_insurance', 'model_prior', 'l_m', 'p_p' , 'true_parms'))
system.time(out <- parLapply(cl, 1:Run, function(k) simulation_abc(k, popSize = 1000,  MC_R = 2000, acc = 1, sp_bounds =c(0.4, 0.7))))
# system.time(out <- lapply(1:4, function(k) simulation_abc(k, popSize = 100,  MC_R = 1000, acc = 1, sp_bounds =c(0.4, 0.7))))
res_sim <- do.call(rbind, out)
stopCluster(cl)

write.csv(res_sim, "../data/res_sim.csv")
print("done")
