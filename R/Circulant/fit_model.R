library(cmdstanr)
library(stdmatern)
library(tidyverse)

dim <- 30
nu <- 2
rho <- 0.1
n_obs <- 10

y <- rmatern(n_obs, dim, rho, nu)

model <- cmdstan_model("Stan/model.stan")


results <- model$sample(
  data = list(dim = dim, nu = nu, n_obs = n_obs, y = y),
  chains = 4, 
  parallel_chains = 4
)


results$summary()
