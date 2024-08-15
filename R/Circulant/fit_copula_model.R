library(cmdstanr)
library(stdmatern)
library(tidyverse)
library(evd)
library(bayesplot)
dim <- 10
nu <- 0
rho <- 0.5
n_obs <- 5

z <- rmatern_copula_circulant(n_obs, dim, rho, nu)
u <- pnorm(z)
y <- qgev(u, loc = 6, scale = 2, shape = 0.1)

gev_fit <- fgev(y)$estimate
inits <- list(
  mu = gev_fit[1],
  sigma = gev_fit[2],
  xi = gev_fit[3],
  rho = 0.5
)

model <- cmdstan_model("Stan/Circulant/copula_model.stan")


results <- model$sample(
  data = list(
    dim = dim,
    n_obs = n_obs,
    y = y,
    nu = nu
  ),
  chains = 4,
  parallel_chains = 4,
  init = list(inits, inits, inits, inits)
)

results$summary()
