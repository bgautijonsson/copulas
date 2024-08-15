library(cmdstanr)
library(stdmatern)
library(tidyverse)
library(evd)
library(bayesplot)
dim <- 40
nu <- 0
rho <- 0.8
n_obs <- 1

z <- rmatern_copula(n_obs, dim, rho, nu)
u <- pnorm(z)
y <- qgev(u, loc = 6, scale = 2, shape = 0.1)

gev_fit <- fgev(y)$estimate
inits <- list(
  mu = gev_fit[1],
  sigma = gev_fit[2],
  xi = pmax(0, gev_fit[3]),
  rho = 0.5
)

model <- cmdstan_model("Stan/Folded/copula_model.stan")


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
