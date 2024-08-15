library(cmdstanr)
library(stdmatern)
library(tidyverse)
library(evd)
library(bayesplot)
dim1 <- 10
dim2 <- 10
nu <- 2
rho1 <- 0.6
rho2 <- 0.2
n_obs <- 20

z <- rmatern_copula(n_obs, c(dim1, dim2), c(rho1, rho2), nu)
u <- pnorm(z)
y <- qgev(u, loc = 6, scale = 2, shape = 0.1)

gev_fit <- fgev(y)$estimate
inits <- list(
  mu = gev_fit[1],
  sigma = gev_fit[2],
  xi = gev_fit[3],
  rho = c(0.5, 0.5)
)

model <- cmdstan_model("Stan/Exact/copula_model.stan")


results <- model$sample(
  data = list(
    dim1 = dim1,
    dim2 = dim2,
    n_obs = n_obs,
    y = y,
    nu = nu
  ),
  chains = 4,
  parallel_chains = 4,
  init = list(inits, inits, inits, inits)
)

results$summary()
