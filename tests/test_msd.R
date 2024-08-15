library(cmdstanr)
library(stdmatern)
library(tidyverse)

dim <- 4
nu <- 0
rho <- 0.3

Q1 <- make_AR_prec_matrix(dim, rho)
E <- eigen(Q1)
A1 <- E$values
V1 <- E$vectors

msd <- marginal_sd_eigen(A1, V1, dim, A1, V1, dim, nu)

model <- cmdstan_model("Stan/test_msd.stan")


results <- model$sample(
  data = list(dim = dim, nu = nu, rho = rho),
  chains = 4, 
  parallel_chains = 4,
  fixed_param = TRUE
)

results$summary("msd") |> 
  select(variable, mean) |> 
  arrange(desc(mean)) |> 
  mutate(
    msd = sort(msd, decreasing = TRUE)
  ) 

A_hat <- results$summary("E:2") |> 
  pull(mean)

V_hat <- matrix(
  results$summary("E:1") |> pull(mean),
  nrow = dim, ncol = dim, byrow = FALSE
)

zapsmall(V_hat %*% diag(A_hat) %*% t(V_hat))

marginal_sd_eigen(A_hat, V_hat, dim, A_hat, V_hat, dim, nu) - msd

results$summary("msd") |> 
  select(variable, mean) |> 
  arrange(desc(mean)) |> 
  mutate(
    msd = sort(msd, decreasing = TRUE)
  )

results$summary("E:2") |> 
  select(variable, mean) |> 
  arrange(desc(mean)) |> 
  mutate(
    A = A1
  )

results$summary("E:1") |> 
  select(variable, mean) |> 
  mutate(
    V = as.numeric(V1),
    diff = mean - V
  ) |> 
  mutate_at(vars(-variable), round, digits = 5)


loglik <- function(pars) {
  -sum(dmatern(y, c(dim_x, dim_y), c(pars, pars), nu))
}





optim(
  c(0.5),
  loglik, 
  lower = c(0), 
  upper = c(1),
  method = "Brent"
)



kron <- function(v1, v2) {
  n1 <- length(v1)
  n2 <- length(v2)
  out <- matrix(
    0,
    nrow = n1,
    ncol = n2
  )

  for (i in 1:n2) {
    out[ , i] = v1[i] * v2;
  }

  as.numeric(out)
}

marg_msd <- function(A, V, dim, nu) {
  N <- dim * dim
  out <- numeric(N)
  for (i in seq_len(dim)) {
    for (j in seq_len(dim)) {
      v <- kron(V[ , j], V[ , i])
      lambda <- (A[i] + A[j])^(nu + 1)
      out <- out + v^2 / lambda
    }
  }
  sqrt(out)
}

dim <- 3
rho <- 0.5
nu <- 0
Q1 <- make_AR_prec_matrix(dim, rho)

E <- eigen(Q1)

marg_msd(E$values, E$vectors, dim, nu) |> sort()

marginal_sd_eigen(E$values, E$vectors, dim, E$values, E$vectors, dim, nu) |> sort()
