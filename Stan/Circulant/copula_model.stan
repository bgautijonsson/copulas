functions {
  matrix create_base_matrix(int dim, real rho) {
    matrix[dim, dim] c = rep_matrix(0, dim, dim);
    real scale = 1.0 / (1.0 + square(rho));
  
    // Set the first row
    c[1, 1] = (2 + 2 * square(rho)) * scale;
    c[1, 2] = -rho * scale;
    c[1, dim] = -rho * scale;

    // Set the second and last row
    c[2, 1] = -rho * scale;
    c[dim, 1] = -rho * scale;

    return c;
  }
  
  complex_matrix compute_and_rescale_eigenvalues(matrix c, int nu) {
    int dim = rows(c);
    complex_matrix[dim, dim] eigs = fft2(c);
    
    // Apply nu and compute inverse
    complex_matrix[dim, dim] inv_eigs = rep_matrix(0, dim, dim);
    for (i in 1:dim) {
      for (j in 1:dim) {
        eigs[i, j] = pow(eigs[i, j], nu + 1);
        inv_eigs[i, j] = 1.0 / eigs[i, j];
      }
    }
    
    // Compute marginal variance
    real mvar = get_real(inv_fft2(inv_eigs)[1, 1]);
    
    // Scale eigenvalues
    for (i in 1:dim) {
      for (j in 1:dim) {
        eigs[i, j] = mvar * eigs[i, j];
      }
    }
    
    return eigs;
  }
  
  vector matvec_prod(complex_matrix eigs, vector v) {
    int dim = rows(eigs);
    complex_matrix[dim, dim] v_mat = to_matrix(v, dim, dim);
    
    complex_matrix[dim, dim] fft_v = fft2(v_mat);
    complex_matrix[dim, dim] prod = eigs .* fft_v;
    complex_matrix[dim, dim] result_complex = inv_fft2(prod);
    
    return to_vector(get_real(result_complex));
  }

  /*
  This function computes the CDF of the Generalized Extreme Value distribution.
  */
  real gev_cdf(real y, real mu, real sigma, real xi) {
  
    if (abs(xi) < 1e-10) {
      real z = (y - mu) / sigma;
      return exp(-exp(z));
    } else {
      real z = 1 + xi * (y - mu) / sigma;
      if (z > 0) {
        return exp(-pow(z, -1/xi));
      } else {
        reject("Found incompatible GEV parameter values");
      }
    }
  }
  /*
  This function computes the log-pdf of the Generalized Extreme Value distribution.
  */
  real gev_lpdf(real y, real mu, real sigma, real xi) {
    if (abs(xi) < 1e-10) {
      real z = (y - mu) / sigma;
      return -log(sigma) - z - exp(-z);
    } else {
      real z = 1 + xi * (y - mu) / sigma;
      if (z > 0) {
        return -log(sigma) - (1 + 1/xi) * log(z) - pow(z, -1/xi);
      } else {
        reject("Found incompatible GEV parameter values");
      }
    }
  }
}

data {
  int<lower=1> dim;
  int<lower=0> nu;  
  int<lower=1> n_obs;
  matrix[dim*dim, n_obs] y;
}

transformed data {
  int D = dim * dim;
  real min_y = min(y);
}

parameters {
  real<lower=0, upper=1> rho;
  real<lower = 0> sigma;
  real<lower = 0> xi;
  real<lower = 0, upper = min_y + sigma/xi> mu;
}

model {
  matrix[dim, dim] c = create_base_matrix(dim, rho);
  complex_matrix[dim, dim] eigs = compute_and_rescale_eigenvalues(c, nu);
  matrix[D, n_obs] u;

  for (i in 1:D) {
    for (j in 1:n_obs) {
      u[i, j] = gev_cdf(y[i, j] | mu, sigma, xi);
      target+= gev_lpdf(y[i, j] | mu, sigma, xi);
    }
  }
  
  // Likelihood
  matrix[D, n_obs] Z = inv_Phi(u);
  real log_det = sum(log(get_real(eigs)));
  for (i in 1:n_obs) {
    vector[D] Qx = matvec_prod(eigs, Z[,i]);
    target += -0.5 * (dot_product(Z[,i], Qx) - log_det - sum(square(Z[,i])));
  }

  // Priors
  target += beta_lpdf(rho | 1, 1);
  target += std_normal_lpdf(xi);
  target += exponential_lpdf(sigma | 1);
}

