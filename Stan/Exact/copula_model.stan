functions {
  /*
  This function computes the Kronecker product of two vectors a and b.
  */
  vector kronecker(vector a, vector b) {
    int n_a = rows(a);
    int n_b = rows(b);
    vector[n_a * n_b] result;
    
    for (i in 1:n_a) {
      for (j in 1:n_b) {
        result[(i-1)*n_b + j] = a[i] * b[j];
      }
    }
    
    return result;
  }
  
  /*
  This function generates the precision matrix for an AR1 process with a given correlation parameter rho.
  */
  matrix ar1_precision(int n, real rho) {
    matrix[n, n] Q;
    real scaling = 1.0 / (1.0 - rho * rho);
    real off_diag = -rho * scaling;
  
    Q = rep_matrix(0, n, n);
    for (i in 1:n) {
      Q[i, i] = (i == 1 || i == n) ? scaling : (1.0 + rho * rho) * scaling;
      if (i > 1) Q[i, i-1] = off_diag;
      if (i < n) Q[i, i+1] = off_diag;
    }
    return Q;
  }
  
  /*
  This function takes as input the eigendecomposition of two matrices Q1 and Q2 and 
  returns the marginal standard deviations of the matrix Q, defined as the
  kronecker sum of Q1 and Q2 with smoothness parameter nu.
  */
  vector marginal_sd(vector A1, matrix V1, int dim1, vector A2, matrix V2, int dim2, real nu) {
      vector[dim1 * dim2] marginal_sds = rep_vector(0.0, dim1 * dim2);
      for (i in 1:dim1) {
        for (j in 1:dim2) {
          vector[dim1 * dim2] v = kronecker(V2[, j], V1[, i]);
          real lambda = nu == 0 ? A1[i] + A2[j] : pow(A1[i] + A2[j], nu + 1);
          marginal_sds += square(v) / lambda;
        }
      }
      return sqrt(marginal_sds);
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
  int dim1;
  int dim2;
  real nu;
  int n_obs;
  matrix[dim1 * dim2, n_obs] y;
}

transformed data {
  int D = dim1 * dim2;
  real min_y = min(y);
}

parameters {
  vector<lower = 0, upper = 1>[2] rho;
  real<lower = 0> sigma;
  real<lower = 0> xi;
  real<lower = 0, upper = min_y + sigma/xi> mu;
}

model {
  matrix[dim1, dim1] Q1 = ar1_precision(dim1, rho[1]);
  matrix[dim2, dim2] Q2 = ar1_precision(dim2, rho[2]);
  tuple(matrix[dim1, dim1], vector[dim1]) E1 = eigendecompose_sym(Q1);
  tuple(matrix[dim2, dim2], vector[dim2]) E2 = eigendecompose_sym(Q2);

  matrix[D, n_obs] u;

  for (i in 1:D) {
    for (j in 1:n_obs) {
      u[i, j] = gev_cdf(y[i, j] | mu, sigma, xi);
      target+= gev_lpdf(y[i, j] | mu, sigma, xi);
    }
  }

  vector[D] marginal_sds = marginal_sd(E1.2, E1.1, dim1, E2.2, E2.1, dim2, nu);
  matrix[D, n_obs] Z = inv_Phi(u);
  real log_det = 0;
  real quadform_sum = 0;

  for (i in 1:dim1) {
    for (j in 1:dim2) {
      vector[D] v = kronecker(E2.1[, j], E1.1[, i]);
      v = v .* marginal_sds;  
      real norm_v = sqrt(sum(square(v)));
      v /= norm_v;  
      
      real lambda = pow(E1.2[i] + E2.2[j], nu + 1) * square(norm_v);
      log_det += log(lambda);
      
      row_vector[n_obs] q = v' * Z;  
      quadform_sum += dot_self(q) * lambda;
    }
  }

  real z_squared = sum(columns_dot_self(Z));
  target += -0.5 * (quadform_sum - n_obs * log_det - z_squared);

  // Priors
  target += beta_lpdf(rho | 1, 1);
  target += std_normal_lpdf(xi);
  target += exponential_lpdf(sigma | 1);
  
}
