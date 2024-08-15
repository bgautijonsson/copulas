functions {
    /*
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
  This function generated a precision matrix for an AR1 process with a given correlation parameter rho.
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
  vector marginal_sd(vector A, matrix V, int dim, int nu) {
      vector[dim * dim] marginal_sds = rep_vector(0.0, dim * dim);
      for (i in 1:dim) {
        for (j in 1:dim) {
          vector[dim * dim] v = kronecker(V[, j], V[, i]);
          real lambda = nu == 0 ? A[i] + A[j] : pow(A[i] + A[j], nu + 1);
          marginal_sds += square(v) / lambda;
        }
      }
      return sqrt(marginal_sds);
  }

  real matern_lpdf(matrix X, int dim, real rho, int nu) {
    int N = cols(X);
    int D = dim * dim;
    real C = D * log(2 * pi());
    vector[N] log_densities;
    matrix[dim, dim] Q1 = ar1_precision(dim, rho);
    
    // Eigendecompositions
    tuple(matrix[dim, dim], vector[dim]) E = eigendecompose_sym(Q1);
    
    
    // Compute log-density
    {
      real log_det = 0;
      vector[N] quadform_sums = rep_vector(0.0, N);
      
      for (i in 1:dim) {
        for (j in 1:dim) {
          vector[D] v = kronecker(E.1[, j], E.1[, i]);
          
          real lambda = pow(E.2[i] + E.2[j], nu + 1);
          
          log_det += log(lambda);
          
          for (n in 1:N) {
            real u = dot_product(v, X[, n]);
            quadform_sums[n] += square(u) * lambda;
          }
        }
      }
      log_densities = -0.5 * (C - log_det + quadform_sums);
    }
    
    return sum(log_densities);
  }
}

data {
  int dim;
  int nu;
  int n_obs;
  matrix[dim * dim, n_obs] y;
}

parameters {
  real<lower = 0, upper = 1> rho;
}

model {
  target += matern_lpdf(y | dim, rho, nu);
}

generated quantities {
}
