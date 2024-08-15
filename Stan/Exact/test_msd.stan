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
          real lambda = pow(A[i] + A[j], nu + 1);
          marginal_sds += square(v) / lambda;
        }
      }
      return sqrt(marginal_sds);
  }
}

data {
  int dim;
  int nu;
  real rho;
}

parameters {
}

model {
}

generated quantities {
  matrix[dim, dim] Q1 = ar1_precision(dim, rho);
  tuple(matrix[dim, dim], vector[dim]) E = eigendecompose_sym(Q1);
  vector[dim * dim] msd = marginal_sd(E.2, E.1, dim, nu);
}
