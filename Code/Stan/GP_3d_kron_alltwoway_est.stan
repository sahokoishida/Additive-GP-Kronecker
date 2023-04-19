functions {
#include GP_helpers.stan
}
data {
  int<lower=1> N1;
  int<lower=1> N2;
  int<lower=1> N3;
  int<lower=1> NN1;
  int<lower=1> NN2;
  int<lower=1> NN3;
  int<lower=1> D1;
  int<lower=1> D2;
  int<lower=1> D3;
  matrix[N1,D1] X1;
  matrix[N2,D2] X2;
  matrix[N3,D3] X3;
  vector[N1*N2*N3] y;
  real<lower=0,upper=1> Hurst1;
  real<lower=0,upper=1> Hurst2;
  real<lower=0,upper=1> Hurst3;
}
transformed data {
  int<lower =1> N = N1 * N2 * N3 ;
  vector[N1] l1;
  vector[N2] l2;
  vector[N3] l3;
  vector[N] m;
  {
    matrix[N1, N1] Q1;
    matrix[N2, N2] Q2;
    matrix[N3, N3] Q3;
    {
      matrix[N1, N1]  K1 = sq_cen_kernel_mat(X1, N1, Hurst1);
      l1 = eigenvalues_sym(K1);
      {
        int k = n_zero_eval(l1,N1);
        if (k > 0){
          Q1 = eigenvectors_sym(K1);
          l1 = eval_zero(l1, k, N1);
          if (k>1){
              Q1 = GS_complete(Q1, k, N1);
          } else {
            Q1[,1] = rep_vector(1/sqrt(N1),N1);
          }
        } else {
          reject("k must be positive; found k=", k);
        }
      }
    }
    {
      matrix[N2, N2]  K2 = sq_cen_kernel_mat(X2, N2, Hurst2);
      l2 = eigenvalues_sym(K2);
      {
        int k = n_zero_eval(l2,N2);
        if (k > 0){
          Q2 = eigenvectors_sym(K2);
          l2 = eval_zero(l2, k, N2);
          if (k>1){
              Q2 = GS_complete(Q2, k, N2);
          } else {
            Q2[,1] = rep_vector(1/sqrt(N2),N2);
          }
        } else {
          reject("k must be positive; found k=", k);
        }
      }
    }
    {
      matrix[N3, N3]  K3 = sq_cen_kernel_mat(X3, N3, Hurst3);
      l3 = eigenvalues_sym(K3);
      {
        int k = n_zero_eval(l3,N3);
        if (k > 0){
          Q3 = eigenvectors_sym(K3);
          l3 = eval_zero(l3, k, N3);
          if (k>1){
              Q3 = GS_complete(Q3, k, N3);
          } else {
            Q3[,1] = rep_vector(1/sqrt(N3),N3);
          }
        } else {
          reject("k must be positive; found k=", k);
        }
      }
    }
    {
      vector[N] s = y;
      s = mat_vec_prod(Q3', s, N3, NN3);
      s = mat_vec_prod(Q2', s, N2, NN2);
      s = mat_vec_prod(Q1', s, N1, NN1);
      m = s;
    }
  }
}
parameters {
  real<lower = 0> alpha0;
  real<lower = 0> alpha1;
  real<lower = 0> alpha2;
  real<lower = 0> alpha3;
  real<lower=0> sigma;
}
model {
  vector[N] eval;
  {
    vector[N1] e1 = square(alpha1) * l1;
    vector[N2] e2 = square(alpha2) * l2;
    vector[N3] e3 = square(alpha3) * l3;
    vector[N1] d1 = rep_vector(0,N1);
    vector[N2] d2 = rep_vector(0,N2);
    vector[N3] d3 = rep_vector(0,N3);
    d1[1] = N1;
    d2[1] = N2;
    d3[1] = N3;
    {
      vector[N] t0 = square(alpha0)*to_vector(to_vector(d3*d2')*d1');
      vector[N] t1 = to_vector(to_vector(d3*d2')*e1');
      vector[N] t2 = to_vector(to_vector(d3*e2')*d1');
      vector[N] t3 = to_vector(to_vector(e3*d2')*d1');
      vector[N] t12 = to_vector(to_vector(d3*e2')*e1');
      vector[N] t13 = to_vector(to_vector(e3*d2')*e1');
      vector[N] t23 = to_vector(to_vector(e3*e2')*d1');
      eval = t0 + t1 + t2 + t3 + t12 + t13 + t23 + square(sigma)*rep_vector(1,N);
     }
  }
  //prior
  target += std_normal_lpdf(alpha0);
  target += std_normal_lpdf(alpha1);
  target += std_normal_lpdf(alpha2);
  target += std_normal_lpdf(alpha3);
  target += std_normal_lpdf(sigma);
  //likelihood
  target += -0.5 * sum(square(m)./eval) - 0.5 * sum(log(eval));
}
