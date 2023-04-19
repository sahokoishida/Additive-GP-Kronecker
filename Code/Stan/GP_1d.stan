functions {
#include GP_helpers.stan
}
data {
  int<lower=1> N1;
  int<lower=1> N2;
  vector[N1] t1;
  vector[N1] y1;
  vector[N2] t2;
  real<lower=0, upper=1> Hurst;
}
transformed data {
  real delta = 1e-9;
  vector[N1] l;
  matrix[N1, N1] Q;
  vector[N1] m;
  vector[N1] h;
  vector[N1] g[N2];
  {
    matrix[N1, N1] K_fbm ;
    matrix[N1, N1] K_sqfbm;
    vector[N1] Browsum ;
    real Bsum ;
    matrix[N1, N1] E;
    matrix[N1, N1] B = rep_matrix(0, N1, N1);
    matrix[N1, N1]tcp = t1 * t1' ;
    vector[N1] dvec = diagonal(tcp);
    vector[N1] d;
    matrix[N1, N1] A = diag_matrix(rep_vector(1,N1)) - (1.0/N1)*rep_matrix(1, N1, N1);
    for (i in 1:(N1-1)){
      d[i] = pow(fabs(dvec[i]), Hurst);
      K_fbm[i,i] = 1;
      for (j in (i+1):N1){
        real r = fabs(t1[i]-t1[j]);
        B[i,j] = pow(square(r), Hurst);
        B[j,i] = B[i,j];
      }
    }
    d[N1] = pow(fabs(dvec[N1]), Hurst);
    E = rep_matrix(d, N1);
    Browsum = B *rep_vector(1,N1);
    Bsum = sum(Browsum);
    K_fbm = 0.5 * (E + E' - B);
    K_fbm = A * K_fbm * A;
    K_sqfbm = K_fbm*K_fbm;
    l = eigenvalues_sym(K_sqfbm);
    {
      int k = n_zero_eval(l,N1);
      if (k > 0){
        Q = eigenvectors_sym(K_sqfbm);
        l = eval_zero(l, k, N1);
        if (k>1){
            Q = GS_complete(Q, k, N1);
        } else {
          Q[,1] = rep_vector(1/sqrt(N1),N1);
        }
      } else {
        reject("k must be positive; found k=", k);
      }
    }
    {
      for (i in 1:N2){
        vector[N1] r2vec;
        vector[N1] k_fbm;
        for (j in 1:N1){
          real r = fabs(t2[i]-t1[j]);
          r2vec[j] = pow(square(r), Hurst);
        }
        k_fbm = 0.5* (-r2vec + (sum(r2vec)/N1)*rep_vector(1,N1) + (1.0/N1)*Browsum - (1.0/square(N1))*Bsum*rep_vector(1,N1));
        k_fbm = K_fbm*k_fbm;
        g[i] = Q' * k_fbm;
      }
    }
    h = Q' * rep_vector(1, N1);
    m = Q' * y1 ;
  }
}
parameters {
  real<lower=0> sigma;
  real<lower=0> alpha;
  real<lower=0> alpha_0;
}
model{
  vector[N1] eval;
  {
    vector[N1] lambda;
    lambda = square(alpha) * l ;
    lambda[1] = square(alpha_0)* N1;
    eval = lambda  + square(sigma) * rep_vector(1,N1);
  }
  //prior
  target += std_normal_lpdf(alpha_0);
  target += std_normal_lpdf(alpha);
  target += std_normal_lpdf(sigma);
  //likelihood
  target += -0.5 * sum(square(m)./eval) -0.5 *sum (log(eval));
}
generated quantities {
  vector[N2] y2 ;
  for (i in 1:N2){
    vector[N1] hg = square(alpha_0) * h + square(alpha) * g[i];
    vector[N1] eval ;
    {
      vector[N1] lambda;
      lambda = square(alpha) * l ;
      lambda[1] = square(alpha_0)* N1;
      eval = lambda  + square(sigma) * rep_vector(1,N1);
    }
    y2[i] = sum((hg .* m)./eval);
  }
}
