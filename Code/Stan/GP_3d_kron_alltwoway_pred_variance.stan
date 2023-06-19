functions {
#include GP_helpers_pred.stan
}
data {
  //train
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
  // test
  int<lower=1> N1_new;
  int<lower=1> N2_new;
  int<lower=1> N3_new;
  matrix[N1_new,D1] X1_new;
  matrix[N2_new,D2] X2_new;
  matrix[N3_new,D3] X3_new;

  // parameters
  real<lower = 0> alpha0;
  real<lower = 0> alpha1;
  real<lower = 0> alpha2;
  real<lower = 0> alpha3;
  real<lower = 0> sigma;
}
transformed data {
  int<lower=1> N = N1 * N2 * N3  ;
  int<lower=1> N_new = N1_new * N2_new * N3_new ;
  //int <lower=0> n_term = 1 + 3 + 3 + 1;
  //vector[N_new] mu_fnew[n_term];
  vector[N_new] var_fnew;
  //real mloglik;
  {
    vector[N] eval;
    //vector[N] m_div_eval ;
    matrix[N1, N1] Q1;
    matrix[N2, N2] Q2;
    matrix[N3, N3] Q3;
    {
      vector[N] m ;
      matrix[N1, N1] K1;
      matrix[N2, N2] K2;
      matrix[N3, N3] K3;
      vector[N1] Browsum1;
      vector[N2] Browsum2;
      vector[N3] Browsum3;
      vector[N1] l1;
      vector[N2] l2;
      vector[N3] l3;
      {
        matrix[N1, (2*N1+1)]  M = sq_cen_kernel_mat(X1, N1, Hurst1);
        matrix[N1, N1] Ksq = M[1:N1, (N1+1):(2*N1)];
        K1  = M[1:N1, 1:N1];
        Browsum1 = to_vector(M[1:N1, (2*N1+1)]);
        l1 = eigenvalues_sym(Ksq);
        {
          int k = n_zero_eval(l1,N1);
          if (k > 0){
            Q1 = eigenvectors_sym(Ksq);
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
        matrix[N2, (2*N2+1)]  M = sq_cen_kernel_mat(X2, N2, Hurst2);
        matrix[N2, N2] Ksq = M[1:N2, (N2+1):(2*N2)];
        K2  = M[1:N2, 1:N2];
        Browsum2 = to_vector(M[1:N2, (2*N2+1)]);
        l2 = eigenvalues_sym(Ksq);
        {
          int k = n_zero_eval(l2,N2);
          if (k > 0){
            Q2 = eigenvectors_sym(Ksq);
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
        matrix[N3, (2*N3+1)]  M = sq_cen_kernel_mat(X3, N3, Hurst3);
        matrix[N3, N3] Ksq = M[1:N3, (N3+1):(2*N3)];
        K3  = M[1:N3, 1:N3];
        Browsum3 = to_vector(M[1:N3, (2*N3+1)]);
        l3 = eigenvalues_sym(Ksq);
        {
          int k = n_zero_eval(l3,N3);
          if (k > 0){
            Q3 = eigenvectors_sym(Ksq);
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
          vector[N] t0 = to_vector(to_vector(d3*d2')*d1');
          vector[N] t1 = to_vector(to_vector(d3*d2')*e1');
          vector[N] t2 = to_vector(to_vector(d3*e2')*d1');
          vector[N] t3 = to_vector(to_vector(e3*d2')*d1');
          vector[N] t12 = to_vector(to_vector(d3*e2')*e1');
          vector[N] t13 = to_vector(to_vector(e3*d2')*e1');
          vector[N] t23 = to_vector(to_vector(e3*e2')*d1');
          eval = square(alpha0)*(t0 + t1 + t2 + t3 + t12 + t13 + t23) + square(sigma)*rep_vector(1,N);
         }
      }

      //mean and variance
      int pos[2] ;
      int ix ;
      real v_star[7];
      vector[N] k_vec[7];
      vector[N1] one1 =  rep_vector(1,N1);
      vector[N2] one2 =  rep_vector(1,N2);
      vector[N3] one3 =  rep_vector(1,N3);
      v_star[1] = 1;
      k_vec[1] =  rep_vector(1,N);
      for (n1 in 1:N1_new){
        real v_all;
        vector[N] z;
        vector[N1] k1;
        {
          vector[N1+1] ktmp = sq_cen_kernel_vec(K1, Browsum1, X1, X1_new[n1,], N1, Hurst1, alpha1);
          k1 = ktmp[1:N1] ;
          v_star[1+1]= ktmp[N1+1];
          k_vec[1+1] = to_vector(to_vector(one3*one2')*k1');
        }
        pos[1] = (n1-1)*(N3_new*N2_new) ;
        for (n2 in 1:N2_new){
          vector[N2] k2;
          {
            vector[N2+1] ktmp = sq_cen_kernel_vec(K2, Browsum2, X2, X2_new[n2,], N2, Hurst2, alpha2);
            k2 = ktmp[1:N2];
            v_star[1+2] = ktmp[N2+1];
            k_vec[1+2] = to_vector(to_vector(one3*k2')*one1');
            //v12 = inteaction between 1 and 2
            v_star[3+1+1] = v_star[1+1] * v_star[2+1];
            k_vec[3+1+1] = to_vector(to_vector(one3*k2')*k1');
          }
          pos[2] = (n2-1)*(N3_new);
          for (n3 in 1:N3_new){
            vector[N3] k3;
            {
              vector[N3+1] ktmp = sq_cen_kernel_vec(K3, Browsum3, X3, X3_new[n3,], N3, Hurst3, alpha3);
              k3 = ktmp[1:N3];
              v_star[1+3] = ktmp[N3+1];
              k_vec[1+3] = to_vector(to_vector(k3*one2')*one1');
              //v13 =inteaction between 1 and 3
              v_star[3+1+2] = v_star[1+1] * v_star[1+3];
              k_vec[3+1+2] = to_vector(to_vector(k3*one2')*k1');
              // interaction between 2 and 3
              v_star[3+1+3] = v_star[1+2] * v_star[1+3];
              k_vec[3+1+3] = to_vector(to_vector(k3*k2')*one1');
            }
            ix = pos[1] + pos[2] + n3 ;
            v_all = square(alpha0)*(v_star[1] + v_star[2] + v_star[3] + v_star[4] +  v_star[5] +  v_star[6] + v_star[7]);
            {
              vector[N] k_vec_all = square(alpha0)*(k_vec[1] +  k_vec[2] + k_vec[3] + k_vec[4] + k_vec[5] + k_vec[6] + k_vec[7]);
              // computing Q^\top*k_vec ;
              vector[N] s = k_vec_all;
              s = mat_vec_prod(Q3', s, N3, NN3);
              s = mat_vec_prod(Q2', s, N2, NN2);
              s = mat_vec_prod(Q1', s, N1, NN1);
              z = s;
            }
            var_fnew[ix] =  v_all - sum(square(z)./eval);
          }
        }
      }
    }
  }
}
model {
}
generated quantities {
  vector[N_new] v_new = to_vector(var_fnew);
}
