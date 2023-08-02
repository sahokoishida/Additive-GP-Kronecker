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
  real f_new0[2];
  vector[N1_new] f_new1[2] ;
  vector[N2_new] f_new2[2] ;
  vector[N3_new] f_new3[2];
  vector[N1_new*N2_new] f_new12[2] ;
  vector[N1_new*N3_new] f_new13[2] ;
  vector[N2_new*N3_new] f_new23[2] ;
  vector[N1_new*N2_new*N3_new] f_new123[2] ;
  real mloglik;
  {
    vector[N] eval;
    vector[N] m_div_eval ;
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
                Q1 = GS_complete(Q2, k, N2);
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
          vector[N] t123 = to_vector(to_vector(e3*e2')*e1');
          eval = square(alpha0)*(t0 + t1 + t2 + t3 + t12 + t13 + t23 + t123) + square(sigma)*rep_vector(1,N);
          m_div_eval = m ./ eval ;
         }
      }
      mloglik = -0.5 * sum(square(m)./eval)  -0.5 * sum(log(eval)) - 0.5*N*log(2*pi());

      //mean and variance
      vector[N1] Qv1 = Q1' * rep_vector(1,N1);
      vector[N2] Qv2 = Q2' * rep_vector(1,N2);
      vector[N3] Qv3 = Q3' * rep_vector(1,N3);
      real sq_alpha0 = square(alpha0);
      {
        vector[N] z0 = sq_alpha0 * to_vector(to_vector(Qv3*Qv2')*Qv1');
        f_new0[1] = z0' * m_div_eval;
        f_new0[2] = sq_alpha0 - sum(square(z0)./eval);
      }
      for (n1 in 1:N1_new){
        vector[N1] Qk1;
        real v1;
        {
          vector[N] z1;
          vector[N1+1] ktmp = sq_cen_kernel_vec(K1, Browsum1, X1, X1_new[n1,], N1, Hurst1, alpha1);
          Qk1 = Q1' * ktmp[1:N1] ;
          v1 = ktmp[N1+1];
          z1 = sq_alpha0*to_vector(to_vector(Qv3*Qv2')*Qk1');
          f_new1[1, n1] = z1' * m_div_eval;
          f_new1[2, n1] = sq_alpha0*v1 - sum(square(z1)./eval);
        }
        //pos[1] = (n1-1)*(N3_new*N2_new) ;
        for (n2 in 1:N2_new){
          vector[N2] Qk2;
          real v2;
          real v12;
          {
            vector[N] z2 ;
            vector[N] z12;
            vector[N2+1] ktmp = sq_cen_kernel_vec(K2, Browsum2, X2, X2_new[n2,], N2, Hurst2, alpha2);
            Qk2 = Q2' * ktmp[1:N2];
            z2 = sq_alpha0*to_vector(to_vector(Qv3*Qk2')*Qv1');
            z12 = sq_alpha0*to_vector(to_vector(Qv3*Qk2')*Qk1');
            v2= ktmp[N2+1];
            v12 = v1 * v2 ; //v12 = inteaction between 1 and 2
            f_new2[1, n2] = z2' * m_div_eval;
            f_new2[2, n2] = sq_alpha0*v2 - sum(square(z2)./eval);
            f_new12[1, (n1-1)*N2_new + n2] = z12' * m_div_eval;
            f_new12[2, (n1-1)*N2_new + n2] = sq_alpha0*v12 - sum(square(z12)./eval);
          }
          //pos[2] = (n2-1)*(N3_new);
          for (n3 in 1:N3_new){
            vector[N3] Qk3;
            real v3;
            real v13;
            real v23;
            real v123;
            vector[N] z3;
            vector[N] z13;
            vector[N] z23;
            vector[N] z123;
            vector[N3+1] ktmp = sq_cen_kernel_vec(K3, Browsum3, X3, X3_new[n3,], N3, Hurst3, alpha3);
            Qk3 = Q3' * ktmp[1:N3];
            z3 = sq_alpha0*to_vector(to_vector(Qk3*Qv2')*Qv1');
            z13 = sq_alpha0*to_vector(to_vector(Qk3*Qv2')*Qk1');
            z23 = sq_alpha0*to_vector(to_vector(Qk3*Qk2')*Qv1');
            z123 = sq_alpha0*to_vector(to_vector(Qk3*Qk2')*Qk1');
            v3 = ktmp[N3+1];
            v13 = v1 * v3; //v13 = inteaction between 1 and 3
            v23 = v2 * v3; //v23 =inteaction between 2 and 3
            v123 = v1 * v2 * v3 ;
            f_new3[1, n3] = z3' * m_div_eval;
            f_new3[2, n3] = sq_alpha0*v3 - sum(square(z3)./eval);
            f_new13[1, (n1-1)*N3_new + n3] = z13' * m_div_eval;
            f_new13[2, (n1-1)*N3_new + n3] = sq_alpha0*v13 - sum(square(z13)./eval);
            f_new23[1, (n2-1)*N3_new + n3] = z23' * m_div_eval;
            f_new23[2, (n2-1)*N3_new + n3] = sq_alpha0*v23 - sum(square(z23)./eval);
            f_new123[1,(n1-1)*(N2_new * N3_new) + (n2-1)*N3_new + n3 ] = z123' * m_div_eval;
            f_new123[2,(n1-1)*(N2_new * N3_new) + (n2-1)*N3_new + n3 ] = sq_alpha0*v123 - sum(square(z123)./eval);
          }
        }
      }
    }
  }
}
model {
}
generated quantities {
  // mean
  real mu0 = f_new0[1];
  vector[N1_new] mu1 = f_new1[1];
  vector[N2_new] mu2 = f_new2[1];
  vector[N3_new] mu3 = f_new3[1];
  vector[N1_new*N2_new] mu12 = f_new12[1] ;
  vector[N1_new*N3_new] mu13 = f_new13[1] ;
  vector[N2_new*N3_new] mu23 = f_new23[1] ;
  vector[N_new] mu123 = f_new123[1] ;
  // variance
  real v0 = f_new0[2];
  vector[N1_new] v1 = f_new1[2];
  vector[N2_new] v2 = f_new2[2];
  vector[N3_new] v3 = f_new3[2];
  vector[N1_new*N2_new] v12 = f_new12[2] ;
  vector[N1_new*N3_new] v13 = f_new13[2] ;
  vector[N2_new*N3_new] v23 = f_new23[2] ;
  vector[N_new] v123 = f_new123[2] ;
  // marginal likelihood
  real mllik = mloglik;
}
