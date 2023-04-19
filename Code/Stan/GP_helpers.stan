int n_zero_eval(vector eval, int n){
  // count number of zero eigen values (smaller than threshold)
  // Input: eval - a vector of eigenvalues
  //        n ---- a size of a matrix
  // Output: number of zero eigen values
  real d = eval[1];
  int i = 1 ;
  while (d < 1e-9){
      i += 1 ;
      d = eval[i];
  }
  return i - 1 ;
}

vector eval_zero(vector eval, int k, int n){
  // replace eigenvalues with zero
  vector[n] evalz = eval;
  for (i in 1:k)  evalz[i] = 0.0;
  return evalz ;
}

matrix GS_complete(matrix Evec, int k, int n){
  // Gram-Schmidt process
  // Input: Evec - original sets of eigenvaectors
  //        k, n --the number of zero eigen values, a size of a matrix
  // Output: sets of eigenvaectors after Gram-Schumidt process
  matrix[n,k] Q;
  matrix[n,k] V ;
  matrix[n,k] X = Evec[,1:k];
  X[,1] =  rep_vector(1/sqrt(n),n) ;
  V = X ;
  Q[,1] =  V[,1];
  for (i in 2:k){
    for (j in 1:(i-1)){
      V[,i] = V[,i] - ((Q[,j]')*X[,i])*Q[,j];
    }
    Q[,i] = V[,i]/sqrt(sum(V[,i].*V[,i]));
  }
  return append_col(Q,Evec[,(k+1):n]);
}

matrix sq_cen_kernel_mat(matrix X, int N, real Hurst){
  // Output: Gram matrix with square centered fractional brownian motion
  // Input: X - covariates
  //        N - Number of observations (nrow of X)
  //        Hurst - Hurst coefficeint for fBM kernel
  matrix[N, N] K;
  matrix[N, N] E ;
  matrix[N, N] B = rep_matrix(0, N, N);;
  vector[N] d;
  matrix[N, N] A = diag_matrix(rep_vector(1,N)) - (1.0/N)*rep_matrix(1, N, N);
  matrix[N, N]  Xcp = X * X' ;
  vector[N] dvec = diagonal(Xcp);
  for (i in 1:(N-1)){
    d[i] = pow(fabs(dvec[i]), Hurst);
    for (j in (i+1):N){
      B[i,j] = pow(fabs(dvec[i] + dvec[j] - 2 * Xcp[i,j]), Hurst);
      B[j,i] = B[i,j];
    }
  }
  d[N] = pow(fabs(dvec[N]), Hurst);
  E = rep_matrix(d, N);
  K = 0.5 * (E + E' - B);
  K = A * K * A;
  return K * K ;
}

vector mat_vec_prod(matrix A, vector s, int Nd, int NNd ){
  //matrix-vector product involving Kronecker product
  matrix[Nd, NNd] S = to_matrix(s, Nd, NNd);
  matrix[NNd, Nd] Z = (A * S)';
  return to_vector(Z);
}
