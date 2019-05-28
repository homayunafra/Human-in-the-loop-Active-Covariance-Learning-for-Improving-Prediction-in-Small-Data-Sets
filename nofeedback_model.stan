data{
  int N; // the number of observations
  int M; // the number of features
  int d; // the length of the feature descriptor
  vector[N] y; // the response
  matrix[N,M] X; // the observed data matrix
  matrix[M,d] FD; // the matrix of feature descriptors
  real a_sigma; // the hyper-prior of the residual noise
  real b_sigma; // the hyper-prior of the residual noise
  real a_tau; // the hyper-prior of the residual noise
  real b_tau; // the hyper-prior of the residual noise
}
parameters{
  vector<lower=0>[d] gamma; // the diagonal elements of the metric matrix A
  real<lower=0> tau;
}
transformed parameters{
  cov_matrix[M] C; // the prior correlation
 
  C = quad_form_sym(diag_matrix(gamma), FD');
  C = exp(-0.5 * (rep_matrix(diagonal(C), M) + rep_matrix(diagonal(C)', M)) + C);
  for(m in 1:M) C[m, m] = C[m, m] + 1e-06;
}
model{
  tau ~ inv_gamma(a_tau, b_tau);
  
  gamma ~ normal(1.0,0.5); // normal(mean,std) and multi_normal(mean,(co)variance): in normal we give std and in multi-normal we give (co)variance
  
  {
    matrix[N,N] S;
    S = tau * quad_form_sym(C, X');
    for(n in 1:N) S[n, n] = S[n, n] + 1.0;
    y ~ multi_student_t(2.0 * a_sigma, rep_vector(0, N), (b_sigma / a_sigma) * S);
  }
}
