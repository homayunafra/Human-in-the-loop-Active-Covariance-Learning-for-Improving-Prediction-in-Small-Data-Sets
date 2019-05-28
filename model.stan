data{
  int N; // the number of observations
  int M; // the number of features
  int N_f; // the number of pairs for which the user gave feedback until now
  int d; // the length of the feature descriptor
  vector[N] y; // the response
  matrix[N,M] X; // the observed data matrix
  matrix[M,d] FD; // the matrix of feature descriptors
  real a_sigma; // the hyper-prior of the residual noise
  real b_sigma; // the hyper-prior of the residual noise
  real a_tau; // the hyper-prior of the scaling parameter
  real b_tau; // the hyper-prior of the scaling parameter
  int<lower = 0, upper = 1> f[N_f]; // the feedback vector with values 1/0 where 1 stands for similarity and 0 stands for dissimilarity
  matrix[N_f,d] W_f; // the matrix containing distance vectors of the pairs for which the user gave feedback
  real threshold_mean;
}
parameters{
  vector<lower=0>[d] gamma; // the diagonal elements of the metric matrix A
  real<lower=0> threshold;
  real<lower=0> tau; // the scale parameters
}
transformed parameters{
  cov_matrix[M] C; // the prior correlation
  vector[N_f] q;
  
  C = quad_form_sym(diag_matrix(gamma), FD');
  C = exp(-0.5 * (rep_matrix(diagonal(C), M) + rep_matrix(diagonal(C)', M)) + C);
  for(m in 1:M) C[m, m] = C[m, m] + 1e-06;
  
  q = 1.0 ./ (1.0 + exp(W_f*gamma - threshold));
}
model{
  tau ~ inv_gamma(a_tau, b_tau);
  
  threshold ~ normal(threshold_mean, threshold_mean/2);
  
  gamma ~ normal(1.0,0.5); // normal(mean,std) and multi_normal(mean,(co)variance): in normal we give std and in multi-normal we give (co)variance
  
  f ~ bernoulli(q);

  {
    matrix[N,N] S;
    S = tau * quad_form_sym(C, X');
    for(n in 1:N) S[n, n] = S[n, n] + 1.0;
    y ~ multi_student_t(2.0 * a_sigma, rep_vector(0, N), (b_sigma / a_sigma) * S);
  }
}
