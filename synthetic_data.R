library(MASS)
library(Metrics)
library(corpcor)

simulated_data <- function(run,n,P)
{
  nte = 1000
  d_g1 = 10
  d_g2 = 5
  d_g3 = 5
  d_g4 = 5
  gamma_opt = as.vector(c(1,1,0,1,0,1,0,0,0))
  sigma_noise = 5
  
  fd_1 = as.vector(c(runif(d_g1,0.0,0.25),runif((P - d_g1),2,2.25)))
  fd_2 = as.vector(c(runif(d_g1,2,2.25),runif(d_g2,0.0,0.25),runif((d_g3+d_g4),2,2.25)))
  fd_3 = as.vector(c(runif(d_g1,1,10),runif((P - d_g1),2,2.25)))
  fd_4 = as.vector(c(runif((d_g1+d_g2),2,2.25),runif(d_g3,0.0,0.25),runif(d_g4,2,2.25)))
  fd_5 = as.vector(c(runif(d_g1,2,2.25),runif(d_g2,1,10),runif((d_g3+d_g4),2,2.25)))
  fd_6 = as.vector(c(runif((P-d_g4),2,2.25),runif(d_g4,0.0,0.25)))
  fd_7 = as.vector(c(runif((d_g1+d_g2),2,2.25),runif(d_g3,1,10),runif(d_g4,2,2.25)))
  fd_8 = as.vector(c(runif((P-d_g4),2,2.25),runif(d_g4,1,10)))
  fd_9 = as.vector(c(runif(P,1,10)))
  
  feature_descriptor = t(rbind(fd_1,fd_2,fd_3,fd_4,fd_5,fd_6,fd_7,fd_8,fd_9))
  
  C = diag(P)
  for(r in 1:(P-1)){
    for(c in (r+1):P){
      C[r,c] = t(feature_descriptor[r,] - feature_descriptor[c,])%*%diag(gamma_opt)%*%(feature_descriptor[r,] - feature_descriptor[c,])
      C[r,c] = exp(-0.5*C[r,c])
    }
  }
  C = C + t(C)
  diag(C) = diag(C) - 1
  
  beta_opt = mvrnorm(1,rep(0,P),sigma_noise*C)
  
  x_tr = matrix(rnorm(n * P), n, P)
  x_te = matrix(rnorm(nte * P), nte, P)
  y_tr = as.vector(x_tr %*% beta_opt + sqrt(sigma_noise)*rnorm(n))
  y_te = as.vector(x_te %*% beta_opt + sqrt(sigma_noise)*rnorm(nte))
  
  mse_ground_truth = mse(x_te%*%beta_opt, y_te)
  
  return(list(modelx = x_tr, modely = y_tr, modelxte = x_te, modelyte = y_te, beta_star = beta_opt, c_star = C, fd = feature_descriptor, mse_gt = mse_ground_truth))
}
