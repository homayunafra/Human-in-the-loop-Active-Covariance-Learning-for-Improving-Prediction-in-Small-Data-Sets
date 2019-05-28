library(MASS)

calWmat <- function(data, points_vec, num_of_comp)
{
  #Calculate the w-vectors needed in distance metric learning from Yang et al. (2007)
  #
  #Args:
  # data: datamatrix, samples in rows, variables in columns
  # points_vec: vector which includes the imseMatrix_Prior_posterior_meand's of points given as feedback,
  # num_of_comp: the number of eigenvectors we use in approximation of the distance metric
  #Returns:
  # The matrix which includes each of the w-vectors in its rows
  if(length(points_vec) == 0)
  {
    return(list(w_mat = NULL, n = 0))
  }
  points_mat <- matrix(points_vec, ncol = 2, byrow = TRUE)
  n_pairs <- nrow(points_mat)
  
  w_mat <- matrix(-1, nrow = n_pairs, ncol = num_of_comp)
  
  for(i in 1:n_pairs)
  {
    p1 <- points_mat[i,1]
    p2 <- points_mat[i,2]
    w_mat[i,1:num_of_comp] <- (data[p1,] - data[p2,])^2
  }
  return(list(w_mat = w_mat, n = n_pairs))
}

feedback_generator <- function(pairIndex, feedbackTable)
{
  sim_vec = c()
  dis_vec = c()
  f = c()
  for(i in 1:ncol(pairIndex)){
    if(feedbackTable[pairIndex[1,i],pairIndex[2,i]] == 1)
    {
      sim_vec = cbind(sim_vec,as.numeric(pairIndex[,i]))
      f = c(f,1)
    }
    else
    {
      dis_vec = cbind(dis_vec,as.numeric(pairIndex[,i]))
      f = c(f,0)
    }
  }
  # feedbackTable[pairIndex] <- 0
  # return(sim = sim, dis = dis, feedbackTable = feedbackTable)
  return(list(sim = sim_vec, dis = dis_vec, feedback = f))
}

compute_utilities <- function(fitted_model, feedback_list, feature_descriptors, N, M, d, y_tr, x_tr, a_sigma, b_sigma, a_tau, b_tau, muStar, vStar, a_sigma_star, b_sigma_star, Wf_total, feedback_total, q, threshold_0)
{
  num_feedbacks = ncol(feedback_list)
  Utility = vector(mode = "numeric",length = num_feedbacks)
  
  xtx = t(x_tr)%*%x_tr
  xty = t(x_tr)%*%y_tr
  yty = t(y_tr)%*%y_tr
  
  cores = detectCores()
  c1 = makeCluster(cores[1]-1)
  clusterExport(c1,c("calWmat","optimizing"))
  registerDoParallel(c1)
  
  # compute the expected utility of each single feedback: prob(feed = 0)*KL(feed = 0) + prob(feed = 1)*KL(feed = 1)
  parOutput = foreach(i = 1:num_feedbacks, .export = c("calWmat","optimizing"), .packages = "rstan") %dopar% {
  # for(i in 1:num_feedbacks){
    query_vec = as.vector(feedback_list[,i])
    feedback_1 = as.vector(c(feedback_total,1))
    
    Wf = calWmat(feature_descriptors,query_vec,d)$w_mat
    Wf_all = rbind(Wf_total, Wf)
    
    model_output.sim = optimizing(fitted_model, data = list(N = N, M = M, N_f = length(feedback_1), d = d, y = as.vector(y_tr), X = x_tr, FD = feature_descriptors, a_sigma = a_sigma, b_sigma = b_sigma, a_tau = a_tau, b_tau = b_tau, f = array(feedback_1,dim = length(feedback_1)), W_f = Wf_all, threshold_mean = threshold_0), as_vector=FALSE, iter = 10000)
    
    tau2C.sim = model_output.sim$par$tau * model_output.sim$par$C
    tau2Cinv.sim = chol(tau2C.sim)
    tau2Cinv.sim = chol2inv(tau2Cinv.sim)
    vStarInv.sim = tau2Cinv.sim + xtx
    vStar.sim = chol(vStarInv.sim)
    vStar.sim = chol2inv(vStar.sim)
    muStar.sim = vStar.sim%*%xty
    
    b_sigma_star_sim = b_sigma + 0.5*(yty - t(muStar.sim)%*%vStarInv.sim%*%muStar.sim)
    
    xvStarx = apply(t(x_tr) * (vStar%*%t(x_tr)),2,sum)    # equivalent to sum_{i=1}^{N} [t(x_tr[i,])%*%vStar%*%x_tr[i,]]
    xvStarx.sim = apply(t(x_tr) * (vStar.sim%*%t(x_tr)),2,sum)  # equivalent to sum_{i=1}^{N} [t(x_tr[i,])%*%vStar.sim%*%x_tr[i,]]
    KL_sim = ((x_tr%*%muStar.sim - x_tr%*%muStar)^2)*(a_sigma_star - 1) / (2*as.numeric(b_sigma_star)*(1+xvStarx))
    KL_sim = KL_sim + ((as.numeric(b_sigma_star_sim/(2*b_sigma_star)) * (1 + xvStarx.sim)) / (1 + xvStarx))
    KL_sim = KL_sim + log(sqrt( as.numeric(b_sigma_star / b_sigma_star_sim) * ( (1 + xvStarx) / (1 + xvStarx.sim) ) ))
    KL_sim = sum(KL_sim - 0.5)
    
    feedback_0 = as.vector(c(feedback_total,0))
    
    model_output.dis = optimizing(fitted_model, data = list(N = N, M = M, N_f = length(feedback_0), d = d, y = as.vector(y_tr), X = x_tr, FD = feature_descriptors, a_sigma = a_sigma, b_sigma = b_sigma, a_tau = a_tau, b_tau = b_tau, f = array(feedback_0,dim = length(feedback_0)), W_f = Wf_all, threshold_mean = threshold_0), as_vector=FALSE, iter = 10000)
    
    tau2C.dis = model_output.dis$par$tau * model_output.dis$par$C
    tau2Cinv.dis = chol(tau2C.dis)
    tau2Cinv.dis = chol2inv(tau2Cinv.dis)
    vStarInv.dis = tau2Cinv.dis + xtx
    vStar.dis = chol(vStarInv.dis)
    vStar.dis = chol2inv(vStar.dis)
    muStar.dis = vStar.dis%*%xty
    
    b_sigma_star_dis = b_sigma + 0.5*(yty - t(muStar.dis)%*%vStarInv.dis%*%muStar.dis)
    
    xvStarx.dis = apply(t(x_tr) * (vStar.dis%*%t(x_tr)),2,sum)  # equivalent to sum_{i=1}^{N} [t(x_tr[i,])%*%vStar.dis%*%x_tr[i,]]
    KL_dis = ((x_tr%*%muStar.dis - x_tr%*%muStar)^2)*(a_sigma_star - 1) / (2*as.numeric(b_sigma_star)*(1+xvStarx))
    KL_dis = KL_dis + ((as.numeric(b_sigma_star_dis/(2*b_sigma_star)) * (1 + xvStarx.dis)) / (1 + xvStarx))
    KL_dis = KL_dis + log(sqrt( as.numeric(b_sigma_star / b_sigma_star_dis) * ( (1 + xvStarx) / (1 + xvStarx.dis) ) ))
    KL_dis = sum(KL_dis - 0.5)
    
    if(length(q) > 1){
      # Utility[i] = q[i] * as.vector(KL_sim) + (1-q[i]) * as.vector(KL_dis)
      Utility = q[i] * as.vector(KL_sim) + (1-q[i]) * as.vector(KL_dis)
    }else{
      # Utility[i] = q * as.vector(KL_sim) + (1-q) * as.vector(KL_dis)
      Utility = q * as.vector(KL_sim) + (1-q) * as.vector(KL_dis)
    }
  }
  
  stopCluster(c1)
  
  return(unlist(parOutput))
}