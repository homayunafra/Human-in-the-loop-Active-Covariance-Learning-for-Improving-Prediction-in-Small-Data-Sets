setwd(getwd())

library(doParallel)
library(foreach)
library(MASS)
library(Metrics)
library(rstan)
library(Rcpp)
library(glmnet)

source("functions.R")
source("synthetic_data.R")
source("simulated_user.R")

feedback_per_round = 1
num_methods = 3     # 1: random, 2: non-sequential info. gain, 3: sequential info. gain
N = 15              # number of samples
M = 25              # number of dimensions
d = 9               # dimensionality of the feature descriptor
num_iterations = 101
num_feedbacks = 100

# hyperprior parameters
a_sigma = 2
b_sigma = 7
a_tau = 2
b_tau = 4

mseMatrix = matrix(0,num_iterations,num_methods)

set.seed(1234)
for(numRun in 1:5){
  data_address = paste("./Data/data_",numRun,".RData",sep="")
  load(data_address)
  x_tr                  = model_data$x_tr
  y_tr                  = model_data$y_tr
  x_te                  = model_data$x_te
  y_te                  = model_data$y_te
  feature_descriptors   = model_data$fd
  # syn_data                = simulated_data(numRun,N,M)
  # x_tr                    = syn_data$modelx
  # y_tr                    = syn_data$modely
  # x_te                    = syn_data$modelxte
  # y_te                    = syn_data$modelyte
  # feature_descriptors     = syn_data$fd
  c_star                  = model_data$c_star
  beta_star               = model_data$beta_star
  mse_gt                  = model_data$mse_gt
  # model_data              = list(x_tr = x_tr, x_te = x_te, y_tr = y_tr, y_te = y_te, c_star = c_star, beta_star = beta_star, mse_gt = mse_gt, fd = feature_descriptors)
  # save(model_data, file = paste("./Data/data_",numRun, ".RData",sep = ""))
  
  xtx           = t(x_tr)%*%x_tr
  xty           = t(x_tr)%*%y_tr
  yty           = t(y_tr)%*%y_tr

  ###################### Generate simulated feedback ######################

  sim_usr_output = simulated_user(beta_star)
  feedback_table = sim_usr_output$feedback_table
  feedback_list  = sim_usr_output$feedback_list

  ###################### Performance of the model with no feedback ######################
  # Computing parameters of the joint posterior distribution of \beta and \sigma^2, i.e. NIG(muStar, vStar, aStar, bStar), without any feedback
  nofeedback_fitted_model = stan_model(file = "nofeedback_model.stan")
  nofeedback_model_output = optimizing(nofeedback_fitted_model, data = list(N = N, M = M, d = d, y = as.vector(y_tr), X = x_tr, FD = feature_descriptors, a_sigma = a_sigma, b_sigma = b_sigma, a_tau = a_tau, b_tau = b_tau), as_vector=FALSE, iter = 50000)

  C_0 = nofeedback_model_output$par$C

  gamma_0 = nofeedback_model_output$par$gamma
  random_gamma = gamma_0
  nonseq_gamma = gamma_0
  seq_gamma = gamma_0

  tau_0 = nofeedback_model_output$par$tau
  random_tau = tau_0
  nonseq_tau = tau_0
  seq_tau = tau_0

  tau2C = tau_0 * C_0
  tau2Cinv = chol(tau2C)
  tau2Cinv = chol2inv(tau2Cinv)
  vStarInv = tau2Cinv + xtx
  vStar = chol(vStarInv)
  vStar = chol2inv(vStar) # covariance matrix of the joint posterior of \beta and \sigma^2
  muStar = vStar%*%xty    # mean vector of the joint posterior of \beta and \sigma^2

  # update the parameters of the inverse gamma distribution, i.e. b_sigma and a_sigma
  b_sigma_star = b_sigma + (1/2)*(yty - t(muStar)%*%vStarInv%*%muStar)
  a_sigma_star = a_sigma + (N/2)

  # compute the expected utility of each single feedback: prob(feed = 0)*KL(feed = 0) + prob(feed = 1)*KL(feed = 1)
  feedbacks = feedback_generator(feedback_list, feedback_table)
  sim_vecs = feedbacks$sim
  dis_vecs = feedbacks$dis
  sim_w = calWmat(feature_descriptors, sim_vecs, d)$w_mat
  dis_w = calWmat(feature_descriptors, dis_vecs, d)$w_mat
  w_all = rbind(sim_w,dis_w)
  threshold_0 = mean(w_all%*%gamma_0)
  random_threshold = threshold_0
  nonseq_threshold = threshold_0
  seq_threshold = threshold_0

  q = 1 / (1 + exp(w_all%*%gamma_0 - threshold_0))

  fitted_model = stan_model(file = "model.stan")
  nonseq_utility = compute_utilities(fitted_model, feedback_list, feature_descriptors, N, M, d, y_tr, x_tr, a_sigma, b_sigma, a_tau, b_tau, muStar, vStar, a_sigma_star, b_sigma_star, c(), c(), q, threshold_0)
  seq_utility = nonseq_utility

  random_feedback_list = nonseq_feedback_list = seq_feedback_list = feedback_list
  random_feedback_table = nonseq_feedback_table = seq_feedback_table = feedback_table

  random_Wf = nonseq_Wf = seq_Wf = c()

  random_Wf_All = nonseq_Wf_All = seq_Wf_All = w_all

  random_feedback_total = nonseq_feedback_total = seq_feedback_total = c()

  for(iteration in 1:num_iterations){
    if(iteration == 1){
      mseMatrix[iteration,1:3] = mse(x_te%*%muStar, y_te)
    }else{
      ###################### Random Query Selection ######################
      
      randomly_selected_index = sample(ncol(random_feedback_list), feedback_per_round)
      random_selected_pairs = as.matrix(random_feedback_list[,randomly_selected_index])
      random_feedback_list = as.matrix(random_feedback_list[,-randomly_selected_index])

      if(ncol(random_feedback_list)>0)
      {
        random_Wf = rbind(random_Wf, random_Wf_All[randomly_selected_index,])
        random_Wf_All = random_Wf_All[-randomly_selected_index,]
      }else{
        random_Wf = rbind(random_Wf, random_Wf_All)
      }

      random_given_feedback = feedback_generator(random_selected_pairs, random_feedback_table)
      random_feedback = random_given_feedback$feedback
      random_feedback_total = c(random_feedback_total, random_feedback)

      random_model_output = optimizing(fitted_model, data = list(N = N, M = M, N_f = feedback_per_round*(iteration-1), d = d, y = as.vector(y_tr), X = x_tr, FD = feature_descriptors, a_sigma = a_sigma, b_sigma = b_sigma, a_tau = a_tau, b_tau = b_tau, f =array(random_feedback_total,dim = feedback_per_round*(iteration-1)), W_f = random_Wf, threshold_mean = threshold_0), init = list(gamma = random_gamma, tau = random_tau),as_vector=FALSE, iter = 50000)

      random_C = random_model_output$par$C

	    random_gamma = random_model_output$par$gamma

      random_tau = random_model_output$par$tau

      random_tau2C = random_tau * random_C
      random_tau2Cinv = chol(random_tau2C)
      random_tau2Cinv = chol2inv(random_tau2Cinv)
      random_vStar = random_tau2Cinv + xtx
      random_vStar = chol(random_vStar)
      random_vStar = chol2inv(random_vStar)
      random_beta = random_vStar %*% xty

      mseMatrix[iteration,1] = mse(x_te%*%random_beta, y_te)

      ##################### Nonsequential ######################

      nonseq_utilities_sorted = sort.int(nonseq_utility,decreasing=TRUE,index.return=TRUE)
      nonseq_selected_index = nonseq_utilities_sorted$ix[1:feedback_per_round]
      nonseq_selected_pairs = as.matrix(nonseq_feedback_list[,nonseq_selected_index])
      nonseq_feedback_list = as.matrix(nonseq_feedback_list[,-nonseq_selected_index])
      nonseq_utility = nonseq_utility[-nonseq_selected_index]
      if(ncol(nonseq_feedback_list)>0){
        nonseq_Wf = rbind(nonseq_Wf, nonseq_Wf_All[nonseq_selected_index,])
        nonseq_Wf_All = nonseq_Wf_All[-nonseq_selected_index,]
      }else{
        nonseq_Wf = rbind(nonseq_Wf, nonseq_Wf_All)
      }

      nonseq_given_feedback = feedback_generator(nonseq_selected_pairs, nonseq_feedback_table)
      nonseq_feedback = nonseq_given_feedback$feedback
      nonseq_feedback_total = c(nonseq_feedback_total, nonseq_feedback)

      nonseq_model_output = optimizing(fitted_model, data = list(N = N, M = M, N_f = feedback_per_round*(iteration-1), d = d, y = as.vector(y_tr), X = x_tr, FD = feature_descriptors, a_sigma = a_sigma, b_sigma = b_sigma, a_tau = a_tau, b_tau = b_tau, f = array(nonseq_feedback_total, dim = feedback_per_round*(iteration-1)), W_f = nonseq_Wf, threshold_mean = threshold_0), init = list(gamma = nonseq_gamma, tau = nonseq_tau), as_vector=FALSE, iter = 50000)

      nonseq_C = nonseq_model_output$par$C

      nonseq_tau = nonseq_model_output$par$tau

      nonseq_gamma = nonseq_model_output$par$gamma

      nonseq_tau2C = nonseq_tau * nonseq_C
      nonseq_tau2Cinv = chol(nonseq_tau2C)
      nonseq_tau2Cinv = chol2inv(nonseq_tau2Cinv)
      nonseq_vStar = nonseq_tau2Cinv + xtx
      nonseq_vStar = chol(nonseq_vStar)
      nonseq_vStar = chol2inv(nonseq_vStar)
      nonseq_beta = nonseq_vStar %*% xty

      mseMatrix[iteration,2] = mse(x_te%*%nonseq_beta, y_te)

      ##################### Infogain - Sequential ######################

      if(iteration > 2){
        q = 1 / (1 + exp(seq_Wf_All%*%seq_gamma - seq_threshold))
        seq_utility = compute_utilities(fitted_model, seq_feedback_list, feature_descriptors, N, M, d, y_tr, x_tr, a_sigma, b_sigma, a_tau, b_tau, seq_beta, seq_vStar, seq_a_sigma_star, seq_b_sigma_star, seq_Wf, seq_feedback_total,q,threshold_0)
      }
      seq_utilities_sorted = sort.int(seq_utility,decreasing=TRUE,index.return=TRUE)
      seq_selected_index = seq_utilities_sorted$ix[1:feedback_per_round]
      seq_selected_pairs = as.matrix(seq_feedback_list[,seq_selected_index])
      seq_feedback_list = as.matrix(seq_feedback_list[,-seq_selected_index])

      if(ncol(seq_feedback_list)>0)
      {
        seq_Wf = rbind(seq_Wf, seq_Wf_All[seq_selected_index,])
        seq_Wf_All = seq_Wf_All[-seq_selected_index,]
      }else{
        seq_Wf = rbind(seq_Wf, seq_Wf_All)
      }

      seq_given_feedback = feedback_generator(seq_selected_pairs, seq_feedback_table)
      seq_feedback = seq_given_feedback$feedback
      seq_feedback_total = c(seq_feedback_total, seq_feedback)

      seq_model_output = optimizing(fitted_model, data = list(N = N, M = M, N_f = feedback_per_round*(iteration-1), d = d, y = as.vector(y_tr), X = x_tr, FD = feature_descriptors, a_sigma = a_sigma, b_sigma = b_sigma, a_tau = a_tau, b_tau = b_tau, f = array(seq_feedback_total,dim = feedback_per_round*(iteration-1)), W_f = seq_Wf, threshold_mean = threshold_0), init = list(gamma = seq_gamma, tau = seq_tau), as_vector=FALSE, iter = 50000)

      seq_C = seq_model_output$par$C

	    seq_gamma = seq_model_output$par$gamma

      seq_tau = seq_model_output$par$tau

      seq_tau2C = seq_tau * seq_C
      seq_tau2Cinv = chol(seq_tau2C)
      seq_tau2Cinv = chol2inv(seq_tau2Cinv)
      seq_vStarInv = seq_tau2Cinv + xtx
      seq_vStar = chol(seq_vStarInv)
      seq_vStar = chol2inv(seq_vStar)
      seq_beta = seq_vStar %*% xty

      # update the parameters of the inverse gamma distribution, i.e. b_sigma and a_sigma
      seq_b_sigma_star = b_sigma + (1/2)*(yty - t(seq_beta)%*%seq_vStarInv%*%seq_beta)
      seq_a_sigma_star = a_sigma + (N/2)

      mseMatrix[iteration,3] = mse(x_te%*%seq_beta, y_te)
    }
    save(mseMatrix, file = paste("./Results/mseMatrix_", numRun, ".RData",sep=""))
  }
}

