simulated_user = function(beta)
{
  num_feature = length(beta)
  indicesG1 = as.matrix(c(1:10))
  indicesG2 = as.matrix(c(11:15))
  indicesG3 = as.matrix(c(16:20))
  indicesG4 = as.matrix(c(21:25))
  
  similarPairsG1 = as.matrix(combn(indicesG1,2))
  similarPairsG2 = as.matrix(combn(indicesG2,2))
  similarPairsG3 = as.matrix(combn(indicesG3,2))
  similarPairsG4 = as.matrix(combn(indicesG4,2))
  similarPairs = cbind(similarPairsG1, similarPairsG2, similarPairsG3, similarPairsG4)
  
  dissimilarPairsFirstRow_1 = t(as.matrix(rep(indicesG1, each = length(as.vector(c(indicesG2,indicesG3,indicesG4))))))
  dissimilarPairsSecondRow_1 = t(as.matrix(rep(as.vector(c(indicesG2,indicesG3,indicesG4)), length(indicesG1))))
  dissimilarPairs_1 = as.matrix(rbind(dissimilarPairsFirstRow_1, dissimilarPairsSecondRow_1))
  
  dissimilarPairsFirstRow_2 = t(as.matrix(rep(indicesG2, each = length(as.vector(c(indicesG3,indicesG4))))))
  dissimilarPairsSecondRow_2 = t(as.matrix(rep(as.vector(c(indicesG3,indicesG4)), length(indicesG2))))
  dissimilarPairs_2 = as.matrix(rbind(dissimilarPairsFirstRow_2, dissimilarPairsSecondRow_2))
  
  dissimilarPairsFirstRow_3 = t(as.matrix(rep(indicesG3, each = length(as.vector(c(indicesG4))))))
  dissimilarPairsSecondRow_3 = t(as.matrix(rep(as.vector(c(indicesG4)), length(indicesG3))))
  dissimilarPairs_3 = as.matrix(rbind(dissimilarPairsFirstRow_3, dissimilarPairsSecondRow_3))
  
  dissimilarPairs = as.matrix(cbind(dissimilarPairs_1,dissimilarPairs_2,dissimilarPairs_3))
  
  feedback_table = matrix(0,num_feature, num_feature)
  for(i in 1:ncol(similarPairs)) {
    feedback_table[similarPairs[1,i], similarPairs[2,i]] = 1
  }
  
  dissimilarPairs = apply(dissimilarPairs,2,sort,decreasing=F)
  for(i in 1:ncol(dissimilarPairs)){
    feedback_table[dissimilarPairs[1,i], dissimilarPairs[2,i]] = -1
  }
  
  feedback_list = cbind(similarPairs, dissimilarPairs)
  return(list(feedback_table = feedback_table, feedback_list = feedback_list))
}