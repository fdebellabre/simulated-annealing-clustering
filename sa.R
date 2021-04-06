
f_normalize <- function(Z) {
  # Output: normalized data matrix
  # Input Z: data matrix
  
  for (col in 1:ncol(Z)) {
    Z[,col] <- (Z[,col] - min(Z[,col])) / (max(Z[,col]) - min(Z[,col]))
  }
  return(Z)
}


f_bary <- function(W, Z, K) {
  ### Output: list of barycentres (vectors)
  ### Inputs:
  # W --> vector of assignment into clusters
  # Z --> data matrix
  # K --> number of clusters
  
  bary_list <- list()
  for (k in 1:K) {
    Zk <- Z[which(W==k),]
    bary_list[[k]] <- apply(Zk, 2, FUN=mean)
  }
  return(bary_list)
}


f_kmeans <- function(W, Z, K) {
  ### Output: k-means objective function value
  ### Inputs:
  # W --> vector of assignment into clusters
  # Z --> data matrix
  # K --> number of clusters
  
  Bary <- f_bary(W, Z, K)
  dist <- 0
  for (k in 1:K) {
    Zk <- Z[which(W==k),]
    dist <- dist + sum(sapply(1:nrow(Zk), FUN=function(x) sum((Zk[x,]-Bary[[k]])^2)))
  }
  return(dist)
}


f_reassign <- function(W, Z, K, Temp) {
  ### Output: new vector of assignment into clusters
  ### Inputs:
  # W --> vector of assignment into clusters
  # Z --> data matrix
  # K --> number of clusters
  # Temp --> temperature parameter
  
  Cl <- 1:K
  Bary <- f_bary(W, Z, K)
  Dist <- list()
  for (c in Cl) {
    Dist[[c]] <- sapply(1:nrow(Z), function(x) sum((Z[x,] - Bary[[c]])^2))
  }
  # for each input, we consider reassignment
  for (xi in 1:nrow(Z)) {
    Cx <- W[xi]
    # we consider it in every other cluster
    for (clust in sample(Cl[Cl!=Cx])) {
      if (runif(1) < exp(-max(0, Dist[[clust]][xi] - Dist[[Cx]][xi])/Temp)) {
        W[xi] <- clust
        break
      }
    }
  }
  return(W)
}


f_simanneal <- function(Z, K, FUN=f_kmeans, N=1, prob=0.3, stop=20, mu=0.95, out=TRUE) {
  ### Minimizes function FUN using simulated annealing.
  ### Output: assignments of points into clusters
  ### Input:
  # Z --> data matrix
  # K --> number of clusters
  # FUN --> criterion to minimize
  # N --> amount of generations for one temperature (how many reassignments before lowering the temperature)
  # stop --> stop parameter: amount of iterations without seeing any improvement
  # prob --> initial probability of accepting a cluster change
  # mu --> temperature multiplier (0<mu<1)
  
  # Abbreviations:
  # _b -> best
  # _t -> trial
  # _c -> current
  
  # Data normalization
  Z <- f_normalize(Z)
  
  # Random assignment into clusters
  W_b <- W_t <- W_c <- ceiling(K*runif(nrow(Z)))
  # Value of the objective function
  J_b <- J_t <- J_c <- FUN(W_c, Z, K)
  
  j_test <- c()
  for (test in 1:50) {
    new_j <- J_c-1
    while(new_j<J_c) {
      new_j <- FUN(ceiling(K*runif(nrow(Z))), Z, K)
    }
    j_test <- c(j_test, new_j)
  }
  Temp <- mean((J_c-j_test)/log(prob))
  
  if (out) message("\nIt\tBest\tCurrent\tTrial\tTemp")
  if (out) message(sprintf("0\t%.2f\t%.2f\t%.2f\t%.2f", J_b, J_c, J_t, Temp))
  Iter <- 0
  no_improvement <- 0
  while (no_improvement < stop) {
    for (k in 1:N) {
      # Getting a test assignment and evaluating its performance
      W_t <- f_reassign(W_c, Z, K, Temp)
      J_t <- FUN(W_t, Z, K)
      
      # We accept the test assignment if it is an improvement over the current one
      if (J_t < J_c) {
        no_improvement <- 0
        W_c <- W_t
        J_c <- J_t
        if (J_t < J_b) {
          W_b <- W_t
          J_b <- J_t
        }
      } else {
        # Otherwise, we still accept it with some probability
        if (runif(1) < exp(-(J_t-J_c)/Temp)) {
          W_c <- W_t
          J_c <- J_t
        }
        no_improvement <- no_improvement + 1
      }
    }
    # Lowering the temperature
    Temp <- mu*Temp
    Iter <- Iter+1
    if (Iter%%10==0 & out) message(sprintf("%i\t%.2f\t%.2f\t%.2f\t%.2f\t", Iter, J_b, J_c, J_t, Temp))
  }
  return(W_b)
}
