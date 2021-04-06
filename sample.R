rm(list=ls())
library(Rcpp)
library(RcppArmadillo)
source("sa.R")
sourceCpp("sa.cpp")

### ------------------------------ Sample code on IRIS data
### -----------------------------------------------------------------

##### -------------------- IRIS

X_iris <- as.matrix(iris[,1:4])
Y_iris <- as.integer(iris[,5])

Yhat_iris <- cpp_simanneal(X_iris, K=3)

PCA <- prcomp(X_iris, rank=2, scale=TRUE)

plot_iris <- data.frame(cbind(predict(PCA, X_iris), Yhat_iris, Y_iris))
table(predicted=Yhat_iris, true=iris$Species)


### ------------------------------ Precision benchmark (lower is better)
### ---------------------------------------------------------------------------

perf_KM <- function(Z, K, n=100) {
  ### Assess performance of Kmeans algorithm, n reps.
  # Returns a vector of the n values of the objective function
  value <- c()
  for (k in 1:n) {
    YKM <- kmeans(Z, K)$cluster
    value <- c(value, cpp_kmeans(YKM, cpp_normalize(Z), K=K))
  }
  return(c(mean=mean(value), min=min(value), max=max(value)))
}

perf_SA <- function(Z, K, n=20) {
  ### Assess performance of SA algorithm, n reps.
  # Returns a vector of the n values of the objective function
  value <- c()
  for (k in 1:n) {
    YSA <- cpp_simanneal(Z, K)
    value <- c(value, cpp_kmeans(YSA, cpp_normalize(Z), K=K))
  }
  return(c(mean=mean(value), min=min(value), max=max(value)))
}

SA_iris <- perf_SA(X_iris, K=3)
KM_iris <- perf_KM(X_iris, K=3)


### ------------------------------ Speed benchmark (lower is better)
### -------------------------------------------------------------------------

library(microbenchmark)
mbm <- microbenchmark(
  "R" = f_simanneal(X_iris, K=3, out=F),
  "C++" = cpp_simanneal(X_iris, K=3),
  "K-M" = kmeans(X_iris, 3),
  times=20, unit="ms"
)
mbm
