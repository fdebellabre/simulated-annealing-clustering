// [[Rcpp::plugins(cpp11)]]
#include <random>
//#include <Rcpp.h>
#include <math.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <chrono>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::plugins(cpp11)]]
using namespace arma;
using namespace Rcpp;
using std::pow;
using std::sqrt;

// [[Rcpp::export]]
NumericMatrix cpp_bary(NumericVector W, NumericMatrix Z, int K) {
  int nrow=Z.nrow();
  int ncol=Z.ncol();
  NumericMatrix output(K, ncol);

  for (int k=0; k<K; k++) {
    for (int col=0; col<ncol; col++) {
      double rowsum_k=0;
      int nrow_k=0;
      for (int row=0; row<nrow; row++) {
        if (W[row]==k+1) {
          rowsum_k += Z(row, col);
          nrow_k++;
        }
      }
      output(k,col) = (rowsum_k/nrow_k);
    }
  }
  return(output);
}

// [[Rcpp::export]]
double cpp_kmeans(NumericVector W, NumericMatrix Z, int K) {
  int nrow=Z.nrow();
  int ncol=Z.ncol();
  NumericMatrix Bary = cpp_bary(W, Z, K);
  double output=0;

  for (int k=0; k<K; k++) {
    for (int col=0; col<ncol; col++) {
      for (int row=0; row<nrow; row++) {
        if (W[row]==k+1) {
          output += pow(Z(row,col)-Bary(k,col), 2);
        }
      }
    }
  }
  return(output);
}

// [[Rcpp::export]]
NumericVector cpp_reassign(NumericVector W, NumericMatrix Z, int K, double Temp) {
  int nrow=Z.nrow();
  int ncol=Z.ncol();
  double dist_clust;
  double dist_k;
  double y;
  double thres;
  std::vector<int> clusters;
  for (int k=0; k<K; k++) clusters.push_back(k);

  NumericMatrix Bary = cpp_bary(W, Z, K);

  for (int row=0; row<nrow; row++) {
    // calcul de la distance à son propre cluster
    dist_clust=0;
    for (int k=0; k<K; k++) {
      if (k+1==W[row]) {
        for (int col=0; col<ncol; col++) {
          dist_clust += pow(Z(row,col)-Bary(k,col), 2);
        }
      }
    }
    // test reassignment to other clusters
    for (int& k: Rcpp::RcppArmadillo::sample(clusters, K, FALSE)) {
      if (k+1!=W[row]) {
        dist_k=0;
        for (int col=0; col<ncol; col++) {
          dist_k += pow(Z(row,col)-Bary(k,col), 2);
        }
        y = runif(1)[0];
        thres = exp(-(dist_k-dist_clust)/Temp);
        if (y < thres) {
          W[row] = k+1;
          break;
        }
      }
    }
  }
  return(W);
}

// [[Rcpp::export]]
NumericMatrix cpp_normalize(NumericMatrix Z) {
  int ncol = Z.ncol();
  int nrow = Z.nrow();
  double max;
  double min;

  for (int col=0; col<ncol; col++) {
    max=Z(0,col);
    min=max;
    for (int row=0; row<nrow; row++) {
      if (Z(row,col)<min) min=Z(row,col);
      else if (Z(row,col)>max) max=Z(row,col);
    }
    for (int row=0; row<nrow; row++) {
      Z(row,col) = (Z(row,col)-min)/(max-min);
    }
  }
  return(Z);
}

// [[Rcpp::export]]
NumericVector cpp_simanneal(NumericMatrix Z, int K, int N=1, double prob=0.3, int stop=20, double mu=0.95) {
  Z = cpp_normalize(Z);
  int nrow=Z.nrow();

  // Random assignment of points into clusters
  NumericVector W_c(nrow);
  for (int i=0; i<nrow; i++) {
    W_c[i] = trunc(1+K*runif(1)[0]);
  }
  NumericVector W_b = W_c;
  NumericVector W_t = W_c;

  // Objective function value
  double J_c = cpp_kmeans(W_c, Z, K);
  double J_b = J_c;
  double J_t = J_c;
  double y;

  // Initial temperature
  double temp_J;
  NumericVector temp_W(nrow);
  double Temp = 50*J_c;
  for (int size=0; size<50; size++) {
    temp_J = J_c-1;
    while (temp_J < J_c) {
      for (int i=0; i<nrow; i++) {
        temp_W[i] = trunc(1+K*runif(1)[0]);
      }
      temp_J = cpp_kmeans(temp_W, Z, K);
    }
    Temp -= temp_J;
  }
  Temp = Temp/log(prob);

  int not_improve = 0;
  while (not_improve < stop) {
    for (int i=0; i<N; i++) {
      // Obtain a trial assignment
      W_t = cpp_reassign(W_c, Z, K, Temp);

      // Assess reassignment performance :
      // if trial assignment is better, replace current assignment with trial
      // otherwise accept trial assignment with some probability
      J_t = cpp_kmeans(W_t, Z, K);
      if (J_t < J_c) {
        not_improve = 0;
        W_c = W_t;
        J_c = J_t;
        if (J_t < J_b) {
          W_b = W_t;
          J_b = J_t;
        }
      } else {
        not_improve++;
        y = runif(1)[0];
        if (y < exp(-(J_t-J_c)/Temp)) {
          W_c = W_t;
          J_c = J_t;
        }
      }
    }
    Temp = mu*Temp;
  }
  return(W_b);
}
