// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//
// [[Rcpp::export]]
Rcpp::List opnmfR_opnmfRcpp(const arma::mat & X, const arma::mat & W0, double tol=0.00001, int maxiter=10000, double eps=1e-16, bool memsave=1) {
  arma::mat W = W0;
  arma::mat Wold = W;
  //arma::mat XXW = X * (X.t()*W);
  double diffW = 9999999999.9;
  
  //Rcout << "The value of maxiter : " << maxiter << "\n";
  //Rcout << "The value of tol : " << tol << "\n";
  
  arma::mat XX;
  if(memsave==0) {
    XX = X * X.t();
  }
  
  int i;
  for (i = 1; i <= maxiter; i++) {
    if(memsave==0) {
      W = W % (XX*W) / (W*(W.t()*(XX*W)));
      //XXW = XX * W;
    } else {
      W = W % (X*(X.t()*W)) / (W*((W.t()*X)*(X.t()*W)));
      //XXW = X * (X.t()*W);
    }
    //W = W % XXW / (W  * (W.t() * XXW));
    
    
    arma::uvec idx = find(W < eps);
    W.elem(idx).fill(eps);
    W = W / norm(W,2);
    diffW = norm(Wold-W, "fro") / norm(Wold, "fro");
    if(diffW < tol) {
      break;
    } else {
      Wold = W;
    }
    
    if(i % 10 == 0) {
      Rcpp::checkUserInterrupt();
    }
    
  }
  return Rcpp::List::create(Rcpp::Named("W")=W,
                            Rcpp::Named("iter")=i,
                            Rcpp::Named("diffW")=diffW);
}


