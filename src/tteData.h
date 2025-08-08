

#ifndef TTEDATA_H
#define TTEDATA_H

#include <RcppArmadillo.h>

// Data class definition
class tteData {
public:  
  arma::vec y;
  arma::uvec event;
  arma::mat X;
  int n;
  int p;
  arma::uvec cenindx;
  int ncen;
  
  tteData(
    const arma::vec& y_, const arma::uvec& event_, const arma::mat& X_
  ) {
    y        = y_;
    event    = event_;
    X        = X_;
    n        = y.n_elem;
    p        = X.n_cols;
    cenindx  = arma::find(event == 0);
    ncen     = cenindx.n_elem;
  } 
};


// Historical data class definition
class tteHistData {
public:  
  arma::vec y;
  arma::uvec event;
  arma::mat X;
  arma::mat Xunexch;
  int n;
  int p;
  int punexch;
  arma::uvec cenindx;
  int ncen;
  
  tteHistData(
    const arma::vec& y_, const arma::uvec& event_, const arma::mat& X_, const arma::mat& Xunexch_
  ) {
    y        = y_;
    event    = event_;
    X        = X_;
    Xunexch  = Xunexch_;
    n        = y.n_elem;
    p        = X.n_cols;
    punexch  = Xunexch.n_cols;
    cenindx  = arma::find(event == 0);
    ncen     = cenindx.n_elem;
  } 
};

#endif // TTEDATA_H



