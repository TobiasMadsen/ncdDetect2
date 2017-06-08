#include <Rcpp.h>
using namespace Rcpp;

double phi(double t, NumericVector & p, NumericVector & s, int idxStart, int idxStop){
  double ret = 0;
  for(int i = idxStart; i < idxStop; ++i){
    ret += p[i]*exp(t*s[i]);
  }
  return ret;
}

double phiD1(double t, NumericVector & p, NumericVector & s, int idxStart, int idxStop){
  double ret = 0;
  for(int i = idxStart; i < idxStop; ++i){
    ret += p[i]*s[i]*exp(t*s[i]);
  }
  return ret;
}

double phiD2(double t, NumericVector & p, NumericVector & s, int idxStart, int idxStop){
  double ret = 0;
  for(int i = idxStart; i < idxStop; ++i){
    ret += p[i]*s[i]*s[i]*exp(t*s[i]);
  }
  return ret;
}

// [[Rcpp::export]]
NumericVector cumulantDerivatives(double t, IntegerVector x, NumericVector p, NumericVector s) {
  // Assert in R that all vectors have the same length
  // Return cumulant and first and second derivative
  NumericVector ret(3);

  int idxStart = 0; 
  // Loop over vector
  for(int i = 1; i <= x.length(); ++i){
    if(i == x.length() || x(i-1) != x(i)){
      double phi_   = phi(t, p, s, idxStart, i);
      double phiD1_ = phiD1(t, p, s, idxStart, i);
      double phiD2_ = phiD2(t, p, s, idxStart, i);

      ret(0) += log( phi_ );
      ret(1) += phiD1_ / phi_;
      ret(2) += (phi_ * phiD2_ - phiD1_ * phiD1_) / phi_ / phi_;

      idxStart = i;
    }
  }

  return ret;
}
