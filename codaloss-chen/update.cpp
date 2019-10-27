/* solution to the minimize problems in iterations */

#include <Rcpp.h>
using namespace Rcpp;



// [[Rcpp::export]]
double sgn(double x){
  if(x > 0)
    return 1.0;
  else if(x < 0)
    return -1.0;
  else return 0.0;
}

// [[Rcpp::export]]
NumericVector S_sl(NumericVector A, NumericVector Wlambda, int p){
  /* solution to the problem
  S(A,B) = argmin{1/2*<X^2, I> - <X, A> + lambda* Wo||X||_{1,off}
  subject to X = t(X), where A is symmetric, W is weight on the penalty */
  NumericVector B(A.size());
  
  int l = 0, i, j;
  for(i = 0; i < p; i++){
    B[l] = A[l];
    l++;
    for(j = (i+1); j < p; j++){
      if(fabs(A[l]) <= Wlambda[l])
        B[l] = 0;
      else 
        B[l] = A[l] - sgn(A[l]) * Wlambda[l];
      l++;
    }
  }
  
  return B;
  }




// [[Rcpp::export]]
NumericVector T_sl_exact(NumericMatrix invm, NumericVector A, 
                         IntegerMatrix loc, int p){
  /* the exact solution to problem T*/
  NumericVector B(p*p), C(A.size());
  int l = 0, i, j, k;
  
  for(i = 0; i < p; i++){
    for(j = 0; j < p; j++){
      B[l++] = A[loc(i,j)];
    }
  }
  
  double s;
  for(i = 0; i < p; i++){
    for(j = i; j < p; j++){
      s = 0;
      for(k = 0; k < (p*p); k++)
        s += invm(i*p+j,k) * B[k];
      C[loc(i,j)] = s;
    }
  }
  
  return C;
}
