/* some basic operation of symmetric matrix */

#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector symmat2vec(NumericMatrix A){
  /* convert a symmatrix to a vector, row major */
  int p = A.nrow();
  NumericVector a(p*(p+1)/2);
  
  int l = 0, i, j;
  for(i = 0; i < p; i++){
    for(j = i; j < p; j++){
      a[l++] = A(i,j);
    }
  }
  
  return a;
}



// [[Rcpp::export]]
NumericMatrix vec2symmat(NumericVector a, int p){
  /* convert a vector to a symmatrix, row major */
  NumericMatrix A(p,p);
  
  int l = 0, i, j;
  for(i = 0; i < p; i++){
    A(i,i) = a[l++];
    for(j = (i+1); j < p; j++){
      A(i,j) = a[l++];
      A(j,i) = A(i,j);
    }
  }
  
  return A;
}



// [[Rcpp::export]]
NumericMatrix addij(NumericVector d){
  /* create a symmatrix with ij-th element is (d_i + d_j)/2 */
  int p = d.size();
  NumericMatrix D(p,p);
  
  int i, j;
  for(i = 0; i < p; i++){
    D(i,i) = d[i];
    for(j = (i+1); j < p; j++){
      D(i,j) = (d[i] + d[j]) / 2;
      D(j,i) = D(i,j);
    }
  }
  
  return D;
}



// [[Rcpp::export]]
NumericMatrix multipij(NumericVector d){
  /* create a symmatrix with ij-th element is d_i * d_j */
  int p = d.size();
  NumericMatrix D(p,p);
  
  int i, j;
  for(i = 0; i < p; i++){
    D(i,i) = d[i]*d[i];
    for(j = (i+1); j < p; j++){
      D(i,j) = d[i] * d[j];
      D(j,i) = D(i,j);
    }
  }
  
  return D;
}



// [[Rcpp::export]]
IntegerMatrix getlocmat(int p){
  /* get the location of the element in the symmetric matrix */
  IntegerMatrix loc(p,p);
  
  int l = 0, i, j;
  for(i = 0; i < p; i++){
    loc(i,i) = l++;
    for(j = (i+1); j < p; j++){
      loc(i,j) = l++;
      loc(j,i) = loc(i,j);
    }
  }
  
  return loc;
}




// [[Rcpp::export]]
NumericVector multip_symmat(NumericVector A,
                            IntegerMatrix loc, int p){
  /* multiplication of symmetric matrix saved by vector, A%*%A */
  
  NumericVector C(A.size());
  int l = 0, i, j, k;
  double s;
  for(i = 0; i < p; i++){
    for(j = i; j < p; j++){
      s = 0;
      for(k = 0; k < p; k++){
        s += A[loc(i,k)] * A[loc(j,k)];
      }
      C[l++] = s;
    }
  }
  
  return C;
}