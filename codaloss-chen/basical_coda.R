####basic solutions to some minimize problems#####



##solution to the problem
##G(A,B) = argmin{1/2*<X^2, A> - <X, B>} subject to X = t(X),
##where A is positive definite, B is symmetric
G_sl = function(A, B, p){
  A = vec2symmat(A, p)
  B = vec2symmat(B, p)
  
  decomp = eigen(A, symmetric = TRUE, only.values = FALSE)
  
  if(any(decomp$values <= 0))
    stop('The matrix A may be not positive definite!')
  
  U = decomp$vectors
  D = addij(decomp$values)
  
  ret = tcrossprod(U %*% ((crossprod(U, B) %*% U) / D), U)
  
  return(symmat2vec(ret))
}



##solution to the problem
##R(A,B) = argmin{1/2*<X^2, I> - <X, A> subject to X >> eps*I,
##where A is symmetric
R_sl = function(A, p, eps = 1e-8){
  A = vec2symmat(A, p)
  decomp = eigen(A, symmetric = TRUE, only.values = FALSE)
  
  D = diag(p)
  diag(D) = pmax(decomp$values, eps)
  U = decomp$vectors
  
  return(symmat2vec(tcrossprod(U %*% D, U)))
}




##solution to the problem
##T(A,B) = argmin{eta*<X^2, A> + kappa*<XF, FX> - <X, A> 
##                - mu*sum(t(Ij) t(X) Ij t(Ij) X Ij)} 
##subject to X = t(X),
##where A is symmetric and Ij is the j_th column of identity matrix
Fm_decomp = function(Fm, p, eta, kappa, fasteta = NULL){
  decomp = eigen(Fm, symmetric = TRUE, only.values = FALSE)
  DF = kappa * multipij(decomp$values)
  
  if(is.null(fasteta)){
    fastDF = NULL
  }else{
    fastDF = fasteta + DF
  }
  
  DF = DF + eta
  
  return(list(DF = DF,
              fastDF = fastDF,
              UF = decomp$vectors))
}

T_sl = function(A, kappa, mu, UF, DF, p, last_X, diagloc){
  A = A / 2 
  A[diagloc] = A[diagloc] + mu * last_X[diagloc]
  A = vec2symmat(A, p)
  A = tcrossprod(UF %*% ((crossprod(UF, A) %*% UF) / DF), UF)
  return(symmat2vec(A))
}



pre_exact = function(p, eta, kappa, mu, Fm, fasteta = NULL){
  I = diag(p)
  
  ind = sapply(1:p, function(j){
    (j-1) * p + j
  })
  
  Iii = rep(0, times = p^2)
  Iii[ind] = 1
  Iii = diag(x = Iii)
  
  kII = kronecker(I, I)
  A =  2 * kappa * kronecker(Fm, Fm) - 2 * mu * Iii
  if(is.null(fasteta)){
    fastinvm = NULL
  }else{
    fastinvm = solve(2 * fasteta * kII + A)
  }
  
  invm = solve(2 * eta * kII  + A)
  return(list(invm = invm,
              fastinvm = fastinvm))
}
