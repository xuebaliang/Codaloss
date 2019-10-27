
####some functions for convenience
checklen1num = function(x){
  stopifnot(is.numeric(x), length(x) == 1, x >= 0)
}
checklen1pnum = function(x){
  stopifnot(is.numeric(x), length(x) == 1, x > 0)
}
checklen1bool = function(x){
  stopifnot(is.logical(x), length(x) == 1)
}

checkparm = function(data, kappa, lambda, pweight, mu, rho, eps, method, 
                     tol, itermax, fast, fastitermax, fasttol,
                     init, validation, runlog, onlyval, show){
  stopifnot(is.matrix(data))
  stopifnot(all(is.numeric(data)))
  if(!is.null(pweight)){
    stopifnot(is.matrix(pweight))
    stopifnot(all(is.numeric(pweight)))
    stopifnot(all(pweight >= 0))
    stopifnot(dim(pweight) == rep(dim(data)[2], 2))
  }
  
  checklen1num(kappa); checklen1num(lambda); checklen1num(mu)
  checklen1pnum(rho); checklen1pnum(eps); checklen1pnum(tol); checklen1pnum(fasttol)
  checklen1num(itermax); checklen1num(fastitermax)
  checklen1bool(fast); checklen1bool(runlog)
  checklen1bool(onlyval); checklen1bool(show)

  
  if(!is.null(validation)){
    stopifnot(is.matrix(validation))
    stopifnot(dim(data)[2] == dim(validation)[2])
  }
  
  if(method != 'approximate' & method != "exact"){
    stop("method must be approximate or exact!")
  }
  
  if(!is.null(init$Theta)){
    stopifnot(is.matrix(init$Theta))
  }
  
  if(!is.null(init$Sigma)){
    stopifnot(is.matrix(init$Sigma))
  }
}



checkparmcv = function(data, folds, foldid, kappaseq, lambdaratio, maxnlambda,
                       pweight, museq, seed, rho, eps, method, tol, itermax, 
                       fast, fastitermax, fasttol, init, runlog, onlyval, show){
  stopifnot(is.matrix(data))
  stopifnot(all(is.numeric(data)))
  if(!is.null(pweight)){
    stopifnot(is.matrix(pweight))
    stopifnot(all(is.numeric(pweight)))
    stopifnot(all(pweight >= 0))
    stopifnot(dim(pweight) == rep(dim(data)[2], 2))
  }
  
  
  stopifnot(is.numeric(folds), folds <= dim(data)[1], folds >= 1)
  if(!is.null(foldid)){
    stopifnot(length(foldid) == dim(data)[1], setequal(unique(foldid), 1:folds))
  }
  
  stopifnot(is.numeric(kappaseq), is.numeric(museq)) 
  checklen1pnum(lambdaratio); checklen1pnum(maxnlambda)
  checklen1pnum(rho); checklen1pnum(eps); checklen1pnum(tol); checklen1pnum(fasttol)
  checklen1num(itermax); checklen1num(fastitermax); checklen1num(seed)
  checklen1bool(fast); checklen1bool(runlog); checklen1bool(onlyval); checklen1bool(show)
  
  if(method != 'approximate' & method != "exact"){
    stop("method must be approximate or exact!")
  }
  
  if(!is.null(init$Theta)){
    stopifnot(is.matrix(init$Theta))
  }
  
  if(!is.null(init$Sigma)){
    stopifnot(is.matrix(init$Sigma))
  }
}



symmatsum = function(x, diagloc){
  sum(x[diagloc]) + 2 * sum(x[-diagloc])
}



evaluate = function(Theta, Theta1, Theta2, Sigma, Sigma1, Sigma2, 
                    lambda11, lambda12, lambda21, lambda22, rho,
                    kappa, Wlambda, mu, p, Fm, diagloc, loc, kFSlnxF, fast){
  Theta1_dif = Theta - Theta1
  Sigma1_dif = Sigma - Sigma1
  Sigmamat = vec2symmat(Sigma, p)
  lagrange = 1/2 * symmatsum(multip_symmat(Theta, loc = loc, p = p) * 
                               multip_symmat(Sigma1, loc = loc, p = p), 
                             diagloc = diagloc) - 
    symmatsum(Theta * Sigma1, diagloc = diagloc) + 
    kappa / 2 * sum((Sigmamat %*% Fm) * (Fm %*% Sigmamat)) - 
    symmatsum(Sigma * kFSlnxF, diagloc = diagloc) + 
    2 * sum(abs(Wlambda[-diagloc] * Theta1[-diagloc])) + 
    mu * 2 * sum((Sigma[-diagloc])^2) + 
    symmatsum(lambda11 * Theta1_dif + lambda21 * Sigma1_dif, 
              diagloc = diagloc) + 
    rho/2 * symmatsum(Theta1_dif^2 + Sigma1_dif^2, diagloc = diagloc)
  
  if(!fast){
    Theta2_dif = Theta - Theta2
    Sigma2_dif = Sigma - Sigma2
    lagrange = lagrange + 
      symmatsum(lambda12 * Theta2_dif + lambda22 * Sigma2_dif,
                diagloc = diagloc) + 
      rho / 2 * symmatsum(Theta2_dif^2 + Sigma2_dif^2, diagloc = diagloc)
  }
  
  Theta = Theta1; Sigma = Sigma1
  Sigmamat = vec2symmat(Sigma, p)
  objvalue = 1/2 * symmatsum(multip_symmat(Theta, loc = loc, p = p) * 
                               multip_symmat(Sigma, loc = loc, p = p), 
                             diagloc = diagloc) - 
    symmatsum(Theta * Sigma, diagloc = diagloc) + 
    kappa / 2 * sum((Sigmamat %*% Fm) * (Fm %*% Sigmamat)) - 
    symmatsum(Sigma * kFSlnxF, diagloc = diagloc) + 
    2 * sum(abs(Wlambda[-diagloc] * Theta[-diagloc])) + 
    mu * 2 * sum((Sigma[-diagloc])^2)
  
  return(data.frame(lagrange = lagrange, objvalue = objvalue, fast = 1*fast))
}



converge_cond = function(lmat, nmat, diagloc){
  sqrt(symmatsum((lmat - nmat)^2, diagloc = diagloc)) /
    sqrt(max(1, symmatsum(lmat^2, diagloc = diagloc),
             symmatsum(nmat^2, diagloc = diagloc)))
}


ispdef = function(A, eps = 1e-8){
  return(all(eigen(A, symmetric = TRUE, only.values = TRUE)$values >= eps))
}


initifelse = function(init, name, mat, p){
  if(is.null(init[[name]]))return(mat) else
    return(symmat2vec(init[[name]]))
}