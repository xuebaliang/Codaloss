#solve 1/2<Theta^2, A> - <Theta, B> + lambda ||Theta||1,off

##solution to the problem
##G(A,B) = argmin{1/2*<X^2, A> - <X, B>} subject to X = t(X),
##where A is positive definite, B is symmetric
Haddij = function(x,p){
  ret = matrix(NA,p,p)
  for(i in 1:p)
    for(j in i:p)
      ret[j,i] = ret[i,j] = 2/(x[i] + x[j])
    return(ret)
}

HG_sl = function(A, B, p){
  decomp = eigen(A, only.values = FALSE)
  
  U = decomp$vectors
  D = Haddij(decomp$values, p)
  
  ret = tcrossprod(U %*% ((crossprod(U, B) %*% U) * D), U)
  return((ret + t(ret))/2)
}



#solution to 1/2<X^2,I> - <X,A> + lambda||X||1,off subject to X = t(X)
HS_sl = function(A, lambda, p){
  for(i in 1:(p-1))
    for(j in (i+1):p){
      x = abs(A[i,j])
      A[i,j] = ifelse(x <= lambda, 0, sign(A[i,j])*(x - lambda))
      A[j,i] = A[i,j]
    }
  return(A)
}


#convergence condition
Hcond = function(X, lastX){
  return(norm(X-lastX, type = "F") / max(c(1, norm(X,type = 'F'), norm(lastX, type = 'F'))))
}


##solution to the problem
##R(A,B) = argmin{1/2*<X^2, I> - <X, A> subject to X >> eps*I,
##where A is symmetric
HR_sl = function(A, p, eps = 1e-8){
  decomp = eigen(A, only.values = FALSE)
  D = diag(p)
  diag(D) = pmax(decomp$values, eps)
  U = decomp$vectors
  ret = tcrossprod(U %*% D, U)
  return((ret+t(ret))/2)
}



Hcdtrace = function(A,B, lambda = 0, rho = 10, eps = 1e-8, 
                   tol = 1e-4, itermax = 1000, fastitermax = 1000, fasttol = 1e-4){
  
  ####define some constant matrix for convenient
  p = nrow(A)
  rhoI = rho * diag(p)
  tworhoI = 2 * rhoI
  
  
  ####initial
  Theta = Theta0 = diag(1/diag(A))
  Lambda = diag(p)
  lastTheta = Theta
  lastTheta0 = Theta0
  fastiter = 0
  iter = 0
  
  ##fast iterations
  convergence = FALSE
  while(!convergence & fastiter < fastitermax){
    Theta = HG_sl(A+rhoI, B+rho*Theta0-Lambda, p)
    Theta0 = HS_sl(Theta + Lambda/rho, lambda/rho, p)
    Lambda = Lambda + rho*(Theta - Theta0)
    convergence = Hcond(Theta, lastTheta) < fasttol & Hcond(Theta0, lastTheta0) < fasttol
    
    lastTheta = Theta; lastTheta0 = Theta0
    fastiter = fastiter + 1
  }
  
  ##iterations
  if(min(eigen(Theta0, only.values = TRUE)$values) < eps | !convergence){
    Theta1 = Theta; lastTheta1 = Theta1; Lambda1 = Lambda
    while(!convergence & iter < itermax){
      Theta = HG_sl(A + tworhoI, B + rho*Theta0 + rho*Theta1 - Lambda - Lambda1, p)
      Theta0 = HS_sl(Theta + Lambda/rho, lambda/rho, p)
      Theta1 = HR_sl(Theta + Lambda1/rho, p, eps)
      Lambda = Lambda + rho*(Theta - Theta0)
      Lambda1 = Lambda1 + rho*(Theta - Theta1)
      convergence = Hcond(Theta, lastTheta) < fasttol & Hcond(Theta0, lastTheta0) < tol & Hcond(Theta1, lastTheta1) < tol
      
      lastTheta = Theta; lastTheta0 = Theta0; lastTheta1 = Theta1
      iter = iter + 1
    }
  }
  
  Theta = Theta0
  return(list(Theta = Theta, fastiter, iter))
}


cdtrace_path = function(data, exact = TRUE, rho = 10, eps = 1e-8, lambda.min.ratio = 1e-4, nlambda = 30,
                        tol = 1e-4, itermax = 1000, fastitermax = 1000, fasttol = 1e-4){
  p = ncol(data); n = nrow(data)
  G = diag(p) - 1/p*matrix(1,p,p)
  A = G %*% (var(log(data))*(1-1/n)) %*% G
  if(exact){
    B = G
  }else{
    B = diag(p)
  }
  
  lambda.max <- max(max(A - diag(p)), -min(A - diag(p)))
  lambda.min <- lambda.min.ratio * lambda.max
  lambda <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
  
  icovpath = list()
  IC = c()
  num = 0
  for(l in lambda){
    res = Hcdtrace(A,B, lambda = l, rho = rho, eps = eps, tol = tol, itermax = itermax,
                  fastitermax = fastitermax, fasttol = fasttol)
    num = num + 1
    icovpath[[num]] = res$Theta
    # icovbic = 2*(1/2*sum(diag(Theta%*%Theta%*%A)) - sum(diag(Theta%*%B))) +
    #   log(n)*sum(res$Theta[upper.tri(res$Theta, diag = FALSE)]!=0)/n
    
    icovIC = n*norm((A%*%res$Theta + res$Theta%*%A)/2 - B, type = 'O') +
      log(n)*sum(res$Theta[upper.tri(res$Theta, diag = FALSE)]!=0)
    #cat(c(log(n)*sum(res$Theta[upper.tri(res$Theta, diag = FALSE)]!=0)/n,icovbic),"\n")
    IC = c(IC, icovIC)
  }
  icov= icovpath[[which.min(IC)]]
  # idx = which.min(IC)
  # tmpicov = icov
  # diag(tmpicov) = 0
  # while(all(tmpicov == 0)){
  #   idx = idx + 1
  #   tmpicov = icovpath[[idx]]
  #   diag(tmpicov) = 0
  # }
  
  return(list(lambda = lambda, icov = icov, icovpath = icovpath, IC = IC))
}



# test
# setwd("C:\\Users\\asus\\Desktop\\CDtrace\\CDtrace")
# source("basicsl.R")
# require(MASS)
# set.seed(1)
# p = 50
# n = 100
# mu = runif(p, -0.5, 0.5)
# Theta = diag(p)
# for(i in 1:(p-1))Theta[i,i+1] = 0.5*2
# for(i in 1:(p-2))Theta[i,i+2] = 0.2*2
# 
# Theta = (Theta + t(Theta))/2
# Sigma = solve(Theta)
# data = exp(mvrnorm(n = n, mu = mu, Sigma = Sigma))
# data = data / rowSums(data)
# # A = var(data)
# # B = diag(p)
# B = diag(p) - 1/p*matrix(1,p,p)
# A = B %*% var(log(data)) %*% B
# 
# lambda = 0
# rho = 10
# eps = 1e-8
# tol = 1e-4
# itermax = 1000
# fastitermax = 1000
# fasttol = 1e-4
# 
# aa = cdtrace(A,B, lambda = 0.03, rho = 10, eps = 1e-8, tol = 1e-4,
#         itermax = 1000, fastitermax = 1000, fasttol = 1e-4)
# aa$Theta %*% Sigma
# bb = cdtrace_path(data, rho = 10, eps = 1e-8, lambda.min.ratio = 1e-3, nlambda = 20,
#                         tol = 1e-4, itermax = 1000, fastitermax = 1000, fasttol = 1e-4)
# bb$icov
# plot(bb$bic)
# which.min(bb$bic)
# 
# 
# path = lapply(bb$icovpath, function(x){
#   x = 1*(x!=0)
#   #diag(x) = 0
#   return(x)
# })
# require(huge)
# #Theta_off = 1*(Theta!=0); diag(Theta_off) = 0
# Theta_off = Theta; diag(Theta_off) = 0
# huge.roc(path, Theta_off)
