###solution to codaloss

codaloss = function(data, kappa = 1, lambda = 0, pweight = NULL, mu = 0, 
                    rho = 10, eps = 1e-8, method = 'approximate', 
                    tol = 1e-4, itermax = 1000, 
                    fast = TRUE, fastitermax = 1000, fasttol = 1e-4,
                    init = list(Theta = NULL, Sigma = NULL), 
                    validation = NULL, runlog = TRUE, 
                    onlyval = FALSE, show = TRUE){
  
  ####check parameters
  checkparm(data = data, kappa = kappa, lambda = lambda, 
            pweight = pweight, mu = mu, rho = rho, eps = eps, 
            method = method, tol = tol, itermax = itermax, 
            fast = fast, fastitermax = fastitermax, 
            fasttol = fasttol, init = init, validation = validation, 
            runlog = runlog, onlyval = onlyval, show = show)
  
  
  
  ####define some constant matrix for convenient
  p = dim(data)[2]
  Fm = diag(p) - 1/p
  
  Slnx = var(log(data))
  kFSlnxF = kappa * symmat2vec(Fm %*% Slnx %*% Fm)
  Slnx = symmat2vec(Slnx)
  
  loc = getlocmat(p = p)
  diagloc = diag(loc) + 1
  
  I = numeric(length = p*(p+1)/2)
  I[diagloc] = 1
  rhoI = rho * I
  tworhoI = 2 * rhoI
  
  if(is.null(pweight)){
    pweight = matrix(1, nrow = p, ncol = p)
  }
  Wlambda = lambda * symmat2vec(pweight)
  
  
  if(fast)fasteta = rho / 2 + mu else fasteta = NULL
  if(method == 'approximate'){
    Fm_dcp = Fm_decomp(Fm = Fm, p = p, eta = rho + mu, kappa = kappa / 2,
                       fasteta = fasteta)
    UF = Fm_dcp$UF; DF = Fm_dcp$DF
    fastDF = Fm_dcp$fastDF
    invm = fastinvm = NULL
    Fm_dep = NULL
  }else 
    if(method == 'exact'){
      prex = pre_exact(p = p, eta = rho + mu, kappa = kappa / 2, mu = mu, Fm = Fm,
                       fasteta = fasteta)
      invm = prex$invm
      fastinvm = prex$fastinvm
      UF = DF = fastDF = NULL
      prex = NULL
    }
  
  
  ####initial
  Theta1 = Theta2 = Theta = initifelse(init, 'Theta', I, p)
  Sigma1 = Sigma2 = Sigma  = initifelse(init, 'Sigma', Slnx, p)
  lambda11 = lambda12 = lambda21 = lambda22 = I
  
  
  ####main loops
  converge = fastconverge = FALSE
  iter = fastiter = 0
  condlog = NULL
  if(runlog){
    evallog = evaluate(Theta = Theta, Theta1 = Theta1, Theta2 = Theta2, 
                       Sigma = Sigma, Sigma1 = Sigma1, Sigma2 = Sigma2, 
                       lambda11 = lambda11, lambda12 = lambda12, 
                       lambda21 = lambda21, lambda22 = lambda22, 
                       rho = rho, kappa = kappa, Wlambda = Wlambda, mu = mu, 
                       p = p, Fm = Fm, kFSlnxF = kFSlnxF, diagloc = diagloc,
                       loc = loc, fast = fast)
  }else{
    evallog = NULL
  }
  
  Theta1ispdef = FALSE
  if(fast){
    while(!fastconverge & fastiter < fastitermax){
      lTheta = Theta; lTheta1 = Theta1; lSigma = Sigma; lSigma1 = Sigma1
      llambda11 = lambda11; llambda21 = lambda21
      
      update = codaloss_optim_fast(lTheta = lTheta, lTheta1 = lTheta1,
                                   lSigma = lSigma, lSigma1 = lSigma1,
                                   llambda11 = llambda11, llambda21 = llambda21, 
                                   kFSlnxF = kFSlnxF, p = p, rhoI = rhoI,
                                   rho = rho, mu = mu, kappa = kappa, Wlambda = Wlambda, 
                                   eps = eps, loc = loc, diagloc = diagloc,
                                   method = method, UF = UF, fastDF = fastDF, fastinvm = fastinvm)
      
      Theta = update$Theta; Theta1 = update$Theta1
      Sigma = update$Sigma; Sigma1 = update$Sigma1
      lambda11 = update$lambda11; lambda21 = update$lambda21
      
      cond = c(converge_cond(lTheta, Theta, diagloc = diagloc), 
               converge_cond(lTheta1, Theta1, diagloc = diagloc), NA, 
               converge_cond(lSigma, Sigma, diagloc = diagloc),
               converge_cond(lSigma1, Sigma1, diagloc = diagloc), NA)
      fastconverge = max(cond[c(-3,-6)]) < fasttol
      
      if(runlog){
        condlog = rbind(condlog, cond)
        evallog = rbind(evallog, evaluate(Theta = Theta, Theta1 = Theta1, Theta2 = Theta2, 
                                          Sigma = Sigma, Sigma1 = Sigma1, Sigma2 = Sigma2, 
                                          lambda11 = lambda11, lambda12 = lambda12, 
                                          lambda21 = lambda21, lambda22 = lambda22, 
                                          rho = rho, kappa = kappa, Wlambda = Wlambda, 
                                          mu = mu, p = p, Fm = Fm, kFSlnxF = kFSlnxF, 
                                          diagloc = diagloc, loc = loc, fast = TRUE))
      }
      
      fastiter = fastiter + 1
    }
    
    if(fastiter == fastitermax)warning('Reach fastitermax!')
    
    Sigma1mat = vec2symmat(Sigma1, p = p)
    Theta1mat = vec2symmat(Theta1, p = p)
    Sigma1ispdef = ispdef(Sigma1mat, eps = eps)
    Theta1ispdef = ispdef(Theta1mat, eps = eps)
    if(!Sigma1ispdef)warning('Sigma may NOT be positive definite after fast iterations!')
    if(!Theta1ispdef)warning('Theta may NOT be positive definite after fast iterations!')
    
    if(show){
      cat('Fast iterations done and', 
          ifelse(fastconverge, 'converge', 'may not converge'),'after', fastiter, 'iterations.\n')
    }
    
    if(Theta1ispdef){
      if(show)cat('Theta is positive definite after fast iterations.\n')
      if(fastconverge){
        Theta = Theta1mat; Sigma = Sigma1mat
      }else{
        if(show)cat('Turn into iterations.\n')
        Theta2 = Theta1; Sigma2 = Sigma1
        lambda12 = lambda11; lambda22 = lambda21
      }
    }else{
      if(show)cat('Theta is not positive definite after fast iterations and
          will turn into iterations that guarantee positive definite.\n')
      Theta2 = Theta1; Sigma2 = Sigma1
      lambda12 = lambda11; lambda22 = lambda21
    }
  }
  
  if(!(fast & fastconverge & Theta1ispdef)){
    while(!converge & iter < itermax){
      lTheta = Theta; lTheta1 = Theta1; lTheta2 = Theta2
      lSigma = Sigma; lSigma1 = Sigma1; lSigma2 = Sigma2
      llambda11 = lambda11; llambda12 = lambda12
      llambda21 = lambda21; llambda22 = lambda22
      
      update = codaloss_optim(lTheta = lTheta, lTheta1 = lTheta1, lTheta2 = lTheta2, 
                              lSigma = lSigma, lSigma1 = lSigma1, lSigma2 = lSigma2, 
                              llambda11 = llambda11, llambda12 = llambda12, 
                              llambda21 = llambda21, llambda22 = llambda22,
                              kFSlnxF = kFSlnxF, p = p, rhoI = rhoI, tworhoI = tworhoI,
                              rho = rho, mu = mu, kappa = kappa, Wlambda = Wlambda, eps = eps,
                              loc = loc, diagloc = diagloc,
                              method = method, UF = UF, DF = DF, invm = invm)
      
      Theta = update$Theta; Theta1 = update$Theta1; Theta2 = update$Theta2
      Sigma = update$Sigma; Sigma1 = update$Sigma1; Sigma2 = update$Sigma2
      lambda11 = update$lambda11; lambda12 = update$lambda12
      lambda21 = update$lambda21; lambda22 = update$lambda22
      
      cond = c(converge_cond(lTheta, Theta, diagloc = diagloc), 
               converge_cond(lTheta1, Theta1, diagloc = diagloc),
               converge_cond(lTheta2, Theta2, diagloc = diagloc), 
               converge_cond(lSigma, Sigma, diagloc = diagloc),
               converge_cond(lSigma1, Sigma1, diagloc = diagloc), 
               converge_cond(lSigma2, Sigma2, diagloc = diagloc))
      converge = max(cond) < tol
      
      if(runlog){
        condlog = rbind(condlog, cond)
        evallog = rbind(evallog, evaluate(Theta = Theta, Theta1 = Theta1, Theta2 = Theta2, 
                                          Sigma = Sigma, Sigma1 = Sigma1, Sigma2 = Sigma2, 
                                          lambda11 = lambda11, lambda12 = lambda12, 
                                          lambda21 = lambda21, lambda22 = lambda22, 
                                          rho = rho, kappa = kappa, Wlambda = Wlambda, 
                                          mu = mu, Fm = Fm, p = p, kFSlnxF = kFSlnxF,
                                          diagloc = diagloc, loc = loc, fast = FALSE))
      }
      
      iter = iter + 1
    }
    
    if(iter == itermax)warning('Reach itermax and may NOT converge!')
    
    if(show){
      cat('Iterations done and', 
        ifelse(converge, 'converge', 'may not converge'), 'after', iter, 'iterations.\n')
    }
    
    ##check positive definite and symmetric
    Theta = vec2symmat(Theta1, p = p)
    Sigma = vec2symmat(Sigma1, p = p)
    Sigma1ispdef = ispdef(Sigma, eps = eps)
    Theta1ispdef = ispdef(Theta, eps = eps)
    if(!Sigma1ispdef)warning('Sigma may NOT be positive definite!')
    if(!Theta1ispdef)warning('Theta may NOT be positive definite!')
    if(show)cat('Theta', ifelse(Theta1ispdef, 'is', 'is not'), 'positive definite after iterations.\n')
  }
  
  
  ####validation value
  if(!is.null(validation)){
    validation_Slnx = var(log(validation))
    #validation = sum((Theta %*% Sigma - diag(p))^2) + kappa * sum((Fm %*% (validation_Slnx - Sigma) %*% Fm)^2)
    validation = norm(Theta %*% Sigma - diag(p), type = "1") + kappa * sum((Fm %*% (validation_Slnx - Sigma) %*% Fm)^2)
  }
  
  nonzeronum = (sum(Theta != 0) + p) / 2
  invbias = (norm(Theta %*% Sigma - diag(p), type = "F") + 
    norm(Sigma %*% Theta - diag(p), type = "F"))/2
  compbias = norm(Fm %*% (Sigma - vec2symmat(Slnx, p)) %*% Fm, type = "F")
  
  if(onlyval){
    return(list(nonzeronum = nonzeronum,
                invbias = invbias, compbias = compbias,
                converge = converge, iter = iter, 
                fastconverge = fastconverge, fastiter = fastiter, 
                validation = validation, 
                Sigmaispdef = Sigma1ispdef, Thetaispdef = Theta1ispdef))
  }else{
    return(list(Sigma = Sigma, Theta = Theta,
                condlog = condlog, evallog = evallog,
                nonzeronum = nonzeronum,
                invbias = invbias, compbias = compbias,
                converge = converge, iter = iter, 
                fastconverge = fastconverge, fastiter = fastiter, 
                validation = validation, 
                Sigmaispdef = Sigma1ispdef, Thetaispdef = Theta1ispdef))
  }
}

Wcodaloss = function(data, kappa = 1, lambda = NULL, nlambda = 20,  mu = 0, 
                     rho = 10, eps = 1e-8, method = 'approximate', 
                     tol = 1e-4, itermax = 1000, lambda.min.ratio = 1e-3, lambda.max.ratio = 10,
                     fast = TRUE, fastitermax = 1000, fasttol = 1e-4,
                     init = list(Theta = NULL, Sigma = NULL), 
                     validation = NULL, runlog = TRUE, 
                     onlyval = FALSE, show = TRUE){
  p = dim(data)[2]
  n = dim(data)[1]
  G = diag(p) - 1/p*matrix(1,p,p)
  A = G %*% (var(log(data))*(1-1/n)) %*% G
  if(!is.null(lambda)){
    nlambda = length(lambda)
  }
  if(is.null(lambda)){
    lambda.max <- max(max(A - diag(p)), -min(A - diag(p)))
    lambda.min <- lambda.min.ratio * lambda.max
    lambda <- exp(seq(log(lambda.max), log(lambda.min), length = nlambda))
    
  }
  cov_ret = list()
  icov_ret = list()
  path = list()
  bic = list()
  for(i in 1:nlambda){
    result = codaloss(data, kappa = kappa, lambda = lambda[i], pweight = NULL, mu = mu, 
                      rho = rho, eps = eps, method = method, 
                      tol = tol, itermax = itermax, 
                      fast = fast, fastitermax = fastitermax, fasttol = fasttol,
                      init = init, 
                      validation = validation, runlog = runlog, 
                      onlyval = onlyval, show = show)
    cov_ret[[i]] = result$Sigma
    icov_ret[[i]] = result$Theta
    path[[i]] = 1 * (icov_ret[[i]]!=0)
    bic = c(bic, n * norm(result$Sigma %*% result$Theta %*% result$Sigma -result$Sigma, type = '1') + log(n) * sum(path[[i]][upper.tri(path[[i]], diag = FALSE)] != 0))
    print(i)
  }
  return(list(Theta = icov_ret, Sigma = cov_ret, path = path, lambda = lambda, nlambda = lambda, bic = bic))
}













