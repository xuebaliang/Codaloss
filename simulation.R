starttime = Sys.time()
library(gtools)
library(MASS)
library(Matrix)
library(dplyr)
library(Rcpp)
library(ROCR)
library(gplots)
library(huge)
library(clime)


path0 = getwd()
setwd(path0)
source(paste0(path0, "\\gCoda-master\\R\\gcoda.R"))
path1 = paste0(path0, "\\SpiecEasi-master\\R")
setwd(path1)
for(f in list.files(path1)){
  source(f)
}
path2 = paste0(path0, "\\cdtrace-yuan")
setwd(path2)
for(f in list.files(path2)){
  source(f)
}
path3 = paste0(path0, "\\cdtrace-he")
setwd(path3)
for(f in list.files(path3)){
  source(f)
}

require(Rcpp)
path4 = paste0(path0, "\\codaloss-chen")
setwd(path4)
sourceCpp("symmatrix.cpp")
sourceCpp("update.cpp")
source("basical_coda.R")
source("updatefunc.R")
source("convenientfunc.R")
source("codalosscv.R")
source("codaloss.R")
setwd(path0)


geteval = function(Theta_real, Sigma_real, Theta_est, Sigma_est, p, data){
  
  Thetabias = norm(Theta_real - Theta_est, type = "F")
  Sigmabias = norm(Sigma_real - Sigma_est, type = "F")
  I = diag(p)
  Fm = I - 1/p
  varln = var(log(data))
  invbias = 1/2 * (norm(Theta_est %*% Sigma_est - I, type = "F") + 
                     norm(Sigma_est %*% Theta_est - I, type = "F"))
  codabias = norm(Fm %*% (varln - Sigma_est) %*% Fm, type = "F")
  ind = upper.tri(Fm,  diag = FALSE)
  require(ROCR)
  pred = prediction(predictions = abs(Theta_est[ind]), labels = 1*(Theta_real[ind]!=0))
  tpr = performance(pred, "tpr")@y.values[[1]]
  fpr = performance(pred, "fpr")@y.values[[1]]
  precision = performance(pred, "prec")@y.values[[1]]
  recall = performance(pred, "rec")@y.values[[1]]
  res = data.frame(tpr, fpr, precision, recall)
  res[is.na(res)] = 1
  AUC = performance(pred, "auc")@y.values[[1]]
  return(list(AUC = AUC, res = res, Thetabias = Thetabias, Sigmabias = Sigmabias, invbias = invbias, codabias = codabias))
}

###
simnetwork = function(p = 50, model = "random", seed = 1){
  set.seed(seed)
  require(MASS)
  Theta = diag(p)
  checkmat = diag(p)
  
  if(model == "random"){
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        if(rbinom(n = 1, size = 1, prob = 0.25)){
          val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * runif(1, 0.7, 0.9)
          Theta[i,j] = val; Theta[j,i] = val
        }
        checkmat[i,j] = 1
      }
    }
  }else if(model == "neighbor"){
    locmat = matrix(0, p, p)
    d  = runif(n = p, 0, 1)
    for(i in 1:p){
      x = abs(d[i] - d)
      loc = order(x)[2:14]
      locmat[i, loc] = 1; locmat[loc, i] = 1
    }
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        if(locmat[i,j] == 1){
          val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * runif(1, 0.7, 0.9)
          Theta[i, j] = val; Theta[j, i] = val
        }
        checkmat[i,j] = 1
      }
    }
  }else if(model == "band"){
    for(i in 1:(p-1)){
      val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * 
        runif(1, 0.7, 0.9)
      Theta[i,i+1] = val; Theta[i+1,i] = val
    }
    for(i in 1:(p-2)){
      val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * 
        runif(1, 0.7, 0.9)
      Theta[i,i+2] = val; Theta[i+2,i] = val
    }
    for(i in 1:(p-3)){
      val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * 
        runif(1, 0.7, 0.9)
      Theta[i,i+3] = val; Theta[i+3,i] = val
    }
    for(i in 1:(p-4)){
      val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * 
        runif(1, 0.7, 0.9)
      Theta[i,i+4] = val; Theta[i+4,i] = val
    }
    for(i in 1:(p-5)){
      val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * 
        runif(1, 0.7, 0.9)
      Theta[i,i+5] = val; Theta[i+5,i] = val
    }
    for(i in 1:(p-6)){
      val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) *
        runif(1, 0.7, 0.9)
      Theta[i,i+6] = val; Theta[i+6,i] = val
    }
    for(i in 1:(p-7)){
      val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) *
        runif(1, 0.7, 0.9)
      Theta[i,i+7] = val; Theta[i+7,i] = val
    }
    for(i in 1:(p-8)){
      val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) *
        runif(1, 0.7, 0.9)
      Theta[i,i+8] = val; Theta[i+8,i] = val
    }
  }else if(model == "hub"){
    hub = sort(sample(p, size = 3, replace = FALSE))
    nothub = (1:p)[-hub]
    
    for(i in hub){
      for(j in nothub){
        if(rbinom(n = 1, size = 1, prob = 0.8)){
          val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) *
            runif(n = 1, 0.7, 0.9)
          Theta[i,j] = val
          Theta[j,i] = val
        }
        checkmat[min(i,j),max(i,j)] = 1
      }
    }
    
    for(i in nothub){
      for(j in nothub[nothub > i]){
        if(rbinom(n = 1, size = 1, prob = 0.2)){
          val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) *
            runif(n = 1, 0.7, 0.9)
          Theta[i,j] = val
          Theta[j,i] = val
        }
        checkmat[i,j] = 1
      }
    }
  }else if(model == "block"){
    block = sample(p, size = p, replace = FALSE)
    blocksize = p / 5
    block = lapply(1:5, function(i){
      sort(block[(blocksize*(i-1)+1):(blocksize*i)])
    })
    
    for(i in 1:5){
      for(j in block[[i]]){
        for(k in block[[i]][block[[i]] > j]){
          if(rbinom(n = 1, size = 1, prob = 0.4)){
            val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) *
              runif(n = 1, 0.7, 0.9)
            Theta[k,j] = val
            Theta[j,k] = val
          }
          checkmat[j,k] = 1
        }
      }
    }
    
    for(i in 1:4){
      for(j in block[[i]]){
        for(k in unlist(lapply(i:5, function(l)block[[l]]))){
          if(rbinom(n = 1, size = 1, prob = 0.2)){
            val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) *
              runif(n = 1, 0.7, 0.9)
            Theta[k,j] = val
            Theta[j,k] = val
          }
          checkmat[min(j,k),max(j,k)] = 1
        }
      }
    }
  }else if(model == "scalefree"){
    start = sample(1:p, size = 12, replace = FALSE)
    for(i in start){
      for(j in start[start > i]){
        val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) *
          runif(n = 1, 0.7, 0.9)
        Theta[i,j] = val
        Theta[j,i] = val
        checkmat[i,j] = 1
      }
    }
    
    left = (1:p)[-start]
    while(length(start) < p){
      if(length(left)>1){
        i = sample(left, size = 1, replace = FALSE)
      }else{
        i = left
      }
      
      choose = sample(start, size = 6, replace = FALSE, 
                      prob = apply(Theta[start,]!=0, 1, sum)-1)
      
      for(j in choose){
        val = (2 * rbinom(n = 1, size = 1, prob = 0.5) - 1) * 
          runif(n = 1, 0.7, 0.9)
        Theta[i,j] = val
        Theta[j,i] = val
        checkmat[min(j,i),max(j,i)] = 1
      }
      
      start = c(start, i)
      left = (1:p)[-start]
    }
  }
  toadd = 0.01
  eps = 0
  eval = eigen(Theta, only.values = TRUE)$values
  if(any(eval <= 0)){
    eps = toadd + abs(min(eval))
    Theta = Theta + diag(eps, p)
    Theta = Theta / (eps + 1)
  }
  
  Sigma = solve(Theta)
  mu = c(rep(0, p))
  return(list(Theta = Theta, Sigma = Sigma, mu = mu,
              eps = eps, checkmat = checkmat))
}


onesim = function(Sigma, Theta, mu, p, n, seed, model, simple = TRUE){
  set.seed(seed)
  data = exp(mvrnorm(n = n, mu = mu, Sigma = Sigma))
  data = data / rowSums(data)
  clr.data = log(data) %*% (diag(p) - 1/p * matrix(1, p, p))
  
  ##se(gl)
  gl.res = sparseiCov(clr.data, method = 'glasso', lambda.min.ratio = 0.01, nlambda = 30)
  gl.res = icov.select(gl.res, criterion = 'ebic')
  gl.icov = as.matrix(gl.res$opt.icov)
  glval = geteval(Theta_real = Theta, Sigma_real = Sigma,
                  Theta_est = as.matrix(gl.icov), Sigma_est = solve(as.matrix(gl.icov)),
                  p = p, data = data)
  print("end the gl train")
  
  ###result of fang
  fang = gcoda(x = data, counts = F, lambda.min.ratio = 1e-3,
               nlambda = 30, ebic.gamma = 0.5)

  fangval = geteval(Theta_real = Theta, Sigma_real = Sigma,
                    Theta_est = fang$opt.icov, Sigma_est = solve(fang$opt.icov),
                    p = p, data = data)
  print("end the fang train ")
  
  ##result of yuan
  yuan = CompDtrace(clr.data = clr.data, lambda = NULL, rho = 0.8,
                    nlambda = 30, lambda.min.ratio = 0.0001,
                    lambda.max.ratio = 10)
  yuanval = geteval(Theta_real = Theta, Sigma_real = Sigma,
                    Theta_est = yuan$Theta, Sigma_est = solve(yuan$Theta),
                    p = p, data = data)
  print("end the yuan train ")

  ###result of he
  he = cdtrace_path(data = data, exact = TRUE, rho = 10, eps = 1e-8, lambda.min.ratio = 1e-4, nlambda = 30,
                    tol = 1e-4, itermax = 1000, fastitermax = 1000, fasttol = 1e-4)
  heval = geteval(Theta_real = Theta, Sigma_real = Sigma,
                  Theta_est = he$icov, Sigma_est = solve(he$icov),
                  p = p, data = data)
  print("end the he train ")

  ###my method
  #get weight
  nopenalty = codaloss(data = data, kappa = 1, lambda = 0, pweight = NULL,
                       mu = 0, rho = 10, eps = 1e-08, method = "approximate",
                       tol = 1e-03, itermax = 1000, fast = TRUE,
                       fastitermax = 1000, fasttol = 1e-03,
                       init = list(Theta = NULL, Sigma = NULL),
                       validation = NULL, runlog = TRUE,
                       onlyval = FALSE, show = FALSE)
  
  #choose lambda for weighted penalty
  pweight = 1 / abs(nopenalty$Theta)
  pweight = (p * p / sum(pweight)) * pweight * 100
  weightcv = codaloss_cv(data = data, folds = 5, foldid = NULL,
                         kappaseq = c(1), lambdaratio = 1e-05, maxnlambda = 30,
                         pweight = pweight, museq = c(0), seed = 1,
                         rho = 10, eps = 1e-08, method = "approximate", tol = 1e-03,
                         itermax = 1000, fast = TRUE, fastitermax = 1000, fasttol = 1e-03,
                         init = list(Theta = NULL, Sigma = NULL),
                         runlog = FALSE, onlyval = FALSE, show = FALSE)
  weightoptind = which.min(weightcv$cv_mat_summary$cvm)
  weightoptlambda = weightcv$cv_mat_summary$lambda[weightoptind]

  #weight result
  weightoptresult = codaloss(data = data, kappa = 1, lambda = weightoptlambda,
                             pweight = pweight, mu = 0, rho = 10, eps = 1e-08,
                             method = "approximate", tol = 1e-03, itermax = 1000,
                             fast = TRUE, fastitermax = 1000, fasttol = 1e-03,
                             init = list(Theta = NULL, Sigma = NULL),
                             validation = NULL, runlog = TRUE,
                             onlyval = FALSE, show = FALSE)
  weightval = geteval(Theta_real = Theta, Sigma_real = Sigma,
                      Theta_est = weightoptresult$Theta,
                      Sigma_est = weightoptresult$Sigma,
                      p = p, data = data)
  weight.icov = weightoptresult$Theta
  print("end the weight train")

  if(simple){
    Theta_gl = gl.icov
    Theta_fang = fang$opt.icov
    Theta_yuan = yuan$Theta
    Theta_he = he$icov
    Theta_weight = weightoptresult$Theta
    Sigma_weight = weightoptresult$Sigma
    metrictab = list(gl = glval, gcoda = fangval, cdtrace_Y = yuanval, cdtrace_H = heval, newmethod = weightval)
    AUC_vector = c(glval$AUC, fangval$AUC, yuanval$AUC, heval$AUC, weightval$AUC)
    thetabias_vector = c(glval$Thetabias, fangval$Thetabias, yuanval$Thetabias, heval$Thetabias, weightval$Thetabias)
    sigmabias_vector = c(fangval$Sigmabias, yuanval$Sigmabias, heval$Sigmabias, weightval$Sigmabias)
    invbias_vector = c(glval$invbias, fangval$invbias, yuanval$invbias, heval$invbias, weightval$invbias)
    codabias_vector = c(fangval$codabias, yuanval$codabias, heval$codabias, weightval$codabias)

    save(data, Theta, Sigma, metrictab, AUC_vector, thetabias_vector, sigmabias_vector, invbias_vector, codabias_vector, 
         Theta_gl, Theta_fang, Theta_yuan, Theta_he, Theta_weight, Sigma_weight,
         file = paste0("ROCR", model, "-p", p, "-n", n, "-seed", seed, ".rdata"))

  }
  return(list(AUC_vector = AUC_vector, thetabias_vector = thetabias_vector, sigmabias_vector = sigmabias_vector, invbias_vector = invbias_vector, codabias_vector = codabias_vector))
}

glassosim = function(Sigma, Theta, mu, p, n, seed, model, simple = TRUE){
  set.seed(seed)
  data = exp(mvrnorm(n = n, mu = mu, Sigma = Sigma))
  data = data / rowSums(data)
  glasso = huge(data, method = "glasso")
  glasso = huge.select(glasso, criterion = "ric")
  glasso = glasso$opt.icov
  glassoval = geteval(Theta_real = Theta, Sigma_real = Sigma,
                      Theta_est = glasso, Sigma_est = solve(glasso),
                      p = p, data = data)
  AUC_vector = c(glassoval$AUC)
  thetabias_vector = c(glassoval$Thetabias)
  save(data, Theta, Sigma, AUC_vector, thetabias_vector, glasso, 
       file = paste0("ROCR", model, "-glasso", "-p", p, "-n", n, "-seed", seed, ".rdata"))
  return(list(AUC_vector = AUC_vector, thetabias_vector = thetabias_vector))
}

climesim = function(Sigma, Theta, mu, p, n, seed, model, simple = TRUE){
  set.seed(seed)
  data = exp(mvrnorm(n = n, mu = mu, Sigma = Sigma))
  data = data / rowSums(data)
  re.clime = clime(data, standardize = FALSE)
  re.cv = cv.clime(re.clime)
  re.clime.opt = clime(data, standardize = FALSE)
  clime.icov = re.clime.opt$Omegalist[[1]]
  climeval = geteval(Theta_real = Theta, Sigma_real = Sigma,
                     Theta_est = clime.icov, Sigma_est = solve(clime.icov),
                     p = p, data = data)
  AUC_vector = c(climeval$AUC)
  thetabias_vector = c(climeval$Thetabias)
  save(data, Theta, Sigma, AUC_vector, thetabias_vector, clime.icov,
       file = paste0("ROCR", model, "-clime", "-p", p, "-n", n, "-seed", seed, ".rdata"))
  return(list(AUC_vector = AUC_vector, thetabias_vector = thetabias_vector))
}


n = 200
p = 50
AUC_result = c()
Thetabias_result = c()
for(model in c("random", "neighbor", "band", "hub", "block", "scalefree")){
  net = simnetwork(p = p, model = model, seed = 1)
  a = glassosim(Sigma = net$Sigma, Theta = net$Theta, mu = net$mu, p = p, n = n, seed = 1, model = model, simple = TRUE)
  b = climesim(Sigma = net$Sigma, Theta = net$Theta, mu = net$mu, p = p, n = n, seed = 1, model = model, simple = TRUE)
  c = onesim(Sigma = net$Sigma, Theta = net$Theta, mu = net$mu, p = p, n = n, seed = 1, model = model, simple = TRUE)
  AUC_result = rbind(AUC_result, c(model, n, p, round(a$AUC_vector, 4), round(b$AUC_vector, 4), round(c$AUC_vector, 4)))
  Thetabias_result = rbind(Thetabias_result, c(model, n, p, round(a$thetabias_vector, 4), round(b$thetabias_vector, 4), round(c$thetabias_vector, 4)))
  print(AUC_result)
  print(Thetabias_result)
}
colnames(AUC_result)= c("network", "n", "p", "glasso", "clime", "SPIEC-EASI(gl)", "gCoda", "CDtrace", "CDtr", "codaloss")
colnames(Thetabias_result) = c("network", "n", "p", "glasso", "clime", "SPIEC-EASI(gl)", "gCoda", "CDtrace", "CDtr", "codaloss")
print(AUC_result)
print(Thetabias_result)

