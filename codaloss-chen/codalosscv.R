codaloss_cv = function(data, folds = 5, foldid = NULL, kappaseq = c(1), 
                       lambdaratio = 1e-4, maxnlambda = 40, pweight = NULL,
                       museq = c(0), seed = 1, rho = 10, eps = 1e-8, 
                       method = 'approximate', tol = 1e-4, itermax = 1000, 
                       fast = TRUE, fastitermax = 1000, fasttol = 1e-4,
                       init = list(Theta = NULL, Sigma = NULL),
                       runlog = FALSE, onlyval = TRUE, show = FALSE){
  checkparmcv(data = data, folds = folds, foldid = foldid, 
              kappaseq = kappaseq, lambdaratio = lambdaratio,
              maxnlambda = maxnlambda, pweight = pweight,
              museq = museq, seed = seed, rho = rho, eps = eps,
              method = method, tol = tol, itermax = itermax, 
              fast = fast, fastitermax = fastitermax, fasttol = fasttol, 
              init = init, runlog = runlog, onlyval = onlyval, show = show)
  require(dplyr)
  set.seed(seed)
  
  ####generate fold id
  n = dim(data)[1]
  p = dim(data)[2]
  
  if(is.null(foldid)){
    a = n %% folds; b = n %/% folds
    foldid = rep(NA, times = n)
    for(i in 1:folds){
      foldid[sample(which(is.na(foldid)), size = b, replace = FALSE)] = i
    }
    if(a > 0){
      foldid[which(is.na(foldid))] = sample(1:folds, size = a, replace = FALSE)
    }
  }
  
  
  ####generate parameter list
  lambdastart = log(max(abs(var(log(data)) - diag(p))) * lambdaratio)
  lambdagap = - log(lambdaratio) / 20
  lambdaseq = c(0, exp(seq(from = lambdastart, by = lambdagap,
                           length = maxnlambda)))
  parmls = list()
  num = 1
  for(kappa in kappaseq)
    for(mu in museq)
      for(lambda in lambdaseq)
        for(fold in 1:folds){
          parmls[[num]] = list(kappa = kappa, mu = mu, lambda = lambda, fold = fold)
          num = num + 1
        }
  
  result = list()
  num = 0
  goon = TRUE
  for(parm in parmls){
    if(goon){
      if(show){
        cat("-------------------------------------------------------------------\n")
        cat("kappa:", parm$kappa, "  lambda:", parm$lambda, 
            "  mu:", parm$mu, "  fold:", parm$fold, "\n")
      }
      
      ret = codaloss(data = data[foldid != parm$fold,], 
                     kappa = parm$kappa, lambda = parm$lambda, pweight = pweight,
                     mu = parm$mu, rho = rho, eps = eps, method = method, 
                     tol = tol, itermax = itermax, init = init, 
                     fast = fast, fastitermax = fastitermax, fasttol = fasttol,
                     validation = data[foldid == parm$fold,],
                     runlog = runlog, onlyval = onlyval, show = show)
      ret$parm = parm
      if(show)cat("-------------------------------------------------------------------\n")
      num = num + 1; result[[num]] = ret
      
      if(parm$fold == folds & ret$nonzeronum <= p)goon = FALSE
    }
  }

  
  cv_mat = t(sapply(1:num, function(i)unlist(parmls[[i]])))
  cv_mat = data.frame(cv_mat, validation = sapply(result, function(x)x$validation),
                      nonzeronum = sapply(result, function(x)x$nonzeronum))
  rownames(cv_mat) = NULL
  
  cv_mat_summary = group_by(cv_mat, kappa, mu, lambda)
  cv_mat_summary = data.frame(summarise(cv_mat_summary, cvm = mean(validation), 
                                        cvsd = sd(validation), nonzeronum = mean(nonzeronum)))
  cv_mat_summary$cvu = cv_mat_summary$cvm + cv_mat_summary$cvsd
  cv_mat_summary$cvl = cv_mat_summary$cvm - cv_mat_summary$cvsd
  
  return(list(cv_mat_summary = cv_mat_summary, 
              cv_mat = cv_mat,
              result = result))
}