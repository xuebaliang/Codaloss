#####main update function
codaloss_optim = function(lTheta, lTheta1, lTheta2, 
                          lSigma, lSigma1, lSigma2, 
                          llambda11, llambda12, llambda21, llambda22,
                          kFSlnxF, p, rhoI, tworhoI, 
                          loc, diagloc, rho, mu, kappa, 
                          Wlambda, eps = 1e-8,
                          method, UF, DF, invm){
  #update Theta
  nTheta = G_sl(A = multip_symmat(lSigma1, loc = loc, p = p) + tworhoI, 
                B = lSigma1 + rho * (lTheta1 + lTheta2) - llambda11 - llambda12,
                p = p)
  
  #update Theta1
  nTheta1 = S_sl(A = nTheta + llambda11 / rho, Wlambda = Wlambda / rho, p = p)
  
  #update Theta2
  nTheta2 = R_sl(A = nTheta + llambda12 / rho, p = p, eps = eps)
  
  
  #update Sigma
  if(method == 'approximate'){
    nSigma = T_sl(A = kFSlnxF + rho * (lSigma1 + lSigma2) - llambda21 - llambda22, 
                  kappa = kappa / 2, mu = mu, UF = UF, DF = DF, 
                  p = p, last_X = lSigma, diagloc = diagloc)
  }else if(method == 'exact'){
    nSigma = T_sl_exact(invm = invm, 
                        A = kFSlnxF + rho * (lSigma1 + lSigma2) - llambda21 - llambda22, 
                        loc = loc, p = p)
  } 
  
  #update Sigma1
  nSigma1 = G_sl(A = multip_symmat(nTheta, loc = loc, p = p) + rhoI, 
                 B = nTheta + rho * nSigma + llambda21, p = p)
  
  #update Sigma2
  nSigma2 = R_sl(A = nSigma + llambda22 / rho, p = p, eps = eps)
  
  
  #update lambda11,lambda12,lambda21,lambda22
  nlambda11 = llambda11 + rho * (nTheta - nTheta1)
  nlambda12 = llambda12 + rho * (nTheta - nTheta2)
  nlambda21 = llambda21 + rho * (nSigma - nSigma1)
  nlambda22 = llambda22 + rho * (nSigma - nSigma2)
  
  return(list(Theta = nTheta, Theta1 = nTheta1, Theta2 = nTheta2,
              Sigma = nSigma, Sigma1 = nSigma1, Sigma2 = nSigma2,
              lambda11 = nlambda11, lambda12 = nlambda12,
              lambda21 = nlambda21, lambda22 = nlambda22))
}

codaloss_optim_fast = function(lTheta, lTheta1,  
                               lSigma, lSigma1, 
                               llambda11, llambda21,
                               kFSlnxF, p, rhoI, 
                               loc, diagloc, rho, mu, kappa,
                               Wlambda, eps = 1e-8,
                               method, UF, fastDF, fastinvm){
  #update Theta
  nTheta = G_sl(A = multip_symmat(lSigma1, loc = loc, p = p) + rhoI, 
                B = lSigma1 + rho * lTheta1 - llambda11,
                p = p)
  
  #update Theta1
  nTheta1 = S_sl(A = nTheta + llambda11 / rho, Wlambda = Wlambda / rho, p = p)
  
  
  #update Sigma
  if(method == 'approximate'){
    nSigma = T_sl(A = kFSlnxF + rho * lSigma1 - llambda21, 
                  kappa = kappa / 2, mu = mu, UF = UF, DF = fastDF, p = p,
                  last_X = lSigma, diagloc = diagloc)
  }else if(method == 'exact'){
    nSigma = T_sl_exact(invm = fastinvm,
                        A = kFSlnxF + rho * lSigma1 - llambda21, 
                        loc = loc, p = p)
  } 
  
  #update Sigma1
  nSigma1 = G_sl(A = multip_symmat(nTheta, loc = loc, p = p) + rhoI, 
                 B = nTheta + rho * nSigma + llambda21, p = p)
  
  
  #update lambda11,lambda12,lambda21,lambda22
  nlambda11 = llambda11 + rho * (nTheta - nTheta1)
  nlambda21 = llambda21 + rho * (nSigma - nSigma1)
  
  return(list(Theta = nTheta, Theta1 = nTheta1,
              Sigma = nSigma, Sigma1 = nSigma1, 
              lambda11 = nlambda11,
              lambda21 = nlambda21))
}