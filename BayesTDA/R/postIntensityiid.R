#' Posterior Intensity for persistence diagrams modeled as IID cluster Point Process
#' @description This function use the Bayesian inference obtained by characterizing persistence diagrams (PDs) as IID cluster point process. 
#' Posterior intensity can be computed from a prior intensity and a set of observed persistence diagrams. 
#' The inputs consists of two density functions estimated by using 'Wedge_Gaussian_Mixture' function: one is the prior intensity and another one is denoted as clutter.
#' The prior density will be depending on the knowledge we have about the underlying truth. The observed PDs can exhibit two features.
#' One which is believed to be associated to underlying truth and one is not anticipated by underlying truth or produce due to the noise. We call the later one surprise features.
#' Usually these densities are considered to have mean near the birth axis as they can be formed due to noise.  
#' @name postIntensityPoisson
#' @usage postIntensityPoisson(prior, clutter, Dy, sigma.y, alpha)
#' @param prior: prior density object. It is defined as a mixed Gaussian and use the class 'Wedge_Gaussian_Mixture'. 
#' The arguments (weights, means and sigmas) for the 'Wedge_Gaussian_Mixture' need to be predefined according to the prior knowledge.
#' see example for more details.
#' @param clutter: clutter density object. It is defined as a mixed Gaussian and use the class 'Wedge_Gaussian_Mixture'. 
#' These are the surprise features we encounter in the observed persistence diagrams. see example for more details. 
#' @param Dy: list of two element numerics representing points observed in a persistence diagram. Here we consider a fixed homological feature.
#' The coordinates (birth, death)/(birth, persistence) need to be a list of two element vectors. 
#' @param sigma.y: positive constant. variance coefficient of the likelihood density.This represents the degree of faith on the observed PDs representing the underlying truth.
#' @param alpha: The probablity of a feature in the prior will be detected in the observation.
#' @return The function 'Posterior' returns posterior intensity given prior and set of observed PDs using Bayesian framework, where PDs are characterized by poisson point process.
#' @details Required packages are TDA and mvtnorm.
#' @references Bayesian Inference for Persistent Homology, V Maroulas, F Nasrin, C Oballe, \url{https://arxiv.org/abs/1901.02034}


el_symmetric = function(values,k){
  if(length(values) == 0){
    return(1)
  } else{
    M = length(values)
    polynomial = as.numeric(poly.calc(values))
    return(((-1)^k)*polynomial[M+1]*polynomial[M-k+1])
  }
}

N_0 = 150 #size chosen must be larger than any expected cardinality

#Gammas in eqn (16)
#Inputs:  x, the point to evalues the posterior at
# Ys, list of K birth and death pairs from D_Y
# prob_dys, a numeric for the value of the probability for surprise D_ys points
# weights_Dx, the c_i^Dx weights for the Gaussian mixture for f_Dx as vector
# alpha, the probability of a point appearing, here taken as constant
# mu_Dx, list of N means for prior
# sig_Dx, vector of N std devs
# sig_Dyo, numeric for sigma*I for D_Y_o
# weights_Dys, list of weights for D_ys
# mu_Dys, list of means for D_ys
# sig_Dys, vector of values for covariance 
# prob_Dx, a numeric for the value of the probability for prior Dx points
# b, numeric equal to 0 or 1
Gamma = function(Ys,prob_dys,weights_Dx,alpha,mu_Dx,sig_Dx,sig_Dyo,weights_Dys,mu_Dys,sig_Dys,prob_Dx,b){
  
  K = length(Ys) 
  pdys = unlist(lapply(0:K,dbinom,size = N_0,prob=prob_dys)) 
  qy = lapply(Ys,function(y){mapply(function(mu,sig1,sig2){dmvnorm(y,mean=mu,
                                                                   sigma=(sig1+sig2)*diag(2))},mu_Dx,sig_Dx,MoreArgs = list(sig2 = sig_Dyo))}) 
  tilde_c = function(means,sigmas,weights){pmvnorm(lower = c(0,0),upper = c(Inf,Inf),mean=means,
                                                   sigma=sigmas*diag(2))[1]*weights} 
  tilde_c_Dx = (1-alpha)*sum(mapply(tilde_c,mu_Dx,sig_Dx,weights_Dx))
  f_Dys = unlist(lapply(Ys,wedge_gaussian_mixture,c_i = weights_Dys,means = mu_Dys,sigmas = sig_Dys))
  input_el_sym = mapply(function(q,fDys){alpha*t(weights_Dx)%*%q*(1/fDys)},qy,f_Dys)
  el_sym = c(1,unlist(lapply(1:K,el_symmetric,values = input_el_sym)))
  
  perm_coeff = function(n,j){
    if(n>=j){
      return(factorial(n)/factorial(n-j))
    }else{return(0)}
  } 
  
  bottom = integer(N_0+1)
  for(tau in 0:N_0){
    bot_sum = integer(min(K,tau)+1)
    for(k in 0:min(K,tau)){
      bot_sum[k+1] =  factorial(K-k)*perm_coeff(tau,k+b)*pdys[K-k+1]*
        (tilde_c_Dx^(tau-k-b))*el_sym[k+1]
    }
    bottom[tau+1] = dbinom(tau,N_0,prob_Dx)*sum(bot_sum)
  }
  return(sum(bottom))
}


posterior = function(x,Ys,prob_dys,weights_Dx,alpha,mu_Dx,sig_Dx,sig_Dyo,weights_Dys,mu_Dys,sig_Dys,prob_Dx){
  #######maybe later set some defaults for the input values if not specified
  first_term = (1-alpha)*wedge_gaussian_mixture(x,weights_Dx,mu_Dx,sig_Dx)*
    Gamma(Ys,prob_dys,weights_Dx,alpha,mu_Dx,sig_Dx,sig_Dyo,weights_Dys,mu_Dys,sig_Dys,prob_Dx,1)*
    (1/Gamma(Ys,prob_dys,weights_Dx,alpha,mu_Dx,sig_Dx,sig_Dyo,weights_Dys,mu_Dys,sig_Dys,prob_Dx,0))
  qy = lapply(Ys,function(y){mapply(function(mu,sig1,sig2){dmvnorm(y,mean=mu,
                                                                   sigma=(sig1+sig2)*diag(2))},mu_Dx,sig_Dx,MoreArgs = list(sig2 = sig_Dyo))}) 
  K = length(Ys) 
  to_sum = matrix(0,nrow = K,length(weights_Dx)) 
  for(j in 1:K){
    for(i in 1:length(weights_Dx)){
      to_sum[j,i]=Gamma(Ys[-j],prob_dys,weights_Dx,alpha,mu_Dx,sig_Dx,sig_Dyo,weights_Dys,mu_Dys,sig_Dys,prob_Dx,1)*
        alpha*qy[[j]][i]*weights_Dx[i]*(1/wedge_gaussian_mixture(Ys[[j]],c_i = weights_Dys,means = mu_Dys,sigmas = sig_Dys))*
        (1/Gamma(Ys,prob_dys,weights_Dx,alpha,mu_Dx,sig_Dx,sig_Dyo,weights_Dys,mu_Dys,sig_Dys,prob_Dx,0))*
        dmvnorm(x,mean = (sig_Dx[i]*Ys[[j]]+sig_Dyo*mu_Dx[[i]])/(sig_Dx[[i]]+sig_Dyo),
                sigma = (sig_Dyo*sig_Dx[[i]]/(sig_Dx[[i]]+sig_Dyo))*diag(2))
    }
  }
  return(first_term+sum(to_sum))
}