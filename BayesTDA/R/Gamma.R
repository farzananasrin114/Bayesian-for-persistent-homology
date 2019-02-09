#'  Gamma
#' @name Gamma
#' @param x: two element vector. The posterior intensity is computed at x.
#' @param Dy: list of two element vectors representing points observed in a persistence diagram. Here we consider a fixed homological feature.
#' The coordinates (birth, death)/(birth, persistence) need to be a list of two element vectors. 
#' @param alpha: The probability of a feature in the prior will be detected in the observation.
#' @param prob.prior: The prior cardinality is defined as a binomial and prob.prior is the probability term, i.e., prior cardinality ~ B(Nmax,prob.prior). 
#' @param weight.prior: a list of mixture weights for the prior density object. This parameter will be an input for 'Wedge_Gaussian_Mixture' function to estimate prior density.
#' @param mean.prior: a list of two element vector, means of the prior density object. This parameter will be an input for 'Wedge_Gaussian_Mixture' function to estimate prior density.
#' @param sigma.prior: a list of positive constants, sigmas in covariance matrices of the prior density object. This parameter will be an input for 'Wedge_Gaussian_Mixture' function to estimate prior density.  
#' @param sigma.Dyo: positive constant. variance coefficient of the likelihood density. This represents the degree of faith on the observed PDs representing the underlying truth.
#' @param prob.clutter: The clutter cardinality is defined as a binomial and prob.clutter is the probability term, i.e., clutter cardinality ~ B(Nmax,clutter). 
#' @param weights.clutter: a list of mixture weights for the clutter density object. This parameter will be an input for 'Wedge_Gaussian_Mixture' function to estimate clutter density.
#' @param mean.clutter: a list of two element vector, means of the clutter density object. This parameter will be an input for 'Wedge_Gaussian_Mixture' function to estimate clutter density.
#' @param sigma.clutter: a list of positive constants, sigmas in covariance matrices of the clutter density object. This parameter will be an input for 'Wedge_Gaussian_Mixture' function to estimate clutter density.  
#' @param Nmax: The maximum number of points at which the posterior cardinality will be truncated, i.e., {p_n} is computed for n=0 to n-Nmax. Also, this Nmax will be used to define prior and clutter cardinality. 
#' So, it should be large enough so that all Binomials involved in the computation make sense
#' @param b: 0 or 1
#' @keywords internal
#' @return numeric

Gamma = function(Dy,alpha,prob.prior,weight.prior,mean.prior,sigma.prior,sigma.Dyo,prob.clutter,weights.clutter,mean.clutter,sigma.clutter,Nmax,b){
  
  if( alpha ==  1)
    stop('alpha can not be 1')
  else
    
    
  K = length(Dy) 
  pdys = unlist(lapply(0:K,dbinom,size = Nmax,prob=prob.clutter)) 
  qy = lapply(Dy,function(y){mapply(function(mu,sig1,sig2){dmvnorm(y,mean=mu,
                                                                   sigma=(sig1+sig2)*diag(2))},mean.prior,sigma.prior,MoreArgs = list(sig2 = sigma.Dyo))}) 
  tilde_c = function(means,sigmas,weights){pmvnorm(lower = c(0,0),upper = c(Inf,Inf),mean=means,
                                                   sigma=sigmas*diag(2))[1]*weights} 
  tilde_c_Dx = (1-alpha)*sum(mapply(tilde_c,mean.prior,sigma.prior,weight.prior))
  f_Dys = unlist(lapply(Dy,Wedge_Gaussian_Mixture,weights = weights.clutter,means = mean.clutter,sigmas = sigma.clutter))
  input_el_sym = mapply(function(q,fDys){alpha*t(weight.prior)%*%q*(1/fDys)},qy,f_Dys)
  el_sym = c(1,unlist(lapply(1:K,el_symmetric,values = input_el_sym)))
  
  perm_coeff = function(n,j){
    if(n>=j){
      return(factorial(n)/factorial(n-j))
    }else{return(0)}
  } 
  
  bottom = integer(Nmax+1)
  for(tau in 0:Nmax){
    bot_sum = integer(min(K,tau)+1)
    for(k in 0:min(K,tau)){
      bot_sum[k+1] =  factorial(K-k)*perm_coeff(tau,k+b)*pdys[K-k+1]*
        (tilde_c_Dx^(tau-k-b))*el_sym[k+1]
    }
    bottom[tau+1] = dbinom(tau,Nmax,prob.prior)*sum(bot_sum)
  }
  return(sum(bottom))
}
