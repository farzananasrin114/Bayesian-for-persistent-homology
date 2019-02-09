#' Posterior Intensity for persistence diagrams modeled as Poisson Point Process
#' @description This function use the Bayesian inference obtained by characterizing persistence diagrams (PDs) as Poisson point process. 
#' Posterior intensity can be computed from a prior intensity and a set of observed persistence diagrams. 
#' The inputs consists of two density functions estimated by using \code{Wedge_Gaussian_Mixture} function: one is the prior intensity and another one is denoted as clutter.
#' The prior density will be depending on the knowledge we have about the underlying truth. The observed PDs can exhibit two features.
#' One which is believed to be associated to underlying truth and one is not anticipated by underlying truth or produce due to the noise. We call the later one surprise features.
#' Usually these densities are considered to have mean near the birth axis as they can be formed due to noise.  
#' @name postIntensityPoisson
#' @usage postIntensityPoisson(x,Dy,alpha,weight.prior,mean.prior,sigma.prior,sigma.y,weights.clutter,mean.clutter,sigma.clutter)
#' @param x: two element vector. The posterior intensity is computed at x.
#' @param Dy: list of two element vectors representing points observed in a persistence diagram. Here we consider a fixed homological feature.
#' The coordinates (birth, death)/(birth, persistence) need to be a list of two element vectors. 
#' @param weight.prior: a vector of mixture weights for the prior density object. This parameter will be an input for \code{Wedge_Gaussian_Mixture} function to estimate prior density.
#' @param mean.prior: a list of two element vector, means of the prior density object. This parameter will be an input for \code{Wedge_Gaussian_Mixture} function to estimate prior density.
#' @param sigma.prior: a vector of positive constants, sigmas in covariance matrices of the prior density object. This parameter will be an input for \code{Wedge_Gaussian_Mixture} function to estimate prior density.  
#' @param weights.clutter: a vector of mixture weights for the clutter density object. This parameter will be an input for \code{Wedge_Gaussian_Mixture} function to estimate clutter density.
#' @param mean.clutter: a list of two element vector, means of the clutter density object. This parameter will be an input for \code{Wedge_Gaussian_Mixture} function to estimate clutter density.
#' @param sigma.clutter: a vector of positive constants, sigmas in covariance matrices of the clutter density object. This parameter will be an input for \code{Wedge_Gaussian_Mixture} function to estimate clutter density.  
#' @param sigma.y: positive constant. variance coefficient of the likelihood density.This represents the degree of faith on the observed PDs representing the underlying truth.
#' @param alpha: The probablity of a feature in the prior will be detected in the observation.
#' @return The function \code{postIntensityPoisson} returns posterior intensity given prior and set of observed PDs using Bayesian framework, where PDs are characterized by poisson point process.
#' @details Required packages are \code{TDA} and \code{mvtnorm}.
#' @references Bayesian Inference for Persistent Homology, V Maroulas, F Nasrin, C Oballe, \url{https://arxiv.org/abs/1901.02034}

#' @examples # sample data created from a unit circle to define prior
#' t = seq(from=0, to = 1, by = 0.01)
#' x = cos(2*pi*t)
#' y = sin(2*pi*t)
#' coord.mat = cbind(x,y)
#' 
#' # persistence diagram using ripsDiag of TDA package.  
#' circle.samp = coord.mat[sample(1:nrow(coord.mat),50),]
#' pd.circle = ripsDiag(circle.samp,maxscale=3,maxdimension = 1)
#' circle.diag = matrix(pd.circle$diagram,ncol=3)
#' circle.holes = matrix(circle.diag[which(circle.diag[,1]==1),],ncol=3)
#' # this is required if we are working on a tilted wedge
#' circle.tilted = cbind(circle.holes[,2],circle.holes[,3]-circle.holes[,2])
#' circle.df = data.frame(Birth = circle.tilted[,1], Persistence = circle.tilted[,2]) ## create the dataframe consisting of pd coordinates.
#' 
#' alpha = 1 
#' # these parameters are required to define the object clutter
#' clut.mean = list(c(0.5,0))
#' clut.noise = 1
#' clut.weight = 1
#' # these parameters are required to define the object prior
#' inf.mean = list(c(0.5,1.2))
#' inf.noise = 0.2
#' inf.weight = 1
#'  
#'  ## next we define the likelihhod from the observed pds. Here we consider a dataset from the same unit circle that is perturbed by a gaussian noise.
#'  
#'    sy = 0.1 #the variance of the likelihood density function
#'    circle.noise = 0.01 ## the variance coefficient of noise in this dataset 
#'    noise = rmvnorm(nrow(coord.mat),mean = c(0,0), sigma = circle.noise*diag(2))
#'    coord.mat = coord.mat + noise ## gives us the dataset which we will consider as observation
#'    
#'    ## following section creates observed PD
#'    circle.samp = coord.mat[sample(1:nrow(coord.mat),50),]
#'    pd.circle = ripsDiag(circle.samp,maxscale=3,maxdimension = 1)
#'    circle.diag = matrix(pd.circle$diagram,ncol=3)
#'    circle.holes = matrix(circle.diag[which(circle.diag[,1]==1),],ncol=3)
#'    circle.tilted = cbind(circle.holes[,2],circle.holes[,3]-circle.holes[,2])
#'    circle.df = data.frame(Birth = circle.tilted[,1], Persistence = circle.tilted[,2]) 
#'    dy = lapply(1:nrow(circle.df),function(x){unlist(c(circle.df[x,][1],circle.df[x,][2]))}) ## creates the parameter Dy
#'    
#'    ## next we compute the posterior over a grid
#'    values = seq(from = 0, to = 2, by = 0.05)
#'    grid = expand.grid(values,values)
#'    intens = apply (grid,1,postIntensityPoisson,dy,alpha,inf.weight,inf.mean,inf.noise,sy,clut.weight,clut.mean,clut.noise)
#'    post.df = data.frame(Birth = grid[,1],Persistence = grid[,2], Intensity = intens)
#'    
#'    ##plot the posterior intensity
#'    ## require(gplot2)
#'    g.post = ggplot(data = post.df, aes(x = Birth, y = Persistence)) + geom_tile(aes(fill=Intensity)) + coord_equal()
#'    g.post = g.post + geom_point(data = circle.df, aes(x = Birth, y= Persistence), pch = 19, cex=3, color = "Green")
#'    g.post = g.post+ labs(title='Posterior',x='Birth',y='Persistence', fill = "") + theme_grey(base_size=14) +  theme(plot.title = element_text(hjust = 0.5))

#'    print(g.post)

postIntensityPoisson = function(x,Dy,alpha,weight.prior,mean.prior,sigma.prior,sigma.y,weights.clutter,mean.clutter,sigma.clutter){
  
  first_term = (1-alpha)*Wedge_Gaussian_Mixture(x,weight.prior,mean.prior,sigma.prior)
  qy = lapply(Dy,function(y){mapply(function(mu,sig1,sig2){dmvnorm(y,mean=mu,
                                                                   sigma=(sig1+sig2)*diag(2))},mean.prior,sigma.prior,MoreArgs = list(sig2 = sigma.y))}) 
  
  w = mapply(function(q,fDys){t(weight.prior)%*%q},qy)
  
  tilde_c = function(means,sigmas){pmvnorm(lower = c(0,0),upper = c(Inf,Inf),mean=means,
                                           sigma=sigmas*diag(2))}
  Q = mapply(tilde_c,mean.prior,sigma.prior)
  f_Dys = unlist(lapply(Dy,Wedge_Gaussian_Mixture, weights= weights.clutter,means = mean.clutter,sigmas = sigma.clutter))
  
  K = length(Dy)
  to_sum = matrix(0,nrow = K,length(weight.prior)) 
  for(j in 1:K){
    for(i in 1:length(weight.prior)){
      to_sum[j,i]=w[[j]][i]* (1/(f_Dys[[j]]+alpha*sum(t(w)*Q)))*
        dmvnorm(x,mean = (sigma.prior[i]*Dy[[j]]+sigma.y*mean.prior[[i]])/(sigma.prior[[i]]+sigma.y),
                sigma = (sigma.y*sigma.prior[[i]]/(sigma.prior[[i]]+sigma.y))*diag(2))
    }
  }
  
  return(first_term+sum(to_sum))
}
