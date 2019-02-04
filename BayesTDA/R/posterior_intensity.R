#' Posterior Intensity from Baysian Framework
#' @description Posterior intensity can be computed from a prior intensity and a set of observed persistence diagrams
#' @name Posterior_Intensity
#' @usage Posterior_Intensity(prior,clutter,Dy,sigma.y,alpha)
#' @param prior: prior density object. It is defined as a mixed Gaussian and use the class 'Wedge_Gaussian_Mixture'. 
#' The arguments (weights, means and sigmas) for the 'Wedge_Gaussian_Mixture' need to be predefined according to the prior knowledge.
#' see example for more details.
#' @param clutter: clutter density object. It is defined as a mixed Gaussian and use the class 'Wedge_Gaussian_Mixture'. 
#' These are the surprise features we encounter in the observed persistence diagrams. see example for more details. 
#' @param Dy: list of two element numerics representing point observed in a persistence diagram. Here we consider a fixed homological feature.
#' The coordinates (birth, death)/(birth, persistence) need to converted as a list of two element vectors. 
#' @param sigma.y: variance coefficient of the likelihhod function. This is a posistive constant. 
#' @param alpha: The probablity of a feature in the prior will be detected in the observation.
#' @return The function 'Posterior' returns posterior intensity given prior and set of observed pd using Bayesian framework by characterizing pd by poisson point process
#' @details Require packages are TDA and mvtnorm.
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
#' clut.mean = c(0.5,0)
#' clut.noise = 1
#' clut.weight = 1
#' # these parameters are required to define the object prior
#' inf.mean = c(0.5,1.2)
#' inf.noise = 0.2
#' inf.weight = 1
#' 
#'  values = seq(from = 0, to = 2, by = 0.05)
#'  grid = expand.grid(values,values)
#'  
#'  ## next we utilize the 'Wedge_Gaussian_Mixture' class to define prior and clutter over the grid
#'  inf.prior = Wedge_Gaussian_Mixture(means = list(inf.mean),weights=list(inf.weight),sigmas = list(inf.noise)) # prior
#'  clut = Wedge_Gaussian_Mixture(means = list(clut.mean),weights = list(clut.weight),sigmas = list(clut.noise))
#'  inf.dens = apply(grid,1,inf.prior$evaluate)
#'  
#'  ## next we define the likelihhod from the observed pds. Here we consider a dataset from the same unit circle that is perturbed by a gaussian noise.
#'  
#'    pri = inf.prior ## defining the parameter prior in the function
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
#'    ## next we compute the posterior over the same grid
#'    posterior = Posterior$new(prior = pri, Dy = dy, clutter = clut, sigma.y = sy, alpha = alpha)
#'    values = seq(from = 0, to = 2, by = 0.05)
#'    grid = expand.grid(values,values)
#'    intens = apply(grid,1,posterior$evaluate)
#'    intens = (1-alpha)*inf.dens + alpha*intens ## compute posterior intensity
#'    post.df = data.frame(Birth = grid[,1],Persistence = grid[,2], Intensity = normalized_intens)
#'    
#'    ##plot the posterior intensity
#'    ## require(gplot2)
#'    g.post = ggplot(data = post.df, aes(x = Birth, y = Persistence)) + geom_tile(aes(fill=Intensity)) + coord_equal()
#'    g.post = g.post + scale_fill_gradient2(low = "blue", high = "red", mid = "yellow",midpoint = 0.5, limit = c(0,1))
#'    g.post = g.post + geom_point(data = circle.df, aes(x = Birth, y= Persistence), pch = 19, cex=3, color = "Green")
#'    g.post = g.post+ labs(title='Posterior',x='Birth',y='Persistence', fill = "") + theme_grey(base_size=14) +  theme(plot.title = element_text(hjust = 0.5))

#'    print(g.post)

Posterior_Intensity <- setRefClass('Posterior_Intensity',
                         fields = list(prior = 'Wedge_Gaussian_Mixture',clutter='Wedge_Gaussian_Mixture',Dy='list',sigma.y='numeric',alpha='numeric',weight.scalars='matrix',posterior.means='array',posterior.sigmas='matrix',Qs = 'matrix', fast.eval = 'logical'),
                         # prior: prior density object
                         # clutter: clutter density object (C_Y density in Oballe Maroulas 2018)
                         # Dy: list of two element numerics representing point observed in a pd
                         # sigma.y: stochastic kernel sigma
                         # weight.scalars: components of w^y_i(x) that do not depend on x
                         # Qs: Q values in Oballe, Maroulas 2018
                         # posterior.means: posterior means
                         # posterior.sigmas: posterior sigmas
                         # i,j entries of parameter matrices (or [i,,j] entry in case of posterior.means) are parameter values corresponding to means[[i]], y[[j]]
                         # fast.eval: a flag indicating whether to carry out integration when computing posterior weights or approximate them with 1

                         methods = list(
                           initialize = function(prior,clutter,Dy,sigma.y,alpha){
                             .self$prior = prior
                             .self$clutter = clutter
                             .self$Dy = Dy
                             .self$sigma.y = sigma.y
                             .self$alpha = alpha
                             .self$fast.eval = FALSE

                             weight.scalars = matrix(nrow=length(prior$means),ncol=length(Dy))
                             Qs = matrix(nrow=length(prior$means),ncol=length(Dy))
                             posterior.means = array(NA,dim=c(length(prior$means),2,length(Dy)))
                             posterior.sigmas = matrix(nrow=length(prior$means),ncol=length(Dy))

                             for(i in 1:length(prior$means)){
                               for(j in 1:length(Dy)){
                                 weight.scalars[i,j] = prior$weights[[i]]*dmvnorm(Dy[[j]],mean=prior$means[[i]],sigma=(sigma.y+prior$sigmas[[i]])*diag(2))
                                 posterior.means[i,,j] = (prior$sigmas[[i]]*Dy[[j]]+sigma.y*prior$means[[i]])/(prior$sigmas[[i]]+sigma.y)
                                 posterior.sigmas[i,j] = (sigma.y*prior$sigmas[[i]])/(prior$sigmas[[i]]+sigma.y)

                                 if(!fast.eval){Qs[i,j] = pmvnorm(lower=c(0,0),upper=c(Inf,Inf),mean=posterior.means[i,,j],sigma=posterior.sigmas[i,j]*diag(2))
                                 }
                                 else{Qs[i,j] = 1}
                               }
                             }

                             .self$weight.scalars = weight.scalars
                             .self$Qs = Qs
                             .self$posterior.means = posterior.means
                             .self$posterior.sigmas = posterior.sigmas
                           }
                           ,

                           evaluate = function(x){
                             # evaluate posterior intensity density at x
                             # x: two element numeric representing element in TW
                             # returns a numeric

                             # compute w^y_i

                             w.mat = .self$weight.scalars

                             # compute denominators for C^y_i
                             noise.vals = rep(0,length(.self$Dy))
                             weight.sums = rep(0,length(.self$Dy))

                             for(j in 1:length(.self$Dy)){
                               noise.vals[j] = .self$clutter$evaluate(.self$Dy[[j]])
                               ws = w.mat[,j]
                               Qs = .self$Qs[,j]
                               weight.sums[j] = alpha*sum(ws*Qs)
                             }

                             # compute value of density in a loop
                             prob = 0
                             for(i in 1:length(.self$prior$means)){
                               for(j in 1:length(.self$Dy)){
                                 C = w.mat[i,j]/(noise.vals[j]+weight.sums[j])
                                 N = dmvnorm(x,mean=.self$posterior.means[i,,j],sigma=.self$posterior.sigmas[i,j]*diag(2))
                                 val = C*N
                                 prob = prob + val
                               }
                             }

                             return(prob)
                           }
                         )

)
