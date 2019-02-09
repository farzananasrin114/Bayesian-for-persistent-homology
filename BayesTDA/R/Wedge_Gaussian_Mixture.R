#' Restricted Gaussian Mixure density
#' @description Computes mixed Gaussian density restricted to wedge W. Require \code{mvtnorm} package.
#' @usage Wedge_Gaussian_Mixture(x, weights, means, sigmas)
#' @name Wedge_Gaussian_Mixture
#' @param x: two element vector. These are the points to evaluate density at.
#' @param weights: a vector of the mixture weights
#' @param means: a list of two element vector where the vectors are the means of the restricted Gaussian
#' @param sigmas: a vector of positive constants, sigmas in covariance matrices of the restricted Gaussian
#' @return The function \code{Wedge_Gaussian_Mixture} returns numerical value of the mixed Gaussian density restricted to wedge W.
#' @examples # input the mean and constant sigma of covariance matrix
#' w <- Wedge_Gaussian_Mixture(x = c(1,2), weights = 1, means = list(c(1,1)), sigmas= 0.1)
#' @export

Wedge_Gaussian_Mixture = function(x, weights, means, sigmas){
  to_sum = weights*mapply(Wedge_Gaussian,means,sigmas,MoreArgs = list(z = x))
  return(sum(to_sum))
}



# Wedge_Gaussian_Mixture <- setRefClass('Wedge_Gaussian_Mixture',
# 
#                                       fields=list(weights='list',means='list',sigmas='list'),
#                                       # weights: list of c_i
#                                       # means: list of u_i
#                                       # sigmas: list of sigma_i
#                                       # matching parameters should be at same position in lists
# 
#                                       methods=list(
#                                         initialize=function(weights=list(),means=list(),sigmas=list()){
#                                           .self$weights=weights
#                                           .self$means=means
#                                           .self$sigmas=sigmas
#                                         },
# 
#                                         evaluate = function(x){
#                                           # x: two element numeric
#                                           # returns value of wedge gaussian mixture density at x
# 
#                                           dens.list = mapply(function(u,s){Wedge_Gaussian$new(mean=u,sigma=s)},.self$means,.self$sigmas,SIMPLIFY = FALSE)
#                                           dens.values = lapply(dens.list,function(y){y$evaluate(x)})
#                                           prod = mapply(function(c,v){c*v},.self$weights,dens.values,SIMPLIFY=FALSE)
#                                           return(Reduce('+',prod))
#                                         }
#                                       )
# )
