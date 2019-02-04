#' Restricted Gaussian Mixure density
#' @description Computes mixed Gaussian density restricted to wedge W. Require 'mvtnorm'package.
#' @usage Wedge_Gaussian_Mixture(weights,mean,sigma)
#' @name Wedge_Gaussian_Mixture
#' @param weights: a list of mixture weights
#' @param mean: a list of two element vector, means of the restricted Gaussian
#' @param sigma: a list of positive constants, sigmas in covariance matrices of the restricted Gaussian
#' @return The function 'Wedge_Gaussian_Mixture' returns numerical value of the mixed Gaussian density restricted to wedge W.
#' @examples # input the mean and constant sigma of covariance matrix
#' w <- Wedge_Gaussian_Mixture(weights = list(1), means = list(c(1,1)), sigmas= list(0.001))
#' w$evaluate(c(0.5,0.5))

#-----------
#IMPORTS
#-----------
require(mvtnorm) # for multivariate normal pdf and cdf
#-----------

#-----------
# CLASSES
#-----------

Wedge_Gaussian_Mixture <- setRefClass('Wedge_Gaussian_Mixture',

                                      fields=list(weights='list',means='list',sigmas='list'),
                                      # weights: list of c_i
                                      # means: list of u_i
                                      # sigmas: list of sigma_i
                                      # matching parameters should be at same position in lists

                                      methods=list(
                                        initialize=function(weights=list(),means=list(),sigmas=list()){
                                          .self$weights=weights
                                          .self$means=means
                                          .self$sigmas=sigmas
                                        },

                                        evaluate = function(x){
                                          # x: two element numeric
                                          # returns value of wedge gaussian mixture density at x

                                          dens.list = mapply(function(u,s){Wedge_Gaussian$new(mean=u,sigma=s)},.self$means,.self$sigmas,SIMPLIFY = FALSE)
                                          dens.values = lapply(dens.list,function(y){y$evaluate(x)})
                                          prod = mapply(function(c,v){c*v},.self$weights,dens.values,SIMPLIFY=FALSE)
                                          return(Reduce('+',prod))
                                        }
                                      )
)
