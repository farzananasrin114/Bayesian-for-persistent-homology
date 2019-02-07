#' Restricted Gaussian density
#' The function Wedge_Gaussian computes Gaussian density restricted to wedge W.
#' @name Wedge_Gaussian
#' @usage Wedge_Gaussian(mean, sigma)
#' @param mean: two element vector, mean of the restricted Gaussian
#' @param sigma: positive constant, sigma in covariance matrix of the restricted Gaussian
#' @details The function 'Wedge_Gaussian' computes Gaussian density restricted to wedge W. 
#' For a death vs. birth persistence diagram wedge is defined as the two dimensional coordinate (b,d) such that d>b>0. 
#' For a tilted represetation, i.e., persistence vs. birth persistence diagram wedge is defined as T(b,d) = (b, d-b).
#' @return The function 'Wedge_Gaussian' returns numerical value of the  Gaussian density restricted to wedge W.
#' @description The function 'Wedge_Gaussian' computes Gaussian density restricted to wedge W. Require 'mvtnorm' package for computing multivariate normal pdf and cdf.
#' @examples # input the mean and constant sigma of covariance matrix
#' w <- Wedge_Gaussian(mean = c(1,1), sigma = 0.001)
#' w$evaluate(c(0.5,0.5))

#-----------
#IMPORTS
#-----------
require(mvtnorm) # for multivariate normal pdf and cdf
#-----------

#-----------
# CLASSES
#-----------

Wedge_Gaussian <- setRefClass('Wedge_Gaussian',

                              fields=list(mean='numeric',sigma='numeric'),
                              # mean: two element vector, mean of the modified Gaussian from Oballe, Maroulas 2018
                              # sigma: positive constant, sigma in covariance matrix from Oballe, Maroulas 2018


                              methods = list(
                                initialize = function(mean,sigma){
                                  .self$mean = mean
                                  .self$sigma = sigma
                                },

                                evaluate = function(x){
                                  # x: two element numeric
                                  # return value of wedge gaussian density at x

                                  return(dmvnorm(x,mean=.self$mean,sigma=.self$sigma*diag(2)))
                                }
                              )
)
