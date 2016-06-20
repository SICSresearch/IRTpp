#######################################################################
#' @name IRTpp-package
#' @docType package
#' @title IRTpp : Item Response Theory Implemented in R and Cpp
#' @description IRTpp supports the IRT and MIRT methodologies by implementing 
#' many IRT relevant algorithms for model exploration, calibration and fitting, IRTpp
#' aims to be the fastest using Cpp where it is needed, focusing in extending
#' IRT not only to classical applications but also in new applications which sometimes 
#' need highly performant implementations.
#' @details
#' \tabular{ll}{
#'Package: \tab IRTpp\cr
#'Type: \tab Package\cr
#'Version: \tab 0.2.5\cr
#'Date: \tab 2016-05-14\cr
#'License: \tab MIT + file LICENSE \cr
#'}
#'@author SICS Research Team
#'@keywords IRT MIRT Psychometry 
#'@useDynLib IRTpp
#'@importFrom Rcpp sourceCpp
#'@importFrom FactoMineR PCA
#'@importFrom numDeriv hessian
#'@importFrom mvtnorm rmvnorm
#'@importFrom MASS mvrnorm
#'@section Getting Started:
#'Get started with the IRTpp package browsing the index of this documentation
#'if you need help the vignettes should be helpful.
#'@section Getting Started:
#'The IRTpp package allows you to use the IRT methodology for simulating, analyzing and scoring tests \cr
#'You can browse the package vignettes to get started.
#'
NULL


#######################################################################
#' Probability function for all models.
#' @param model The model to calculate the probability to
#' @return A list of probability functions that are used 
#' according to the required model
#' @keywords internal
 
irtpp.p<- function(model){
  model = irtpp.model(model)
  if(model=="3PL") probability.3pl
  if(model=="2PL") probability.2pl
  if(model=="1PL") probability.1pl
  if(model=="Rasch") probability.2pl
}


####Fast probability functions

#######################################################################
#'@name prob.3pl
#'@title 3PL fast probability function
#'@description The probability function in the 3PL models.
#'@param z List containing the item parameters a, d and c.
#'@param theta a Vector that contains the latent traits of the individual.
#' @seealso
#' \code{\link{probability.3pl}}, \code{\link{prob.2pl}}, 
#' \code{\link{probability.2pl}}, \code{\link{probability.1pl}}
#'@examples 
#'## Simulate the test
#'# data <- simulateTest(model = "3PL", items = 20, individuals = 200)
#'## item parameters
#'# zita <-data$itempars
#'# zita <- list(zita$a, zita$b, zita$c)
#'## transformation of the parameter b
#'# zita <- model.transform(zita, model = "3PL", src = "b", target = "d")
#'# zita <- list(a = zita[,1], d = zita[,1], c = zita[,1])
#'## Latent traits
#'# thetas <- data$latentTraits
#'## Probability model
#'# prob.3pl(zita, thetas)
#'@export
prob.3pl<- function(z, theta){
  prob <- z$c + (1 - z$c)/(1 + exp(-(sum(theta*z$a)+z$d)))
  return(prob)
}

#######################################################################
#'@name prob.2pl
#'@title 2PL fast probability function
#'@description The probability function in the 2PL model.
#'@param z List containing the item parameters a and d
#'@param theta Vector, contains the latent traits of the individual
#'@return The value of the probability for the 2PL model.
#' @seealso
#' \code{\link{probability.2pl}}, \code{\link{prob.3pl}}, 
#' \code{\link{probability.3pl}}, \code{\link{probability.1pl}}
#'@examples 
#'## Simulate the test
#'# data <- simulateTest(model = "2PL", items = 20, individuals = 200)
#'## item parameters
#'# zita <-data$itempars
#'# zita <- list(zita$a, zita$b, zita$c)
#'## transformation of the parameter b
#'# zita <- model.transform(zita, model = "3PL", src = "b", target = "d")
#'# zita <- list(a = zita[,1], d = zita[,1], c = zita[,1])
#'## Latent traits
#'# thetas <- data$latentTraits
#'## Probability model
#'# prob.2pl(zita, thetas)
#'@export
prob.2pl <- function(z,theta){
  prob <- (1)/(1 + exp(-(sum(theta*z$a)+z$d)))
  return(prob)
}

#######################################################################
#'@name probability.3pl
#'@title 3PL probability function
#'@description The probability function in the 3PL model.
#'@param z Optional. A list with the parameters a b and c specified by keys.
#'@param a The discrimination parameter
#'@param b The difficulty parameter. (Optional if d is specified)
#'@param c The guessing parameter
#'@param theta The subject's latent trait.
#'@param d Optional. Overrides the b parameter, it is equal to -a*b. Used in some functions.
#'@param cp Optional. Overrides the c parameter, it is logit(c)
#'@return The value of the probability for the 3pl model.
#' @seealso
#' \code{\link{probability.1pl}}, \code{\link{probability.2pl}}, 
#' \code{\link{prob.3pl}}, \code{\link{prob.2pl}} 
#'@examples
#'## Simulate the test
#'# data <- simulateTest(model = "3PL", items = 20, individuals = 200)
#'## item parameters
#'# zita <-data$itempars
#'## Latent traits
#'# thetas <- data$latentTraits
#'## Probability model
#'# probability.3pl(zita, theta = thetas)
#'@export
probability.3pl <- function (z = NULL, a = z$a, b = z$b, c = NULL, d = z$d, cp = z$cp, 
                           theta) 
{
  if("c" %in% names(z) && is.null(c))
  {
    c = z$c;
  }
  if( length(a)>1 && length(a) == length(theta)){
    ##Multidimensional case.
    if( is.null(c) && !is.null(cp)){
      c = plogis(cp)
    }
    if( is.null(c) && is.null(cp)){
      c = 0;
    }
    if( is.null(d)){
      stop("d parameter must be specified in multidim case.")
    }
    prob <- c + (1 - c)/(1 + exp(-(sum(theta*a)+d)))
    return(prob)
  }
  else{
  if (is.null(d)) {
    d = -a * b
  }
  if (is.null(b)) {
    b = -d/a
  }
  if (is.null(cp)) {
    c + ((1 - c)/(1 + exp(-a * (theta - b))))
  }
  else {
    exp(cp)/(1 + exp(cp)) + (1 - (exp(cp)/(1 + exp(cp)))) * 
      (1 + exp(-(a * theta + d)))^(-1)
    }
  }
}

#######################################################################
#'@name probability.2pl
#'@title 2PL probability function
#'@description The probability function in the 2PL model.
#'@param z Optional. A list with the parameters a and b specified by keys.
#'@param a The discrimination parameter
#'@param b The difficulty parameter. (Optional if d is specified)
#'@param theta The subject's latent trait.
#'@param d Optional. Overrides the b parameter, it is equal to -a*b. Used in some functions.
#'@return Tue value of the probability for the model 2pl.
#' @seealso
#' \code{\link{probability.1pl}}, \code{\link{probability.3pl}}, 
#' \code{\link{prob.3pl}}, \code{\link{prob.2pl}} 
#'@examples
#'## Simulate the test
#'# data <- simulateTest(model = "2PL", items = 20, individuals = 200)
#'## item parameters
#'# zita <-data$itempars
#'## Latent traits
#'# thetas <- data$latentTraits
#'## Probability model
#'# probability.2pl(zita, theta = thetas)
#'@export
probability.2pl <- function(z,a=z$a,b=z$b,theta, d=-a*b)((1)/(1+exp(-a*(theta-b))))

#######################################################################
#'@name probability.1pl
#'@title 1PL probability function
#'@description The probability function in the 1PL model.
#'@param z Optional. A list with the parameter b specified by keys.
#'@param b The difficulty parameter. (Optional if d is specified)
#'@param theta The subject's latent trait.
#'@return Tue value of the probability for the model 1pl.
#' @seealso
#' \code{\link{probability.2pl}}, \code{\link{probability.3pl}}, 
#' \code{\link{prob.3pl}}, \code{\link{prob.2pl}} 
#'@examples 
#'## Simulate the test
#'# data <- simulateTest(model = "1PL", items = 20, individuals = 200)
#'## item parameters
#'# zita <-data$itempars
#'## Latent traits
#'# thetas <- data$latentTraits
#'## Probability model
#'# probability.1pl(zita, theta = thetas)
#'@export
probability.1pl <- function(z,b=z$b,theta)((1)/(1+exp(-(theta-b))))

#######################################################################
#'@name loglik
#'@title LogLikelihood of a IRT model
#'@description LogLikelihood of a IRT model for UIRT
#'@param test A matrix of 0's and 1's
#'@param traits A vector with each individual parameter, or list of 
#'@param z A list with the parameters a b and c specified by keys.
#' Each key must contain a vector of the item parameters for each parameter
#' @return The value of the loglikelihood.
#' @seealso
#' \code{\link{loglik.f}}
#' @keywords internal
loglik<-function (test, traits, z) 
{
  pm = lapply(traits, function(x) probability.3pl(z = z, theta = x));
  Reduce("+", mapply(function(x, y) {
    ifelse(y, log(x), log(1 - x))
  }, unlist(pm), c(t(test)), SIMPLIFY = F))
}

#######################################################################
#'@name loglik.f
#'@title LogLikelihood of a IRT model
#'@description LogLikelihood of a IRT model for UIRT (Fast version)
#'@param test A data frame or a matrix with the responses of the test.
#'@param traits A vector with each individual parameter.
#'@param z A list with the parameters a, b and c specified by keys.
#' Each key must contain a vector of the item parameters for each parameter.
#'@return The value of the loglikelihood.
#' @seealso
#' \code{\link{loglik}}
#' @keywords internal
loglik.f <- function(test, traits, z){
  pm = lapply(traits, function(x) prob.3pl(z = z, theta = x));
  Reduce("+", mapply(function(x, y) {
    ifelse(y, log(x), log(1 - x))
  }, unlist(pm), c(t(test)), SIMPLIFY = F))
}


