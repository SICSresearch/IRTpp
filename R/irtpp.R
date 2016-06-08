#######################################################################
#' @name irtpp
#' @title Estimate a test item parameters.
#' @description Estimate a test item parameters according to Item Response Theory.
#' @param dataset The matrix with the responses from the individuals.
#' @param model The model used to calibrate the parameters.
#' @param dims The dimensions to use on the estimation, remember to use the initial parameters if you want highquality estimation.
#' @param initialvalues The matrix with the initial values for the optimization process.
#' @param filename Optional argument specifying a CSV file to read instead of a dataset in memory.
#' @param output Optional.  Additonal arguments that need. 
#' @param loglikflag Optional. Show the loglikelihood at the end of estimation procedure. Also shows AIC and BIC statistic.
#' @return A list containing the estimates of the model parameters, 
#' the number of iterations, the loglikelihood final, the final values
#' of the estimation procedure EM.
#' @examples 
#' ## Simulation data for the model "1pl"
#' # data <- simulateTest(model = "1PL", items = 10, individuals = 500)
#' ## Estimation of the parameters
#' # irtpp(data$test, model = "1PL")
#'
#'## Simulation data for the model "2pl"
#'# data <- simulateTest(model = "2PL", items = 20, individuals = 800)
#'## Estimation of the parameters
#'# irtpp(data$test, model = "2PL")
#'
#'## Simulation data for the model "3pl"
#'# data <- simulateTest(model = "3PL", items = 100, individuals = 1000)
#'## Estimation of the parameters
#'# irtpp(data$test, model = "3PL") 
#'@export

irtpp <- function(dataset          = NULL,
                  model,
                  dims             = 1   ,
                  initialvalues    = NULL,
                  filename         = NULL,
                  output           = NULL,
                  loglikflag       = F,
				  convergenceEpsilon = 0.001)
{
  if(dims > 1) {
	writeLines("Multidimensional analysis not yet implemented.\nPlease wait for the next release of the IRTpp package")
  }
  else
  {
	dataset = data.matrix(dataset);
    mod = irtpp.model(model,asnumber=T);
    ret = uirtestimate(dataset,mod,convergenceEpsilon)

    ret$z = ret$z[1:(mod*ncol(dataset))]
    
    ret$z = matrix(ret$z,ncol=mod,byrow=T)
    
    if(mod == 2)
    {
      ret$z = cbind(ret$z,rep(0,ncol(dataset)))
    }
    if(mod == 1)
    {
      ret$z = cbind(rep(1,ncol(dataset)),ret$z,rep(0,ncol(dataset)))
    }

    ret$z = parameter.list(ret$z)

    if(loglikflag)
    {
      z       = ret$z
      theta   = ret$theta
      weights = ret$weights

      thsum = NULL
      idx   = 1;
      for (th in theta)
      {
        thsum[[idx]] = probability.3pl(z=z,theta=th)
        idx = idx +1;
      }
      thsum
      thmat   = t(matrix(unlist(thsum),ncol = length(theta)))
      i       = 1
      logliks = 0;
      idx     = 1;
      pfrq    = pattern.freqs(data = dataset)
      pfrqs   = pfrq[,ncol(dataset)+1]
      pat     = pfrq[,1:ncol(dataset)]
      pat     = exp(as.matrix(pat)%*%t(log(thmat))+(1-as.matrix(pat))%*%t(log(1-thmat)))

      pwei    = pat%*%weights
      LL      = -sum(log(rep(pwei,pfrqs)))
      print(paste("Loglikelihood : ",LL));
      AIC     = -2*(-LL)+ 2* ncol(dataset)*3;
      BIC     = -2*(-LL)+ log(nrow(dataset)) * ncol(dataset)*3;
      print(paste("AIC : ",AIC));
      print(paste("BIC : ",BIC));
    }
  }
  ret
}

#######################################################################
#' @name individual.traits
#' @title Estimate the latent traits. 
#' @description Estimate the latent traits of the individuals in a test with some given item parameters.
#' @param model The model used to calibrate the parameters.
#' @param itempars The item parameters for the model in matrix form.
#' @param method The method for estimating the traits, may be "EAP" or "MAP".
#' @param dataset The matrix with the responses from the individuals.
#' @param probability_matrix The probability matrix in case it does not need to be recalculated.
#' @return A list with the patterns and the estimated latent traits.
#' @examples 
#' ## Simulation data for the model "3pl"
#' # data <- simulateTest(model = "3PL", items = 100, individuals = 1000)
#' ## Estimation of the parameters
#' # est <- irtpp(data$test, model = "3PL")
#' ## Parameter Matrix 
#' # z <- parameter.matrix(est$z)
#' #trait <- individual.traits(model = "3PL", itempars = z, dataset = data$test,
#' #                  method = "EAP", probability_matrix = est$prob_mat)
#' @export
individual.traits<-function(model,
                            itempars,
                            method,
                            dataset             = NULL,
                            probability_matrix  = NULL)
{
  model = irtpp.model(model,asnumber=T)
  cuads = as.matrix(read.table(system.file("extdata","Cuads.csv",package="IRTpp"),sep=",",header=T))

  if(is.null(dataset)) { stop("Please provide a dataset or filename") }
  est.method = ifelse(method == "EAP", eapinterface, mapinterface)

  est = est.method(zita_par     = itempars,
                   dat          = dataset,
                   e_model      = model,
                   matrix_flag  = !is.null(probability_matrix),
                   prob_matrix  = probability_matrix
                   )

  est   = individual.traits.aux(dataset=dataset, est=est)
  freqs = pattern.freqs(dataset);
  freqs = freqs[,ncol(freqs)];
  trait = est$trait;
  est   = cbind.data.frame(est$patterns,freqs,trait);
  est
}

individual.traits.aux <- function(dataset, est)
{
  est = list(matrix(est[[1]],ncol=dim(dataset)[[2]],byrow=T),est[[2]])
  names(est) <- c("patterns","trait")
  est
}
