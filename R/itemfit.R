#######################################################################
#' @name x2_itemf
#' @title X2 Statistic
#' @description calculates the values of statistical X2 from Steven P. Reise.
#' @usage x2_itemf(model,z,patterns,pval.sim,G,FUN,B=NULL)
#' @param model	type of model ( "1PL", 2PL", "3PL" )
#' @param z	list of estimations of the parameters of the items (discrimination, difficulty, guessing).
#' @param patterns matrix of patterns response, the frequency of each pattern and the latent traits.
#' @param pval.sim if TRUE simulates a p-value using bootstrap.
#' @param G	number of groups.
#' @param FUN function to calculate the frecuency expected in the each group (mean, median,..)
#' @param B number of iterations bootstrap. If "NULL" B=100.
#' @return X2 Statistic and p-value.
#' @author SICS Research, National University of Colombia \email{ammontenegrod@@unal.edu.co}
#' @export
#'
#' @seealso
#' \code{\link{z3_personf}}, \code{\link{plotenvelope}}
#'
#' @references
#'
#' Steven P. Reise (1990). A Comparison of item-and person-fit Methods of Assesing Model-Data fit in IRT \emph{University of Minnesota}.
#'
#' @examples
#'
#' #Simulates a test and returns a list:
#' test <- simulateTest()
#'
#' #the simulated data:
#' datos <- test$test
#'
#' #model:
#' mod <- irtpp(dataset = datos,model = "3PL")
#'
#' #latent trait:
#' zz <- parameter.matrix(mod$z)
#' p_mat <- mod$prob_mat
#' traits <- individual.traits(model="3PL",method = "EAP",dataset = datos,itempars = zz,
#' probability_matrix=p_mat)
#'
#' #X2-Statistic
#' x2_itemf(model = "3PL",z = mod$z,patterns = traits,pval.sim = TRUE,G = 10,FUN = median,B=4)
#'

# <- function(model,z,patterns,pval.sim,G,FUN,B=NULL){

#  check.model(model)
#  if(is.null(B)){B=100}

#  x2obs <- x2(model=model,z=z,patterns=patterns,G=G,FUN=FUN)

#    if (pval.sim==F) {
#      df <- G-irtpp.model(model,asnumber=T);
#      pvals <- pchisq(x2obs, df = df, lower.tail = FALSE)
#    }

#  else {
#    T.boot <- matrix(0, ncol(patterns[,-c(ncol(patterns)-1,ncol(patterns))]), B)
#      for (b in 1:B) {
#      X.new <- simulateTest(model = model,ncol(patterns)-2,individuals = sum(patterns[,ncol(patterns)-1]),itempars = z)
#      X.new <- X.new$test
      
#      object.new <- irtpp(X.new,model = model)
#      print("Iteration:")
#      print(b)
#      z <- parameter.matrix(object.new$z)
#      z <- model.transform(z,model=model,"b","d")
#      z  <- parameter.list(z,dpar = T)
#      z_mat <- parameter.matrix(object.new$z)
#      p_mat <- object.new$prob_mat
#      pat.new <- individual.traits(model=model,z_mat,method="EAP",dataset = X.new,probability_matrix=p_mat)
#      freqs <- cbind(pattern.freqs(X.new),pat.new$trait)
      
#      T.boot[, b] <- x2(model=model,z=z,patterns=freqs,G=G,FUN=FUN)
#      }
#    pvals <- (rowSums(T.boot >= rep(x2obs, B), na.rm = TRUE) +
#                1)/(B + 1) 
#  }

#  return(cbind(x2obs,pvals))

#}

#######################################################################
#' @name z3_personf
#' @title Z3 Person fit statistic
#' @description Calculates the values of statistical Z3 for individuals.
#' @param data a data frame or a matrix with the test.
#' @param zita a list of estimations of the parameters of the items (discrimination,difficulty, guessing).
#' @param patterns matrix of patterns response, the frequency of each pattern and the latent traits.
#' @author SICS Research, National University of Colombia \email{ammontenegrod@@unal.edu.co}
#' @export
#'
#' @seealso
#' \code{\link{z3_itemf}}, \code{\link{orlando_itemf}}
#'
#' @references
#'
#' Fritz Drasgow, Michael V. Levine and Esther A. Williams (1985). Appropiateness measurement with polychotomous item response models and standarized indices.
#'
#' @examples
#'
#' #Simulates a test and returns a list:
#' test <- simulateTest()
#'
#' #the simulated data:
#' data <- test$test
#'
#' #model:
#' mod <- irtpp(dataset = data,model = "3PL")
#'
#' #latent trait:
#' zz <- parameter.matrix(mod$z)
#' p_mat <- mod$prob_mat
#' traits <- individual.traits(model="3PL",method = "EAP",dataset = data,itempars = zz,
#' probability_matrix=p_mat)
#'
#' #Z3 PERSONFIT-Statistic
#' z3_personf(data = data,zita = mod$z,patterns = traits)
#'

z3_personf <- function(data,zita,patterns){
  scores <- patterns[,-(ncol(patterns)-1)]
  nitems <- ncol(data) 
  nscores <- nrow(patterns) 
  ninds <- nrow(data) 

  index <- indexPat(data,patterns[,-ncol(patterns)]) 
  scoresTot <- numeric(nrow(data)) 
  for(mm in 1:nrow(scores)){
    scoresTot[index[[mm]]] <- scores[mm,ncol(data) +1]  
  }

  P <- lapply(scoresTot,FUN=function(x){probability.3pl(theta=x,z=zita)}) 
  P <- matrix(unlist(P),ncol=nitems,byrow=T)

  LL <- matrix(0,ncol = ncol(P),nrow = nrow(P)) 
  LL[data == 1] <- P[data == 1] 
  LL[data == 0] <- 1 - P[data == 0] 
  LL <- rowSums(log(LL)) 

  mu <- sigmaCuad <- rep(0,ninds)
  for( i in 1:nitems){
    Pi <- cbind(P[,i],1 - P[,i])
    logPi <- log(Pi)
    mu <- mu+rowSums(Pi * logPi)
    sigmaCuad <- sigmaCuad + (Pi[,1] * Pi[,2] * (log(Pi[,1]/Pi[,2])^2))
  }

  Z3 <- (LL - mu) / sqrt(sigmaCuad)
  Z3

}

#######################################################################
#' @name z3_itemf
#' @title Z3 statistic.
#' @description Calculates the values of statistical Z3 for items.
#' @param data a data frame or a matrix with the test.
#' @param zita a list of estimations of the parameters of the items (discrimination,difficulty, guessing).
#' @param patterns matrix of patterns response, the frequency of each pattern and the latent traits.
#' @author SICS Research, National University of Colombia \email{ammontenegrod@@unal.edu.co}
#' @export
#'
#' @seealso
#' \code{\link{z3_personf}}
#'
#' @references
#'
#' Fritz Drasgow, Michael V. Levine and Esther A. Williams (1985). Appropiateness measurement with polychotomous item response models and standarized indices.
#'
#' @examples
#'
#' #Simulates a test and returns a list:
#' test <- simulateTest()
#'
#' #the simulated data:
#' data <- test$test
#'
#' #model:
#' mod <- irtpp(dataset = data,model = "3PL")
#'
#' #latent trait:
#' zz <- parameter.matrix(mod$z)
#' p_mat <- mod$prob_mat
#' traits <- individual.traits(model="3PL",method = "EAP",dataset = data,itempars = zz,
#' probability_matrix=p_mat)
#'
#' #Z3 PERSONFIT-Statistic
#' z3_itemf(data = data,zita = mod$z,patterns = traits)
#'

z3_itemf <- function(data,zita,patterns){

  scores <- patterns[,-(ncol(patterns)-1)]
  nitems <- ncol(data) 
  nscores <- nrow(patterns) 
  ninds <- nrow(data) 
  
  index <- indexPat(data,patterns[,-ncol(patterns)]) 
  scoresTot <- numeric(nrow(data)) 
  for(mm in 1:nrow(scores)){
    scoresTot[index[[mm]]] <- scores[mm,ncol(data) +1]  
  }

  P <- lapply(scoresTot,FUN=function(x){probability.3pl(theta=x,z=zita)}) #p_i (probabilidad de contestar correctamente al item i)
  P <- matrix(unlist(P),ncol=nitems,byrow=T)

  LL <- matrix(0,ncol = ncol(P),nrow = nrow(P)) 
  LL[data == 1] <- P[data == 1] 
  LL[data == 0] <- 1 - P[data == 0] 
  LL <- colSums(log(LL)) 

  mu <- sigmaCuad <- rep(0,nitems)  
  for( i in 1:nitems){
    Pi <- cbind(P[,i],1 - P[,i]) 
    logPi <- log(Pi) 
    mu[i] <- sum(Pi * logPi)    
    sigmaCuad[i] <- sum(Pi[,1] * Pi[,2] * (log(Pi[,1]/Pi[,2])^2)) 

  }
  Z3 <- (LL - mu) / sqrt(sigmaCuad) 
  Z3
}

#######################################################################
#' @name orlando_itemf
#' @title Orlando's statistic
#' @description Calculate the values of the statistics S_x2 from
#' Maria Orlando and David Thisen.
#' @param patterns matrix of patterns response, the frequency of each pattern and the latent traits.
#' @param G number of quadrature points.
#' @param zita matrix of estimations of the parameters of the items (discrimination, difficulty, guessing).
#' @param model type of model ( "1PL", 2PL", "3PL" ).
#' @return Orlando's statistic, degrees of freedom and p-value for each item.
#' @export
#'
#' @seealso
#' \code{\link{z3_itemf}}
#' 
#' @author 
#' SICS Research, National University of Colombia \email{ammontenegrod@@unal.edu.co}
#' 
#' @references
#' Orlando, M. & Thissen, D. (2000). Likelihood-based item fit indices for dichotomous item
#' response theory models. \emph{Applied Psychological Measurement, 24}, 50-64.
#'
#' @examples
#'
#' #Simulates a test and returns a list:
#' test <- simulateTest()
#'
#' #the simulated data:
#' data <- test$test
#'
#' #model:
#' mod=irtpp(dataset = data,model = "3PL")
#' 
#' #Convert parameters to a matrix
#' zz <- parameter.matrix(mod$z,byrow = FALSE)
#' 
#' #Estimating Latent Traits
#' p_mat <- mod$prob_mat
#' trait <- individual.traits(model="3PL", itempars = zz,method = "EAP",dataset = data,
#' probability_matrix = p_mat)
#' 
#' #Z3 PERSONFIT-Statistic
#' orlando_itemf(patterns = as.matrix(trait),G = 61,zita = mod$z,model = "3PL")
#'

orlando_itemf <- function(patterns,G,zita,model){

  if(model=="3PL"){mo=3}
  if(model=="2PL"){mo=2}
  if(model=="1PL"){mo=1}

  pats <- as.matrix(patterns[,-ncol(patterns)])
  frec <- patterns[,ncol(patterns)-1]  
  patsSinFrec <- pats[,-ncol(pats)]  
  nitems <- ncol(patsSinFrec)

  seq <- seq(-6,6,length=G)
  pesos <- dnorm(seq)/sum(dnorm(seq))
  Cuad <- matrix(c(seq,pesos),byrow=F,ncol=2)

  theta <- Cuad
  w.cuad <- theta[,2] 
  thetaG <- theta[,1] 

  pr <- lapply(thetaG,FUN=function(x){probability.3pl(theta=x,z=zita)}) #probabilidad para cada punto de cuadratura
  pr <- matrix(unlist(pr),ncol=nitems,byrow=T)

  score <- rowSums(patsSinFrec )  
  Nk <- NULL
  for(i in 1:nitems - 1){
    inds <- which(score == i) 
    patsCoin <- pats[inds,]   
    if(class(patsCoin) == "matrix"){
      if(dim(patsCoin)[1] != 0){     
        Nk[i] <- sum(patsCoin[,ncol(patsCoin)])
      }else{
        Nk[i] <- 0
      }
    }else{
      if(class(patsCoin) == "numeric"){
        Nk[i] <- patsCoin[length(patsCoin)]
      }
    }
  }

  O <- list()
  Oik <- matrix(0,ncol = nitems -1 ,nrow = nitems)  
  for(i in 1:nitems - 1){  
    inds <-which(score == i); 
    for(j in 1:(nitems)){  
      patsCoin <- pats[inds,]   
      if(class(patsCoin) == "matrix"){
        if(dim(patsCoin)[1] != 0){     
          Oik[j,i] <- sum(apply(X = patsCoin,MARGIN = 1,FUN = function(x){ifelse(x[j] == 1,yes = x[nitems + 1],0)}))
        }else{
          Oik[j,i] <- 0
        }
      }else{
        if(class(patsCoin) == "numeric"){
          Oik[j,i] <- ifelse(patsCoin[j] == 1,yes = patsCoin[nitems + 1],0)
        }
      }
      O[[j]] <- cbind(Nk-Oik[j,],Oik[j,])
    }
  }

  sact <- s_ss(pr,nitems=nitems,G=G)
  Denom <- colSums(matrix(rep(w.cuad,nitems -1 ),ncol = nitems - 1) * sact[,-c(1,ncol(sact))])

  nitems <- nitems - 1
  smono <- list()

  for(p in 1:(nitems+1)){
    smono[[p]] <- s_ss(pr[,-p],nitems=nitems,G=G)
  }

  nitems <- ncol(patsSinFrec)
  E <- list()
  for(i in 1:length(smono)){
    E[[i]] <- cbind(1-colSums(smono[[i]][,-ncol(smono[[i]])]*(pr[,i]*w.cuad))/Denom,colSums(smono[[i]][,-ncol(smono[[i]])]*(pr[,i]*w.cuad))/Denom)
  }

  S_X2 <- NULL
  df.S_X2 <- NULL
  for (i in 1:nitems) E[[i]] <- E[[i]] * Nk
  coll <- collapseCells(O, E, mincell = 1)
  O <- coll$O
  E <- coll$E
  for (i in 1:length(O)) {
    S_X2[i] <- sum((O[[i]] - E[[i]])^2/E[[i]], na.rm = TRUE)
    df.S_X2[i] <- sum(!is.na(E[[i]])) - nrow(E[[i]]) - mo
  }

  pval <- pchisq(S_X2,df.S_X2,lower.tail = F)

  lista <- cbind("S_X2"=S_X2,"df.SX2"=df.S_X2,"p.val.S_X2"=pval)

  return(lista)

}




#######################################################################
#' @name plotevenlope
#' @title Confidence bands
#' @description Graphic bands of confidence of an item, to evaluate 
#' the goodness of fit of the model.
#' @param item a number indicating the item that you want evaluate.
#' @param numboot number of iterations bootstrap, used to plot the envelopes.
#' @param alpha level of significance to plot the envelopes.
#' @param model object irtpp() type.
#' @param data the data frame or a matrix with the test data.
#' @return plot with the envelopes and the caracteristic curve of the item.
#' @author SICS Research, National University of Colombia \email{ammontenegrod@@unal.edu.co}
#' @export
#'
#' @seealso
#' \code{\link{orlando_itemf}}, \code{\link{z3_itemf}}
#'
#' @references
#'
#' David Thissen, Howard Wainer  D. (1990). Confidence Envelopes for Item Response Theory. \emph{Journal of Educational Statistics, Vol 15, No 2}, 113-128.
#' @examples
#'
#' #Simulates a test and returns a list:
#' test <- simulateTest()
#'
#' #the simulated data:
#' t_data <- test$test
#'
#' #model:
#' mod <- irtpp(dataset = t_data,model = "3PL")
#'
#' #Envelopes:
#' item <- 8
#' numboot <- 100
#' alpha <- 0.05
#'
#' #call the function:
#' plotenvelope(item=item,numboot=numboot,alpha=alpha,model=mod,data=t_data)

plotenvelope=function(item, numboot, alpha, model, data){


  r <- matrix(model$r,ncol=ncol(data),byrow=F)
  f <- matrix(rep(model$f,ncol(data)),ncol=ncol(data),byrow=F)
  theta <- model$theta
  params <- model.transform(model$z,model = "3PL","b","d")

  nitems <- ncol(r) 
  x <- seq(from = -6,to = 6,by=0.1)
  y <- sapply(X = x,FUN = gg,a = params[item,1],d = params[item,2],cp = qlogis(params[item,3])) 
  inf <- solve(-hessian(func=llikm,x=as.vector(params[item,]),args=list(theta,r,f)))
  media <- params[item,] 
  boot <- rmvnorm(numboot,mean = media,sigma = inf)
  boot[,3] <- ifelse(boot[,3] >= 1,1 - 1e-6,boot[,3])
  boot[,3] <- ifelse(boot[,3] <= 0,sqrt(2.2e-16),boot[,3])
  boot[,1] <- ifelse(boot[,1] <= 0,sqrt(2.2e-16),boot[,1])

  envelop <- matrix(0,nrow = length(x),ncol = 2)
  for(i in 1:length(x)){
    trace <- apply(X = boot, MARGIN = 1,
                  FUN = function(puntoBoot){
                    gg(a = puntoBoot[1], d = puntoBoot[2],cp = qlogis(puntoBoot[3]),theta = x[i])
                  })
    trace <- sort(trace) 
    lower <- trace[floor(length(trace) * alpha / 2)] 
    upper <- trace[ceiling(length(trace) * (1 - (alpha / 2)))] 
    envelop[i,] = c(lower,upper)

  }

  plot(x,y,xlim=c(-6,6),ylim=c(0,1),type="l",main=paste("Item Characteristic Curve",sep=" "),
       xlab = expression(theta),ylab = expression(P(theta)))
  lines(x,envelop[,1],col="red",lty=2)
  lines(x,envelop[,2],col="red",lty=2)

}
