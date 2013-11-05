### question 3 -- simulation study
library(MASS) # for multivariate normal random variables
library(glmnet) # for regression routines
library(ggplot2) # for plotting
library(reshape2)

.ar1CovMatrix <- function(n,rho) {
  ### there's got to be a slick way to do this.
  rows <- 1:n
  exp.mat <- matrix(data=sapply(X=rows,FUN=function(i){abs(i-rows)}),ncol=n,byrow=T)
  return(rho ^ exp.mat)
}

### simulation parameters
J <- 10 # number of covariates
beta <- c((-J/2 + 1):0,0:(J/2-1)) / 10 # true coefficients
sigma <- 2 # variance for Y | X
rho <- 0.5 # variance for AR(1) structure

n.learning.set <- 100
n.training.set <- 1000

.generateDataSet <- function(n.obs) {
  ### make X according to N(0,Gamma), Gamma follows AR(1)
  gamma <- .ar1CovMatrix(n=J, rho=rho)
  X.learningSet <- mvrnorm(n=n.obs, mu=rep(0,times=J), Sigma=gamma)
  
  ### make Y according to Y | X ~ N(Xbeta,sigma^2)
  Y.learningSet <- rnorm(n=n.obs, mean=X.learningSet %*% beta, sd=sigma)

  ### format output
  dataSet <- cbind(X.learningSet,Y.learningSet)
  colnames(dataSet) <- c(paste("X", 1:J,sep=""),"Y")
  return(dataSet)
}

learning.set <- .generateDataSet(n.learning.set)
training.set <- .generateDataSet(n.training.set)

### visualize learning set
melted.learning.set <- melt(data=as.data.frame(learning.set),id.vars=c("Y"))

### make scatterplot grid of Y vs. X_j
scatterplot.grid <- ggplot(data=melted.learning.set,aes(y=Y,x=value,group=variable)) + geom_point() + facet_wrap(~ variable,ncol=5)

### look at individual distributions of X_j's, Y, and X*beta
viz.data <- cbind(learning.set,learning.set[,-11]%*%beta)
colnames(viz.data) <- c(colnames(learning.set),"X*beta")

covariate.boxplot <- ggplot(data=melt(data=as.data.frame(viz.data))) + stat_boxplot(aes(x=variable,y=value))

### examine properties of elastic net
### glm-net parametrization has the following:
###   alpha = 0 == ridge
###   alpha = 1 == LASSO
###   alpha = 1/3 == elastic net
# 
# regularization <- "elasticnet"
.obtainCoefficients <- function(regularization="ridge",standardize=T,intercept=T) {
  # set glmnet's alpha and lambda to coincide with classical parametrization
  switch(regularization,
         ridge = {alpha <-0; lambda <- 0:100 / n.learning.set},
         lasso = {alpha <- 1; lambda <- (0:100)/(2*n.learning.set)},
         elasticnet = {alpha <- 1/3; lambda <- 3*(0:100)/(2*n.learning.set)}
         )
  
  glmnet(x=learning.set[,1:J],y=learning.set[,J+1],family="gaussian",alpha=alpha,lambda=lambda,standardize=standardize,intercept=intercept)
}

.calcRidgeDF <- function(lambda) {
  df.seq <- sapply(X=lambda,FUN=function(l){
    design.mat.svd <- svd(learning.set[,1:J])
    singular.values <- design.mat.svd$d
    numerator <- singular.values^2
    denominator <- (numerator + l)
    sum(numerator / denominator)
  })
  return(df.seq)
}
.calculateMSE <- function(glmnet.object, new.data, true.values) {
  Y.hat <- predict(object=glmnet.object,newx=new.data,type="response")
  apply(X=Y.hat, MARGIN=2,
        FUN=function(predicted.values){mean((true.values - predicted.values)^2)}
  )
}

.makePlots <- function(glmnet.object,regularization) {
  require(ggplot2)
  require(reshape2)
  
  ### extract lambda seq, transform it back to original scale
  ### and calculate df
  switch(regularization,
         ridge={
           lambda <- glmnet.object$lambda * n.learning.set
           df <- .calcRidgeDF(lambda)
         },
         lasso={
           lambda <- glmnet.object$lambda * (2*n.learning.set)
           df <- glmnet.object$df
         },
         elasticnet={
           lambda <- glmnet.object$lambda * (2*n.learning.set/3)
           df <- glmnet.object$df
           }
         )

  ### extract coefficients, and package everything
  coef.df <- data.frame(cbind(lambda,df,t(as.matrix(glmnet.object$beta))))
  colnames(coef.df)[1:2] <- c("lambda","df")
  
  ### melt down for visualization
  viz.data1 <- melt(data=coef.df,id.vars=c("lambda","df"))
  
  ### check out beta vs. lambda
  beta.vs.lambda <- ggplot(data=viz.data1,aes(x=lambda,y=value,colour=variable,shape=variable)) + geom_hline(yintercept=beta,alpha=0.25,lty=2) + geom_line() + geom_point() + labs(x=expression(lambda), y=expression(beta))
  
  show(beta.vs.lambda)
  
  ### check out DF vs. lambda:
  df.vs.lambda <- switch(regularization,
         ridge={
           ggplot(data=coef.df) + geom_line(aes(x=lambda,y=df))
         },
         lasso={
           ggplot(data=coef.df) + geom_step(aes(x=lambda,y=df)) + geom_hline(yintercept=0,lty=2,alpha=0.5)
         },
         elasticnet={
           ggplot(data=coef.df) + geom_step(aes(x=lambda,y=df)) + geom_hline(yintercept=0,lty=2,alpha=0.5)
         }
  )
  df.vs.lambda <- df.vs.lambda + labs(x=expression(lambda),y="Effective Degrees of Freedom")
  show(df.vs.lambda)
  
  ### look at beta vs. df
  beta.vs.df <- ggplot(data=viz.data1,aes(y=value,x=df,colour=variable,shape=variable)) + geom_line() + labs(x="Effective Degrees of Freedom", y=expression(beta))
  show(beta.vs.df)
  
  ### check out MSE's
  learning.set.MSE <- .calculateMSE(glmnet.object, new.data=learning.set[,1:J], true.values=learning.set[,J+1])
  training.set.MSE <- .calculateMSE(glmnet.object, new.data=training.set[,1:J], true.values=training.set[,J+1])
  
  ls.MSE.plot <- ggplot(data=data.frame(lambda=lambda, MSE=learning.set.MSE),aes(x=lambda)) + geom_line(aes(y=MSE)) + labs(x=expression(lambda),y="MSE",title=expression(paste("Learning set MSE vs. ",lambda)))
  
  ts.MSE.plot <- ggplot(data=data.frame(lambda=lambda, MSE=training.set.MSE),aes(x=lambda)) + geom_line(aes(y=MSE)) + labs(x=expression(lambda),y="MSE",title=expression(paste("Training set MSE vs. ",lambda)))
  
  show(ls.MSE.plot)
  show(ts.MSE.plot)
  return(list(p1=beta.vs.lambda,p2=df.vs.lambda,p3=beta.vs.df,p4=ls.MSE.plot,p5=ts.MSE.plot))
}


### .calculateMSE




glmnet.object <- .obtainCoefficients("elasticnet")
.makePlots(glmnet.object,"elasticnet")

