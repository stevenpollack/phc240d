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
lambda.values <- 0:250

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
  colnames(dataSet) <- c(paste("X[", 1:J, "]",sep=""),"Y")
  return(dataSet)
}

learning.set <- .generateDataSet(n.learning.set)
training.set <- .generateDataSet(n.training.set)

### visualize learning set
melted.learning.set <- melt(data=as.data.frame(learning.set),id.vars=c("Y"))

### make scatterplot grid of Y vs. X_j
scatterplot.grid <- ggplot(data=melted.learning.set,aes(y=Y,x=value,group=variable)) + geom_point() + facet_wrap(~ variable)

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
.makeAndAnalyzeGlmnet <- function(regularization="ridge",standardize=T,intercept=T) {
  # set glmnet's alpha and lambda to coincide with classical parametrization
  switch(regularization,
         ridge = {alpha <-0; lambda <- lambda.values / n.learning.set},
         lasso = {alpha <- 1; lambda <- lambda.values/(2*n.learning.set)},
         elasticnet = {alpha <- 1/3; lambda <- 3*lambda.values/(2*n.learning.set)}
         )
  
  glmnet.obj <- glmnet(x=learning.set[,1:J],y=learning.set[,J+1],family="gaussian",alpha=alpha,lambda=lambda,standardize=standardize,intercept=intercept)
  
  ### place lambda back on original scale and calculate df
  switch(regularization,
         ridge = {lambda <- glmnet.obj$lambda * n.learning.set;
                  df <- .calcRidgeDF(lambda)},
         lasso = {lambda <-  2 * n.learning.set * glmnet.obj$lambda;
                  df <- glmnet.obj$df},
         elasticnet = {lambda <- 2*n.learning.set * glmnet.obj$lambda / 3 ;
                       df <- glmnet.obj$df}
  )
  ### extract coefficients and bundle relevant data into a df
  glmnet.df <- data.frame(regularization, lambda, df, t(as.matrix(glmnet.obj$beta)),check.names=F ) 
  colnames(glmnet.df)[1:3] <- c("regularization", "lambda", "df")
  
  ### check out MSE's
  glmnet.df <- within(data=glmnet.df,expr={
    learning.set.MSE <- .calculateMSE(glmnet.obj, new.data=learning.set[,1:J], true.values=learning.set[,J+1])
    training.set.MSE <- .calculateMSE(glmnet.obj, new.data=training.set[,1:J], true.values=training.set[,J+1])
  })
   
  return(list(glmnet.object=glmnet.obj,data.frame=glmnet.df))
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

.makePlots <- function(glmnet.df,display.plots=F) {

  ### melt down for visualization; exclude MSE's
  viz.data1 <- melt(data=glmnet.df[,-c(J+4,J+5)],id.vars=c("regularization","lambda","df"))
  
  ### check out beta vs. lambda
  beta.vs.lambda <- ggplot(data=viz.data1,aes(x=lambda,y=value,colour=variable,shape=variable)) + geom_hline(yintercept=beta,alpha=0.25,lty=2) + geom_line() + labs(x=expression(lambda), y=expression(beta)) # + geom_point() 
  
  
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
  
  ### look at beta vs. df
  beta.vs.df <- ggplot(data=viz.data1,aes(y=value,x=df,colour=variable,shape=variable)) + geom_line() + labs(x="Effective Degrees of Freedom", y=expression(beta))
  
  ### check out MSE's
  learning.set.MSE <- .calculateMSE(glmnet.object, new.data=learning.set[,1:J], true.values=learning.set[,J+1])
  training.set.MSE <- .calculateMSE(glmnet.object, new.data=training.set[,1:J], true.values=training.set[,J+1])
  
  ls.MSE.plot <- ggplot(data=data.frame(lambda=lambda, MSE=learning.set.MSE),aes(x=lambda)) + geom_line(aes(y=MSE)) + labs(x=expression(lambda),y="MSE",title=expression(paste("Learning set MSE vs. ",lambda)))
  
  ts.MSE.plot <- ggplot(data=data.frame(lambda=lambda, MSE=training.set.MSE),aes(x=lambda)) + geom_line(aes(y=MSE)) + labs(x=expression(lambda),y="MSE",title=expression(paste("Training set MSE vs. ",lambda)))
  
  if (display.plots) {
    show(beta.vs.lambda)
    show(df.vs.lambda)
    show(beta.vs.df)
    show(ls.MSE.plot)
    show(ts.MSE.plot)
  } 
  
  return(list(p1=beta.vs.lambda,p2=df.vs.lambda,p3=beta.vs.df,p4=ls.MSE.plot,p5=ts.MSE.plot))
}


set.seed(1234)

.calcOLSCoefs <- function(X,Y) {
  t.X <- t(X)
  solve(t.X %*% X) %*% t.X %*% Y
}

simulation.study1 <- .makeAndAnalyzeGlmnet()
simulation.study3 <- .makeAndAnalyzeGlmnet(regularization="lasso")
simulation.study2 <- .makeAndAnalyzeGlmnet(regularization="elasticnet")

combined.df <- rbind(simulation.study1$data.frame,simulation.study2$data.frame,simulation.study3$data.frame)

ridge.t.MSE.min.indx <- which.min(simulation.study1$data.frame$training.set.MSE)
enet.t.MSE.min.indx <- which.min(simulation.study2$data.frame$training.set.MSE)
lasso.t.MSE.min.indx <- which.min(simulation.study3$data.frame$training.set.MSE)

location.of.MSE.mins <- rbind(simulation.study1$data.frame[ridge.t.MSE.min.indx,],simulation.study2$data.frame[enet.t.MSE.min.indx,],simulation.study3$data.frame[lasso.t.MSE.min.indx,])

melted.df <- melt(data=combined.df,id.vars=c("lambda","df","regularization")); head(melted.df)
mse.indx <- which(melted.df$variable %in% c("learning.set.MSE","training.set.MSE"))

### look at MSE vs. lambda for training and learning sets
ggplot(data=melted.df[mse.indx,],aes(x=lambda,y=value,color=regularization,lty=variable)) + geom_line() + geom_vline(data=location.of.MSE.mins,aes(xintercept=lambda,color=regularization),lty=8,alpha=0.75) + labs(x=expression(lambda),y="MSE")

### look at how coefficients vary with lambda
ggplot(data=melted.df[-mse.indx,],aes(x=lambda,y=value,color=variable)) + geom_line(show_guide=F) + geom_hline(yintercept=beta,alpha=0.5,lty=3) +  geom_vline(data=location.of.MSE.mins,aes(xintercept=lambda),lty=8,alpha=0.5) + geom_text(data=subset(x=melted.df[-mse.indx,],subset={lambda==0}),aes(label=variable,x=0,y=value,color=variable,group=regularization),show_guide=F,hjust=1,parse=T) + facet_grid(regularization~.) + labs(x=expression(lambda),y=expression(beta))

### compare "optimal" coefficients between all 4 methods
beta.OLS <- .calcOLSCoefs(X=learning.set[,1:J],Y=learning.set[,J+1])
extra.info <- rbind(data.frame(regularization="OLS",lambda=0,df=J,t(beta.OLS),learning.set.MSE=0,training.set.MSE=0,check.names=F),location.of.MSE.mins)
ggplot(data=melt(extra.info[,-(14:15)],id.vars=c("lambda","df","regularization")),aes(x=regularization,y=value,color=variable,group=variable)) + geom_point(show_guide=F) + geom_line(show_guide=F) + geom_hline(yintercept=beta,alpha=0.5,lty=3) + geom_text(data=melt(extra.info[1,-(14:15)],id.vars=c("lambda","df","regularization")),aes(label=variable,x=as.factor("OLS"),y=value),hjust=1.25,vjust=0,show_guide=F,parse=T) + labs(x="method",y=expression(beta))

### look at beta vs. df
ggplot(data=melted.df[-mse.indx,],aes(x=df,y=value,color=variable)) + geom_line(show_guide=F) + geom_point(data=melt(location.of.MSE.mins[,-(14:15)],id.vars=c("lambda","df","regularization")),aes(x=df,y=value),size=3,shape=2,show_guide=F) + geom_text(data=subset(x=melted.df[-mse.indx,],subset={df==10 & lambda==0}),aes(label=variable,x=10,y=value,color=variable,group=regularization),show_guide=F,hjust=-.05,parse=T) +  facet_grid(regularization~.) + labs(x="Effective Degrees of Freedom",y=expression(beta))

### check out df vs. lambda
ggplot(data=melted.df[-mse.indx,],aes(x=lambda,y=df,color=regularization)) + geom_line() + geom_vline(data=location.of.MSE.mins,aes(xintercept=lambda,color=regularization),lty=8,alpha=0.75)  + labs(x=expression(lambda),y="Effective Degrees of Freedom",title=expression(X[i]))

