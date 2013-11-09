### question 3 -- simulation study

### ------------------------
### simulation libraries
### ------------------------
set.seed(1234)
library(MASS) # for multivariate normal random variables
library(glmnet) # for regression routines
library(ggplot2) # for plotting
library(reshape2)

### -----------------------
### helper functions
### -----------------------
.ar1CovMatrix <- function(n,rho) {
  ### there's got to be a slick way to do this.
  rows <- 1:n
  exp.mat <- matrix(data=sapply(X=rows,FUN=function(i){abs(i-rows)}),ncol=n,byrow=T)
  return(rho ^ exp.mat)
}
.centerMat <- function(mat) {
  apply(X=mat,MARGIN=2,FUN=function(col){col-mean(col)})
}
.generateDataSet <- function(n.obs,J,rho) {
  ### make X according to N(0,Gamma), Gamma follows AR(1)
  ### make sure X is column-centered
  gamma <- .ar1CovMatrix(n=J, rho=rho)
  X.learningSet <- .centerMat(mvrnorm(n=n.obs, mu=rep(0,times=J), Sigma=gamma))
  
  ### make Y according to Y | X ~ N(Xbeta,sigma^2)
  Y.learningSet <- rnorm(n=n.obs, mean=X.learningSet %*% beta, sd=sigma)
  
  ### format output
  dataSet <- cbind(X.learningSet,Y.learningSet)
  colnames(dataSet) <- c(paste("X[", 1:J, "]",sep=""),"Y")
  return(dataSet)
}
### routine for graphical portion of EDA of learning set
.visualizeLearningSet <- function(learning.set) {
  melted.learning.set <- melt(data=as.data.frame(learning.set),id.vars=c("Y"))
  
  ### make scatterplot grid of Y vs. X_j
  scatterplot.grid <- ggplot(data=melted.learning.set,aes(y=Y,x=value,group=variable)) + geom_point() + facet_wrap(~ variable)
  
  ### look at individual distributions of X_j's, Y, and X*beta
  viz.data <- cbind(learning.set,learning.set[,-11]%*%beta)
  colnames(viz.data) <- c(colnames(learning.set),"X*beta")
  
  covariate.boxplot <- ggplot(data=melt(data=as.data.frame(viz.data))) + stat_boxplot(aes(x=variable,y=value))
}
### helper routines to do calculations for analysis
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
### routine to calculate beta's and overal MSE for various regularization schemes
.makeAndAnalyzeGlmnet <- function(regularization="ridge",standardize=F,intercept=T) {
  # note: analysis does not require columns of X to be anything other
  # than column-centered, which is taken care of during generation.
  # so perform glmnet without standardization
  
  # set glmnet's alpha and lambda to coincide with classical parametrization
  switch(regularization,
         ridge = {alpha <-0; lambda <- lambda.values / n.learning.set},
         lasso = {alpha <- 1; lambda <- lambda.values/(2*n.learning.set)},
         elasticnet = {alpha <- 1/3; lambda <- 3*lambda.values/(2*n.learning.set)}
  )
  
  # actually do the regression
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
### ridge-specific routines:
.calculateRidgeCoefs <- Vectorize(vectorize.args="lambda",FUN=function(X,Y,lambda){
  J <- dim(X)[2] # number of features
  X.t <- t(X); M <- X.t %*% X
  A <- solve(M + lambda * diag(J))
  A %*% X.t %*% Y
})
.calculateRidgeMSE <- Vectorize(vectorize.args="lambda",FUN=function(X,lambda,beta,sigma) {
  J <- dim(X)[2] # number of features
  X.t <- t(X); M <- X.t %*% X
  A <- solve(M + lambda * diag(J))
  bias <- (A %*% M - diag(J)) %*% beta
  var <- sigma^2 * diag(A %*% M %*% A)
  return(list(bias=bias,var=var,MSE=bias^2 + var))
} )
### OLS-specific routine:
.calcOLSCoefs <- function(X,Y) {
  t.X <- t(X)
  solve(t.X %*% X) %*% t.X %*% Y
}


### -----------------------
### simulation parameters
### -----------------------
J <- 10 # number of covariates
beta <- c((-J/2 + 1):0,0:(J/2-1)) / 10 # true coefficients
sigma <- 2 # variance for Y | X
rho <- 0.5 # variance for AR(1) structure
lambda.values <- 0:250 # domain of hyper-parameter space to optimize over
n.learning.set <- 100 # number of obs in learning set
n.training.set <- 1000 # number of obs in test set


### ------------------------
### simulation study
### ------------------------

### 1) make learning and test set
learning.set <- .generateDataSet(n.learning.set,J,rho)
training.set <- .generateDataSet(n.training.set,J,rho)

### 2) perform EDA on learning set

### 3) examine properties of regularization

### make coefficients for each regularization scheme
ridge.study <- .makeAndAnalyzeGlmnet(); ridge.study.df <- ridge.study$data.frame
lasso.study <- .makeAndAnalyzeGlmnet(regularization="lasso"); lasso.study.df <- lasso.study$data.frame
enet.study <- .makeAndAnalyzeGlmnet(regularization="elasticnet"); enet.study.df <- enet.study$data.frame


### find lambda associated to minimal traning risk
ridge.t.MSE.min.indx <- which.min(ridge.study.df$training.set.MSE)
enet.t.MSE.min.indx <- which.min(enet.study.df$training.set.MSE)
lasso.t.MSE.min.indx <- which.min(lasso.study.df$training.set.MSE)

location.of.MSE.mins <- rbind(ridge.study.df[ridge.t.MSE.min.indx,],enet.study.df[enet.t.MSE.min.indx,],lasso.study.df[lasso.t.MSE.min.indx,])

### reshape data for visualization
combined.df <- rbind(ridge.study.df,enet.study.df,lasso.study.df)
melted.df <- melt(data=combined.df,id.vars=c("lambda","df","regularization")); head(melted.df)
mse.indx <- which(melted.df$variable %in% c("learning.set.MSE","training.set.MSE"))

### look at MSE vs. lambda for training and learning sets
mse.vs.lambda <- ggplot(data=melted.df[mse.indx,],aes(x=lambda,y=value,color=regularization,lty=variable)) + geom_line() + geom_vline(data=location.of.MSE.mins,aes(xintercept=lambda,color=regularization),lty=8,alpha=0.75) + labs(x=expression(lambda),y="MSE")
show(mse.vs.lambda)

### look at how coefficients vary with lambda
beta.vs.lambda <- ggplot(data=melted.df[-mse.indx,],aes(x=lambda,y=value,color=variable)) + geom_line(show_guide=F) + geom_hline(yintercept=beta,alpha=0.5,lty=3) +  geom_vline(data=location.of.MSE.mins,aes(xintercept=lambda),lty=8,alpha=0.5) + geom_text(data=subset(x=melted.df[-mse.indx,],subset={lambda==0}),aes(label=variable,x=0,y=value,color=variable,group=regularization),show_guide=F,hjust=1,parse=T) + facet_grid(regularization~.) + labs(x=expression(lambda),y=expression(beta))
show(beta.vs.lambda)

### compare "optimal" coefficients between all 4 methods
beta.OLS <- .calcOLSCoefs(X=learning.set[,1:J],Y=learning.set[,J+1])
extra.info <- rbind(data.frame(regularization="OLS",lambda=0,df=J,t(beta.OLS),learning.set.MSE=0,training.set.MSE=0,check.names=F),location.of.MSE.mins)

coefs.vs.method <- ggplot(data=melt(extra.info[,-(14:15)],id.vars=c("lambda","df","regularization")),aes(x=regularization,y=value,color=variable,group=variable)) + geom_point(show_guide=F) + geom_line(show_guide=F) + geom_hline(yintercept=beta,alpha=0.5,lty=3) + geom_text(data=melt(extra.info[1,-(14:15)],id.vars=c("lambda","df","regularization")),aes(label=variable,x=as.factor("OLS"),y=value),hjust=1.25,vjust=0,show_guide=F,parse=T) + labs(x="method",y=expression(beta))
show(coefs.vs.method)

### look at beta vs. df
beta.vs.df <- ggplot(data=melted.df[-mse.indx,],aes(x=df,y=value,color=variable)) + geom_line(show_guide=F) + geom_point(data=melt(location.of.MSE.mins[,-(14:15)],id.vars=c("lambda","df","regularization")),aes(x=df,y=value),size=3,shape=2,show_guide=F) + geom_text(data=subset(x=melted.df[-mse.indx,],subset={df==10 & lambda==0}),aes(label=variable,x=10,y=value,color=variable,group=regularization),show_guide=F,hjust=-.05,parse=T) +  facet_grid(regularization~.) + labs(x="Effective Degrees of Freedom",y=expression(beta))
show(beta.vs.df)

### check out df vs. lambda
df.vs.lambda <- ggplot(data=melted.df[-mse.indx,],aes(x=lambda,y=df,color=regularization)) + geom_line() + geom_vline(data=location.of.MSE.mins,aes(xintercept=lambda,color=regularization),lty=8,alpha=0.75)  + labs(x=expression(lambda),y="Effective Degrees of Freedom",title=expression(X[i]))
show(df.vs.lambda)

### investigate properties of Ridge estimator;
# For the simulation model of a), provide and comment on graphical displays of
# the bias, variance, and MSE of the ridge estimators based on the learning set.
# For each coefficient, provide the value of the shrinkage parameter λ minimizing
# the MSE and the corresponding estimate.

### use closed form expression to investigate bias, variance, and MSE of beta_ridge
library(plyr)
ridge.sampling.stats.list.mat <- .calculateRidgeMSE(learning.set[,1:J],lambda.values,beta,sigma)
ridge.stats <- adply(.data=ridge.sampling.stats.list.mat, .margins=1,
                     .fun=function(row){cbind(lambda=lambda.values,ldply(.data=row, .fun='t'))}
                     )
colnames(ridge.stats)[1] <- "stat"

### find lambda that minimizes MSE for individual coefficient
ridge.indv.study.df <- adply(.data=subset(ridge.stats,subset={stat=="MSE"})[,-(1:2)],.margins=2,.fun=function(column){
  col <- unlist(column)
  col.min.indx <- which.min(col)
  lambda.min <- lambda.values[col.min.indx]
  MSE.min <- min(col)
  coef.at.min <-ridge.study.df[col.min.indx,colnames(column)]
  c(coef=coef.at.min,lambda=lambda.min,MSE=MSE.min)
  })

### calculate MSE of this frankenstein estimator
frank.MSE.lset <- mean( (learning.set[,J+1] - learning.set[,1:J] %*% ridge.indv.study.df$coef)^2 )
frank.MSE.tset <- mean( (training.set[,J+1] - training.set[,1:J] %*% ridge.indv.study.df$coef)^2 )

### compare this to minimum MSE found above
library(xtable)
mse.comparison <- rbind(c(Learning=frank.MSE.lset,Test=frank.MSE.tset),c(Learning=location.of.MSE.mins$learning.set.MSE[1],Test=location.of.MSE.mins$training.set.MSE[1]))
mse.comparison <- data.frame(method=c("aggregated","argmin"),mse.comparison)
mse.comparison.plot <- ggplot(data=melt(mse.comparison),aes(x=method)) + geom_bar(aes(weight=value,fill=variable),position="dodge") + geom_text(aes(x=method,y=value,label=round(value,digits=3)),hjust=c(1.5,-0.5,1.5,-0.5),vjust=-0.2,color='black',position="dodge") + labs(y="MSE") + scale_fill_discrete(name="")
show(mse.comparison.plot)

### compare frankenstein coefficients to argmin coefs
coef.comparison <- data.frame(ridge.indv.study.df[,1:2],t(location.of.MSE.mins[1,4:13]))
colnames(coef.comparison) <- c("covariate","aggregated","argmin")
coef.comparison.plot <- ggplot(data=melt(coef.comparison),aes(x=covariate,fill=variable)) + geom_bar(aes(weight=value),position="dodge") + labs(x="",y=expression(beta)) + scale_fill_discrete(name="method")

# fix the x-axis labels
cmd2 <- paste("coef.comparison.plot <- coef.comparison.plot + scale_x_discrete(labels=c(",
      paste("\"",ridge.indv.study.df$X1,"\" = ",llply(ridge.indv.study.df$X1,.fun=function(cov){parse(text=cov)}),sep="",collapse=",")
      ,"))",sep="")
eval(parse(text=cmd2))
show(coef.comparison.plot)

### check out bias, variance, and MSE as functions of lambda
melted.df <- melt(data=ridge.stats,id.vars=c("stat","lambda")); head(melted.df)
### dig through this lapply to figure out what's up.
lapply(X=1:J,FUN=function(i){
  ### double check that this is returning the min.
  coef <- paste("beta[",i,"]",sep="")
  relev.subset <- subset(x=melted.df,subset={variable==coef})
  lambda.min <- with(data=subset(x=relev.subset,subset={stat=="MSE"}),expr={
    lambda.values[which.min(value)]
  })
  ### calculate beta.ridge corresponding to lambda.min
  beta.ridge <- .calculateRidgeCoefs(X=learning.set[,1:J],Y=learning.set[,J+1],lambda=lambda.min)[i]
  plot.title <- substitute(atop(paste("bias, variance, and MSE for ", hat(beta)[i]), paste("minimal MSE at ", group("(",list(lambda,hat(beta)[i]),")") == group("(",list(lambda.min,beta.ridge),")"))),list(i=i,lambda.min=lambda.min,beta.ridge=beta.ridge))
  
  output <- ggplot(data=relev.subset,aes(x=lambda,y=value,group=stat:variable,lty=stat)) + geom_line() + geom_vline(xintercept=lambda.min,color='red') + labs(x=expression(lambda),y="",title=plot.title)
  show(output)
  return(list(plot=output,lambda.min=lambda.min,beta.ridge=beta.ridge))
})

ggplot(data=melted.df,aes(x=lambda,y=value,color=variable)) + geom_line(aes(group=stat:variable,lty=stat),show_guide=T) #+ facet_wrap(facets=~variable,nrow=4)#,labeller=label_parsed) 

### -------------------
### code hell
### -------------------

# ridge.sampling.stats <- lapply(X=1:dim(ridge.sampling.stats.list.mat)[1],FUN=function(row.i){
#   row <- ridge.sampling.stats.list.mat[row.i,]
#   stat.name <- rownames(ridge.sampling.stats.list.mat)[row.i]
#   mat <- matrix(unlist(row),ncol=J,byrow=T,dimnames=list(lambda.values,paste("beta[",1:J,"]",sep="")))
#   return(data.frame(stat=stat.name,lambda=lambda.values,mat,check.names=F))
# })
# ridge.bias.df <- ridge.sampling.stats[[1]]
# ridge.var.df <- ridge.sampling.stats[[2]]
# ridge.mse.df <- ridge.sampling.stats[[3]]
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
