
## ----globalParameters,echo=FALSE,cache=FALSE,message=FALSE---------------
opts_chunk$set(comment="#", cache=T, echo=F, tidy=F, warning=FALSE, message=FALSE, highlight=T, autodep=T)
library(doMC)
registerDoMC(cores=detectCores()) # this could break someone's shit.


## ----stayInside,echo=FALSE-----------------------------------------------
  options(width=60)
  listing <- function(x, options) {
     paste("\\begin{lstlisting}[basicstyle=\\ttfamily,breaklines=true]\n", x,"\\end{lstlisting}\n", sep = "")
  }
  
  knit_hooks$set(source=listing, output=listing)


## ----simulationSetup-----------------------------------------------------
### question 3 -- simulation study

### ------------------------
### simulation libraries
### ------------------------
set.seed(1234)
library(MASS) # for multivariate normal random variables
library(glmnet) # for regression routines
library(ggplot2) # for plotting
library(GGally) # ggplot2 equiv of pairs
library(grid)
library(reshape2)
library(plyr)
library(xtable)


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
.makeAndAnalyzeGlmnet <- function(learning.set, training.set, regularization="ridge",standardize=F,intercept=T) {
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

### make learning and test set
learning.set <- .generateDataSet(n.learning.set,J,rho)
training.set <- .generateDataSet(n.training.set,J,rho)


## ----prob3aEDA-----------------------------------------------------------
pairs.plot <- ggpairs(data.frame(learning.set[,1:J]),upper=list(continuous='cor'),alpha=0.5)


## ----prob3aPairPlot,out.width="0.85\\textwidth",fig.align='center',fig.cap="Paired scatterplot for learning set covariates.",dev='pdf'----
show(pairs.plot)


## ----studyCov,results='asis'---------------------------------------------
sim.cov <- xtable(cov2cor(.ar1CovMatrix(J,rho)),align=rep("",J+1),
                  digits=4,
                  label="eq:simCov")
print(sim.cov,tabular.environment="pmatrix",include.rownames=F,include.colnames=F,hline.after=NULL,latex.environments=c("equation"))


## ----prob3aCorPlot,out.width="0.49\\textwidth",fig.align='left',fig.cap="Heatmap for correlation structure of learning set (sample) juxtapose to theoretical correlation structure.",dev='pdf',fig.keep='all', fig.show='hold'----
makeCorrelationPlot <- function(cor.mat,show.guide=T,plot.title="") {
  require(grid)
  require(reshape2)
  melted.cor <- melt(cor.mat)
  colnames(melted.cor)[1:2] <- c("Var1", "Var2")
  output <- ggplot(data=melted.cor, aes(x=Var1,y=Var2,fill=value)) + geom_tile() + labs(x="",y="",title=plot.title) + theme_bw() + theme(legend.position="bottom") + theme(legend.key.width=unit(2,"cm")) 
  
  output <- if (show.guide) {
    output + scale_fill_continuous(name="Correlation",limits=c(-1,1),breaks=seq(-1,1,length.out=5))
  } else {
    output + scale_fill_continuous(guide = show.guide, name="Correlation",limits=c(-1,1),breaks=seq(-1,1,length.out=5))
  }
  
  # fix the axis labels, variable names are nonsense
  tmp <- paste("\"",colnames(cor.mat),"\" = ",llply(colnames(cor.mat),.fun=function(cov){parse(text=cov)}),sep="",collapse=",")
  cmd2 <- paste("output <- output + scale_x_discrete(limits = colnames(cor.mat), labels=c(", tmp, ")) + scale_y_discrete(limits=rev(colnames(cor.mat)), labels=c(", tmp, "))",sep="")
  eval(parse(text=cmd2))
  return(output)
}

cov.X <- .ar1CovMatrix(J,rho); cov.Y <- cov.X %*% beta
theor.cov.mat <- rbind(cbind(cov.X,cov.Y),c(cov.Y,sigma^2))
theor.cor.mat <- cov2cor(theor.cov.mat)
colnames(theor.cor.mat) <- colnames(learning.set)
theor.cor.plot <- makeCorrelationPlot(theor.cor.mat,T,"Theoretical correlation")
sample.cor.plot <- makeCorrelationPlot(cor(learning.set),T,"Sample correlation")
show(sample.cor.plot)
show(theor.cor.plot)


## ----prob3aJointDensity--------------------------------------------------
require(grid)
sim.df <- data.frame(rbind(cbind(learning.set,set="learning"),cbind(training.set,set="test")),check.names=F)
melted.sim.df <- melt(sim.df,id.vars=c("Y","set"))
melted.sim.df <- transform(melted.sim.df,Y=as.numeric(as.character(Y)),value=as.numeric(as.character(value)))
joint.density.plots <- ggplot(data=melted.sim.df,aes(y=Y,x=value,group=set:variable)) + geom_point(aes(shape=set),alpha=0.25) + stat_smooth(method="loess",se=F,n=100,aes(color=set)) + facet_wrap(facets=~variable,ncol=5,scales="free_x") + theme(legend.position="bottom") + theme(legend.key.width=unit(2,"cm")) + labs(x=expression(paste("covariate value, ", X[i],sep="")), y=expression(paste("Response, ", Y, sep=""))) + scale_shape_discrete(name="",label=c("learning"="Learning set", "test"="Test set")) + scale_color_discrete(name="", label=c("learning"="Learning set", "test"="Test set"))


## ----prob3aJointDensityPlot,out.width="0.85\\textwidth",fig.align='center',fig.cap="Joint density plots of $Y$ against $X_{i}$ with loess applied to both learning and test sets.",dev='pdf',dependson="prob3aJointDensity"----
show(joint.density.plots)


## ----prob3aBoxPlots,out.width="0.85\\textwidth",fig.align='center',fig.cap="Box plots for learning and test set variables.",dev='pdf',dependson="prob3aJointDensity"----
require(grid)
melted.sim.df <- melt(sim.df, id.vars=c("set"))
melted.sim.df <- transform(melted.sim.df, value=as.numeric(as.character(value)))
boxplots <- ggplot(data=melted.sim.df, aes(x=variable,y=value,group=set:variable,color=set)) + geom_boxplot() + scale_color_discrete(name="", label=c("learning"="Learning set", "test"="Test set")) + theme(legend.position="bottom") + theme(legend.key.width=unit(2,"cm")) + labs(x="", y="")

parseAxisLabels <- function(raw.axis.labels, plot.obj, axis="x") {
  output <- plot.obj
  label.string <- paste("\"",raw.axis.labels,"\" = ",
                          llply(raw.axis.labels,.fun=function(tick.mark){parse(text=tick.mark)}),sep="",collapse=",")
  switch(axis,
         x = {
             output.update <- paste("output <- output + scale_x_discrete(limits = raw.axis.labels, labels=c(", label.string, "))",sep="")             
         },
         y = {
              output.update <- paste("output <- output + scale_y_discrete(limits = raw.axis.labels, labels=c(", label.string, "))",sep="")  
         })
  eval(parse(text=output.update))
  return(output)
}
parseAxisLabels(levels(melted.sim.df$variable),boxplots)


## ----prob3bDataGeneration------------------------------------------------
### make coefficients for each regularization scheme
ridge.study <- .makeAndAnalyzeGlmnet(learning.set, training.set)
ridge.study.df <- ridge.study$data.frame

lasso.study <- .makeAndAnalyzeGlmnet(learning.set, training.set, regularization="lasso")
lasso.study.df <- lasso.study$data.frame

enet.study <- .makeAndAnalyzeGlmnet(learning.set, training.set, regularization="elasticnet")
enet.study.df <- enet.study$data.frame


### find lambda associated to minimal traning risk
ridge.t.MSE.min.indx <- which.min(ridge.study.df$training.set.MSE)
enet.t.MSE.min.indx <- which.min(enet.study.df$training.set.MSE)
lasso.t.MSE.min.indx <- which.min(lasso.study.df$training.set.MSE)

location.of.MSE.mins <- rbind(ridge.study.df[ridge.t.MSE.min.indx,],enet.study.df[enet.t.MSE.min.indx,],lasso.study.df[lasso.t.MSE.min.indx,])


## ----prob3b-cPlots-------------------------------------------------------
### reshape data for visualization
combined.df <- rbind(ridge.study.df,enet.study.df,lasso.study.df)
melted.df <- melt(data=combined.df,id.vars=c("lambda","df","regularization"))
mse.indx <- which(melted.df$variable %in% c("learning.set.MSE","training.set.MSE"))

mse.vs.lambda <- ggplot(data=melted.df[mse.indx,],aes(x=lambda,y=value,color=regularization,lty=variable)) + geom_line() + geom_vline(data=location.of.MSE.mins,aes(xintercept=lambda,color=regularization),lty=8,alpha=0.75) + labs(x=expression(lambda),y="MSE")

### look at how coefficients vary with lambda
beta.vs.lambda <- ggplot(data=melted.df[-mse.indx,],aes(x=lambda,y=value,color=variable)) + geom_line(show_guide=F) + geom_hline(yintercept=beta,alpha=0.5,lty=3) +  geom_vline(data=location.of.MSE.mins,aes(xintercept=lambda),lty=8,alpha=0.5) + geom_text(data=subset(x=melted.df[-mse.indx,],subset={lambda==0}),aes(label=variable,x=0,y=value,color=variable,group=regularization),show_guide=F,hjust=1,parse=T) + facet_grid(regularization~.) + labs(x=expression(lambda),y=expression(beta))

### compare "optimal" coefficients between all 4 methods
beta.OLS <- .calcOLSCoefs(X=learning.set[,1:J],Y=learning.set[,J+1])
extra.info <- rbind(data.frame(regularization="OLS",lambda=0,df=J,t(beta.OLS),learning.set.MSE=0,training.set.MSE=0,check.names=F),location.of.MSE.mins)

coefs.vs.method <- ggplot(data=melt(extra.info[,-(14:15)],id.vars=c("lambda","df","regularization")),aes(x=regularization,y=value,color=variable,group=variable)) + geom_point(show_guide=F) + geom_line(show_guide=F) + geom_hline(yintercept=beta,alpha=0.5,lty=3) + geom_text(data=melt(extra.info[1,-(14:15)],id.vars=c("lambda","df","regularization")),aes(label=variable,x=as.factor("OLS"),y=value),hjust=1.25,vjust=0,show_guide=F,parse=T) + labs(x="method",y=expression(beta))

### look at beta vs. df
beta.vs.df <- ggplot(data=melted.df[-mse.indx,],aes(x=df,y=value,color=variable)) + geom_line(show_guide=F) + geom_point(data=melt(location.of.MSE.mins[,-(14:15)],id.vars=c("lambda","df","regularization")),aes(x=df,y=value),size=3,shape=2,show_guide=F) + geom_text(data=subset(x=melted.df[-mse.indx,],subset={df==10 & lambda==0}),aes(label=variable,x=10,y=value,color=variable,group=regularization),show_guide=F,hjust=-.05,parse=T) +  facet_grid(regularization~.) + labs(x="Effective Degrees of Freedom",y=expression(beta))

### check out df vs. lambda
df.vs.lambda <- ggplot(data=melted.df[-mse.indx,],aes(x=lambda,y=df,color=regularization)) + geom_line() + geom_vline(data=location.of.MSE.mins,aes(xintercept=lambda,color=regularization),lty=8,alpha=0.75)  + labs(x=expression(lambda),y="Effective Degrees of Freedom")


## ----prob3plot5,out.width='0.85\\textwidth',dev='pdf',fig.align='center',fig.cap="Effective degrees of Freedom against $\\lambda$. Vertical lines indicate location of $\\lambda$ which minimize test set MSE."----
show(df.vs.lambda)


## ----prob3plot2,out.width='0.85\\textwidth',dev='pdf',fig.align='center',fig.cap="Behavior of coefficient estimates across $\\lambda$ for each method. Horizontal lines indicate true coefficient values, and vertical lines indicate $\\lambda$ which minimze test set MSE."----
show(beta.vs.lambda)


## ----prob3plot1,out.width='0.85\\textwidth',dev='pdf',fig.align='center',fig.cap="MSE vs. $\\lambda$ for training and learning sets. Vertical lines indicate location of minima."----
### look at MSE vs. lambda for training and learning sets
show(mse.vs.lambda)


## ----prob3plot3,out.width='0.85\\textwidth',dev='pdf',fig.align='center',fig.cap="Optimal coefficient estimates for each method. Horizontal lines indicate true coefficient values."----
show(coefs.vs.method)


## ----prob3plot4,out.width='0.85\\textwidth',dev='pdf',fig.align='center',fig.cap="Value of coefficients against Effective degrees of Freedom. Triangle points indicate optimal coefficient values."----

show(beta.vs.df)


## ----demo,eval=FALSE,echo=TRUE,tidy=TRUE,highlight=TRUE------------------
## ### coefs.MSE.df is a length(lambda) x J data.frame
## ### coefs.MSE.df[i,j] has the MSE for the j-th
## ### coefficient estimate at lambda[i]
## min.MSE.indices <- apply(coefs.MSE.df, MARGIN=2, 'which.min')
## ### coefs.df is a length(lambda) x J data.frame
## ### coefs.df[i,j] has the j-th coef estimate at
## ### lambda[i]
## aggregated.estimator <- coefs.df[cbind(min.MSE.indices,1:J)]


## ----prob3dplots---------------------------------------------------------
### investigate properties of Ridge estimator;
# For the simulation model of a), provide and comment on graphical displays of
# the bias, variance, and MSE of the ridge estimators based on the learning set.
# For each coefficient, provide the value of the shrinkage parameter λ minimizing
# the MSE and the corresponding estimate.

### use closed form expression to investigate bias, variance, and MSE of beta_ridge
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
mse.comparison <- rbind(c(Learning=frank.MSE.lset,Test=frank.MSE.tset),c(Learning=location.of.MSE.mins$learning.set.MSE[1],Test=location.of.MSE.mins$training.set.MSE[1]))
mse.comparison <- data.frame(method=c("aggregated","argmin"),mse.comparison)

mse.comparison.plot <- ggplot(data=melt(mse.comparison),aes(x=method)) + geom_bar(aes(weight=value,fill=variable),position="dodge") + geom_text(aes(x=method,y=value,label=round(value,digits=3)),hjust=c(1.5,-0.5,1.5,-0.5),vjust=-0.2,color='black',position="dodge") + labs(y="MSE") + scale_fill_discrete(name="")

### compare frankenstein coefficients to argmin coefs
coef.comparison <- data.frame(ridge.indv.study.df[,1:2],t(location.of.MSE.mins[1,4:13]),beta)
colnames(coef.comparison) <- c("covariate","aggregated","argmin","truth")

coef.comparison.plot <- ggplot(data=melt(coef.comparison),aes(x=covariate,fill=variable)) + geom_bar(aes(weight=value),position="dodge",show_guide=F) + labs(x="",y=expression(beta)) + scale_fill_discrete(name="method")

# fix the x-axis labels
cmd2 <- paste("coef.comparison.plot <- coef.comparison.plot + scale_x_discrete(limits = coef.comparison$covariate, labels=c(",
      paste("\"",ridge.indv.study.df$X1,"\" = ",llply(ridge.indv.study.df$X1,.fun=function(cov){parse(text=cov)}),sep="",collapse=",")
      ,"))",sep="")
eval(parse(text=cmd2))

### check out bias, variance, and MSE for each coefficient as functions of lambda
melted.ridge.stats <- melt(data=ridge.stats,id.vars=c("stat","lambda"))
coef.plots <- dlply(.data=melted.ridge.stats,.variables=.(variable),.fun=function(subset){
  ### subset has the viz data, minima has the annotation meta-data
  minima <- subset(ridge.indv.study.df, {X1 == subset$variable[1]})
  i <- as.numeric(gsub(pattern="X\\[([1]{0,1}[0-9])\\]",replacement="\\1",x=minima$X1))
  plot.title <- substitute(paste("Minimal MSE = ", MSE, " at ", lambda, " = ", lambda.val, " with ", hat(beta)[i], " = ", beta.val),list(MSE=round(minima$MSE,4), lambda.val=minima$lambda, i=i, beta.val=round(minima$coef,3)))
  output <- ggplot() + geom_vline(xintercept=minima$lambda,color='black',lty=3,alpha=0.65) + geom_hline(yintercept=minima$MSE,color='black',lty=3,alpha=0.65) + geom_line(data=subset,aes(x=lambda,y=value,group=stat:variable,color=stat),show_guide=F)  + geom_point(data=minima,aes(x=lambda,y=MSE),shape=4,color='red',size=2) + labs(y="",x=expression(lambda), title=plot.title)
#   show(output)
  return(output)
})


## ----optimalCoefsOut,tidy=TRUE,results='asis'----------------------------
convert.to.tex <- function(str){gsub(pattern="X\\[([01]{0,1}[0-9]{1})\\]",replacement="$X_{\\1}$",x=str)}
tbl.out <- xtable(ridge.indv.study.df, align="lcccc",
                  label="tab:ridgeOptimal",
                  caption="Locally optimal (ridge) coefficient estimates.")
colnames(tbl.out) <- c("Covariate", "Est. Coef.", "$\\lambda$", "MSE")
digits(tbl.out) <- c(0,0,4,0,4)
print(tbl.out, latex.environments=c("center"),
      sanitize.text.function=convert.to.tex,
      floating=T,
      include.rownames=F)


## ----prob3plot7,out.width='0.85\\textwidth',dev='pdf',fig.align='center',fig.cap="Coefficient comparison for ``aggregated'' (red) and $\\argmin$ (green) estimators, next to true coefficient values (blue)."----
show(coef.comparison.plot)


## ----prob3plot6,out.width='0.5\\textwidth',dev='pdf',fig.align='center',fig.cap="MSE comparison for ``aggregated'' and $\\argmin$ estimators"----
show(mse.comparison.plot)


## ----prob3plot8,out.width='0.3\\textwidth',dev='pdf',fig.align='left',fig.cap="Bias (red), Variance (green), and MSE (blue) for each ridge coefficient estimate. Dashed lines indicate location and value of MSE's minimum.",fig.show='hold',fig.keep='all'----
l_ply(.data=coef.plots,.fun='show')


## ----prob3eDataGen-------------------------------------------------------
### investigate properties of LASSO via bootstrap
###   fix learning set, redraw Y
B <- 1000 # number of redraws
X.fixed <- learning.set[,1:J]

p.time2 <- system.time(lasso.sampling.dist <- ldply(.inform=F,.parallel=T,.data=1:B,.fun=function(i){
  ### remake learning set, and run glmnet
  learning.set[,J+1] <- rnorm(n=n.learning.set, mean=X.fixed %*% beta, sd=sigma)
  out <- .makeAndAnalyzeGlmnet(learning.set, training.set, regularization="lasso")$data.frame
  return(out)
}))


## ----eval=FALSE,echo=TRUE,tidy=TRUE,highlight=TRUE-----------------------
## Y.star <- rnorm(n=n.learning.set, mean = X.n %*% beta, sd = sigma*diag(n.learning.set) )


## ----timeOut,eval=FALSE--------------------------------------------------
## p.time2


## ----prob3eDataProcessing,dependson="prob3eDataGen"----------------------
### this is doing some fucked up shit. Need to make sure lambda is cast as numeric
lasso.sampling.dist$lambda <- as.numeric(as.character(lasso.sampling.dist$lambda))
### for each covariate, find the average coefficient estimate, as well as 
### bias, variance, and MSE
lasso.analysis.coef.stats <- daply(.data=lasso.sampling.dist,.variables=.(lambda),.fun=function(subset){  
  covariates <- subset[,4:13]
  exp.val <- colMeans(covariates) # 4:13 are coefs
  bias <- exp.val - beta
  var.hat <- aaply(.data=covariates,.margins=2,.fun='var')
  mse.hat <- bias^2 + var.hat
  return(rbind(exp.val[1:10],bias,var.hat,mse.hat))
})
### for each covariate, find lambda which corresponds to minimal MSE
lasso.analysis.optimal.coefs <- adply(.data=lasso.analysis.coef.stats,.margins=3,.fun=function(covariate){
  min.lambda.indx <- which.min(covariate[,4])
  return(c(covariate[min.lambda.indx,],lambda=lambda.values[min.lambda.indx]))
})

### make bar plot depicting optimal coefficients versus truth
coef.df <- data.frame(lasso.analysis.optimal.coefs,truth=beta)
melted.df<-melt(coef.df,id=c("X1","lambda"))

lasso.comparison.plot <- ggplot(data=subset(melted.df,{variable %in% c("V1","truth")}), aes(x=X1,weight=value,fill=variable)) + geom_bar(position="dodge",show_guide=F) + scale_fill_discrete(name="",breaks=c("V1","truth"),labels=c(expression(hat(beta)),expression(beta[0]))) + labs(x="", y=expression(beta))

# fix the x-axis labels
cmd2 <- paste("lasso.comparison.plot <- lasso.comparison.plot + scale_x_discrete(limits = lasso.analysis.optimal.coefs$X1, labels=c(",
              paste("\"",lasso.analysis.optimal.coefs$X1,"\" = ",llply(as.character(lasso.analysis.optimal.coefs$X1),.fun=function(cov){parse(text=cov)}),sep="",collapse=",")
              ,"))",sep="")
eval(parse(text=cmd2))


## ----optimalCoefsOut2,echo=FALSE,tidy=TRUE,results='asis'----------------
tbl2.out <- xtable(lasso.analysis.optimal.coefs, align="lcccccc",
                  label="tab:lassoOptimal",
                  caption="Locally optimal (LASSO) coefficient estimates.")
colnames(tbl2.out) <- c("Covariate", "$\\hat{\\beta}_{i}$", "$\\hat{\\text{bias}}$", "$\\widehat{\\sigma}^2$", "$\\hat{\\text{MSE}}$", "$\\lambda$")
digits(tbl2.out) <- c(0,0,4,4,4,4,0)
print(tbl2.out, latex.environments=c("center"),
      sanitize.text.function=convert.to.tex,
      floating=T,
      include.rownames=F)


## ----prob3plot9,out.width='0.3\\textwidth',dev='pdf',fig.align='left',fig.cap="Bias (red), Variance (green), and MSE (blue) for each LASSO coefficient estimate. Dashed lines indicate location and value of MSE's minimum.",fig.show='hold',fig.keep='all'----
### make plots similar to before
lasso.plots <- llply(.data=dimnames(lasso.analysis.coef.stats)[[3]],.fun=function(covariate){
  ### melted.df has the viz data, minima has the annotation meta-data
  minima <- subset(lasso.analysis.optimal.coefs, {X1 == covariate})
  i <- as.numeric(gsub(pattern="X\\[([1]{0,1}[0-9])\\]",replacement="\\1",x=minima$X1))
  plot.title <- substitute(paste("Minimal MSE = ", MSE, " at ", lambda, " = ", lambda.val, " with ", hat(beta)[i], " = ", beta.val),list(MSE=round(minima$mse.hat,4), lambda.val=minima$lambda, i=i, beta.val=round(minima$V1,3)))
  
  df <- data.frame(lasso.analysis.coef.stats[,,covariate],lambda=lambda.values)
  # reshape data
  colnames(df)[1] <- "coefficient"
  melted.df <- melt(df,id.vars="lambda")
  stat.plot <- ggplot() + geom_vline(xintercept=minima$lambda,color='black',lty=3,alpha=0.65) + geom_hline(yintercept=minima$mse.hat,color='black',lty=3,alpha=0.65) + geom_point(data=minima,aes(x=lambda,y=mse.hat),shape=4,color='red',size=2) + geom_line(data=subset(melted.df,{variable != "coefficient"}), aes(x=lambda,y=value,color=variable), show_guide=F) + labs(y="",x=expression(lambda),title=plot.title) 
  show(stat.plot)
  return(stat.plot)
})


## ----prob3plot10,out.width='0.75\\textwidth',dev='pdf',fig.align='center',fig.cap="Optimal LASSO estimates (red) next to true parameter values (green)."----
show(lasso.comparison.plot)


