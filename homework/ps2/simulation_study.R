### question 3 -- simulation study
library(MASS)
library(ggplot2)

ar1CovMatrix <- function(n,rho) {
  ### there's got to be a slick way to do this.
  rows <- 1:n
  exp.mat <- matrix(data=sapply(X=rows,FUN=function(i){abs(i-rows)}),ncol=n,byrow=T)
  return(rho ^ exp.mat)
}
