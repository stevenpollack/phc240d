# In particular, the matrix doxorubicinNCI60Scaled provides Affymetrix HG-U95Av2 expression measures for the training set of 22 cell lines from the NCI60 panel, classified as “Resistant” or “Sensitive” to doxorubicin. The column names are the NCI60 cell line names; their drug sensitivity status can be found in cellLinesUsed (use cellLinesUsed$doxorubicin$listPotti06CorrAug08version). The row names are Affymetrix gene IDs. The data frame doxorubicin07Numbers provides another version of the same data for the training set of 22 cell lines, as well as Affymetrix HG-U133Av2 measures for a test set of 122 acute lymphoblastic leukemia (ALL) patients. The data frame doxorubicin07Info provides the training/test set membership of the samples and their drug sensitivity status, but not the actual names for the 22 cell lines comprising the training set and corresponding to the columns of doxorubicinNCI60Scaled. The row names are Affymetrix gene IDs; LocusLink IDs were used to match probes across the two Affymetrix platforms.
# Note that, as indicated in Baggerly and Coombes (2009): The data were log-transformed, the values for each row (gene) were centered and scaled to have mean zero and unit variance (separately for training and test data), the data were exponentiated to undo the log-transform, and the final results were rounded to two decimal places. As microarray expression measures are based on fluorescence intensities (measured on a 16-bit scale), it may be appropriate to re-log-transform the data.

load("assignments/examiningDoxorubicinInDetail.RData")
library(Biobase)

# Question 1. Examine and prepare datasets. Examine the different objects in the R dataset examiningDoxorubicinInDetail.RData. Store the expression measures and the sample and gene annotation data related to doxorubicinNCI60Scaled and doxorubicin07Numbers in objects of class ExpressionSet (Bioconductor R package Biobase).

### format NCI60Scaled data into ESet
NCI60AssayData <- doxorubicinNCI60Scaled

### reshape sensitivty data to a form
### that can be used as phenoData in ESet
resistanceStatus <- cellLinesUsed$doxorubicin$listPotti06CorrAug08
statusVec <- rep("resistant",dim(NCI60AssayData)[2])
statusVec[ colnames(NCI60AssayData) %in% resistanceStatus$Sensitive ] <- "sensitive"

NCI60PhenoData <- data.frame(status=as.factor(statusVec),row.names=colnames(NCI60AssayData))

NCI60ESet <- ExpressionSet(assayData=NCI60AssayData,annotation="Affymetrix HG-U95Av2",phenoData=new("AnnotatedDataFrame", data=NCI60PhenoData))

### format 07Numbers data into ESet
doxo07AssayData <- as.matrix(doxorubicin07Numbers) # ESet() can't handle data.frame's
doxo07PhenoData <- new("AnnotatedDataFrame", data=doxorubicin07Info)
doxo07ESet <- ExpressionSet(assayData=doxo07AssayData,phenoData=doxo07PhenoData,annotation="Affymetrix HG-U133Av2")

# Question 2. QA/QC for training set. Reconcile the training data in doxorubicinNCI60Scaled and doxorubicin07Numbers, i.e., match genes and samples and compare the expression measures and the sensitivity status assigned to the cell lines in the two datasets.

### First, dim(doxo07ESet)[1] < dim(NCI60ESet)[1],
### hence we check if all genes in doxo07ESet were
### also measured in NCI60ESet
all(rownames(doxo07ESet) %in% rownames(NCI60ESet)) #TRUE

### next we line them up: 
### doxo07ESet[i,] == NCI60ESet[indexOf...[i],]
indexOfRowInNCI60 <- match(rownames(doxo07ESet),rownames(NCI60ESet)) #sapply(X=rownames(doxo07ESet),FUN=function(geneName){which(rownames(NCI60ESet) %in% geneName)})
rowMap <- data.frame(NCI60index=indexOfRowInNCI60,row.names=rownames(doxo07ESet))
featureData(doxo07ESet) <- new("AnnotatedDataFrame",data=rowMap)

### example usage
NCI60and07GeneESet <- NCI60ESet[fData(doxo07ESet)$NCI60index,]

### check out how well the columns line up
diffs <- sapply(X=seq(0,0.06,0.01),FUN=function(diff) {
  sum(abs(exprs(doxo07ESet[,1:22]) - exprs(NCI60and07GeneESet)) > diff)
})
names(diffs) <- seq(0,0.06,0.01)

### report percent of cells which aren't
### equal, by size of difference
diffs / (dim(doxo07ESet)[1]*22) * 100 

### there are no expressions that differ by more than 0.06

### reconcile column names:
set.seed(1234)
random.row <- sample(1:8958,size=1)
exprs(NCI60and07GeneESet[rownames(doxo07ESet)[random.row],]) - exprs(doxo07ESet[random.row,1:22])

### clearly, they line up...
sampleNames(doxo07ESet)[1:22] <- sampleNames(NCI60and07GeneESet)

### investigate resistances;
all(rownames(pData(doxo07ESet)[1:22,]) == rownames(pData(NCI60and07GeneESet))) # T
### row names line up, so check statii
levels(pData(NCI60and07GeneESet)$status) <- levels(pData(doxo07ESet)[1:22,2])
which(pData(doxo07ESet)[1:22,2] != pData(NCI60and07GeneESet)$status)

# all of them are wrong!

#Question 3. QA/QC for test set. Consider now the test data in doxorubicin07Numbers. Is there anything unusual with the samples and their assigned sensitivity stati? Hint: Use dimensionality reduction and clustering methods.

doxo07TestESet <- doxo07ESet[,-(1:22)]

### add 0.1 to data, log it, and transpose it
transTestData <- t(log(exprs(doxo07TestESet)+0.1))

### look at 2-cluster of testData 
library(cluster)
pam.results <- pam(transTestData,2)

### consider the stati of the medoids
doxo07TestESet$status[pam.results$id.med]
### both are resistant!

### table the clusters to get an idea of heterogeneity
table(doxo07TestESet$status[which(pam.results$clustering == 1)])
# 90 resistant, 19 sensitive

table(doxo07TestESet$status[which(pam.results$clustering == 2)])
# 9 resistant, 4 sensitive

### not very compelling evidence to show that genetic signatures => resistance
### consider PCA
pca.results <- prcomp(transTestData)

### kaiser criterion
max(which(pca.results$sdev^2 >= 1)) # keep the first 78

### consider how many until 90% of variation is obtained
min(which(cumsum(pca.results$sdev)/sum(pca.results$sdev)>=0.9)) # keep first 66

library(ggplot2)
ggplot(data=as.data.frame(pca.results$x),aes(x=PC1, y=PC2,color=doxo07TestESet$status)) + geom_text(aes(label=rownames(pca.results$x)),alpha=0.5,size=4)
### you can almost make the case that PC1 does an okay job
### creating a classification criterion, but even then...


### biplot is unuseful. Too many features, and no real separating inside the PC1-PC2 plane
PCbiplot <- function(PC, x="PC1", y="PC2") {
  require(ggplot2)
  # PC being a prcomp object
  scores <- PC$x[,1:6]
  loadings <- PC$rotation[1:10,1:6]
  data <- data.frame(obsnames=row.names(scores), scores)
  plot <- ggplot(data, aes_string(x=x, y=y)) + geom_text(alpha=.4, size=3, aes(label=obsnames))
  plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  datapc <- data.frame(varnames=rownames(loadings), loadings)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_text(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1, color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), alpha=0.75, color="red") # arrow=arrow(length=unit(0.2,"cm")),
  plot
}
PCbiplot(pca.results)
# need first 47 components for 95% of variation


# Question 4. Differential expression analysis for training set. The goal of this question is to derive, as in Potti et al. (2006), a “signature” of sensitivity to doxorubicin, which could be used eventually to predict patient response to the drug. This involves identifying genes that are differentially expressed (DE) between “Resistant” and “Sensitive” cell lines based on the training set of microarray expression measures for the 22 NCI60 cell lines. Use the expression measures stored in doxorubicinNCI60Scaled and the drug sensitivity status given by cellLinesUsed$doxorubicin$listPotti06CorrAug08.

# a) Between-sample normalization. In order to derive expression measures and compare these measures between samples and/or genomic regions of interest (ROI), one first needs to normalize the fluorescence intensities to adjust for a variety of sample-level and gene-level technical effects, such as, library preparation/hybridization/scanning effects and nucleotide composition effects (e.g., GC-content). Normalization is essential to ensure that observed differences in expression measures between samples and/or ROI are truly due to differential expression and not experimental/technical artifacts. For the purpose of this assignment, it is sufficient to perform between-sample normalization using standard methods for Affymetrix (one-channel) microarrays, e.g., loess and full-quantile procedures implemented in Bioconductor R packages such as affy and limma. Provide and comment on numerical and graphical summaries of the data that suggest the need for normalization. Perform between-sample normalization and comment on the results.

### do the same transformation to training data
### but don't transpose it. We want to look at
### inter-column differences...
transTrainData <- log(exprs(NCI60ESet)+0.1)
transTrainData <- as.data.frame(transTrainData)

### look at mean-difference plots for a few lines
mdPlots <- apply(X=combn(x=22,m=2),MARGIN=2,FUN=function(cols){
  cell1 <- transTrainData[,cols[1]]
  cell1.name <- colnames(transTrainData)[cols[1]]
  cell2 <- transTrainData[,cols[2]]
  cell2.name <- colnames(transTrainData)[cols[2]]
  title <- paste("Mean-Difference plot for X=", cell1.name, " and Y=", cell2.name, " cell lines", sep="")
  mdPlot <- ggplot(data=data.frame(X=cell1,Y=cell2),aes(x=(X+Y)/2, y=Y-X)) + geom_point(alpha=0.25) + geom_hline(color='red',yintercept=0) + labs(title=title)
})

plots <- sample(x=1:choose(22,2),size=4,replace=F) # pick 4 random plots
show(mdPlots[plots])

### all plots seem to exhibit a high amount of dispersion

library(reshape2)
melted.df <- melt(transTrainData)
melted.df <- transform(melted.df, value=as.numeric(value))

ggplot(data=melted.df,aes(x=variable,y=value)) + stat_boxplot() + geom_hline(yintercept=0,color='red') + scale_color_discrete(guide=F) + labs(x="cell line", y="log(intensity)",title="non-normalized expression values") 

### the distributions vary wildly between samples

### get affy and do loess normalization
#source("http://bioconductor.org/biocLite.R")
biocLite("affyPLM")

library(affyPLM)
transTrainESet <- ExpressionSet(as.matrix(transTrainData))
normTrainESet <- normalize.ExpressionSet.quantiles(transTrainESet,transfn="none")

### couldn't run loess-based method for some reason.

### look at box plot
normTrainData <- as.data.frame(exprs(normTrainESet))
melted.df2 <- melt(normTrainData)
melted.df2 <- transform(melted.df2, value=as.numeric(value))

ggplot(data=melted.df2,aes(x=variable,y=value)) + stat_boxplot() + geom_hline(yintercept=0,color='red') + scale_color_discrete(guide=F) + labs(x="cell line", y="log(intensity)",title="quantile-based normalized expression values") 

# b) Cluster analysis. Apply dimensionality reduction and clustering (both                                                                hierarchical and partitioning) methods to the training set samples. Comment on the results and, in particular, relate the clustering to drug sensitivity.

pam.results2 <- pam(k=2,x=t(normTrainData))

### consider the stati of the medoids
NCI60ESet$status[pam.results2$id.med]
### great, medoids are of different stati

### table the clusters to get an idea of heterogeneity
table(NCI60ESet$status[which(pam.results2$clustering == 1)])
# not good, 6 res, 8 sensitive

table(NCI60ESet$status[which(pam.results2$clustering == 2)])
# 6 resistant, 2 sensitive

### so even inside the training data,
### we can't cluster on sensitivity

### do PCA
pca.results2 <- prcomp(t(normTrainData))
ggplot(data=as.data.frame(pca.results2$x),aes(x=PC1, y=PC2,color=NCI60ESet$status)) + geom_text(aes(label=rownames(pca.results2$x)),alpha=0.75,size=4) + geom_abline(intercept=0.5,slope=0.4,lty=2)

### seem we may be able to classify based on position in PC1-PC2 plane
hclust.results <- hclust(d=dist(t(normTrainData)))
#plot(hclust.results)

### color dendrogram
reds <<- as.factor(colnames(normTrainData)[ NCI60ESet$status == "resistant"]) 
greens <<- as.factor(colnames(normTrainData)[ NCI60ESet$status == "sensitive"]) 

#define a function for coloring and sizing node elements: 
colLab <- function(n) { 
  if(is.leaf(n)) { 
    a <- attributes(n) 
    if ( length(which(greens == a$label)) == 1 ) { 
      attr(n, "nodePar") <- c(a$nodePar, list(lab.col = "green", lab.cex=.7, col="blue", pch=16 )) 
    } else { 
      attr(n, "nodePar") <- c(a$nodePar, list(lab.col = "red", lab.cex=.7,col="red", pch=16)) 
    }   
  } 
  return(n)
} 

#modfiy dendrogram nodes and re-plot 
dend_colored <- dendrapply(as.dendrogram(hclust.results), colLab) 
plot(dend_colored)
### hierarchiel clustering fails to separate the sensitive and resitant

# c) Differential expression. Identify genes that are differentially expressed between “Resistant” and “Sensitive” cell lines and compare your list with the signature of Potti et al. (2006) stored in doxorubicinGenes. In particular, state which test statistics you are using and, in the case of any probabilistic statement, state and justify the underlying assumptions. Provide and comment on numerical and graphical summaries of the results.

normSensESet <- normTrainESet[,NCI60ESet$status == "sensitive"]
normResESet <- normTrainESet[,NCI60ESet$status == "resistant"]

normSensDataTrans <- t(exprs(normSensESet))
normResDataTrans <- t(exprs(normResESet))

tTestResults <- sapply(X=1:12625,FUN=function(gene){
  t.test(x=normResDataTrans[,gene],y=normSensDataTrans[,gene])$statistic
})

ggplot(data=data.frame(X=tTestResults)) + stat_boxplot(aes(x="t statistic",y=X))

### compare my top 80 with theirs
mostDifferentialGenes <- colnames(normSensDataTrans)[which(tTestResults %in% sort(tTestResults,decreasing=T)[1:80])]

length(which(mostDifferentialGenes %in% doxorubicinGenes)) # 39 match

### check for top 80 absolutely differential expression
mostAbsDiffGenes <- colnames(normSensDataTrans)[which(abs(tTestResults) %in% sort(abs(tTestResults),decreasing=T)[1:80])]

length(which(mostAbsDiffGenes %in% doxorubicinGenes)) # 67 match

### top 80 negatively expressed
mostNegDiffGenes <- colnames(normSensDataTrans)[which(tTestResults %in% sort(tTestResults)[1:80])]

length(which(mostNegDiffGenes %in% doxorubicinGenes)) # 41 match

### personally there's no reason we shouldnt look
### at magnitude of t-stat, so that's my list.

### my own analysis
mostAbsExprGenesESet <- normTrainESet[mostAbsDiffGenes,]
transMostAbsExprGeneData <- t(exprs(mostAbsExprGenesESet))
my.hclust.results <- hclust(dist(transMostAbsExprGeneData))
dend_colored2 <- dendrapply(as.dendrogram(my.hclust.results), colLab) 
plot(dend_colored2) # these genes yield the right kind of clustering

my.pca.results <- prcomp(transMostAbsExprGeneData)
ggplot(data=as.data.frame(my.pca.results$x),aes(x=PC1, y=PC2,color=NCI60ESet$status)) + geom_text(aes(label=rownames(my.pca.results$x)),alpha=1,size=4) # much better separation

### try an SVM on the straight problem
library(e1071)
svm.model <- svm(x=transMostAbsExprGeneData, y=NCI60ESet$status)
table(predict(svm.model,transMostAbsExprGeneData),NCI60ESet$status) # wow!

### train svm to common genes
commonDiffGenes <- mostAbsDiffGenes[mostAbsDiffGenes %in% rownames(doxo07ESet)]
common.svm.model <- svm(x=transMostAbsExprGeneData[,commonDiffGenes], y=NCI60ESet$status)
table(predict(common.svm.model,transMostAbsExprGeneData[,commonDiffGenes]),NCI60ESet$status) # good

### test this out on patients
testESet <- doxo07ESet[commonDiffGenes,-(1:22)]
transTestData <- log(exprs(testESet)+0.1)
normTestESet <- normalize.ExpressionSet.quantiles(ExpressionSet(transTestData),transfn="none")
normTestData <- as.data.frame((exprs(normTestESet)))

### check box plot
melted.df3 <- melt(normTestData)
melted.df3 <- transform(melted.df3, value=as.numeric(value))
ggplot(data=melted.df3,aes(x=variable,y=value)) + stat_boxplot() + geom_hline(yintercept=0,color='red') + scale_color_discrete(guide=F) + labs(x="test patient", y="log(intensity)",title="quantile-based normalized expression values")

### try svm on patients
### refactor status first
levels(doxo07ESet$status) <- levels(NCI60ESet$status)
table(predict(common.svm.model,t(normTestData)),doxo07ESet$status[-(1:22)]) # just okay
### calculate accuracy
sum(predict(common.svm.model,t(normTestData)) == doxo07ESet$status[-(1:22)])/122 # 45% =\

### to see why this fails look at PC1 vs. PC2
common.pca.results <- prcomp(transMostAbsExprGeneData[,commonDiffGenes])
test.principle.coords <- as.data.frame(t(normTestData) %*% common.pca.results$rotation)

ggplot(data=as.data.frame(common.pca.results$x),aes(x=PC1, y=PC2,color=NCI60ESet$status)) + geom_text(aes(label=rownames(common.pca.results$x)),alpha=1,size=4) + geom_text(data=test.principle.coords, aes(x=PC1, y=PC2, label=rownames(test.principle.coords),color=doxo07ESet$status[-(1:22)]),alpha=0.5) + scale_color_discrete(name="")


