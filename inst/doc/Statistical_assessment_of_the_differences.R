## ---- prompt=TRUE--------------------------------------------------------
library("scmamp")
library("ggplot2")
library("Rgraphviz")
data(data_blum_2015)
data(data_gh_2008)
head(data.blum.2015)
head(data.gh.2008)

## ----prompt=TRUE , fig.width=10, fig.height=5----------------------------
plotDensities (data=data.gh.2008, size=1.1)

## ----prompt=TRUE , fig.width=10, fig.height=5----------------------------
qqplot <- qqplotGaussian (data.gh.2008[,"k-NN(k=1)"], size=5 , col="orchid")
qqplot + theme_classic()

## ----prompt=TRUE---------------------------------------------------------
friedmanTest(data.gh.2008)
imanDavenportTest(data.gh.2008)
friedmanAlignedRanksTest(data.gh.2008)
quadeTest(data.gh.2008)

## ----prompt=TRUE---------------------------------------------------------
test <- nemenyiTest (data.gh.2008, alpha=0.05)
test
test$diff.matrix
abs(test$diff.matrix) > test$statistic

## ----prompt=TRUE,fig.width=12 , fig.height=4-----------------------------
plotCD (data.gh.2008, alpha=0.05, cex=1.25)
plotCD (data.gh.2008, alpha=0.01, cex=1.25)

## ----prompt=TRUE---------------------------------------------------------
friedmanPost(data=data.gh.2008, control=NULL)
quadePost(data=data.gh.2008, control=NULL)
pv.matrix <- friedmanAlignedRanksPost(data=data.gh.2008, control=NULL)

## ----prompt=TRUE , warning=FALSE-----------------------------------------
pv.matrix
adjustShaffer(pv.matrix)
pv.adj <- adjustBergmannHommel(pv.matrix)
pv.adj

## ----prompt=TRUE,eval=FALSE----------------------------------------------
#  source("http://www.bioconductor.org/biocLite.R")
#  biocLite("Rgraphviz")

## ----prompt=TRUE,fig.width=10 , fig.height=5-----------------------------
r.means <- colMeans(rankMatrix(data.gh.2008))
drawAlgorithmGraph(pvalue.matrix=pv.adj, mean.value=r.means, alpha=0.05,
                 font.size=10, node.width=3, node.height=1)

## ----prompt=TRUE,fig.width=10 , fig.height=5-----------------------------
r.means <- colMeans (rankMatrix(data.gh.2008))
drawAlgorithmGraph (pvalue.matrix=pv.adj, mean.value=r.means, alpha=0.05, 'fdp',
                    highlight.color="red", node.color="white", font.color="black",
                    font.size=10, node.width=2, node.height=1)

## ----prompt=TRUE,fig.width=10 , fig.height=6, warning=FALSE--------------
plt <- plotPvalues(pvalue.matrix=pv.adj, 
                   alg.order=order(r.means, decreasing=FALSE))
plt + 
  labs(title="Corrected p-values using Bergmann and Hommel procedure") + 
  scale_fill_gradientn("Corrected p-values" , colours = c("skyblue4" , "orange"))

## ----prompt=TRUE---------------------------------------------------------
friedmanAlignedRanksPost(data.gh.2008, control = "NaiveBayes")
pv <- quadePost(data.gh.2008, control = 2)

## ----prompt=TRUE---------------------------------------------------------
adjustHolland(pvalues=pv)
adjustFinner(pvalues=pv)
adjustRom(pvalues=pv, alpha=0.05)
adjustLi(pvalues=pv)

## ----full_process_1, prompt=TRUE-----------------------------------------
alpha <- 0.05
data <- data.gh.2008

friedmanTest(data)

## ----full_process_2, prompt=TRUE-----------------------------------------
multipleComparisonTest(data=data, test="iman")

## ----full_process_3, prompt=TRUE , fig.width=10 , fig.height=5-----------
post.results <- postHocTest(data=data, test="aligned ranks", correct="bergmann", 
                            use.rank=TRUE)
post.results

alg.order <- order(post.results$summary)
plt <- plotPvalues(post.results$corrected.pval, alg.order=alg.order) 
plt + labs(title=paste("Corrected p-values using Bergmann and Hommel procedure",sep=""))
drawAlgorithmGraph(post.results$corrected.pval, mean.value=post.results$summary, 
                   alpha=alpha,  font.size=10)

## ----full_process_4, prompt=TRUE-----------------------------------------
data <- data.blum.2015
group.by <- c("Size","Radius")
multipleComparisonTest(data=data, group.by=group.by, 
                       test="quade", correct="finner")

control <- "FrogCOL"
post.results <- postHocTest(data=data, group.by=group.by, control=control, 
                            test="aligned ranks", correct="rom", use.rank=FALSE)

## ----full_process_5, prompt=TRUE-----------------------------------------
avg.val <- post.results$summary
best <- apply(avg.val, MARGIN=1, 
              FUN=function(x){
                m <- max(x[-(1:2)])
                return(c(FALSE, FALSE, x[-(1:2)]==m))
              })
best <- t(best)
no.diff <- post.results$corrected.pval > alpha
# The size and radius columns set as false
no.diff[,1:2] <- FALSE
no.diff[is.na(no.diff)] <- FALSE
writeTabular(table=avg.val, format='f', bold=best, italic=no.diff, 
             hrule=c(0, 10, 20, 30), vrule=2, digits=c(0, 3, rep(2, 8)), 
             print.row.names = FALSE)

## ----full_process_6, prompt=TRUE, fig.width=10 , fig.height=5------------
control <- NULL
group.by <- "Size"
post.results <- postHocTest(data=data, algorithms=3:10, group.by=group.by, 
                            control=control, test="aligned ranks", correct="holland", 
                            use.rank=TRUE)

# Plot the matrix for the first group
i <- 1
alg.order <- order(post.results$summary[i,-1])
plotPvalues(post.results$corrected.pval[, , i], alg.order=alg.order) 

# Plot the matrix for the second group
i <- 2
alg.order <- order(post.results$summary[i,-1])
plotPvalues(post.results$corrected.pval[, , i], alg.order=alg.order) 

# Plot the matrix for the third group
i <- 3
alg.order <- order(post.results$summary[i,-1])
plotPvalues(post.results$corrected.pval[, , i], alg.order=alg.order) 

