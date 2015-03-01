## ----, prompt=TRUE-------------------------------------------------------
library("scma")
data(data.garcia.herrera)
head(data.garcia.herrera)

## ----,prompt=TRUE , fig.width=10, fig.height=5---------------------------
plot.densities (results.matrix = data.garcia.herrera , size=1.1)

## ----,prompt=TRUE , fig.width=10, fig.height=5---------------------------
qqplot <- gaussian.qqplot (data.garcia.herrera[,"k-NN(k=1)"], size=5 , col="orchid")
qqplot + theme_classic()

## ----,prompt=TRUE--------------------------------------------------------
shapiro.test(data.garcia.herrera[,"k-NN(k=1)"])

## ----,prompt=TRUE--------------------------------------------------------
friedman.test(data.garcia.herrera)
iman.daveport.test(data.garcia.herrera)

## ----,prompt=TRUE--------------------------------------------------------
test <- nemenyi.test (data.garcia.herrera)
test
test$diff.matrix
abs(test$diff.matrix) > test$statistic

## ----,prompt=TRUE,fig.width=12 , fig.height=4----------------------------
critical.difference.plot (data.garcia.herrera , alpha = 0.1 , cex=1.25)
critical.difference.plot (data.garcia.herrera , alpha = 0.01 , cex=1.25)

## ----,prompt=TRUE , warning=FALSE----------------------------------------
pwcomp.holm <- pairwise.test(results.matrix = data.garcia.herrera , 
                             test = "Wilcoxon" , 
                             correction = "holm")
adhoc.test <- function(x,y,...) wilcox.test(x , y , paired=TRUE)$p.value
pwcomp.shaffer <- pairwise.test(results.matrix = data.garcia.herrera , 
                                test = adhoc.test ,  
                                correction = "Shaffer")

## ----,prompt=TRUE--------------------------------------------------------
pwcomp.holm$corrected.pvalues

## ----,prompt=TRUE , warning=FALSE----------------------------------------
pwcomp.bh <- pairwise.test(results.matrix = data.garcia.herrera , 
                             test = "Wilcoxon" , 
                             correction = "Bergmann Hommel")

## ----,prompt=TRUE , warning=FALSE----------------------------------------
pwcomp.bh$corrected.pvalues - pwcomp.shaffer$corrected.pvalues

## ----,prompt=TRUE--------------------------------------------------------
h <- pwcomp.shaffer$corrected.pvalues < 0.05
h
id <- which(h , arr.ind = TRUE)

## ----,prompt=TRUE--------------------------------------------------------
id <- which(h , arr.ind = TRUE)
id <- subset(id , id[,1]<id[,2])    ## Remove the repetitions
names <- colnames(h)
differences <- apply(id , MARGIN = 1 , 
                     FUN = function(x) paste(names[x[1]] , "-" ,  names[x[2]]))
names(differences) <- NULL
differences

## ----,prompt=TRUE,eval=FALSE---------------------------------------------
#  source("http://www.bioconductor.org/biocLite.R")
#  biocLite("Rgraphviz")

## ----,prompt=TRUE,fig.width=10 , fig.height=5----------------------------
rmeans <- colMeans (rank.matrix(data.garcia.herrera))
algorithm.graph (hypothesis.matrix = !h , mean.value = rmeans , 
                 font.size = 10 , node.width = 3 , node.height = 1)

## ----,prompt=TRUE,fig.width=10 , fig.height=5----------------------------
rmeans <- colMeans (rank.matrix(data.garcia.herrera))
rmeans <- colMeans (rank.matrix(data.garcia.herrera))
algorithm.graph (hypothesis.matrix = !h , mean.value = rmeans , 'fdp' , 
                 highlight.color = "red" , node.color = "white" , font.color = "black" ,
                 font.size = 10 , node.width = 2 , node.height = 1)

## ----,prompt=TRUE,fig.width=10 , fig.height=6, warning=FALSE-------------
plot.pvalues(pvalue.matrix = pwcomp.bh$corrected.pvalues , alg.order = order (rmeans,decreasing = FALSE))

