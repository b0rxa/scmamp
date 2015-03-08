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
iman.davenport.test(data.garcia.herrera)

## ----,prompt=TRUE--------------------------------------------------------
test <- nemenyi.test (data.garcia.herrera , alpha = 0.05)
test
test$diff.matrix
abs(test$diff.matrix) > test$statistic

## ----,prompt=TRUE,fig.width=12 , fig.height=4----------------------------
critical.difference.plot (data.garcia.herrera , alpha = 0.05 , cex=1.25)
critical.difference.plot (data.garcia.herrera , alpha = 0.01 , cex=1.25)

## ----,prompt=TRUE , warning=FALSE----------------------------------------
pwcomp.shaffer <- pairwise.test(results.matrix = data.garcia.herrera , 
                             test = "Friedman post" , 
                             correction = "Shaffer")
pwcomp.holm <- pairwise.test(results.matrix = data.garcia.herrera , 
                                test = wilcox.test ,  
                                correction = "holm")

## ----,prompt=TRUE--------------------------------------------------------
pwcomp.holm$corrected.pvalues

## ----,prompt=TRUE , warning=FALSE----------------------------------------
pwcomp.bh <- pairwise.test(results.matrix = data.garcia.herrera , 
                             test = "Friedman post" , 
                             correction = "Bergmann Hommel")

## ----,prompt=TRUE , warning=FALSE----------------------------------------
pwcomp.bh$corrected.pvalues - pwcomp.shaffer$corrected.pvalues

## ----,prompt=TRUE--------------------------------------------------------
rej.h.shaffer <- pwcomp.shaffer$corrected.pvalues < 0.05
rej.h.shaffer

## ----,prompt=TRUE--------------------------------------------------------
id <- which(rej.h.shaffer , arr.ind = TRUE)
id <- subset(id , id[,1]<id[,2])    ## Remove the repetitions
names <- colnames(rej.h.shaffer)
differences <- apply(id , MARGIN = 1 , 
                     FUN = function(x) paste(names[x[1]] , "-" ,  names[x[2]]))
names(differences) <- NULL
differences

## ----,prompt=TRUE,eval=FALSE---------------------------------------------
#  source("http://www.bioconductor.org/biocLite.R")
#  biocLite("Rgraphviz")

## ----,prompt=TRUE,fig.width=10 , fig.height=5----------------------------
rmeans <- colMeans (rank.matrix(data.garcia.herrera))
algorithm.graph (pvalue.matrix = pwcomp.shaffer$corrected.pvalues , alpha = 0.05 , mean.value = rmeans , 
                 font.size = 10 , node.width = 3 , node.height = 1)

## ----,prompt=TRUE,fig.width=10 , fig.height=5----------------------------
rmeans <- colMeans (rank.matrix(data.garcia.herrera))
rmeans <- colMeans (rank.matrix(data.garcia.herrera))
algorithm.graph (pvalue.matrix = pwcomp.bh$corrected.pvalues , alpha = 0.05 , mean.value = rmeans , 'fdp' , 
                 highlight.color = "red" , node.color = "white" , font.color = "black" ,
                 font.size = 10 , node.width = 2 , node.height = 1)

## ----,prompt=TRUE,fig.width=10 , fig.height=6, warning=FALSE-------------
plt <- plot.pvalues(pvalue.matrix = pwcomp.bh$corrected.pvalues , alg.order = order (rmeans,decreasing = FALSE))
plt + 
  labs(title="Corrected p-values using Bergmann and Hommel procedure") + 
  scale_fill_gradientn("Corrected p-values" , colours = c("skyblue4" , "orange"))

## ----,full_process_1, prompt=TRUE----------------------------------------
alpha <- 0.05
data <- data.garcia.herrera

friedman.test(data)

## ----,full_process_2, prompt=TRUE , fig.width=10 , fig.height=5----------
correction <- ifelse(dim(data)[2]<=length(E) , "Bergmann Hommel" , "Shaffer") 
pwcomp <- pairwise.test(data , correction = correction)

mean.rank <- colMeans(rank.matrix(data))
alg.order <- order(mean.rank)
plot.pvalues(pwcomp$corrected.pvalues , alg.order = alg.order) + labs(title=paste("Corrected p-values using ", correction , " procedure",sep=""))
algorithm.graph(pwcomp$corrected.pvalues, mean.rank, alpha = alpha,  font.size = 10)

