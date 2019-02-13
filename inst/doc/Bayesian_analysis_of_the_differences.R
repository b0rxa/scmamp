## ---- prompt=TRUE--------------------------------------------------------
library("scmamp")
library("ggplot2")
data(data_kcv_example)

## ----prompt=TRUE , fig.width=7, fig.height=5, warning=FALSE--------------
algorithms <- names(data.kcv.example)[4:7]
db <- 5
plotDensities (data=data.kcv.example[data.kcv.example$DB==db, algorithms], size=1.1)

## ---- prompt=TRUE--------------------------------------------------------
db <- 5
sample.a <- data.kcv.example[data.kcv.example$DB==db, "AlgA"]
sample.b <- data.kcv.example[data.kcv.example$DB==db, "AlgB"]
results <- bCorrelatedTtest(x=sample.a, y=sample.b, rho=0.1, rope=c(-0.01, 0.01))
results$posterior.probabilities

## ----prompt=TRUE , fig.width=7, fig.height=5-----------------------------
plotPosterior(results, plot.rope=TRUE)

## ---- prompt=TRUE--------------------------------------------------------
db <- 5
sample.a <- data.kcv.example[data.kcv.example$DB==db, "AlgC"]
sample.b <- data.kcv.example[data.kcv.example$DB==db, "AlgD"]
results <- bCorrelatedTtest(x=sample.a, y=sample.b, rho=0.1, rope=c(-0.01, 0.01))
results$posterior.probabilities

## ----prompt=TRUE , fig.width=7, fig.height=5-----------------------------
plotPosterior(results, plot.rope=TRUE)

## ---- prompt=TRUE--------------------------------------------------------
db <- 9
summarized.data <- aggregate(data.kcv.example[, algorithms], 
                             by=data.kcv.example[, 1:2], FUN=mean)
sample.a <- summarized.data[summarized.data$DB==db, "AlgC"]
sample.b <- summarized.data[summarized.data$DB==db, "AlgD"]

## ---- prompt=TRUE, message=FALSE-----------------------------------------
results <- bSignedRankTest(x=sample.a, y=sample.b,rope=c(-0.01, 0.01))
results$posterior.probabilities

## ---- prompt=TRUE, message=FALSE-----------------------------------------
head(results$posterior)

## ---- prompt=TRUE, message=FALSE-----------------------------------------
colMeans(results$posterior)

## ----prompt=TRUE , fig.width=7, fig.height=7-----------------------------
plotSimplex(results, A="Algorithm C", B="Algorithm D")

## ----prompt=TRUE , fig.width=7, fig.height=7-----------------------------
db <- 8
summarized.data <- aggregate(data.kcv.example[, algorithms], 
                             by=data.kcv.example[, 1:2], FUN=mean)
sample.a <- summarized.data[summarized.data$DB==db, "AlgC"]
sample.b <- summarized.data[summarized.data$DB==db, "AlgD"]
results <- bSignedRankTest(x=sample.a, y=sample.b,rope=c(-0.01, 0.01))
results$posterior.probabilities
colMeans(results$posterior)
plotSimplex(results, plot.density=FALSE, A="Algorithm C", B="Algorithm D", 
            posterior.label=TRUE)

## ---- prompt=TRUE--------------------------------------------------------
summarized.data <- aggregate(data.kcv.example[, algorithms], 
                             by=data.frame(DB=data.kcv.example[, 1]), FUN=mean)
sample.a <- summarized.data[, "AlgC"]
sample.b <- summarized.data[, "AlgD"]

## ---- prompt=TRUE, message=FALSE-----------------------------------------
results <- bSignedRankTest(x=sample.a, y=sample.b,rope=c(-0.01, 0.01))
results$posterior.probabilities

## ----prompt=TRUE , fig.width=7, fig.height=7-----------------------------
plotSimplex(results, A="Algorithm C", B="Algorithm D", plot.density=FALSE, alpha=0.5)

## ---- prompt=TRUE--------------------------------------------------------
summarized.data <- aggregate(data.kcv.example[, algorithms], 
                             by=data.frame(DB=data.kcv.example[, 1]), FUN=mean)
sample.a <- summarized.data[, "AlgA"]
sample.b <- summarized.data[, "AlgB"]
results <- bSignedRankTest(x=sample.a, y=sample.b,rope=c(-0.01, 0.01))
results$posterior.probabilities

## ----prompt=TRUE , fig.width=7, fig.height=7-----------------------------
plotSimplex(results, A="Algorithm A", B="Algorithm B", plot.density=FALSE, alpha=0.5)

## ---- prompt=TRUE, message=FALSE-----------------------------------------
colMeans(results$posterior)

## ---- prompt=TRUE, message=FALSE-----------------------------------------
sample.a <- matrix(data.kcv.example$AlgC, byrow=TRUE, nrow=10)
sample.b <- matrix(data.kcv.example$AlgD, byrow=TRUE, nrow=10)

## ---- prompt=TRUE, message=FALSE, warning=FALSE--------------------------
results <- bHierarchicalTest(sample.a, sample.b, rho=0.1, rope=c(-0.01, 0.01), nsim=2000, nchains=5)

## ----prompt=TRUE , fig.width=7, fig.height=7-----------------------------
plotSimplex(results, A="Alg. C", B="Alg. D", posterior.label=TRUE, alpha=0.5)

## ---- prompt=TRUE, message=FALSE-----------------------------------------
results$additional$per.dataset

