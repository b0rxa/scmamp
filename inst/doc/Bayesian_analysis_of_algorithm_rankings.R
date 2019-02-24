## ------------------------------------------------------------------------
library(scmamp)
data("data_blum_2015")
head(data.blum.2015)

## ---- message=FALSE------------------------------------------------------
data.100    <- subset(data.blum.2015, subset=data.blum.2015$Size==100)
results.100 <- bPlackettLuceModel(x.matrix=data.100[,-c(1,2)], min=FALSE, 
                                  nsim=2000, nchains=10,parallel=TRUE)

## ------------------------------------------------------------------------
results.100$expected.win.prob
results.100$expected.mode.rank

## ---- , fig.width=10, fig.height=6, out.width="100%"---------------------
hist(results.100$posterior.weights[,"FrogCOL"], 
     main="", xlab="Prob. FrogCOL being the best")

## ----message=FALSE-------------------------------------------------------
data.1000    <- subset(data.blum.2015, subset=data.blum.2015$Size==1000)
results.1000 <- bPlackettLuceModel(x.matrix=data.1000[,-c(1,2)], 
                                   min=FALSE, nsim=2000, nchains=10,parallel=TRUE)

## ---- fig.width=10, fig.height=6, out.width="100%"-----------------------
boxplot(results.1000$posterior.weights)

## ---- fig.width=8, fig.height=8, out.width="100%"------------------------
weights <- results.1000$posterior.weights[,c(1, 7, 8)]
weights <- weights / rowSums(weights)
plotBarycentric(weights)

## ---- fig.width=8, fig.height=8, out.width="100%"------------------------
fc.better.fm <- weights[, 2] / (weights[, 2] + weights[, 3])
fc.better.ff <- weights[, 2] / (weights[, 2] + weights[, 1])
exp.fc.vs.fm <- mean(fc.better.fm)
exp.fc.vs.ff <- mean(fc.better.ff)

hist(fc.better.fm, main=paste("Expected probability =", 
                               round(exp.fc.vs.fm,3)), 
     xlab="Probability of FrogCOL better than FrogMIS")

hist(fc.better.ff, main=paste("Expected probability =", 
                               round(exp.fc.vs.ff,3)), 
     xlab="Probability of FrogCOL better than FruitFly")

## ------------------------------------------------------------------------
mean(data.1000[,"FrogCOL"]>data.1000[,"FruitFly"])

## ------------------------------------------------------------------------
mean(data.1000[,"FrogCOL"]>data.1000[,"Ikeda"])
mean(data.1000[,"FrogMIS"]>data.1000[,"Ikeda"])
mean(data.1000[,"FruitFly"]>data.1000[,"Ikeda"])

## ----message=FALSE-------------------------------------------------------
data.1000.sub    <- subset(data.blum.2015[, c(3, 9,10)], subset=data.blum.2015$Size==1000)
results.1000.sub <- bPlackettLuceModel(x.matrix=data.1000.sub, 
                                   min=FALSE, nsim=2000, nchains=10,parallel=TRUE)

## ---- fig.width=8, fig.height=8, out.width="100%"------------------------
weights.sub <- results.1000.sub$posterior.weights

plotBarycentric(weights.sub)

fc.better.fm <- weights.sub[, 2] / (weights.sub[, 2] + weights.sub[, 3])
fc.better.ff <- weights.sub[, 2] / (weights.sub[, 2] + weights.sub[, 1])
exp.fc.vs.fm <- mean(fc.better.fm)
exp.fc.vs.ff <- mean(fc.better.ff)

hist(fc.better.fm, main=paste("Expected probability =", 
                               round(exp.fc.vs.fm,3)), 
     xlab="Probability of FrogCOL better than FrogMIS")

hist(fc.better.ff, main=paste("Expected probability =", 
                               round(exp.fc.vs.ff,3)), 
     xlab="Probability of FrogCOL better than FruitFly")

