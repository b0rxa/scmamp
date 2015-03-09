## ----, echo=F------------------------------------------------------------
library(scmamp)
data <- read.csv(paste(system.file("loading_tests",package="scmamp") , "beta_complete_comparison.out" , sep="/"))
data[c(271,301,331),]

## ----, eval=FALSE--------------------------------------------------------
#  data <- read.comparison.file(file = "results.dat" , alg.cols = c('Alg_1', 'Alg_2', 'Alg_3') , skip = 5 , sep = ";")

## ------------------------------------------------------------------------
dir <- paste(system.file("loading_tests",package="scmamp") , "comparison_files" , sep="/")
list.files(dir)

## ------------------------------------------------------------------------
pattern <- "beta_([0-9]*),([0-9]*)_size_([0-9]*)\\.out"

## ------------------------------------------------------------------------
var.names <- c('alpha' , 'beta' , 'size')

## ------------------------------------------------------------------------
alg.names <- c('kakizawa','vitale','boundarykernel','betakernel')

## ------------------------------------------------------------------------
data <- read.comparison.dir (directory = dir , alg.cols = alg.names , col.names = NULL , 
                             names = var.names , fname.pattern = pattern)
dim(data)
head(data)

## ------------------------------------------------------------------------
dir <- system.file("loading_tests",package="scmamp")
file <- paste(dir , "beta_complete_experiment.out" , sep="/")
content <- read.csv(file)
content[c(1,901,181),]

## ----, cache=TRUE--------------------------------------------------------
data <- read.experiment.file (file = file , alg.col = 'algorithm' , value.col = 'error')
head(data)

## ----, cache=TRUE--------------------------------------------------------
dir <- paste(system.file("loading_tests",package="scmamp") , "experiment_files" , sep="/")
list.files(dir)[1]
pattern <- "beta_([0-9]*),([0-9]*)_size_([0-9]*)_([a-z]*).out"
var.names <- c('alpha' , 'beta' , 'size', 'estimator')
data <- read.experiment.dir (directory = dir , names = var.names , 
                             alg.var.name = 'estimator' , value.col = 1 , 
                             fname.pattern = pattern)
head(data)

## ----, echo=-1-----------------------------------------------------------
data <- read.csv(paste(system.file("loading_tests",package="scmamp") , "beta_complete_comparison.out" , sep="/"))
summarize.data(data = data , fun = mean , group.by = c('size') , ignore = c('alpha' , 'beta'))

## ------------------------------------------------------------------------
filtered.data <- filter.data(data = data , condition = "size == 100 & alpha != beta" , remove.cols = 'size')
dim(filtered.data)
dim(data)

## ------------------------------------------------------------------------
summarize.data(filtered.data , group.by = c('alpha' , 'beta'))

## ----, warning=FALSE , cache=TRUE , echo=-1------------------------------
data <- read.csv(paste(system.file("loading_tests",package="scmamp") , "beta_complete_comparison.out" , sep="/"))
test <- wilcox.test
g <- c('size','alpha','beta')
a <- 4:7
result <- all.vs.best.test(data , test = test , group.by = g , 
                        alg.col = a , best = 'min')

## ------------------------------------------------------------------------
summ <- result$summary
pval <- result$adj.pvalues
bold <- is.na(pval)
mark <- pval > 0.05
mark[,(1:3)] <- FALSE
mark[is.na(mark)] <- FALSE
digits <- c(rep(0,3) , rep(3,4))

write.tabular(table = summ , format = 'f' , bold = bold , mark = mark , mark.char = '+' , hrule = c(0,4,8,12) , vrule = c(1,2,3) , digits=digits , print.row.names = FALSE)

