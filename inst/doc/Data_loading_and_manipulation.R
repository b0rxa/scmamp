## ---- prompt=TRUE--------------------------------------------------------
library(scmamp)
data(data_gh_2008)
head(data.gh.2008)
data(data_gh_2010)
head(data.gh.2010)
data(data_blum_2015)
head(data.blum.2015)

## ---- eval=FALSE, prompt=TRUE--------------------------------------------
#  data.raw <- readComparisonFile(file="results.dat", alg.cols=c('Alg_1', 'Alg_2', 'Alg_3'),
#                                 skip=5, sep=";")

## ---- prompt=TRUE--------------------------------------------------------
data.dir <- system.file("loading_tests",package="scmamp")
file.path <- paste(data.dir, "rgg_complete_comparison.out", sep="/")
data.raw <- readComparisonFile(file=file.path, alg.cols=3:10, col.names=NULL)
head(data.raw)

## ---- prompt=TRUE--------------------------------------------------------
dir <- paste(system.file("loading_tests",package="scmamp"), 
             "comparison_files", sep="/")
list.files(dir)

## ---- prompt=TRUE--------------------------------------------------------
fname.pattern <- "rgg_size_([0-9]*)_r_([0-9]*.[0-9]*)\\.out"

## ---- prompt=TRUE--------------------------------------------------------
var.names <- c("Size", "Radius")

## ---- prompt=TRUE--------------------------------------------------------
alg.names <- c("FruitFly", "Shukla", "Ikeda", "Turau", "Rand1", "Rand2", "FrogCOL", "FrogMIS")

## ---- prompt=TRUE--------------------------------------------------------
rm("data.raw")
data.raw <- readComparisonDir (directory=dir, alg.cols=alg.names, col.names=NULL,
                           names=var.names, fname.pattern=fname.pattern)
head(data.raw)

## ---- prompt=TRUE--------------------------------------------------------
dir <- system.file("loading_tests", package="scmamp")
file <- paste(dir, "rgg_complete_experiment.out", sep="/")
content <- read.csv(file)
content[c(1,901,181),]

## ---- prompt=TRUE, cache=TRUE--------------------------------------------
rm("data.raw")
data.raw <- readExperimentFile (file=file, alg.col="Algorithm", value.col="Evaluation")
head(data.raw)

## ---- prompt=TRUE, cache=TRUE---------------------------------------------------------------------------------------
rm("data.raw")
dir <- paste(system.file("loading_tests", package="scmamp"), 
             "experiment_files", sep="/")
list.files(dir)[1:10]
pattern <- "rgg_size_([0-9]*)_r_(0.[0-9]*)_([a-z, A-Z, 1, 2]*).out"
var.names <- c("Size", "Radius", "Algorithm")
data.raw <- readExperimentDir (directory=dir, names=var.names, fname.pattern=pattern,
                               alg.var.name='Algorithm', value.col=1, col.names="Evaluation")
head(data.raw)

## ---- echo=-1-------------------------------------------------------------------------------------------------------
summarizeData(data=data.raw, fun=median, group.by=c("Size"), ignore=c("Radius"))

## -------------------------------------------------------------------------------------------------------------------
data.filtered <- filterData(data=data.raw, 
                            condition="Size == 100 & Rand1 <= Rand2", 
                            remove.cols="Size")
dim(data.filtered)
dim(data.raw)

## -------------------------------------------------------------------------------------------------------------------
summarizeData(data.filtered, group.by=c("Radius"))

## ---- warning=FALSE , cache=TRUE , echo=-1--------------------------------------------------------------------------
test <- "wilcoxon"
group.by <- c("Size","Radius")
alg.cols <- 3:10
result <- postHocTest(data=data.raw, algorithms=alg.cols, group.by=group.by,
                      test=test, control="max", correct="holland")

## -------------------------------------------------------------------------------------------------------------------
summ <- result$summary
pval <- result$corrected.pval
bold <- is.na(pval)
mark <- pval > 0.05
mark[, (1:2)] <- FALSE
mark[is.na(mark)] <- FALSE
digits <- c(0, 3, rep(2, 8))

writeTabular(table=summ, format="f", bold=bold, mark=mark, mark.char="+", 
             hrule=c(0, 10, 20, 30), vrule = c(2, 4), digits=digits, 
             print.row.names=FALSE)

