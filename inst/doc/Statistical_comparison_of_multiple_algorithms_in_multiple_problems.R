## ----prompt=TRUE---------------------------------------------------------
library(scmamp)
head(data.gh.2008)

## ----first_analysis, prompt=TRUE-----------------------------------------
imanDavenportTest(data.gh.2008)

## ----first_analysis_2, prompt=TRUE,fig.keep='all', fig.show='hide', fig.path='./images/', fig.width=8, fig.height=4----
plotCD(results.matrix=data.gh.2008, alpha=0.05)

## ----first_analysis_3, prompt=TRUE,fig.keep='all', fig.show='hide', fig.path='./images/', fig.width=4, fig.height=2----
res <- postHocTest(data=data.gh.2008, test="friedman", correct="bergmann", 
                   use.rank=TRUE)

# corrected p-value matrix
c.pval <- res$corrected.pval
c.pval

# LaTeX formated: No signficances higlighted in bold
bold <- c.pval > 0.05
bold[is.na(bold)] <- FALSE
writeTabular(table=c.pval, format="f", bold=bold, hrule=0, vrule=0)

# Graph including the average ranking and links between algorithms whose 
# are not significant differences
drawAlgorithmGraph(pvalue.matrix=res$corrected.pval, mean.value=res$summary)

## ---- prompt=TRUE--------------------------------------------------------
dir <- paste(system.file('loading_tests', package='scmamp'), 
             "experiment_files", sep="/")
list.files(dir)[1:5]
pattern <- 'rgg_size_([0-9]*)_r_(0.[0-9]*)_([a-z,A-Z,1,2]*).out'
var.names <- c('Size', 'Radius', 'Algorithm')
dataset <- readExperimentDir (directory=dir, names=var.names, fname.pattern=pattern, 
                              alg.var.name='Algorithm', value.col='Evaluation',
                              col.names='Evaluation')
head(dataset)

## ---- prompt=TRUE--------------------------------------------------------
sub.dataset <- filterData(data=dataset, condition='Size==1000', 
                          remove.cols='Size')

## ---- prompt=TRUE--------------------------------------------------------
res <- postHocTest(data=sub.dataset, group.by='Radius', test='wilcoxon',
                   correct='finner', control='FrogCOL')

## ---- prompt=TRUE--------------------------------------------------------
tab <- res$summary
best.res <- t(apply(tab, MARGIN=1, 
                    FUN=function(x) {
                      aux <- rep(FALSE, length(x))
                      aux[x == max(x)] <- TRUE
                      return(aux)
                    }))

no.diff <- res$corrected.pval > 0.05
no.diff[is.na(no.diff)] <- FALSE
writeTabular(table=tab, format='f', bold=best.res, mark=no.diff, 
             hrule=0, vrule=1, print.row.names=FALSE, digits=1)

