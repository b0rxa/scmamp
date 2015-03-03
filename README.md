# scma
This is a simple R package aimed at simplifying the statistical analysis of the results in the comparison of algorithms in different problems.

It is mainly focused on non parametric methods and implements Shaffer static and Bergmann and Hommel dynamic corrections for pairwise tests.

The package also includes some plotting tools, such as the critical difference plots shown in _Demšar, J., 2006_. Indeed, the package is mainly based on the papers:

Demšar, J. (2006) Statistical Comparisons of Classifiers over Multiple Data Sets. _Journal of Machine Learning Research_, 7, 1-30.

García S. and Herrera, F. (2008) An Extension on "Statistical Comparisons of Classifiers over Multiple Data Sets" for all Pairwise Comparisons. _Journal of Machine Learning Research_, 9, 2677-2694.

## Package installation

The last version of the package can be installed running the following commands:

```r
if (!require("devtools"))
  install.packages("devtools")

devtools::install_github("b0rxa/scma")
```
The package can also be installed using the tar.gz files in the root directory. First, download the `scma_*.tar.gz` file and install it running:

```r
install.packages(path.tar.gz.file , reps=NULL)
```

where `path.tar.gz.file` refers to the path of the downloaded file. Note that these files may not be up to date.


## Documentation

The package includes a [vignette](http://htmlpreview.github.io/?https://github.com/b0rxa/scma/blob/master/inst/doc/Statistical%20comparison%20of%20multiple%20algorithms%20in%20multiple%20problems.html) explaining the basic use of the package. It can be accessed typing (once the package is installed):

```r
library("scma")
browseVignettes("scma")
``` 
