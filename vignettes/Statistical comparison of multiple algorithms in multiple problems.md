---
title: "scmampscma: Statistical comparison of multiple algorithms in multiple problems"
author: "Borja Calvo and Guzmán Santafé"
date: "2015-06-01"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Statistical comparison of multiple algorithms in multiple problems}
  %\VignetteEngine{knitr::docco_linear}
  \usepackage[utf8]{inputenc}
---
# Statistical comparison of multiple algorithms in multiple problems

This vignette shows the use of the package `scmamp` to assess the statistical differences between the results obtained by a number of algorithms in different problems. This is a typical task in areas such as Machine Learning or Optimization, where algorithms are typically compared measuring their performance in different instances of problems, datasets, etc. However, a similar procedure may be used in other contexts.

The package and this vignette is based mainly on the paper _García and Herrera (2008)_, which is an extenstion of Demšar's paper (_Demšar, 2006_).

If you are familiar with these papers and want a quick guide, jump to the last section of this document (_Summary_). Then, you can review the rest of the vignette for more details.

A a guiding example, we will use the results included in ther first paper (_García and Herrera ,2008_, Table 2). These data is available in the package under the name `garcia.herrera`. We can load the data by typing:







































