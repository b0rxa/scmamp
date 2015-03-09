---
title: "Data loading and manipulation with scmamp"
author: "Borja Calvo and Guzmán Santafé"
date: "2015-03-09"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Data loading and manipulation with scmamp}
  %\VignetteEngine{knitr::docco_linear}
  \usepackage[utf8]{inputenc}
---
# Data loading and manipulation with scmamp

The main goal of this package is makeing the statistical analysis of emprical comparisions of algorithms easy and fast. For that reason, the package includes functions to load data and manipulate data, as well as to format the results for its further use in publications. This vignettes shows the use of these functions.

## Loading the data

The data matrices required by the package funtions should have one row per problem and a number of columns. The columns can be divided into two subsets, descriptors of the problem and results obtained by the algorithms applied to that problem. This is an example of the type of matrix needed:































