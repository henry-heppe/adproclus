
<!-- README.md is generated from README.Rmd. Please edit that file -->

# adproclus

<!-- badges: start -->
<!-- badges: end -->

\[package under development, will be submitted to CRAN end of 2023\]
This package is an implementation of the additive profile clustering
(ADPROCLUS) method in R. It can be used to obtain overlapping clustering
models for object-by-variable data matrices. It also contains the low
dimensional ADPROCLUS method, which achieves a simultaneous dimension
reduction when searching for overlapping clusters. This can be used when
the object-by-variable data contains a very large number of variables.

## Installation

You can install the development version of ADPROCLUS from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("henry-heppe/adproclus")
```

Or install the latest version from CRAN:

``` r
install.packages("adproclus")
```

## Example

This is a basic example which shows you how to use the regular ADPROCLUS
and the low dimensional ADPROCLUS:

``` r
library(adproclus)
#import data
our_data <- adproclus::CGdata

#perform ADPROCLUS to get an overlapping clustering model
model_full <- adproclus(data = our_data, nclusters = 2)
#> [1] "Clustering completed."
#> [1] "Additive Profile Clustering with 2 overlapping clusters."
#> [1] "Total processing time: 1.4 seconds."
#> [1] "Algorithm starts: 3 random and 3 semi random"

#perform low dimensional ADPROCLUS to get an overlapping clustering model in terms of a smaller number of variables
model_lowdim <- adproclusLD(data = our_data, nclusters = 3, ncomponents = 2)
#> [1] "Clustering completed."
#> [1] "Low dimensional Additive Profile Clustering with  3  overlapping clusters, and  2 components (dimensions)."
#> [1] "Total processing time:  0.5  seconds."
#> [1] "Algorithm starts:  3 random and  3  semi-random."
```

The package also provides functionality to obtain membership matrices,
that the algorithm can start the alternating least squares procedure on.
There are three different possibilities to obtain such matrices: random,
semi-random and rational (see respective function documentation for
details).

``` r
library(adproclus)
#import data
our_data <- adproclus::CGdata
#Obtaining a membership matrix were the entries are randomly assigned values of 0 or 1
start_allocation1 <- getRandom(our_data, 3)
#Obtaining a membership matrix based on a profile matrix consisting of randomly selected rows of the data
start_allocation2 <- getSemiRandom(our_data, 3)
# Obtaining a user-defined rational start profile matrix (here the first 3 rows of the data)
start_allocation3 <- getRational(our_data,our_data[1:3,])$A
```
