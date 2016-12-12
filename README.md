
<!-- README.md is generated from README.Rmd. Please edit README.Rmd -->
[![Travis-CI Build Status](https://travis-ci.org/swihart/stable.svg?branch=master)](https://travis-ci.org/swihart/stable) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/stable)](https://cran.r-project.org/package=stable) ![downloads](http://cranlogs.r-pkg.org/badges/grand-total/stable)

`stable` R package
==================

This package is intended to be the developmental version to the CRAN version of [Jim Lindsey's stable](http://www.commanster.eu/rcode.html). The .zip files listed on his homepage have been listed as version 1.0 since 2005. For the subsequent maintenance on this github and CRAN, we will start at 1.1.

To compare this version with the static v1.0 files on [Jim Lindsey's Homepage](http://www.commanster.eu/rcode.html), it may be useful to use [the compare page for this repo's two branches](https://github.com/swihart/stable/compare/jim-lindsey-homepage-version-1.0...master?diff=split&name=master).

comparisons with `stabledist` R package
=======================================

normal distribution
-------------------

``` r
q <- 3
    stable::pstable(q, tail =2, skew=0, disp =1, loc  =0)
#> [1] 0.9830526
stabledist::pstable(q, alpha=2, beta=0, gamma=1, delta=0)
#> [1] 0.9830526
```

cauchy distribution
-------------------

``` r
q <- 3
    stable::pstable(q, tail =1, skew=0, disp =1, loc  =0)
#> [1] 0.8975836
stabledist::pstable(q, alpha=1, beta=0, gamma=1, delta=0)
#> [1] 0.8975836
```

levy
----

``` r
q <-  0.9
a <-  0.5
b <-  1
c <-  .25
d <-  0.8
    stable::pstable(q, tail =a, skew=b, disp =c, loc  =d)
#> [1] 0.2650645
stabledist::pstable(q, alpha=a, beta=b, gamma=c, delta=d, pm=1)
#> [1] 0.1138462
rmutil::plevy(q, m=d, s=c)
#> [1] 0.1138463
```

play with alpha not 2 and not 1
-------------------------------

``` r
q <- -1.97
a <-  0.8
b <-  0
c <-  1
d <-  0
    stable::pstable(q, tail =a, skew=b, disp =c, loc  =d)
#> [1] 0.1722953
stabledist::pstable(q, alpha=a, beta=b, gamma=c, delta=d)
#> [1] 0.1722945
```

play with skew
--------------

``` r
q <- -1.97
a <-  1.3
b <-  -.4
c <-  2
d <-  0
    stable::pstable(q, tail =a, skew=-b, disp =c^(a), loc  =d, eps=1.0e-12)
#> [1] 0.2193218
stabledist::pstable(q, alpha=a, beta=b, gamma=c, delta=d, pm=1)
#> [1] 0.1844804
```
