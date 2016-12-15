
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

levy cdf (not working as expected...)
-------------------------------------

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

levy pdf (not working as expected...)
-------------------------------------

``` r
q <-  0.9
a <-  0.5
b <-  1
c <-  .25
d <-  0.8
    stable::dstable(q, tail =a, skew=b, disp =c, loc  =d)
#> [1] 2.387435
stabledist::dstable(q, alpha=a, beta=b, gamma=c, delta=d, pm=1)
#> Warning in uniroot(function(th) log(g(th)), lower = l.th, upper = u.th, : -
#> Inf replaced by maximally negative value

#> Warning in uniroot(function(th) log(g(th)), lower = l.th, upper = u.th, : -
#> Inf replaced by maximally negative value
#> Warning in .integrate2(g1, lower = a, upper = b, subdivisions =
#> subdivisions, : roundoff error is detected in the extrapolation table
#> [1] 1.807224
rmutil::dlevy(q, m=d, s=c)
#> [1] 1.807224
```

levy quantile (not working as expected...)
------------------------------------------

``` r
p <-  .3
a <-  0.5
b <-  1
c <-  .25
d <-  0.8
    stable::qstable(p, tail =a, skew=b, disp =c, loc  =d)
#> [1] 0.9156615
stabledist::qstable(p, alpha=a, beta=b, gamma=c, delta=d, pm=1)
#> [1] 1.032735
rmutil::qlevy(p, m=d, s=c)
#> [1] 1.032733
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

Notice when a=1 or a=2, it just uses \[d/p/q/r\]cauchy or \[d/p/q/r\]norm respectively, so they match up. And when b=0 and a != 1 or a!=2, the results seem to match up. Something going on with skew parameterization; consult two pdfs in repo above.

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
