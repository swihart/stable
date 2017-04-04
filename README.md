
<!-- README.md is generated from README.Rmd. Please edit README.Rmd -->
[![Travis-CI Build Status](https://travis-ci.org/swihart/stable.svg?branch=master)](https://travis-ci.org/swihart/stable) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/stable)](https://cran.r-project.org/package=stable) ![downloads](http://cranlogs.r-pkg.org/badges/grand-total/stable)

`stable` R package
==================

This package is intended to be the developmental version to the CRAN version of [Jim Lindsey's stable](http://www.commanster.eu/rcode.html). The .zip files listed on his homepage have been listed as version 1.0 since 2005. For the subsequent maintenance on this github and CRAN, we will start at 1.1.

To compare this version with the static v1.0 files on [Jim Lindsey's Homepage](http://www.commanster.eu/rcode.html), it may be useful to use [the compare page for this repo's two branches](https://github.com/swihart/stable/compare/jim-lindsey-homepage-version-1.0...master?diff=split&name=master).

comparisons with `stabledist` R package
=======================================

In brief, the parameters have different names and are transformations for each other. First, the names:

| stabledist | stable |
|------------|--------|
| alpha      | tail   |
| beta       | skew   |
| gamma      | disp   |
| delta      | loc    |

If you read the Lindsey PDF in this repo, be aware that location is given the greek letter gamma and scale is given the greek letter delta. The Nolan PDF does the opposite and is used for `stabledist`.

For some values for some distributions things match up nicely, as we see with Normal and Cauchy:

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

However, to make `stable` equivalent to `stabledist` in general, some transformations are needed. Please see the following examples. Between `stabledist` and `stable`, the `alpha` is equivalent to `tail` and the `delta` is equivalent to `loc` with no transformation. For the `beta` (`skew`) and `gamma` (`disp`) parameters, a transformation is needed to get equivalent calls. Note differences still may exist to numerical accuracy.

levy cdf
--------

``` r
q <-  0.9

# nolan pm=1 parameters:
a <-  0.5
b <-  1
c <-  .25
d <-  0.8

# lindsey-(3) page 415 conversion:
# tail/alpha and location stay the same
a3 <- a
d3 <- d 
# the others require calcs:
DEL2 <- cos(pi/2 * a)^2 + (-b)^2*sin(pi/2 * a)^2
DEL <- sqrt(DEL2) * sign(1-a)
eta_a <- min(a, 2-a)
# the lindsey-(3) beta:
b3 <- 2/(pi*eta_a)*acos( cos(pi/2 * a) / DEL )
# the lindsey-(3) scale:
c3 <- ( (DEL*c^a) / cos(pi/2 * a) )^(1/a)

    stable::pstable(q, tail =a, skew=b3, disp =c3, loc  =d)
#> [1] 0.1154242
stabledist::pstable(q, alpha=a, beta=b , gamma=c , delta=d, pm=1)
#> [1] 0.1138462
rmutil::plevy(q, m=d, s=c)
#> [1] 0.1138463

# more accuracy!!!!?!
    stable::pstable(q, tail =a, skew=b3, disp =c3, loc  =d, eps = 0.13*1e-7)
#> [1] 0.1138786
```

levy pdf
--------

``` r
q <-  0.9

# nolan pm=1 parameters:
a <-  0.5
b <-  1
c <-  .25
d <-  0.8

# lindsey-(3) page 415 conversion:
# tail/alpha and location stay the same
a3 <- a
d3 <- d 
# the others require calcs:
DEL2 <- cos(pi/2 * a)^2 + (-b)^2*sin(pi/2 * a)^2
DEL <- sqrt(DEL2) * sign(1-a)
eta_a <- min(a, 2-a)
# the lindsey-(3) beta:
b3 <- 2/(pi*eta_a)*acos( cos(pi/2 * a) / DEL )
# the lindsey-(3) scale:
c3 <- ( (DEL*c^a) / cos(pi/2 * a) )^(1/a)

    stable::dstable(q, tail =a, skew=b3, disp =c3, loc  =d)
#> [1] 1.806389
stabledist::dstable(q, alpha=a, beta=b , gamma=c , delta=d, pm=1)
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

levy quantile
-------------

``` r
p <-  .3

# nolan pm=1 parameters:
a <-  0.5
b <-  1
c <-  .25
d <-  0.8

# lindsey-(3) page 415 conversion:
# tail/alpha and location stay the same
a3 <- a
d3 <- d 
# the others require calcs:
DEL2 <- cos(pi/2 * a)^2 + (-b)^2*sin(pi/2 * a)^2
DEL <- sqrt(DEL2) * sign(1-a)
eta_a <- min(a, 2-a)
# the lindsey-(3) beta:
b3 <- 2/(pi*eta_a)*acos( cos(pi/2 * a) / DEL )
# the lindsey-(3) scale:
c3 <- ( (DEL*c^a) / cos(pi/2 * a) )^(1/a)

    stable::qstable(p, tail =a, skew=b3, disp =c3, loc  =d)
#> [1] 1.031301
stabledist::qstable(p, alpha=a, beta=b , gamma=c , delta=d, pm=1)
#> [1] 1.032735
rmutil::qlevy(p, m=d, s=c)
#> [1] 1.032733
```

play with alpha not 2 and not 1
-------------------------------

``` r
q <- -1.97

# nolan pm=1 parameters:
a <-  0.8
b <-  0
c <-  1
d <-  0

# lindsey-(3) page 415 conversion:
# tail/alpha and location stay the same
a3 <- a
d3 <- d 
# the others require calcs:
DEL2 <- cos(pi/2 * a)^2 + (-b)^2*sin(pi/2 * a)^2
DEL <- sqrt(DEL2) * sign(1-a)
eta_a <- min(a, 2-a)
# the lindsey-(3) beta:
b3 <- 2/(pi*eta_a)*acos( cos(pi/2 * a) / DEL )
# the lindsey-(3) scale:
c3 <- ( (DEL*c^a) / cos(pi/2 * a) )^(1/a)

    stable::pstable(q, tail =a, skew=b3, disp =c3, loc  =d)
#> [1] 0.1722953
stabledist::pstable(q, alpha=a, beta=b , gamma=c , delta=d)
#> [1] 0.1722945
```

play with skew
--------------

``` r
q <- -1

# nolan pm=1 parameters:
a <-  1.3
b <-  0.4
c <-  2
d <-  0.75

# lindsey-(3) page 415 conversion:
# tail/alpha and location stay the same
a3 <- a
d3 <- d 
# the others require calcs:
DEL2 <- cos(pi/2 * a)^2 + (-b)^2*sin(pi/2 * a)^2
DEL <- sqrt(DEL2) * sign(1-a)
eta_a <- min(a, 2-a)
# the lindsey-(3) beta:
b3 <- -sign(b)*2/(pi*eta_a)*acos( cos(pi/2 * a) / DEL )
# the lindsey-(3) scale:
c3 <- ( (DEL*c^a) / cos(pi/2 * a) )^(1/a)

    stable::pstable(q, tail =a, skew=b3, disp =c3, loc  =d)
#> [1] 0.4349168
stabledist::pstable(q, alpha=a, beta=b , gamma=c , delta=d, pm=1)
#> [1] 0.4348957

    stable::dstable(q, tail =a, skew=b3, disp =c3, loc  =d)
#> [1] 0.1454112
stabledist::dstable(q, alpha=a, beta=b , gamma=c , delta=d, pm=1)
#> [1] 0.1454111
```

The example above, but using `sd2s` and `s2sd`
----------------------------------------------

``` r
q <- -1
# nolan pm=1 parameters:
a <-  1.3
b <-  -0.4
c <-  2
d <-  0.75
# sd2s takes nolan (stabledist) parameters and returns lindsey (stable)
s <- stable::sd2s(alpha=a, beta=b, gamma=c, delta=d)
stable::pstable(q, tail = s$tail, skew=s$skew, disp = s$disp, loc  = s$loc)
#> [1] 0.196531
stabledist::pstable(q, alpha=a, beta=b , gamma=c , delta=d, pm=1)
#> [1] 0.1965513
# s2sd takes lindsey (stable) parameters and returns nolan (stabledist)
sd <- stable::s2sd(tail = s$tail, skew=s$skew, disp = s$disp, loc  = s$loc)
stabledist::pstable(q, alpha=sd$alpha, beta=sd$beta , gamma=sd$gamma , delta=sd$delta, pm=1)
#> [1] 0.1965513
```

pm1\_to\_pm0
------------

``` r
q <- -1
# nolan pm=1 parameters:
a1 <-  1.3
b1 <-  -0.4
c1 <-  2
d1 <-  0.75
# for a1 != 1
d0 <- d1 + b1*c1*tan(pi*a1/2)

# Calculate d0 by hand or use pm1_to_pm0():
# Convert to nolan pm=0 parameters:
pm0 <- stable::pm1_to_pm0(a1,b1,c1,d1)
a0 <- pm0$a0
b0 <- pm0$b0
c0 <- pm0$c0
d0 <- pm0$d0
# check:
stabledist::pstable(q, alpha=a1, beta=b1 , gamma=c1 , delta=d1, pm=1)
#> [1] 0.1965513
# only change delta=d0 for pm=0
stabledist::pstable(q, alpha=a1, beta=b1 , gamma=c1 , delta=d0, pm=0)
#> [1] 0.1965513
stabledist::pstable(q, alpha=a0, beta=b0 , gamma=c0 , delta=d0, pm=0)
#> [1] 0.1965513


stabledist::dstable(q, alpha=a1, beta=b1 , gamma=c1 , delta=d1, pm=1)
#> [1] 0.0572133
# only change delta=d0 for pm=0
stabledist::dstable(q, alpha=a1, beta=b1 , gamma=c1 , delta=d0, pm=0)
#> [1] 0.0572133
stabledist::dstable(q, alpha=a0, beta=b0 , gamma=c0 , delta=d0, pm=0)
#> [1] 0.0572133
```

mode of a stable distribution
-----------------------------

``` r
q <- -1
# nolan pm=1 parameters:
# a1 <-  1.3
# b1 <-  0.4
# c1 <-  2
# d1 <-  0.75
a1 <-  1.3
b1 <-  .5
c1 <-  1
d1 <-  0
# for a1 != 1
d0 <- d1 + b1*c1*tan(pi*a1/2)


s <- stable::sd2s(alpha=a1, beta=b1, gamma=c1, delta=d1)
stable::stable.mode(tail = s$tail, skew=s$skew, disp = s$disp, loc  = s$loc)$ytilde
#> [1] -1.13224

c1*stabledist::stableMode(alpha=a1, beta=b1)+d0
#> [1] -1.133257

xran <- seq(-2.5,2.6,0.001)
ysd <- stabledist::dstable(xran, alpha=a1, beta=b1, gamma=c1, delta=d1, pm=1)
plot(xran, ysd)

xran[ysd == max(ysd)]
#> [1] -1.133

ys <- stable::dstable(xran, tail = s$tail, skew=s$skew, disp = s$disp, loc  = s$loc)
points(xran, ys, col="blue")
```

![](README-unnamed-chunk-11-1.png)

``` r

xran[ys == max(ys)]
#> [1] -1.133
```
