#' Stable Distribution
#'
#'These functions provide information about the stable distribution
#'with the location, the dispersion, the skewness and the tail thickness
#'respectively modelled by the parameters \code{loc}, \code{disp},
#'\code{skew} and  \code{tail}.
#'
#'\code{dstable}, \code{pstable}, \code{qstable} and \code{hstable}
#'compute the density, the distribution, the quantile and the hazard functions
#'of a stable variate. \code{rstable} generates random deviates with
#'the prescribed stable distribution.
#'
#'\code{loc} is a location parameter in the same way as the mean
#'in the normal distribution: it can take any real value.
#'
#'\code{disp} is a dispersion parameter in the same way as the standard
#'deviation in the normal distribution: it can take any positive value.
#'
#'\code{skew} is a skewness parameter: it can take any value in \eqn{(-1,1)}.
#'The distribution is right-skewed, symmetric and left-skewed
#'when \code{skew} is negative, null or positive respectively.
#'
#'\code{tail} is a tail parameter (often named the characteristic exponent):
#'  it can take any value in \eqn{(0,2)} (with \code{tail=1} and \code{tail=2}
#'                                        yielding the Cauchy and the normal distributions respectively
#'                                        when symmetry holds).
#'
#'If \code{loc}, \code{disp}, \code{skew}, or \code{tail} are not
#'specified they assume the default values of \eqn{0}, \eqn{1/sqrt(2)},
#'\eqn{0} and \eqn{2} respectively. This corresponds to a normal
#'variate with mean\eqn{=0} and variance\eqn{=1/2 disp^2}.
#'
#'The stable characteristic function is given by
#'\deqn{greekphi(t) = i loca t - disp {|t|}^{tail}  [1+i skew sign(t) greekomega(t,tail)]}{phi(t) = i loc t - disp |t|^tail [1+i skew sign(t) omega(t,tail)]}
#'where
#'\deqn{greekomega(t,tail) = \frac{2}{\pi} LOG(ABS(t))}{omega(t,tail) = (2/pi) log|t|}
#'when \code{tail=1}, and
#'\deqn{greekomega(t,tail) = tan(\frac{\pi tail}{2})}{omega(t,tail) = tan(pi alpha / 2)}
#'otherwise.
#'
#'The characteristic function is inverted using Fourier's transform
#'to obtain the corresponding stable density. This inversion requires the
#'numerical evaluation of an integral from \eqn{0} to \eqn{\infty}{infinity}.
#'Two algorithms
#'are proposed for this. The default is Romberg's method
#'(\code{integration}="Romberg") which is used to evaluate the integral
#'with an error bounded by \code{eps}.
#'The alternative method is Simpson's integration
#'(\code{integration}="Simpson"): it approximates the
#'integral from \eqn{0} to \eqn{\infty}{infinity} by an integral
#'from \eqn{0} to \code{up}
#'with \code{npt} points subdividing \eqn{(O, up)}.
#'These three extra arguments -- \code{integration}, \code{up} and
#'\code{npt} -- are only available when using \code{dstable}.
#'The other functions are all based on Romberg's algorithm.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilites.
#' @param n number of observations.
#' @param loc vector of (real) location parameters.
#' @param disp vector of (positive) dispersion parameters.
#' @param skew vector of skewness parameters (in [-1,1]).
#' @param tail vector of parameters (in [0,2]) related to the tail thickness.
#' @param eps scalar giving the required precision in computation.
#' @param npt,up,integration As detailed herein -- only available when using \code{dstable}.
#'
#' @references Lambert, P. and Lindsey, J.K. (1999) Analysing financial returns using
#' regression models based on non-symmetric stable distributions. Applied
#' Statistics, 48, 409-424.
#'
#' @seealso \code{\link[stable]{stablereg}} to fit generalized nonlinear regression models
#' for the stable distribution parameters.
#'
#' @author Philippe Lambert (Catholic University of Louvain, Belgium, \email{phlambert@stat.ucl.ac.be})
#' @author Jim Lindsey
#'
#'
#  stable : A Library of Functions for Stable Distributions
#  Copyright (C) 1998, 1999, 2000, 2001 P. Lambert and J.K. Lindsey
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
#  SYNOPSIS
#
#     dstable(y, loc=0, disp=1/sqrt(2), skew=0, tail=2,
#		npt=501, up=10, eps=1.0e-6, integration="Romberg")
#     dstable2(y, loc, disp, skew, tail, npt=501, up=10,
#		integration="Romberg")
#     pstable(y, loc=0,disp=1/sqrt(2),skew=0,tail=2,eps=1.0e-6)
#     qstable(q, loc=0, disp=1/sqrt(2), skew=0, tail=2, eps=1.0e-6)
#     rstable(n=1,loc=0,disp=1/sqrt(2),skew=0,tail=2,eps=1.0e-6)
#     hstable(y, loc=0,disp=1/sqrt(2),skew=0,tail=2,eps=1.0e-6)
#     stablereg(y=NULL, loc=0, disp=1, skew=0, tail=1.5,
#		oloc=TRUE, odisp=TRUE, oskew=TRUE, otail=TRUE, noopt=FALSE,
#		iloc=NULL, idisp=NULL,iskew=NULL, itail=NULL,
#		loc_h=NULL, disp_h=NULL, skew_h=NULL, tail_h=NULL,
#		weights=1, exact=FALSE, delta=1, envir=parent.frame(),
#		integration="Romberg", eps=1e-6, up=10, npoint=501,
#               hessian=TRUE, llik.output=FALSE, print.level=0, ndigit=10,
#               steptol=0.00001, gradtol=0.00001, fscale=1, typsize=abs(p0),
#               stepmax=sqrt(p0%*%p0), iterlim=100)
#     stable.mode(loc,disp,skew,tail)
#
#  DESCRIPTION
#
#    Functions to obtain information about a stable distribution
#  and to fit regression models using this distribution.

# .First.lib <- function(lib, pkg)
# 	library.dynam("stable", pkg, lib)
# require(rmutil)

###########################################################################
# Density of a stable distribution
#
# This function produces the stable density computed at y.
# Note that it is obtained by integrating a function from 0 to Infinity.
# This integral is approximated by the finite integral from 0 to UP
# using the Simpson's method with npt points or Romberg's integration
#' @describeIn stable density
#' @export
#' @examples
#' par(mfrow=c(2,2))
#' x <- seq(-5,5,by=0.1)
#'
#' # Influence of loc (location)
#' plot(x,dstable(x,loc=-2,disp=1/sqrt(2),skew=-0.8,tail=1.5),
#'   type="l",ylab="",main="Varying LOCation")
#' lines(x,dstable(x,loc=0,disp=1/sqrt(2),skew=-0.8,tail=1.5))
#' lines(x,dstable(x,loc=2,disp=1/sqrt(2),skew=-0.8,tail=1.5))
#'
#' # Influence of disp (dispersion)
#' plot(x,dstable(x,loc=0,disp=0.5,skew=0,tail=1.5),
#'   type="l",ylab="",main="Varying DISPersion")
#' lines(x,dstable(x,loc=0,disp=1/sqrt(2),skew=0,tail=1.5))
#' lines(x,dstable(x,loc=0,disp=0.9,skew=0,tail=1.5))
#'
#' # Influence of skew (skewness)
#' plot(x,dstable(x,loc=0,disp=1/sqrt(2),skew=-0.8,tail=1.5),
#'   type="l",ylab="",main="Varying SKEWness")
#' lines(x,dstable(x,loc=0,disp=1/sqrt(2),skew=0,tail=1.5))
#' lines(x,dstable(x,loc=0,disp=1/sqrt(2),skew=0.8,tail=1.5))
#'
#' # Influence of tail (tail)
#' plot(x,dstable(x,loc=0,disp=1/sqrt(2),skew=0,tail=0.8),
#'   type="l",ylab="",main="Varying TAIL thickness")
#' lines(x,dstable(x,loc=0,disp=1/sqrt(2),skew=0,tail=1.5))
#' lines(x,dstable(x,loc=0,disp=1/sqrt(2),skew=0,tail=2))
#'
#' stabledist::dstable(x=1, 1, 0)
dstable <- function(x, loc=0, disp=1/sqrt(2), skew=0, tail=2,
		npt=501, up=10, eps=1.0e-6, integration="Romberg"){
  y<-x
	ly <- max(length(y),length(loc),length(disp),length(skew),length(tail))
	if(length(y)!=ly){
		if(length(y)==1)y <- rep(y,ly)
		else stop("length of y incorrect")}
	if(any(disp<0))stop("disp must be positive")
	if(any(skew< -1)||any(skew>1))stop("skew must lie in [-1,1]")
	if(length(skew)!=ly){
		if(length(skew)==1)skew <- rep(skew,ly)
		else stop("length of skew incorrect")}
	if(any(tail<=0)||any(tail>2))stop("tail must lie in (0,2]")
	if(length(tail)!=ly){
		if(length(tail)==1)tail <- rep(tail,ly)
		else stop("length of tail incorrect")}
	yy <- (y-loc)/disp
	ctail1 <- tail==1&skew==0
	ctail2 <- tail==2
	res <- vector(mode="numeric",ly)
	if(any(ctail1)){
		res[ctail1] <- dcauchy(yy[ctail1])
		ly <- ly-sum(ctail1)}
	if(any(ctail2)){
		res[ctail2] <- dnorm(yy[ctail2],0,sqrt(2))
		ly <- ly-sum(ctail2)}
	if(ly>0)res[!ctail1&!ctail2] <- .C("stable",
		as.integer(ly),
		as.double(yy[!ctail1&!ctail2]),
		as.double(skew[!ctail1&!ctail2]),
		as.double(tail[!ctail1&!ctail2]),
		as.integer(npt),
		as.double(up),
		as.double(eps),
		as.integer(if(integration=="Simpson")1 else 2),
		err=integer(1),
		ffy=double(ly),
		PACKAGE="stable")$ffy
	res/disp}

###########################################################################
# c.d.f. of a stable distribution
#
#' @describeIn stable cdf
#' @export
pstable <- function(q, loc=0,disp=1/sqrt(2),skew=0,tail=2,eps=1.0e-6){
  y <- q
	yy <- (y-loc)/disp
	ly <- max(length(y),length(loc),length(disp),length(skew),length(tail))
	if(length(y)!=ly){
		if(length(y)==1)y <- rep(y,ly)
		else stop("length of y incorrect")}
	if(any(disp<0))stop("disp must be positive")
	if(any(skew< -1)||any(skew>1))stop("skew must lie in [-1,1]")
	if(length(skew)!=ly){
		if(length(skew)==1)skew <- rep(skew,ly)
		else stop("length of skew incorrect")}
	if(any(tail<=0)||any(tail>2))stop("tail must lie in (0,2]")
	if(length(tail)!=ly){
		if(length(tail)==1)tail <- rep(tail,ly)
		else stop("length of tail incorrect")}
	ctail1 <- tail==1&skew==0
	ctail2 <- tail==2
	res <- vector(mode="numeric",ly)
	if(any(ctail1)){
		res[ctail1] <- pcauchy(yy[ctail1])
		ly <- ly-sum(ctail1)}
	if(any(ctail2)){
		res[ctail2] <- pnorm(yy[ctail2],0,sqrt(2))
		ly <- ly-sum(ctail2)}
	if(ly>0)res[!ctail1&!ctail2] <- .C("pstable_c",
		as.integer(ly),
		as.double(yy[!ctail1&!ctail2]),
		as.double(skew[!ctail1&!ctail2]),
		as.double(tail[!ctail1&!ctail2]),
		as.double(eps),
		err=integer(1),
		ffy=double(ly),
		PACKAGE="stable")$ffy
	res}

###########################################################################
# Quantile function of a stable random variable
#
#' @describeIn stable quantiles
#' @export
qstable <- function(p, loc=0, disp=1/sqrt(2), skew=0, tail=2, eps=1.0e-6){
  q <- p
	h <- function(y).C("pstable_c",
		as.integer(1),
		as.double((y-loc[i])/disp[i]),
		as.double(skew[i]),
		as.double(tail[i]),
		as.double(eps),
		err=integer(1),
		ffy=double(1),
		PACKAGE="stable")$ffy-q[i]
	lq <- max(length(q),length(loc),length(disp),length(skew),length(tail))
	if(length(q)!=lq){
		if(length(q)==1)q <- rep(q,lq)
		else stop("length of q incorrect")}
	if(length(loc)!=lq){
		if(length(loc)==1)loc <- rep(loc,lq)
		else stop("length of loc incorrect")}
	if(any(disp<0))stop("disp must be positive")
	if(length(disp)!=lq){
		if(length(disp)==1)disp <- rep(disp,lq)
		else stop("length of disp incorrect")}
	if(any(skew< -1)||any(skew>1))stop("skew must lie in [-1,1]")
	if(length(skew)!=lq){
		if(length(skew)==1)skew <- rep(skew,lq)
		else stop("length of skew incorrect")}
	if(any(tail<=0)||any(tail>2))stop("tail must lie in (0,2]")
	if(length(tail)!=lq){
		if(length(tail)==1)tail <- rep(tail,lq)
		else stop("length of tail incorrect")}
	tmp <- vector(mode="numeric",lq)
	for (i in (1:lq)){
		if(tail[i]==1&&skew[i]==0)tmp[i] <- qcauchy(q[i],loc[i],disp[i])
		else if(tail[i]==2)tmp[i] <- qnorm(q[i],loc[i],disp[i]*sqrt(2))
		else {
			if(tail[i]<1&&abs(skew[i])==1){
				interval <- sign(skew[i])*c(0.001,10)
				while(TRUE){
					tmp1 <- h(interval[1])
					tmp2 <- h(interval[2])
					if(is.na(tmp1*tmp2)||tmp1*tmp2<0)break
					interval[1] <- 0.5*interval[1]
					interval[2] <- 2*interval[2]}}
			else {
				interval <- loc[i]+disp[i]*c(-5,5)
				while(TRUE){
					tmp1 <- h(interval[1])
					tmp2 <- h(interval[2])
					if(is.na(tmp1*tmp2)||tmp1*tmp2<0)break
					interval <- 2*interval}}
			tmp[i] <- if(is.na(tmp1*tmp2))NA
				  else uniroot(h,interval)$root}}
	tmp}

###########################################################################
# Generation of stable random deviates
#
#' @describeIn stable random deviates
#' @export
rstable <- function(n=1,loc=0,disp=1/sqrt(2),skew=0,tail=2,eps=1.0e-6){
	tmp1 <- qstable(runif(n),loc=loc,disp=disp,skew=skew,tail=tail,eps=eps)
	tmp2 <- tmp1[!is.na(tmp1)]
	n1 <- n
	while(any(is.na(tmp1))){
		n1 <- n1-sum(!is.na(tmp1))
		tmp1 <- qstable(runif(n1),loc=loc,disp=disp,skew=skew,tail=tail,eps=eps)
		tmp2 <- c(tmp2,tmp1[!is.na(tmp1)])}
	tmp2}

###########################################################################
# Stable hazard
#
#' @describeIn stable hazard
#' @export
hstable <- function(x, loc=0,disp=1/sqrt(2),skew=0,tail=2,eps=1.0e-6){
  y <- x
	dstable(y,loc=loc,disp=disp,skew=skew,tail=tail,eps=eps)/
		(1-pstable(y,loc=loc,disp=disp,skew=skew,tail=tail,eps=eps))}

###########################################################################
# Link and inverse link functions for use in stablereg

#' Links
#'
#' Link and inverse functions for use in stablereg
#'
#' @name Links
#' @aliases loc_g loc_h disp_g disp_h skew_g skew_h tail_g tail_h
#' @param x the function argument
#' @rdname Links
#' @export loc_g
loc_g <- function(x) x
#' @rdname Links
#' @export loc_h
loc_h <- function(x) x
#' @rdname Links
#' @export disp_g
disp_g <- function(x) log(x)
#' @rdname Links
#' @export disp_h
disp_h <- function(x) exp(x)
#' @rdname Links
#' @export skew_g
skew_g <- function(x) log((1+x)/(1-x))
#' @rdname Links
#' @export skew_h
skew_h <- function(x) 2/(1+exp(-x))-1
#' @rdname Links
#' @export tail_g
tail_g <- function(x) log(x/(2-x))
#' @rdname Links
#' @export tail_h
tail_h <- function(x) 2/(1+exp(-x))

###########################################################################
# Regression models for the four parameters in the stable distribution
# Note that the returned optimized parameters are on the normalized
# scale, i.e. that they have to be transformed back to the right scale.



#' Stable Generalized Regression Models
#'
#' \code{stablereg} fits user specified generalized linear and nonlinear
#' regression models based on the stable distribution to (uncensored, right
#' and/or left censored) data. This allows the location, the dispersion, the
#' skewness and the tails of the fitted stable distribution to vary with
#' explanatory variables.
#'
#'
#' @aliases stablereg fitted.stable df.residual.stable deviance.stable
#' aic.stable
#' @param y The response vector or a \code{repeated} data object. If the
#' \code{repeated} data object contains more than one response variable, give
#' that object in \code{envir} and give the name of the response variable to be
#' used here.
#'
#' For censored data, two columns with the second being the censoring indicator
#' (1: uncensored, 0: right censored, -1: left censored.)
#' @param loc,loc_h,oloc,iloc Describe the regression model fitted for the
#' location parameter of the stable distribution, perhaps after transformation
#' by the link function \code{loc_g} (set to the identity by default. The
#' inverse link function is denoted by \code{loc_h}. Note that these functions
#' cannot contain unknown parameters).
#'
#' Two specifications are possible:
#'
#' (1) \code{loc} is a linear or nonlinear language expression beginning with ~
#' or an R function, describing the regression function for the location
#' parameter (after transformation by \code{loc_g}, the link function).
#'
#' \code{iloc} is a vector of initial conditions for the parameters in the
#' regression for this parameter.
#'
#' \code{oloc} is a boolean indicating if an optimization of the likelihood has
#' to be carried out on these parameters. If \code{oloc} is set to TRUE, a
#' default zero value is considered for the starting values \code{iloc}. But if
#' no optimization is desired on the location parameters, i.e. when the
#' likelihood has to be evaluated or optimized at a fixed location, then
#' \code{iloc} has to be explicitely specified.
#'
#' (2) \code{loc} is a numeric expression (i.e. a scalar or a vector of the
#' same size as the data vector \code{y}, or \code{y[,1]} when censoring is
#' considered).
#'
#' If \code{oloc} is set to TRUE, i.e. when an optimization of the likelihood
#' has to be carried out on the location parameter, then the location parameter
#' (after transformation by the link function loc_g) is set to an unknown
#' parameter with initial value equal to \code{iloc[1]} or \code{loc[1]} when
#' \code{iloc} is not specified.
#'
#' But when \code{oloc} is set to FALSE, i.e. when the likelihood has to be
#' evaluated or optimized at a fixed location, then the transformed location is
#' assumed to be equal to \code{loc} when it is of the same length as the data
#' vector \code{y} (or \code{y[,1]} when censoring is considered), and to
#' \code{loc[1]} otherwise.
#'
#' Specification (1) is especially useful in ANOVA-like situations where the
#' location is assumed to change with the levels of some factor variable.
#' @param disp,disp_h,odisp,idisp describe the regression model for the
#' dispersion parameter of the fitted stable distribution, after transformation
#' by the link function \code{disp_g} (set to the \code{log} function by
#' default). The inverse link function is denoted by \code{disp_h}. Again these
#' functions cannot contain unknown parameters. The same rules as above apply
#' when specifying the generalized regression model for the dispersion
#' parameter.
#' @param skew,skew_h,oskew,iskew describe the regression model for the
#' skewness parameter of the fitted stable distribution, after transformation
#' by the link function \code{skew_g} (set to \code{log{(1 + [.])/(1 - [.])}}
#' by default). The inverse link function is denoted by \code{skew_h}. Again
#' these functions cannot contain unknown parameters. The same rules as above
#' apply when specifying the generalized regression model for the skewness
#' parameter.
#' @param tail,tail_h,otail,itail describe the regression model considered
#' for the tail parameter of the fitted stable distribution, after
#' transformation by the link function \code{tail_g} (set to \code{log{([.] -
#' 1)/(2 - [.])}} by default. The inverse link function is denoted by
#' \code{tail_h}. Again these functions cannot contain unknown parameters). The
#' same rules as above apply when specifying the generalized regression model
#' for the tail parameter.
#' @param noopt When set to TRUE, it forces \code{oloc}, \code{odisp},
#' \code{oskew} and \code{otail} to FALSE, whatever the user choice for these
#' last three arguments. It is especially useful when looking for appropriate
#' initial values for the regression model parameters, before undertaking the
#' optimization of the likelihood.
#' @param weights Weight vector.
#' @param exact If TRUE, fits the exact likelihood function for continuous data
#' by integration over intervals of observation, i.e. interval censoring.
#' @param delta Scalar or vector giving the unit of measurement for each
#' response value, set to unity by default. For example, if a response is
#' measured to two decimals, \code{delta=0.01}. If the response is transformed,
#' this must be multiplied by the Jacobian.  For example, with a log
#' transformation, \code{delta=1/y}. (The \code{delta} values for the censored
#' response are ignored.) The transformation cannot contain unknown parameters.
#' @param envir Environment in which model formulae are to be interpreted or a
#' data object of class, \code{repeated}, \code{tccov}, or \code{tvcov}; the
#' name of the response variable should be given in \code{y}. If \code{y} has
#' class \code{repeated}, it is used as the environment.
#' @param integration,eps,up,npoint \code{integration} indicates which
#' algorithm must be used to evaluate the stable density when the likelihood is
#' computed with \code{exact} set to FALSE. See the man page on
#' \code{stable} for extra information.
#' @param llik.output is TRUE when the likelihood has to be displayed at each
#' iteration of the optimization.
#' @param print.level Arguments controlling the optimization procedure \code{\link{nlm}}.
#' @param ndigit Arguments controlling the optimization procedure \code{\link{nlm}}.
#' @param steptol Arguments controlling the optimization procedure \code{\link{nlm}}.
#' @param gradtol Arguments controlling the optimization procedure \code{\link{nlm}}.
#' @param fscale Arguments controlling the optimization procedure \code{\link{nlm}}.
#' @param typsize Arguments controlling the optimization procedure \code{\link{nlm}}.
#' @param stepmax Arguments controlling the optimization procedure \code{\link{nlm}}.
#' @param iterlim Arguments controlling the optimization procedure \code{\link{nlm}}.
#' @param hessian Arguments controlling the optimization procedure \code{\link{nlm}}.
#' @return A list of class \code{stable} is returned. The printed output
#' includes the -log-likelihood, the corresponding AIC, the maximum likelihood
#' estimates, standard errors, and correlations. It also include all the
#' relevant information calculated, including error codes.
#' @section Warning: Because of the numerical integrations involved,
#' convergence can be very sensitive to the initial parameter values supplied
#' and to the settings of the arguments controlling \code{\link{nlm}}. If nlm
#' feeds extreme parameter values in the tails of the distribution to the
#' likelihood function, the integration may hang for a long time.
#' @author Philippe Lambert (Catholic University of Louvain, Belgium,
#' \email{phlambert@@stat.ucl.ac.be}) and Jim Lindsey.
#' @seealso \code{\link{lm}}, \code{\link{glm}}, \code{stable}
#' and \code{stable.mode}.
#' @references Lambert, P. and Lindsey, J.K. (1999) Analysing financial returns
#' using regression models based on non-symmetric stable distributions. Applied
#' Statistics 48, 409-424.
#' @keywords models
#' @examples
#'
#' ## Share return over a 50 day period (see reference above)
#' # shares
#' y <- c(296,296,300,302,300,304,303,299,293,294,294,293,295,287,288,297,
#' 305,307,307,304,303,304,304,309,309,309,307,306,304,300,296,301,298,
#' 295,295,293,292,297,294,293,306,303,301,303,308,305,302,301,297,299)
#'
#' # returns
#' ret <- (y[2:50]-y[1:49])/y[1:49]
#' # hist(ret, breaks=seq(-0.035,0.045,0.01))
#'
#' day <- seq(0,0.48,by=0.01) # time measured in days/100
#' x <- seq(1,length(ret))-1
#'
#' # Classic stationary normal model tail=2
#' print(z1 <- stablereg(y = ret, delta = 1/y[1:49],
#' 	loc = ~1, disp= ~1, skew = ~1, tail = tail_g(1.9999999),
#' 	iloc = 0, idisp = -3, iskew = 0, oskew = FALSE, otail = FALSE))
#'
#' # Normal model (tail=2) with dispersion=disp_h(b0+b1*day)
#' print(z2 <- stablereg(y = ret, delta = 1/y[1:49], loc = ~day,
#' 	disp = ~1, skew = ~1, tail = tail_g(1.999999), iloc = c(0.003,0),
#' 	idisp = -4.5, iskew = 0, oskew = FALSE, otail = FALSE))
#'
#' # Stable model with loc(ation)=loc_h(b0+b1*day)
#' print(z3 <- stablereg(y = ret, delta = 1/y[1:49],
#' 	loc = ~day, disp = ~1, skew = ~1, tail = ~1,
#' 	iloc = c(0.001,-0.004), idisp = -4.8, iskew = 0, itail = 0.6))
#'
#' # Stable model with disp(ersion)=disp_h(b0+b1*day)
#' print(z4 <- stablereg(y = ret, delta = 1/y[1:49],
#' 	loc = ~1, disp = ~day, skew = ~1, tail = ~1,
#' 	iloc = 0.003, idisp = c(-4.8,0), iskew = -0.03, itail = 1.6))
#'
#' # Stable model with skew(ness)=skew_h(b0+b1*day)
#' # Evaluation at fixed parameter values (because noopt is set to TRUE)
#' print(z5 <- stablereg(y = ret, delta = 1/y[1:49],
#' 	loc = ~1, disp = ~1, skew = ~day, tail = ~1,
#' 	iloc = 5.557e-04, idisp = -4.957, iskew = c(2.811,-2.158),
#' 	itail = 1.57, noopt=TRUE))
#'
#' # Stable model with tail=tail_h(b0+b1*day)
#' print(z6 <- stablereg(y = ret, delta = 1/y[1:49], loc = ret ~ 1,
#' 	disp = ~1, skew = ~1, tail = ~day, iloc = 0.002,
#' 	idisp = -4.8, iskew = -2, itail = c(2.4,-4), hessian=FALSE))
#'
#' @importFrom stats dcauchy dnorm nlm pcauchy pnorm qcauchy qnorm runif uniroot
#' @import  rmutil
#' @useDynLib stable
#' @export
stablereg <- function(y=NULL, loc=0, disp=1, skew=0, tail=1.5,
	oloc=TRUE, odisp=TRUE, oskew=TRUE, otail=TRUE, noopt=FALSE,
	iloc=NULL, idisp=NULL,iskew=NULL, itail=NULL,
	loc_h=NULL, disp_h=NULL, skew_h=NULL, tail_h=NULL,
	weights=1, exact=FALSE, delta=1, envir=parent.frame(),
	integration="Romberg", eps=1e-6, up=10, npoint=501,
	hessian=TRUE, llik.output=FALSE, print.level=0, ndigit=10,
	steptol=0.00001, gradtol=0.00001, fscale=1, typsize=abs(p0),
	stepmax=sqrt(p0%*%p0), iterlim=100){
call <- sys.call()

if(is.null(loc_h))loc_h <- function(x) x
else if(!is.function(loc_h))stop("loc_h must be a link function")
if(is.null(disp_h))disp_h <- function(x) exp(x)
else if(!is.function(disp_h))stop("disp_h must be a link function")
if(is.null(skew_h))skew_h <- function(x) 2/(1+exp(-x))-1
else if(!is.function(skew_h))stop("skew_h must be a link function")
if(is.null(tail_h))tail_h <- function(x) 2/(1+exp(-x))
else if(!is.function(tail_h))stop("tail_h must be a link function")

# If don't want to optimize at all, just set noopt=TRUE
if(noopt)oloc <- odisp <- oskew <- otail <- FALSE

respenv <- exists(deparse(substitute(y)),envir=parent.frame())&&
	inherits(y,"repeated")&&!inherits(envir,"repeated")
if(respenv){
	if(dim(y$response$y)[2]>1)
		stop("stablereg only handles univariate responses")
	if(!is.null(y$NAs)&&any(y$NAs))
		stop("stablereg does not handle data with NAs")}
envname <- if(respenv)deparse(substitute(y))
	else if(!is.null(class(envir)))deparse(substitute(envir))
	else NULL

nploc <- length(iloc)
npdisp <- length(idisp)
npskew <- length(iskew)
nptail <- length(itail)
loc1 <- disp1 <- skew1 <- tail1 <- loc3 <- disp3 <- skew3 <- tail3 <- NULL
if(respenv||inherits(envir,"repeated")||inherits(envir,"tccov")||inherits(envir,"tvcov")){
	if(is.null(envname))envname <- deparse(substitute(envir))
	if(inherits(loc,"formula")){
		loc3 <- if(respenv)finterp(loc,.envir=y,.name=envname)
			else finterp(loc,.envir=envir,.name=envname)}
	else if(is.function(loc)){
		if(is.null(attr(loc,"model"))){
		        tmp <- parse(text=deparse(loc)[-1])
		        loc <- if(respenv)fnenvir(loc,.envir=y,.name=envname)
		        	else fnenvir(loc,.envir=envir,.name=envname)
		        nploc <- length(attr(loc,"parameters"))
		        loc3 <- loc
			attr(loc3,"model") <- tmp}
		else loc3 <- loc}
	if(inherits(disp,"formula")){
		disp3 <- if(respenv)finterp(disp,.envir=y,.name=envname)
			else finterp(disp,.envir=envir,.name=envname)}
	else if(is.function(disp)){
		if(is.null(attr(disp,"model"))){
		        tmp <- parse(text=deparse(disp)[-1])
		        disp <- if(respenv)fnenvir(disp,.envir=y,.name=envname)
		        	else fnenvir(disp,.envir=envir,.name=envname)
		        npdisp <- length(attr(disp,"parameters"))
		        disp3 <- disp
		        attr(disp3,"model") <- tmp}
		else disp3 <- disp}
	if(inherits(skew,"formula")){
		skew3 <- if(respenv)finterp(skew,.envir=y,.name=envname)
			else finterp(skew,.envir=envir,.name=envname)}
	else if(is.function(skew)){
		if(is.null(attr(skew,"model"))){
		        tmp <- parse(text=deparse(skew)[-1])
		        skew <- if(respenv)fnenvir(skew,.envir=y,.name=envname)
		        	else fnenvir(skew,.envir=envir,.name=envname)
		        npskew <- length(attr(skew,"parameters"))
		        skew3 <- skew
		        attr(skew3,"model") <- tmp}
		else skew3 <- skew}
	if(inherits(tail,"formula")){
		tail3 <- if(respenv)finterp(tail,.envir=y,.name=envname)
			else finterp(tail,.envir=envir,.name=envname)}
	else if(is.function(tail)){
		if(is.null(attr(tail,"model"))){
		        tmp <- parse(text=deparse(tail)[-1])
		        tail <- if(respenv)fnenvir(tail,.envir=y,.name=envname)
		        	else fnenvir(tail,.envir=envir,.name=envname)
		        nptail <- length(attr(tail,"parameters"))
		        tail3 <- tail
		        attr(tail3,"model") <- tmp}
		else tail3 <- tail}}
#
# if data object supplied, find response information in it
#
type <- "unknown"
if(respenv){
	if(dim(y$response$y)[2]>1)
		stop("stablereg only handles univariate responses")
	if(!is.null(y$NAs)&&any(y$NAs))
		stop("stablereg does not handle data with NAs")
	if(inherits(envir,"repeated")&&(length(nobs(y))!=length(nobs(envir))||any(nobs(y)!=nobs(envir))))
		stop("y and envir objects are incompatible")
	if(!is.null(y$response$delta))
		delta <- as.vector(y$response$delta)
	type <- y$response$type
	y <- response(y)
	if(dim(y)[2]==1)y <- as.vector(y)}
else if(inherits(envir,"repeated")){
	if(!is.null(envir$NAs)&&any(envir$NAs))
		stop("stablereg does not handle data with NAs")
	cn <- deparse(substitute(y))
	if(length(grep("\"",cn))>0)cn <- y
	if(length(cn)>1)stop("only one y variable allowed")
	col <- match(cn,colnames(envir$response$y))
	if(is.na(col))stop(paste("response variable",cn,"not found"))
	type <- envir$response$type[col]
	y <- envir$response$y[,col]
	if(!is.null(envir$response$n)&&!all(is.na(envir$response$n[,col])))
		y <- cbind(y,envir$response$n[,col]-y)
	else if(!is.null(envir$response$censor)&&!all(is.na(envir$response$censor[,col])))
		y <- cbind(y,envir$response$censor[,col])
	if(!is.null(envir$response$wt))wt <- as.vector(envir$response$wt)
	if(!is.null(envir$response$delta))
		delta <- as.vector(envir$response$delta[,col])}
else if(inherits(y,"response")){
	if(inherits(y,"multivariate"))
		stop("nordr only handles univariate responses")
	if(!is.null(y$delta))delta <- as.vector(y$delta)
	type <- y$type
	y <- response(y)
	if(dim(y)[2]==1)y <- as.vector(y)}
if(any(is.na(y)))stop("NAs in y - use rmna")
if(type!="unknown"&&type!="continuous"&&type!="duration")
	stop("continuous or duration data required")
if(is.matrix(y)){
	if(length(dim(y))!=2||dim(y)[2]!=2)
		stop("Two column matrix required for response\nTimes and censor indicator")
	else n <- dim(y)[1]
	censor <- TRUE}
else {
	if(!is.vector(y))stop("y must be a vector")
	n <- length(y)
	censor <- FALSE}
if(censor){
	y[,2] <- as.integer(y[,2])
	if(y[,2]!=-1&y[,2]!=0&y[,2]!=1)
		stop("Censor indicator must be -1s, 0s, and 1s")
	cc <- if(y[,2]==1)1 else 0  # observed
	rc <- if(y[,2]==0)1
		else {
			if(y[,2]==-1)-1 else 0}  # right censored
	lc <- if(y[,2]==-1)0 else 1 # left censored
	if(delta<=0&y[,2]==1)
		stop("All deltas for uncensored data must be positive")}

if(censor)y <- y[,1]

if(length(weights)==1)weights <- rep(weights,n)

# LOC
# From the coming lines, we see that if one does not want to optimize over
# loc (ie. oloc=FALSE), then there are 2 alternatives:
# (1) set loc equal to a value identical for all units: use loc=<scalar>
# (2) set loc equal to values varying through units: use loc=<language>
#      and iloc=<corresponding vector of initial conditions>
# Note that we work with the log.g link (identity link by default).
# This means that language description of the loc model and
# initial conditions are understood on that scale!!
if(inherits(loc,"formula")){
	nploc <- length(iloc)
	loc2 <- if(respenv)finterp(loc,.envir=y,.name=envname)
		else finterp(loc,.envir=envir,.name=envname)
	npt1 <- length(attr(loc2,"parameters"))
	if(is.character(attr(loc2,"model"))){
		if(length(attr(loc2,"model"))==1){
			loc1 <- function(p) loc_h(p[1]*rep(1,n))
			attributes(loc1) <- attributes(loc2)
			loc2 <- NULL}}
	else {
		if(nploc!=npt1){
			cat("\nParameters are ")
			cat(attr(loc2,"parameters"),"\n")
			stop(paste("iloc should have",npt1,"estimates"))}
		if(is.list(iloc)){
			if(!is.null(names(iloc))){
				o <- match(attr(loc2,"parameters"),names(iloc))
				iloc <- unlist(iloc)[o]
				if(sum(!is.na(o))!=length(iloc))stop("invalid estimates for loc - probably wrong names")}
			else iloc <- unlist(iloc)}}
	if(!is.null(loc2)){
		if(inherits(envir,"tccov")){
			cv <- covind(y)
			loc1 <- function(p) loc_h(loc2(p)[cv])
			attributes(loc1) <- attributes(loc2)}
		else loc1 <- function(p) loc_h(loc2(p))
		attributes(loc1) <- attributes(loc2)}}
else if(is.function(loc))loc1 <- loc
if(!is.null(loc1)&&is.null(attr(loc1,"parameters"))){
	attributes(loc1) <- if(is.function(loc)){
		if(!inherits(loc,"formulafn")){
			if(respenv)attributes(fnenvir(loc,.envir=y))
			else attributes(fnenvir(loc,.envir=envir))}
		else attributes(loc)}
		else {
			if(respenv)attributes(fnenvir(loc1,.envir=y))
			else attributes(fnenvir(loc1,.envir=envir))}}
nlp <- if(is.function(loc1))length(attr(loc1,"parameters"))
	else if(is.null(loc1))NULL
	else npt1
if(!is.null(nlp)&&nlp!=nploc)
	stop(paste("iloc should have",nlp,"initial estimates"))
if(inherits(loc,"formula")||is.function(loc)){
	if(oloc){
		if(is.numeric(iloc) && length(iloc)!=nploc)
			stop(paste("iloc must be of size ",nploc))
		if(!is.numeric(iloc)) iloc <- rep(0,nploc)
		fnloc <- loc1}
	else {
		if(!is.numeric(iloc))
			stop("Missing initial conditions for loc")
		else if(length(iloc)!=nploc)
			stop(paste("iloc must be of size ",nploc))
		else fnloc <- function(p) loc1(iloc)}}
else if(oloc){
	fnloc <- function(p) loc_h(rep(p[1],n))
	if(!is.numeric(iloc)){
		iloc <- loc[1]
		nploc <- 1}}
else {
					# IMPORTANT
	if(length(loc)==n)fnloc <- function(p) loc_h(loc)
					# IMPORTANT
	else fnloc <- function(p) rep(loc_h(loc[1]),n)
	nploc <- 1}

# DISP
# Note that we work with the disp_g link (log link by default).
# This means that language description of the disp model and
# initial conditions are understood on that scale!!
# Non language description of disp are also understood on the disp_g scale
# and specified using disp=<value>, yielding disp_h(<value>) for the
# parameter disp.

npl2 <- if(odisp)nploc*oloc+1 else 1
if(inherits(disp,"formula")){
	npdisp <- length(idisp)
	disp2 <- if(respenv)finterp(disp,.envir=y,.start=npl2,.name=envname)
		else finterp(disp,.envir=envir,.start=npl2,.name=envname)
	npt2 <- length(attr(disp2,"parameters"))
	if(is.character(attr(disp2,"model"))){
		if(length(attr(disp2,"model"))==1){
			disp1 <- function(p) disp_h(p[npl2]*rep(1,n))
			attributes(disp1) <- attributes(disp2)
			disp2 <- NULL}}
	else {
		if(npdisp!=npt2){
			cat("\nParameters are ")
			cat(attr(disp2,"parameters"),"\n")
			stop(paste("idisp should have",npt2,"estimates"))}
		if(is.list(idisp)){
			if(!is.null(names(idisp))){
				o <- match(attr(disp2,"parameters"),names(idisp))
				idisp <- unlist(idisp)[o]
				if(sum(!is.na(o))!=length(idisp))stop("invalid estimates for disp - probably wrong names")}
			else idisp <- unlist(idisp)}}
	if(!is.null(disp2)){
		if(inherits(envir,"tccov")){
			cv <- covind(y)
			disp1 <- function(p) disp_h(disp2(p)[cv])}
		else disp1 <- function(p) disp_h(disp2(p))
		attributes(disp1) <- attributes(disp2)}}
else if(is.function(disp))disp1 <- function(p) disp(p[npl2:(npl2+npdisp-1)])
if(!is.null(disp1)&&is.null(attr(disp1,"parameters"))){
	attributes(disp1) <- if(is.function(disp)){
		if(!inherits(disp,"formulafn")){
			if(respenv)attributes(fnenvir(disp,.envir=y))
			else attributes(fnenvir(disp,.envir=envir))}
		else attributes(disp)}
		else {
			if(respenv)attributes(fnenvir(disp1,.envir=y))
			else attributes(fnenvir(disp1,.envir=envir))}}
nlp <- if(is.function(disp1))length(attr(disp1,"parameters"))
	else if(is.null(disp1))NULL
	else npt2
if(!is.null(nlp)&&nlp!=npdisp)
	stop(paste("idisp should have",nlp,"initial estimates"))
if(inherits(disp,"formula")||is.function(disp)){
	if(odisp){
		if(is.numeric(idisp) && length(idisp)!=npdisp)
			stop(paste("idisp must be of size ",npdisp))
		else if(!is.numeric(idisp)) idisp <- rep(0,npdisp)
		fndisp <- disp1}
	else {
		if(!is.numeric(idisp))
			stop("Missing initial conditions for disp")
		else if(length(idisp)!=npdisp)
			stop(paste("idisp must be of size ",npdisp))
		else fndisp <- function(p) disp1(idisp)}}
else if(odisp){
	fndisp <- function(p) disp_h(rep(p[npl2],n))
	if(!is.numeric(idisp)){
		idisp <- disp[1]
		npdisp <- 1}}
else {
				# IMPORTANT
	if(length(disp)==n)fndisp <- function(p) disp_h(disp)
				# IMPORTANT
	else fndisp <- function(p) rep(disp_h(disp[1]),n)
	npdisp <- 1}

# SKEW
# Note that, y default, we work with the skew_g([.])=log{(1+[.])/(1-[.])} link.
# This means that language description of the skew model and
# the possible initial conditions are understood on that scale!!
# Non language description of skew are also understood on the skew_g scale
# and specified using skew=<value>, yielding skew_h(<value>) for the
# parameter skew.

npl3 <- if(oskew)nploc*oloc+npdisp*odisp+1 else 1
if(inherits(skew,"formula")){
	npskew <- length(iskew)
	skew2 <- if(respenv)finterp(skew,.envir=y,.start=npl3,.name=envname)
		else finterp(skew,.envir=envir,.start=npl3,.name=envname)
	npt3 <- length(attr(skew2,"parameters"))
	if(is.character(attr(skew2,"model"))){
		if(length(attr(skew2,"model"))==1){
			skew1 <- function(p) skew_h(p[npl3]*rep(1,n))
			attributes(skew1) <- attributes(skew2)
			skew2 <- NULL}}
	else {
		if(npskew!=npt3){
			cat("\nParameters are ")
			cat(attr(skew2,"parameters"),"\n")
			stop(paste("iskew should have",npt3,"estimates"))}
		if(is.list(iskew)){
			if(!is.null(names(iskew))){
				o <- match(attr(skew2,"parameters"),names(iskew))
				iskew <- unlist(iskew)[o]
				if(sum(!is.na(o))!=length(iskew))stop("invalid estimates for skew - probably wrong names")}
			else iskew <- unlist(iskew)}}
	if(!is.null(skew2)){
		if(inherits(envir,"tccov")){
			cv <- covind(y)
			skew1 <- function(p) skew_h(skew2(p)[cv])}
		else skew1 <- function(p) skew_h(skew2(p))
		attributes(skew1) <- attributes(skew2)}}
else if(is.function(skew))skew1 <- function(p) skew(p[npl3:(npl3+npskew-1)])
if(!is.null(skew1)&&is.null(attr(skew1,"parameters"))){
	attributes(skew1) <- if(is.function(skew)){
		if(!inherits(skew,"formulafn")){
			if(respenv)attributes(fnenvir(skew,.envir=y))
			else attributes(fnenvir(skew,.envir=envir))}
		else attributes(skew)}
		else {
			if(respenv)attributes(fnenvir(skew1,.envir=y))
			else attributes(fnenvir(skew1,.envir=envir))}}
nlp <- if(is.function(skew1))length(attr(skew1,"parameters"))
	else if(is.null(skew1))NULL
	else npt3
if(!is.null(nlp)&&nlp!=npskew)
	stop(paste("iskew should have",nlp,"initial estimates"))
if(inherits(skew,"formula")||is.function(skew)){
	if(oskew){
		if(is.numeric(iskew) && length(iskew)!=npskew)
			stop(paste("iskew must be of size ",npskew))
		else if(!is.numeric(iskew)) iskew <- rep(0,npskew)
		fnskew <- skew1}
	else {
		if(!is.numeric(iskew))
			stop("Missing initial conditions for skew")
		else if(length(iskew)!=npskew)
			stop(paste("iskew must be of size ",npskew))
		else fnskew <- function(p) skew1(iskew)}}
else if(oskew){
	fnskew <- function(p) skew_h(rep(p[npl3],n))
	if(!is.numeric(iskew)){
		iskew <- skew[1]
		npskew <- 1}}
else {
	if(length(skew)==n) fnskew <- function(p) skew_h(skew)
	else fnskew <- function(p) rep(skew_h(skew[1]),n) # IMPORTANT
	npskew <- 1}

# TAIL
# Note that we work with the tail_g([.])=log{([.]-1)/(2-[.])} link.
# This means that language description of the tail model and
# the possible initial conditions are understood on that scale!!
# Non language description of tail are also understood on the tail_g scale
# and specified using tail=<value>, yielding tail_h(<value>) for the
# parameter tail.

npl4 <- if(otail)nploc*oloc+npdisp*odisp+npskew*oskew+1 else 1
if(inherits(tail,"formula")){
	nptail <- length(itail)
	tail2 <- if(respenv)finterp(tail,.envir=y,.start=npl4,.name=envname)
		else finterp(tail,.envir=envir,.start=npl4,.name=envname)
	npt4 <- length(attr(tail2,"parameters"))
	if(is.character(attr(tail2,"model"))){
		if(length(attr(tail2,"model"))==1){
			tail1 <- function(p) tail_h(p[npl4]*rep(1,n))
			attributes(tail1) <- attributes(tail2)
			tail2 <- NULL}}
	else {
		if(nptail!=npt4){
			cat("\nParameters are ")
			cat(attr(tail2,"parameters"),"\n")
			stop(paste("itail should have",npt4,"estimates"))}
		if(is.list(itail)){
			if(!is.null(names(itail))){
				o <- match(attr(tail2,"parameters"),names(itail))
				itail <- unlist(itail)[o]
				if(sum(!is.na(o))!=length(itail))stop("invalid estimates for tail - probably wrong names")}
			else itail <- unlist(itail)}}
	if(!is.null(tail2)){
		if(inherits(envir,"tccov")){
			cv <- covind(y)
			tail1 <- function(p) tail_h(tail2(p)[cv])}
		else tail1 <- function(p) tail_h(tail2(p))
		attributes(tail1) <- attributes(tail2)}}
else if(is.function(tail))tail1 <- function(p) tail(p[npl4:(npl4+nptail-1)])
if(!is.null(tail1)&&is.null(attr(tail1,"parameters"))){
	attributes(tail1) <- if(is.function(tail)){
		if(!inherits(tail,"formulafn")){
			if(respenv)attributes(fnenvir(tail,.envir=y))
			else attributes(fnenvir(tail,.envir=envir))}
		else attributes(tail)}
		else {
			if(respenv)attributes(fnenvir(tail1,.envir=y))
			else attributes(fnenvir(tail1,.envir=envir))}}
nlp <- if(is.function(tail1))length(attr(tail1,"parameters"))
	else if(is.null(tail1))NULL
	else npt4
if(!is.null(nlp)&&nlp!=nptail)
	stop(paste("itail should have",nlp,"initial estimates"))
if(inherits(tail,"formula")||is.function(tail)){
	if(otail){
		if(is.numeric(itail) && length(itail)!=nptail)
			stop(paste("itail must be of size ",nptail))
		else if(!is.numeric(itail)) itail <- rep(0,nptail)
		fntail <- tail1}
	else {
		if(!is.numeric(itail))
			stop("Missing initial conditions for tail")
		else if(length(itail)!=nptail)
			stop(paste("itail must be of size ",nptail))
		else fntail <- function(p) tail1(itail)}}
else if(otail){
	fntail <- function(p) tail_h(rep(p[npl4],n))
	if(!is.numeric(itail)){
		itail <- tail[1]
		nptail <- 1}}
else {
	if(length(tail)==n)fntail <- function(p) tail_h(tail)
	else fntail <- function(p) rep(tail_h(tail[1]),n) # IMPORTANT
	nptail <- 1}

# Computation of -log-likelihood
ly <- length(y)
z0 <- rep(0,ly)
llikstable <- function(p){  # ,up=up,integration=integration) {
	loc <- fnloc(p)
	disp <- fndisp(p)
	skew <- fnskew(p)
	tail <- fntail(p)

	if(!censor){
		if(exact){
			tamp <- .C("pstable_c",
				as.integer(ly),
				yy=as.double((y+delta/2-loc)/disp),
				skew=as.double(skew+z0),
				tail=as.double(tail+z0),
				as.double(eps),
				err=integer(1),
				ffy=double(ly))$ffy-
				.C("pstable_c",
				as.integer(ly),
				yy=as.double((y-delta/2-loc)/disp),
				skew=as.double(skew+z0),
				tail=as.double(tail+z0),
				as.double(eps),
				err=integer(1),
				ffy=double(ly),
				PACKAGE="stable")$ffy
			llikcomp <- -(log(tamp))*weights}
		else {
			tamp <- .C("stable",
				as.integer(ly),
				yy=as.double((y-loc)/disp),
				skew=as.double(skew+z0),
				tail=as.double(tail+z0),
				as.integer(npoint),
				as.double(up),
				as.double(eps),
				as.integer(if(integration=="Simpson")1 else 2),
				err=integer(1),
				ffy=double(ly),
				PACKAGE="stable")$ffy/disp
			llikcomp <- -(log(tamp)+log(delta))*weights}}
	else {
		if(exact){
			p1 <- .C("pstable_c",
				as.integer(ly),
				yy=as.double((y+delta/2-loc)/disp),
				skew=as.double(skew+z0),
				tail=as.double(tail+z0),
				as.double(eps),
				err=integer(1),
				ffy=double(ly),
				PACKAGE="stable")$ffy
			ps <- .C("pstable_c",
				as.integer(ly),
				yy=as.double((y-delta/2-loc)/disp),
				skew=as.double(skew+z0),
				tail=as.double(tail+z0),
				as.double(eps),
				err=integer(1),
				ffy=double(ly),
				PACKAGE="stable")$ffy
			llikcomp <- -weights*(cc*log(p1-ps)+log(lc-rc*ps))}
		else {
			tamp <- .C("stable",
				as.integer(ly),
				yy=as.double((y-loc)/disp),
				skew=as.double(skew+z0),
				tail=as.double(tail+z0),
				as.integer(npoint),
				as.double(up),
				as.double(eps),
				as.integer(if(integration=="Simpson")1 else 2),
				err=integer(1),
				ffy=double(ly),
				PACKAGE="stable")$ffy/disp
			llikcomp <- -weights*(cc*(log(tamp)+log(delta))
				+log(lc-rc*
				.C("pstable_c",
				as.integer(ly),
				yy=as.double((y-loc)/disp),
				skew=as.double(skew+z0),
				tail=as.double(tail+z0),
				as.double(eps),
				err=integer(1),
				ffy=double(ly),
				PACKAGE="stable")$ffy))}}
	llik <- sum(llikcomp)
	if(llik.output){
		if(length(p)==0)cat("-LogLik: ",sum(llikcomp),"\n")
		else cat("-LogLik: ",sum(llikcomp)," ",p,"\n")}
	z <- list(
		llik=llik,
		llikcomp=llikcomp,
		loc=loc,
		disp=disp,
		skew=skew,
		tail=tail)
	z}

# This function returns -llik (and NOT the deviance) to get the
# correct s.e's with the hessian returned by nlm (optimization routine).
optstable <- function(p){
	tamp <- llikstable(p)$llik
	if(llik.output) cat("-LogLik: ",tamp," (",p,")","\n")
	if(is.na(tamp))1e20 else tamp}

# initial conditions
p0 <- c()
if(oloc) {
	p0 <- c(iloc)
	names(p0) <- c(rep("iloc",length(iloc)))}
if(odisp) {
	tamp <- names(p0)
	p0 <- c(p0,idisp)
	names(p0) <- c(tamp,rep("idisp",length(idisp)))}
if(oskew) {
	tamp <- names(p0)
	p0 <- c(p0,iskew)
	names(p0) <- c(tamp,rep("iskew",length(iskew)))}
if(otail) {
	tamp <- names(p0)
	p0 <- c(p0,itail)
	names(p0) <- c(tamp,rep("itail",length(itail)))}
if(llik.output) {
	cat("No. of parameters: ",nploc,"",npdisp,"",npskew,"",nptail,"\n")
	if(oloc||odisp||oskew||otail){
		cat("Vector of initial conditions on IR^p:","\n")
		print(p0)}}

# -Log-likelihood at p0
llik0 <- llikstable(p=p0)
np0 <- length(p0)

if(np0>0){
	p.opt <- nlm(optstable, p=p0, hessian=hessian, fscale=fscale,
		typsize=rep(1,length(p0)), print.level=print.level,
		ndigit=ndigit, gradtol=gradtol, steptol=steptol,
		iterlim=iterlim, stepmax=stepmax)
	z <- llikstable(p.opt$estimate)}
else z <- llikstable(p0)

ytilde.tamp <- stable.mode(z$loc,z$disp,z$skew,z$tail)$ytilde
  # corresponding mode
tamp <- dstable(x=ytilde.tamp, loc=z$loc, disp=z$disp, skew=z$skew,
		tail=z$tail, npt=npoint, up=up, integration=integration)
llik.ytilde <- -(log(tamp)+log(delta))*weights

ytilde.tamp <- stable.mode(z$loc,z$disp,z$skew,z$tail)$ytilde

np <- nploc+npdisp+npskew+nptail
llik <- z$llik
aic <- llik+np0
nobs <- sum(as.numeric(weights))
df <- nobs-np0
loc <- as.vector(z$loc)
disp <- as.vector(z$disp)
skew <- as.vector(z$skew)
tail <- as.vector(z$tail)
ytilde <- stable.mode(loc,disp,skew,tail)$ytilde  # corresponding mode
residuals <- as.vector((y-loc)/disp)

if(!is.null(loc3))loc1 <- loc3
if(!is.null(disp3))disp1 <- disp3
if(!is.null(skew3))skew1 <- skew3
if(!is.null(tail3))tail1 <- tail3
if(np0>0){
	cov <- diag(np0)
	if(hessian){
		if(np0==1)cov <- 1/p.opt$hessian
		else {
			a <- if(any(is.na(p.opt$hessian))||any(abs(p.opt$hessian)==Inf))0
			else qr(p.opt$hessian)$rank
			if(a==np0)cov <- solve(p.opt$hessian)
			else cov <- matrix(NA,ncol=np0,nrow=np0)}}
	se <- if(hessian)sqrt(diag(cov)) else NA
	corr <- if(hessian)cov/(se%o%se) else NA
	coefficients <- p.opt$estimate
	gradient <- p.opt$gradient
	error <- p.opt$error
	code <- p.opt$code
	iterations <- p.opt$iter}
else iterations <- coefficients <- se <- cov <- corr <- gradient <-
	hessian <- error <- code <- NULL
z1 <- list(
	call=call,
	y=y,
	loc=loc,
	disp=disp,
	skew=skew,
	tail=tail,
	oloc=oloc,
	odisp=odisp,
	oskew=oskew,
	otail=otail,
	loc1=loc1,
	disp1=disp1,
	skew1=skew1,
	tail1=tail1,
	ytilde=ytilde,
	start=p0,
	weights=weights,
	llik=llik,
	llikcomp=z$llikcomp,
	residuals=residuals,
	aic=aic,
	npar=np,
	nploc=nploc,
	npdisp=npdisp,
	npskew=npskew,
	nptail=nptail,
	npar.opt=np0,
	n=nobs,
	df=df,
	coefficients=coefficients,
	se=se,
	cov=cov,
	corr=corr,
	gradient=gradient,
	hessian=hessian,
	iterations=iterations,
	error=error,
	code=code,
	integration=integration)

class(z1) <- "stable"
z1}

###########################################################################
# Stable class functions
#
fitted.stable <- function(z) z$loc

df.residual.stable <- function(z) z$df

aic <- function (object, ...)
	UseMethod("aic")
aic.stable <- function(z) z$aic

deviance.stable <- function(z) 2*z$llik

print.stable <- function(z,digits=max(4,.Options$digits-3),correlation=TRUE) {
np <- z$npar.opt
if(!z$oloc)nploc <- 0
if(!z$odisp)npdisp <- 0
if(!z$oskew)npskew <- 0
if(!z$otail)nptail <- 0
cat("\nCall:",deparse(z$call),sep="\n")
cat("\n")
if(!is.null(z$code)&&z$code>2)
	cat("Warning: no convergence - error",z$code,"\n\n")
cat("-Log likelihood             ",z$llik,"\n")
cat("No. of obs                  ",z$n,"\n")
cat("No. of estimated parameters ",np,"\n")
cat("No. of parameters           ",z$npar,"\n")
cat("Degrees of freedom          ",z$df,"\n")
cat("AIC                         ",z$aic,"\n")
if(!is.null(z$iterations))
	cat("Iterations                  ",z$iterations,"\n")

if(z$oloc&&z$nploc>0){
	cat("\nLocation parameters\n")
	if(!is.null(attr(z$loc1,"formula")))
		cat(deparse(attr(z$loc1,"formula")),sep="\n")
	else if(!is.null(attr(z$loc1,"model"))){
		t <- deparse(attr(z$loc1,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	cname <- if(is.character(attr(z$loc1,"model")))
		attr(z$loc1,"model")
		else if(!is.null(attr(z$loc1,"parameters")))
		attr(z$loc1,"parameters")
		else paste(1:z$nploc)
	coef.table <- cbind(z$coef[1:z$nploc],z$se[1:z$nploc])
	colname <- c("estimate","se")
	dimnames(coef.table) <- list(cname,colname)
	print.default(coef.table, digits=digits, print.gap=2)}
if(z$odisp&&z$npdisp>0){
	cat("\nDispersion parameters\n")
	if(!is.null(attr(z$disp,"formula")))
		cat(deparse(attr(z$disp,"formula")),sep="\n")
	else if(!is.null(attr(z$disp,"model"))){
		t <- deparse(attr(z$disp,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	cname <- if(is.character(attr(z$disp,"model")))
		attr(z$disp,"model")
		else if(!is.null(attr(z$disp,"parameters")))
		attr(z$disp,"parameters")
		else paste(1:z$npdisp)
	coef.table <- cbind(z$coef[(z$nploc*z$oloc+1):(z$nploc*z$oloc+z$npdisp)],z$se[(z$nploc*z$oloc+1):(z$nploc*z$oloc+z$npdisp)])
	colname <- c("estimate","se")
	dimnames(coef.table) <- list(cname,colname)
	print.default(coef.table, digits=digits, print.gap=2)}
if(z$oskew&&z$npskew>0){
	cat("\nSkew parameters\n")
	if(!is.null(attr(z$skew1,"formula")))
		cat(deparse(attr(z$skew1,"formula")),sep="\n")
	else if(!is.null(attr(z$skew1,"model"))){
		t <- deparse(attr(z$skew1,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	cname <- if(is.character(attr(z$skew1,"model")))
		attr(z$skew1,"model")
		else if(!is.null(attr(z$skew1,"parameters")))
		attr(z$skew1,"parameters")
		else paste(1:z$npskew)
	coef.table <- cbind(z$coef[(z$nploc*z$oloc+z$npdisp*z$odisp+1):(z$nploc*z$oloc+z$npdisp*z$odisp+z$npskew)],z$se[(z$nploc*z$oloc+z$npdisp*z$odisp+1):(z$nploc*z$oloc+z$npdisp*z$odisp+z$npskew)])
	colname <- c("estimate","se")
	dimnames(coef.table) <- list(cname,colname)
	print.default(coef.table, digits=digits, print.gap=2)}
if(z$otail&&z$nptail>0){
	cat("\nTail parameters\n")
	if(!is.null(attr(z$tail1,"formula")))
		cat(deparse(attr(z$tail1,"formula")),sep="\n")
	else if(!is.null(attr(z$tail1,"model"))){
		t <- deparse(attr(z$tail1,"model"))
		t[1] <- sub("expression\\(","",t[1])
		t[length(t)] <- sub("\\)$","",t[length(t)])
		cat(t,sep="\n")}
	cname <- if(is.character(attr(z$tail1,"model")))
		attr(z$tail1,"model")
		else if(!is.null(attr(z$tail1,"parameters")))
		attr(z$tail1,"parameters")
		else paste(1:z$nptail)
	coef.table <- cbind(z$coef[(z$nploc*z$oloc+z$npdisp*z$odisp+z$npskew*z$oskew+1):(z$nploc*z$oloc+z$npdisp*z$odisp+z$npskew*z$oskew+z$nptail)],z$se[(z$nploc*z$oloc+z$npdisp*z$odisp+z$npskew*z$oskew+1):(z$nploc*z$oloc+z$npdisp*z$odisp+z$npskew*z$oskew+z$nptail)])
	colname <- c("estimate","se")
	dimnames(coef.table) <- list(cname,colname)
	print.default(coef.table, digits=digits, print.gap=2)}
if(!is.null(z$hessian)&&z$hessian&&np>1&&correlation){
	cat("\nCorrelations:\n")
	dimnames(z$corr) <- list(seq(1,np),seq(1,np))
	print.default(z$corr, digits=digits)}}

###########################################################################
# Computation of the mode ytilde as a function of (loc,disp,skew,tail)



#' Mode of a Stable Distribution
#'
#' This function gives a reliable approximation to the mode of a stable
#' distribution with location, dispersion, skewness and tail thickness
#' specified by the parameters \code{loc}, \code{disp}, \code{skew} and
#' \code{tail}. \code{tail} must be in (1,2).
#'
#' \code{loc} is a location parameter in the same way as the mean in the normal
#' distribution: it can take any real value.
#'
#' \code{disp} is a dispersion parameter in the same way as the standard
#' deviation in the normal distribution: it can take any positive value.
#'
#' \code{skew} is a skewness parameter: it can take any value in \eqn{(-1,1)}.
#' The distribution is right-skewed, symmetric and left-skewed when \code{skew}
#' is negative, null or positive respectively.
#'
#' \code{tail} is a tail parameter (often named the characteristic exponent):
#' it can take any value in \eqn{(0,2)} (with \code{tail=1} and \code{tail=2}
#' yielding the Cauchy and the normal distributions respectively when symmetry
#' holds).
#'
#' The simplest empirical formula found to give a satisfactory approximation to
#' the mode for values of \code{tail} in \eqn{(1,2)} is \deqn{
#' loc+disp*a*skew*exp(-b*abs(skew))} with \deqn{ a =
#' 1.7665114+1.8417675*tail-2.2954390*tail^2+0.4666749*tail^3} and \deqn{ b =
#' -0.003142967+632.4715*tail*exp(-7.106035*tail)}.
#'
#'
#' @param loc vector of (real) location parameters.
#' @param disp vector of (positive) dispersion parameters.
#' @param skew vector of skewness parameters (in [-1,1]).
#' @param tail vector of parameters (in [1,2]) related to the tail thickness.
#' @return A list of size 3 giving the mode, \eqn{a} and \eqn{b}.
#' @author Philippe Lambert (Catholic University of Louvain, Belgium,
#' \email{phlambert@@stat.ucl.ac.be}) and Jim Lindsey.
#' @seealso \code{stable} for more details on the stable
#' distribution.
#'
#' \code{stablereg} to fit generalized linear models for the
#' stable distribution parameters.
#' @references Lambert, P. and Lindsey, J.K. (1999) Analysing financial returns
#' using regression models based on non-symmetric stable distributions. Applied
#' Statistics, 48, 409-424.
#' @keywords distribution
#' @examples
#'
#' x <- seq(-5,5,by=0.1)
#' plot(x,dstable(x,loc=0,disp=1,skew=-1,tail=1.5),type="l",ylab="f(x)")
#' xhat <- stable.mode(loc=0,disp=1,skew=-1,tail=1.5)$ytilde
#' fxhat <- dstable(xhat,loc=0,disp=1,skew=-1,tail=1.5)
#' lines(c(xhat,xhat),c(0,fxhat),lty="dotted")
#'
#' @export stable.mode
stable.mode <- function(loc, disp, skew, tail){
	if(any(tail<1))
		warning("stable.mode is only reliable for tail in (1,2)")
	coef1 <- 1.7665114+1.8417675*tail-2.2954390*tail^2+0.4666749*tail^3
	coef2 <- -0.003142967+632.4715*exp(-7.106035*tail)*tail
	list(ytilde=loc+disp*coef1*exp(-coef2*abs(skew))*skew,
		c1=coef1,c2=coef2)}
