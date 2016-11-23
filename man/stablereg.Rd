\name{stablereg}
\title{Stable Generalized Regression Models}
\usage{
stablereg(y=NULL, loc=0, disp=1, skew=0, tail=1.5, 
	oloc=TRUE, odisp=TRUE, oskew=TRUE, otail=TRUE, noopt=FALSE,
	iloc=NULL, idisp=NULL,iskew=NULL, itail=NULL,
	loc.h=NULL, disp.h=NULL, skew.h=NULL, tail.h=NULL,
	weights=1, exact=FALSE, delta=1, envir=parent.frame(),
	integration="Romberg", eps=1e-6, up=10, npoint=501,
	hessian=TRUE, llik.output=FALSE, print.level=0, ndigit=10,
	steptol=0.00001, gradtol=0.00001, fscale=1, typsize=abs(p0),
	stepmax=sqrt(p0\%*\%p0), iterlim=100)

fitted(obj.stable)
df.residual(obj.stable)
deviance(obj.stable)
aic(obj.stable)
}
\alias{stablereg}
\alias{fitted.stable}
\alias{df.residual.stable}
\alias{deviance.stable}
\alias{aic.stable}
\description{
\code{stablereg} fits user specified generalized linear and nonlinear
regression models based on the stable distribution to (uncensored,
right and/or left censored) data. This allows the location, the
dispersion, the skewness and the tails of the fitted stable
distribution to vary with explanatory variables.
}
\section{Warning}{Because of the numerical integrations involved, convergence
can be very sensitive to the initial parameter values supplied and to
the settings of the arguments controlling \code{\link{nlm}}. If nlm
feeds extreme parameter values in the tails of the distribution to the
likelihood function, the integration may hang for a long time.}
\arguments{
\item{y}{The response vector or a \code{repeated} data object. If the
\code{repeated} data object contains more than one response variable,
give that object in \code{envir} and give the name of the response
variable to be used here.

For censored data, two columns with the second being the censoring
indicator (1: uncensored, 0: right censored, -1: left censored.)}

\item{loc, loc.h, oloc, and iloc}{Describe the regression model
fitted for the location parameter of the stable distribution, perhaps
after transformation by the link function \code{loc.g} (set to the
identity by default. The inverse link function is denoted by \code{loc.h}.
Note that these functions cannot contain unknown parameters).

Two specifications are possible:

(1) \code{loc} is a linear or nonlinear language expression beginning
with ~ or an R function, describing the regression function for the
location parameter (after transformation by \code{loc.g}, the link function).

\code{iloc} is a vector of initial conditions for the parameters 
in the regression for this parameter.

\code{oloc} is a boolean indicating if
an optimization of the likelihood has to be carried out on these parameters.
If \code{oloc} is set to TRUE, a default zero value is considered for the
starting values \code{iloc}. But if no optimization is desired on the location
parameters, i.e. when the likelihood has to be evaluated or optimized
at a fixed location, then \code{iloc} has to be explicitely specified.

(2) \code{loc} is a numeric expression (i.e. a scalar or a vector of the same 
size as the data vector \code{y}, or \code{y[,1]} when censoring is
considered).

If \code{oloc} is set to TRUE, i.e. when an optimization of the likelihood has
to be carried out on the location parameter, then the location parameter
(after transformation by the link function loc.g) is set to an unknown
parameter with initial value equal to \code{iloc[1]} or \code{loc[1]}
when \code{iloc} is not specified.

But when \code{oloc} is set to FALSE, i.e. when the likelihood has to
be evaluated  or optimized at a fixed location, then the transformed
location is assumed to be equal to \code{loc} when it is of the same
length as the data vector \code{y} (or \code{y[,1]} when censoring is
considered), and to \code{loc[1]} otherwise.

Specification (1) is especially useful in ANOVA-like
situations where the location is assumed to change with the levels
of some factor variable.
}

\item{disp, disp.h, odisp, and idisp}{describe the regression model
for the dispersion parameter of the fitted stable
distribution, after transformation by the link function \code{disp.g}
(set to the \code{log} function by default).
The inverse link function is denoted by \code{disp.h}.
Again these functions cannot contain unknown parameters.
The same rules as above apply when specifying the generalized regression
model for the dispersion parameter.}

\item{skew, skew.h, oskew, and iskew}{describe the regression model
for the skewness parameter of the fitted stable
distribution, after transformation by the link function \code{skew.g}
(set to  \code{log{(1 + [.])/(1 - [.])}} by default). The inverse link
function is denoted by \code{skew.h}. Again these functions cannot
contain unknown parameters.
The same rules as above apply when specifying the generalized regression
model for the skewness parameter.}

\item{tail, tail.h, otail, and itail}{describe the regression model
considered for the tail parameter of the fitted stable distribution,
after transformation by the link function \code{tail.g} (set to 
\code{log{([.] - 1)/(2 - [.])}} by default. The inverse link function
is denoted by \code{tail.h}. Again these functions cannot contain
unknown parameters). The same rules as above apply when specifying
the generalized regression model for the tail parameter.}

\item{noopt}{When set to TRUE, it forces \code{oloc},
\code{odisp}, \code{oskew} and \code{otail} to FALSE, whatever the user choice
for these last three arguments. It is especially useful when looking 
for appropriate initial values for the regression model parameters, before
undertaking the optimization of the likelihood.}

\item{weights}{Weight vector.}

\item{exact}{If TRUE, fits the exact likelihood function for
continuous data by integration over intervals of observation,
i.e. interval censoring.}

\item{delta}{Scalar or vector giving the unit of measurement
for each response value, set to unity by default. For example,
if a response is measured to two decimals, \code{delta=0.01}.
If the response is transformed, this must be multiplied by
the Jacobian.  For example, with a log transformation,
\code{delta=1/y}. (The \code{delta} values for the censored response are
ignored.) The transformation cannot contain unknown
parameters.}

\item{envir}{Environment in which model formulae are to be
interpreted or a data object of class, \code{repeated}, \code{tccov},
or \code{tvcov}; the name of the response variable should be given in
\code{y}. If \code{y} has class \code{repeated}, it is used as the
environment.}

\item{integration, eps, up and npoint}{\code{integration} indicates
which algorithm must be used to evaluate the stable density when
the likelihood is computed with \code{exact} set to FALSE.
See the man page on \code{\link[stable]{stable}} for extra information.}

\item{llik.output}{is TRUE when the likelihood has to be displayed 
at each iteration of the optimization.}

\item{others}{Arguments controlling the optimization procedure
\code{\link{nlm}}.}
}
\value{
A list of class \code{stable} is returned.
The printed output includes the -log-likelihood, the corresponding AIC,
the maximum likelihood estimates, standard errors, and correlations.
It also include  all the relevant information calculated, including
error codes.
}
\references{
Lambert, P. and Lindsey, J.K. (1999) Analysing financial
returns using regression models based on non-symmetric stable
distributions. Applied Statistics 48, 409-424.
}
\author{
Philippe Lambert (Catholic University of Louvain, Belgium,
\email{phlambert@stat.ucl.ac.be}) and Jim Lindsey.
}
\seealso{
\code{\link{lm}}, \code{\link{glm}}, \code{\link[stable]{stable}} and 
\code{\link[stable]{stable.mode}}.
}
\examples{
## Share return over a 50 day period (see reference above)
# shares
y <- c(296,296,300,302,300,304,303,299,293,294,294,293,295,287,288,297,
305,307,307,304,303,304,304,309,309,309,307,306,304,300,296,301,298,
295,295,293,292,297,294,293,306,303,301,303,308,305,302,301,297,299)  

# returns
ret <- (y[2:50]-y[1:49])/y[1:49]
# hist(ret, breaks=seq(-0.035,0.045,0.01))

day <- seq(0,0.48,by=0.01) # time measured in days/100
x <- seq(1,length(ret))-1

# Classic stationary normal model tail=2
print(z1 <- stablereg(y = ret, delta = 1/y[1:49],
	loc = ~1, disp= ~1, skew = ~1, tail = tail.g(1.9999999),
	iloc = 0, idisp = -3, iskew = 0, oskew = FALSE, otail = FALSE))

# Normal model (tail=2) with dispersion=disp.h(b0+b1*day)
print(z2 <- stablereg(y = ret, delta = 1/y[1:49], loc = ~day,
	disp = ~1, skew = ~1, tail = tail.g(1.999999), iloc = c(0.003,0),
	idisp = -4.5, iskew = 0, oskew = FALSE, otail = FALSE))

# Stable model with loc(ation)=loc.h(b0+b1*day)
print(z3 <- stablereg(y = ret, delta = 1/y[1:49],
	loc = ~day, disp = ~1, skew = ~1, tail = ~1,
	iloc = c(0.001,-0.004), idisp = -4.8, iskew = 0, itail = 0.6))

# Stable model with disp(ersion)=disp.h(b0+b1*day)
print(z4 <- stablereg(y = ret, delta = 1/y[1:49],
	loc = ~1, disp = ~day, skew = ~1, tail = ~1,
	iloc = 0.003, idisp = c(-4.8,0), iskew = -0.03, itail = 1.6))

# Stable model with skew(ness)=skew.h(b0+b1*day)
# Evaluation at fixed parameter values (because noopt is set to TRUE)
print(z5 <- stablereg(y = ret, delta = 1/y[1:49],
	loc = ~1, disp = ~1, skew = ~day, tail = ~1,
	iloc = 5.557e-04, idisp = -4.957, iskew = c(2.811,-2.158),
	itail = 1.57, noopt=TRUE))

# Stable model with tail=tail.h(b0+b1*day)
print(z6 <- stablereg(y = ret, delta = 1/y[1:49], loc = ret ~ 1,
	disp = ~1, skew = ~1, tail = ~day, iloc = 0.002,
	idisp = -4.8, iskew = -2, itail = c(2.4,-4), hessian=FALSE))
}
\keyword{models}
