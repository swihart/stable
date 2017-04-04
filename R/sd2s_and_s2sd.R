#' Easy conversion of parameters between stabledist and stable
#'
#' \code{sd2s} has stabledist parameter inputs and returns stable parameters.
#' \code{s2sd} has stable     parameter inputs and returns stabledist parameters.
#' 
#' This is a generic function: methods can be defined for it directly
#' or via the \code{\link{Summary}} group generic. For this to work properly,
#' the arguments \code{...} should be unnamed, and dispatch is on the
#' first argument.
#'
#' @name Parameter_Conversion
#' @aliases sd2s s2sd
#' @param alpha the stabledist 'alpha'
#' @param beta the stabledist 'beta'
#' @param gamma the stabledist 'gamma'
#' @param delta the stabledist 'delta'
#' @param pm    default 1; currently only value supported. the stabledist parameterization 'pm'
#' @param tail  the stable 'tail' analogous to 'alpha'
#' @param skew  the stable 'skew' analogous to 'beta'
#' @param disp  the stable 'disp' analogous to 'gamma'
#' @param loc  the stable 'loc' analogous to 'delta' 
#' 
#' @return What you need.  See examples.  
#' @export
#' @examples
#' q <- -1
#' # nolan pm=1 parameters:
#' a <-  1.3
#' b <-  -0.4
#' c <-  2
#' d <-  0.75
#' s <- sd2s(alpha=a, beta=b, gamma=c, delta=d)
#' stable::pstable(q, tail = s$tail, skew=s$skew, disp = s$disp, loc  = s$loc)
#' stabledist::pstable(q, alpha=a, beta=b , gamma=c , delta=d, pm=1)
#' sd <- s2sd(tail = s$tail, skew=s$skew, disp = s$disp, loc  = s$loc)
#' stabledist::pstable(q, alpha=sd$alpha, beta=sd$beta , gamma=sd$gamma , delta=sd$delta, pm=1)
sd2s <- function(alpha, beta, gamma, delta, pm=1){
  
  if(pm==1){  
    a=alpha
    b=beta
    c=gamma
    d=delta
    
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
  }
  
  return(list(tail=a3, skew=b3, disp=c3, loc=d3 ))
  
}
#' @rdname Parameter_Conversion
#' @export s2sd
s2sd <- function(tail, skew, disp, loc, pm=1){
  
  if(pm==1){  
    #a=alpha
    #b=beta
    #c=gamma
    #d=delta
    
    # lindsey-(3) page 415 conversion:
    # tail/alpha and location stay the same
    #a3 <- a
    #d3 <- d 
    ## the others require calcs:
    #DEL2 <- cos(pi/2 * a)^2 + (-b)^2*sin(pi/2 * a)^2
    #DEL <- sqrt(DEL2) * sign(1-a)
    #eta_a <- min(a, 2-a)
    ## the lindsey-(3) beta:
    #b3 <- -sign(b)*2/(pi*eta_a)*acos( cos(pi/2 * a) / DEL )
    ## the lindsey-(3) scale:
    #c3 <- ( (DEL*c^a) / cos(pi/2 * a) )^(1/a)
    
    ## reverse!
    a3 <- tail
    b3 <- skew
    c3 <- disp
    d3 <- loc
    ## solve the above for b, c
    eta_a <- min(a3, 2-a3)
    DEL <- cos( pi * a3/2) / cos(b3 * pi * eta_a/2)
    DEL2 <- DEL^2
    b <- -sign(b3) * sqrt( (DEL2 - cos(pi*a3/2)^2) / sin(pi*a3/2)^2 )
    c <- c3 * ( cos(pi*a3/2) / DEL )^(1/a3)
  }
  
  return(list(alpha=a3, beta=b, gamma=c, delta=d3))
  
}
