#' Easy conversion of parameters between stabledist (Nolan 1-parameterization) and stable (Lambert and Lindsey 1999)
#'
#' \code{sd2s} has stabledist parameter (Nolan 1-parameterization) inputs and returns stable parameters as put forth in Lambert and Lindsey (1999) and used in this package.
#' \code{s2sd} has stable     parameter (Lambert and Lindsey (1999)) inputs and returns stabledist parameters (Nolan 1-parameterization).
#' See examples and the readme.  There's also more context and references in `?stable::dstable`.
#' 
#' 
#' [Swihart 2022 update:] See the examples and README for how to make equivalent calls
#' to those of 'stabledist' (i.e., Nolan's 1-parameterization 
#' as detailed in Nolan (2020)) using these functions and this package. 
#' See github for Lambert and Lindsey 1999 JRSS-C journal article, 
#' which details the parameterization of the Buck (1995) stable distribution which allowed
#' a Fourier inversion to arrive at a form of the $g_d$ function as detailed in Nolan (2020),
#' The Buck (1995) parameterization most closely resembles the Zolotarev B parameterization
#' outlined in Definition 3.6 on page 93 of Nolan (2020) -- except that Buck (1995) did
#' not allow the scale parameter to multiply with the location parameter.  
#' This explains why the `Zolotarev B` entry in Table 3.1 on page 97 of Nolan (2020) has
#' the location parameter being multiplied by the scale parameter whereas in converting the Lindsey and Lambert (1999)
#' to Nolan 1-parameterization the location parameter stays the same.  
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
#' @references Lambert, P. and Lindsey, J.K. (1999) Analysing financial returns using
#' regression models based on non-symmetric stable distributions. Applied
#' Statistics, 48, 409-424.
#' 
#' @references Nolan, John P. Univariate stable distributions. Berlin/Heidelberg, Germany: Springer, 2020.
#' 
#' @return 
#' \code{sd2s} returns stable parameters as put forth in Lambert and Lindsey (1999) and used in this package.
#' \code{s2sd} returns stabledist parameters (Nolan 1-parameterization).
#' @export
#' @examples
#' \dontrun{
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
#' stabledist::pstable(q, alpha=sd$alpha, beta=sd$beta , gamma=sd$gamma , delta=sd$delta, pm=1)}
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
