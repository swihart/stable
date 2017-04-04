#' Easy conversion of parameters between stabledist and stable
#'
#' \code{pm0_to_pm1} has stabledist parameter inputs for pm=0 and returns pm=1 equivalent parameterization.
#' \code{pm1_to_pm0} has stabledist parameter inputs for pm=1 and returns pm=0 equivalent parameterization.
#' 
#'
#' @name Parameter_Conversion_Nolan_pm1_pm0
#' @aliases pm1_to_pm0 pm0_to_pm1
#' @param a0 the stabledist 'alpha' for pm=0 in 'stabledist'
#' @param b0 the stabledist 'beta' for pm=0 in 'stabledist'
#' @param c0 the stabledist 'gamma' for pm=0 in 'stabledist'
#' @param d0 the stabledist 'delta' for pm=0 in 'stabledist'
#' @param a1 the stabledist 'alpha' for pm=1 in 'stabledist'
#' @param b1 the stabledist 'beta' for pm=1 in 'stabledist'
#' @param c1 the stabledist 'gamma' for pm=1 in 'stabledist'
#' @param d1 the stabledist 'delta' for pm=1 in 'stabledist'
#' 
#' @return What you need.  See examples.  
#' @export
#' @examples
#' q <- -1
#' # nolan pm=1 parameters:
#' a1 <-  1.3
#' b1 <-  -0.4
#' c1 <-  2
#' d1 <-  0.75
#' # Convert to nolan pm=0 parameters:
#' pm0 <- pm1_to_pm0(a1,b1,c1,d1)
#' a0 <- pm0$a0
#' b0 <- pm0$b0
#' c0 <- pm0$c0
#' d0 <- pm0$d0
#' # check:
#' stabledist::pstable(q, alpha=a1, beta=b1 , gamma=c1 , delta=d1, pm=1)
#' #> [1] 0.1965513
#' # only change delta=d0 for pm=0
#' stabledist::pstable(q, alpha=a1, beta=b1 , gamma=c1 , delta=d0, pm=0)
#' stabledist::pstable(q, alpha=a0, beta=b0 , gamma=c0 , delta=d0, pm=0)
#' #> [1] 0.1965513
#' stabledist::dstable(q, alpha=a1, beta=b1 , gamma=c1 , delta=d1, pm=1)
#' #> [1] 0.0572133
#' # only change delta=d0 for pm=0
#' stabledist::dstable(q, alpha=a1, beta=b1 , gamma=c1 , delta=d0, pm=0)
#' stabledist::dstable(q, alpha=a0, beta=b0 , gamma=c0 , delta=d0, pm=0)
#' #> [1] 0.0572133
pm0_to_pm1 <- function(a0, b0, c0, d0){
  
  if(a0!=1){  
    a1 <- a1
    b1 <- b1
    c1 <- c1
    d1 <- d0 - b1*c1*tan(pi*a1/2)
      }
  
  return(list(a1=a1, b1=b1, c1=c1, d1=d1 ))
  
}
#' @rdname Parameter_Conversion_Nolan_pm1_pm0
#' @export pm1_to_pm0
pm1_to_pm0 <- function(a1, b1, c1, d1){
  
  if(a1!=1){  
    a0 <- a1
    b0 <- b1
    c0 <- c1
    d0 <- d1 + b1*c1*tan(pi*a1/2)
      }
  
  return(list(a0=a0, b0=b0, c0=c0, d0=d0))
  
}
