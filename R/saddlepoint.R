#' Saddlepoint Approximation
#'
#' Let \eqn{X_1,X_2,\ldots X_N} be independent Bernouilli trials
#' with \eqn{X_i \sim p_i}.
#' We evaluate the probability
#' \deqn{P(s_1X_1+s_2X_2\ldots s_NX_N > t)}.
#' 
#' @param t   Point where tail probability is evaluated (note this function is not vectorized over t)
#' @param dat A data.table with columns 'x', 'y', and 'probabilities'. Column 'x' is a numeric indicator unique for each random variable to be convoluted together. Column 'y' is the outcome value, and colum 'probabilities' contains the probability corresponding to each value of the 'y' column.
#' @param lattice Lattice size in minimal lattice
#' @param log Return log probability
#' 
#' @examples
#' saddlepoint(300, dat = data.table(x = rep(1:10000, each = 2),
#'                                  y = rep(c(1,0),10000),
#'                                  probability = rep(c(0.01,0.99), 100)))
#' pbinom(299, 10000, 0.01, lower.tail = F, log.p = T)
#' @importFrom Rcpp evalCpp
#' @import data.table
#' @export
saddlepoint <- function(t, dat, lattice = 1L, log = T){
  # Sort data frame after x
  dat <- dat[order(x),]
  
  x <- dat$x
  p <- dat$probability
  s <- dat$y
  
  stopifnot(all(s >= 0))
  
  # Return 0 and 1 if outside range
  min_score <- dat[,min(y), by = x][,sum(V1)]
  if(t == min_score)
    return(ifelse(log, 0, 1))
  max_score <- dat[,max(y), by = x][,sum(V1)]
  if(t > max_score)
    return(ifelse(log, -Inf, 0))
  
  # Compute expected score
  expected_score <- dat[,sum(y*probability)]
  lower_than_expected <- t < expected_score
  t <- ifelse(lower_than_expected, t - lattice, t)
  if(lower_than_expected && (t <= min_score + lattice) ){
    # Case not well-handled by saddlepoint approximation
    # 1 - probability that all variables takes minimum value
    ret <- prod(dat[, .("min_y" = min(y), probability, y), x][y == min_y, .(probability)])
    return( ifelse(log, log1p(-ret), 1-ret) )
  }
  # Giver saddelpunkts approksimationen paa log-skala
  # Use a combination of bisection and Newton raphson
  
  theta_old <- -Inf
  theta <- 0
  
  # 
  theta_largest_non_na <- 0
  theta_smallest_non_na <- 0
  theta_largest_negative <- -Inf # Largest theta that gives a negative value so far
  theta_smallest_positive <- Inf # Smallest theta that gives a positive value so far
  
  iter <- 0
  while(abs(theta - theta_old) > 1e-8){
    cumDer <- cumulantDerivatives(theta, x, p, s)
    
    # Escape too large jumps in Newton-Raphson
    if(any(is.na(cumDer))){
      if(theta > theta_largest_non_na){
        theta <- (theta+theta_largest_non_na) / 4
        next
      } else if (theta < theta_smallest_non_na){
        theta <- (theta+theta_smallest_non_na) / 4
        next
      }
    } else {
      theta_largest_non_na <- max(theta, theta_largest_non_na)
      theta_smallest_non_na <- min(theta, theta_smallest_non_na)
    }
    
    # For convergence properties
    iter <- iter + 1
    theta_old <- theta
    
    # Newton raphson
    val <- cumDer[2] - t
    if(val < 0 && theta > theta_largest_negative)
      theta_largest_negative <- theta
    if(val > 0 && theta < theta_smallest_positive)
      theta_smallest_positive <- theta
    deriv <- cumDer[3]
    step <- - val / deriv 
    
    # Next theta if newton raphson
    theta <- theta + step
    if(theta < theta_largest_negative || theta > theta_smallest_positive){
      # Fall back to bisection
      theta <- (theta_largest_negative+theta_smallest_positive)/2
    }
  }
  
  cumDer <- cumulantDerivatives(theta, x, p, s)
  
  # El approximazione
  v <- cumDer[3]
  
  # Uden lattice-korrektion
  if(lattice == 0){
    
    ret <- ifelse(lower_than_expected,
                  log( - expm1(cumDer[1]-t*theta+theta^2*v/2+pnorm(theta*sqrt(v), lower.tail = T, log.p = T)) ),
                  cumDer[1]-t*theta+theta^2*v/2+pnorm(theta*sqrt(v), lower.tail = F, log.p = T)
    )
  } else{
    ret <- ifelse(lower_than_expected,
                  log( -expm1(cumDer[1]-t*theta+theta^2*v/2+pnorm(theta*sqrt(v), lower.tail = T, log.p = T) +
                             log(abs(theta*lattice))- log((1-exp(-lattice*abs(theta))))) ),
                  cumDer[1]-t*theta+theta^2*v/2+pnorm(theta*sqrt(v), lower.tail = F, log.p = T) +
                    log(abs(theta*lattice))- log((1-exp(-lattice*abs(theta))))
    )
  }
  
  ifelse(log, ret, exp(ret))
}
