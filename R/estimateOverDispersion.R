#' Estimate Overdispersion
#' @param obs    Observed number of mutations per region
#' @param exp    Expected number of mutations per region
#' @param length Total length of bases in region \emph{accumulated} across samples
#'
#' @export
estimateOverDispersion <- function(obs, exp, length){
  N <- length(obs)

  # Probability
  prob <- exp / length

  # Function to be optimized over
  loglikm <- function(logsd){
    sd <- exp(logsd)
    res <- 0
    for(i in 1:N){
      res <- res + VGAM::dbetabinom(obs[i],
                                    length[i],
                                    prob[i],
                                    sd^2*prob[i]/(1-prob[i]), log = T)
    }
    -res
  }

  # Optim
  suppressWarnings({
    opt <- optim(par = log(0.335), fn = loglikm, hessian = T, method = "BFGS")
  })

  # Return o.d. value and confidence intervals
  c("Estimate" = exp(opt$par),
    "CI_5%" = exp(opt$par-1.96/sqrt(opt$hessian)),
    "CI_95%" = exp(opt$par+1.96/sqrt(opt$hessian)))
}
