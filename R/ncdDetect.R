#' Logit function
logit <- function(p){
  log(p) - log1p(-p)
}

#' Logistic function
logistic <- function(x){
  1 / (1 + exp(-x))
}

#' ncdDetect
#'
#' Find regions with an excessive mutation load
#' indicating positive selection
#'
#' @param od The standard deviation of a random effect in logit space. The random effect is not added to the first column. Can be estimated using \code{\link{estimateOverDispersion}}
#'
#' @export
ncdDetect <- function(predictions, scores, observed_score, od, N = 11, log = T){
  ##################################################
  ## Data Handling
  ##################################################

  ## Make initial checks on input --------------------------------------------
  ## - data dimensions
  if (!identical(nrow(predictions), nrow(scores))) {
    stop("the row numbers of predictions and scores must be identical")
  }

  ## - negative scores (cannot be handled)
  if (any(scores < 0)) {
    stop("ncdDetect can currently not handle negative scores")
  }

  ## row sums of predictions matrix
  if (unique(round(rowSums(predictions), 10)) != 1) {
    stop("the row sums of the prediction matrix must sum to one")
  }

  ## Convert data into data.tables -------------------------------------------
  predictions <- as.data.table(predictions)
  scores <- data.table(scores)
  ## Add x-values to data ----------------------------------------------------
  ## (needed in the convolution step in order to know which positions go together)
  scores[, x := 1:.N]

  ## Get scores in long format -------------------------------------------------
  scores_long <- melt(scores, id.vars = c("x"), measure.vars = setdiff(names(scores), "x"),
                      variable.name = "mutation_type", value.name = "y", variable.factor = F)
  ## - set key columns for merging
  setkeyv(scores_long, c("x", "mutation_type"))

  ##################################################
  ## Setup numeric integration
  ##################################################

  val <- seq(-3*od, 3*od, length.out = N)
  if(N == 1){val = 0}
  m <- diff(pnorm(c(-Inf,val,Inf), sd = od)) # N+1
  weight <- (m[1:N]+m[2:(N+1)])/2
  weight[1] <- weight[1] + m[1]/2 # Make constant approximation for remainder of integrand
  weight[N] <- weight[N] + m[N+1]/2 # Same

  ##################################################
  ## Run saddlepoint approximation
  ##################################################
  log_p_values <- list()
  for(i in 1:N){
    random_effect <- val[i]
    pred_cols <- names(predictions)[-1]
    pred_od <- data.table(logistic(logit(as.matrix(predictions[, pred_cols, with = F])) + random_effect))
    ## - add column with probability of no mutation
    no_mut_prob <- data.table(NO = 1-rowSums(pred_od))
    setnames(no_mut_prob, names(predictions)[1])
    pred_od <- cbind(no_mut_prob, pred_od)

    ## - add x-values to predictions
    pred_od[, x := 1:.N]

    ## - get predictions in long format
    predictions_long <- melt(pred_od, id.vars = c("x"), measure.vars = setdiff(names(pred_od), "x"),
                             variable.name = "mutation_type", value.name = "probability", variable.factor = F)


    ## - set key columns for merging
    setkeyv(predictions_long, c("x", "mutation_type"))

    ## - merge predictions and scores; set key columns
    dat <- predictions_long[scores_long]
    setkeyv(dat, c("x", "mutation_type"))

    ## - perform significance evaluation using saddlepoint approximation
    log_p_values[[i]] <- saddlepoint(t = observed_score, dat, log = T)

    ## - clean up
    rm(random_effect, pred_cols, pred_od, no_mut_prob, predictions_long, dat)
  }


  log_prod <- log(weight) + unlist(log_p_values)
  final_p <- log(sum(exp(log_prod - max(log_prod)))) + max(log_prod)

  ## Return result
  return(p = final_p)
}
