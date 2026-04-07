#' balance cohorts (e.g. on ancestry or sex or whatever) 
#' 
#' @param    x      a named vector, c(SUBJECT=FACTOR, ...) (source cohort)
#' @param    y      a target cohort, c(SUBJECT=FACTOR, ...) (target cohort) 
#' @param    size   size of the sample of y? (if NULL, then length(y))
#' @param    exact  use exact sized subsamples instead of (fast) probs? (FALSE)
#'
#' @return   samples of x such that FACTOR is distributed similarly to y
#'
#' @details  Assumes length(x) > length(y). Use set.seed to replicate. If any 
#'           group in y is larger than the corresponding group in x, then 
#'           sampling will be done with replacement, otherwise not.
#'
#' @examples
#'  \dontrun{
#'    joint <- with_ancestry[, c('cohort','eaf.max'), drop=FALSE]
#'    byCohort <- with(joint, split(setNames(eaf.max, rownames(joint)), cohort))
#'
#'    controlApprox <- lapply(byCohort[-1], balance, x=byCohort[[1]])
#'    sapply(controlApprox, function(x) table(attr(x, 'grouping')))
#'
#'    controlExact <- lapply(byCohort[-1], balance, x=byCohort[[1]], exact=TRUE)
#'    sapply(controlExact, function(x) table(attr(x, 'grouping')))
#'  }
#' @export 
#'
balance <- function(x, y, size=NULL, exact=FALSE) { 

  if (length(x) <= length(y)) warning("x should be larger than y!")
  gx <- table(x) 
  gy <- table(y) 
  stopifnot(all(names(gy) %in% names(gx)))
  ixy <- x[x %in% names(gy)] 
  gxy <- table(ixy)[names(gy)]
  repl <- any(gy > gxy)
  if (is.null(size)) size <- length(y)
  
  if (exact) { 
    res <- character()
    for (i in names(gy)) {
      res <- c(res, sample(names(ixy[ixy==i]), size=gy[i], replace=gy[i]>gx[i]))
    }
  } else { 
    pry <- gy/sum(gy)
    prxy <- gxy/sum(gxy)
    adj <- setNames(as.numeric(pry/prxy), names(gy))
    prob <- setNames(adj[ixy], names(ixy))/sum(adj[ixy])
    res <- sample(names(prob), size=size, replace=repl, prob=prob)
  }
  attr(res, 'grouping') <- ixy[res]
  return(res)

}
