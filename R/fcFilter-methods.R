setMethod("show", sig = c("fcFilter"),
          definition = function(object) {
  cat("fcFilter:")
  callNextMethod(object)
})

#' get posteriors from a \code{fcFilter} object
#' 
#' @param x \code{fcFilter}
#' @param y \code{character} or \code{missing} that specifiy which channel to look for
#' 
#' @export 
setGeneric("posteriors", function(x, y, ...) standardGeneric("posteriors"))
setMethod("posteriors", sig = c("fcFilter", "ANY"),
          definition = function(x, y = "missing") {
  x@posteriors
})

setMethod("posteriors", sig = c("fcFilter", "character"),
          definition = function(x, y) {
  post <- posteriors(x)
  
  postNames <- names(post)
  ind <- match(y, postNames)
  if (is.na(ind)) {
    stop("FcFilter not found for:", y)
  } else {
    post[[ind]]
  }
})

#' get priors from a \code{fcFilter} object
#' 
#' @inheritParams posteriors
#' @export
setGeneric("priors", function(x, y, ...) standardGeneric("priors"))
setMethod("priors", sig = c("fcFilter", "ANY"),
          definition = function(x, y = "missing") {
  x@priors
})

setMethod("priors", sig = c("fcFilter", "character"),
          definition = function(x, y) {
  prior <- priors(x)
  
  priorNames <- names(prior)
  ind <- match(y, priorNames)
  if (is.na(ind)) {
    stop("FcFilter not found for:", y)
  } else {
    prior[[ind]]
  }
}) 
