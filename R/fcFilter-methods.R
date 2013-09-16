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
#' @aliases posteriors,fcFilter,ANY-method 
setGeneric("posteriors", function(x, y, ...) standardGeneric("posteriors"))
setMethod("posteriors", sig = c("fcFilter", "ANY"),
          definition = function(x, y = "missing") {
  x@posteriors
})
#' @aliases posteriors,fcFilter,character-method
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
#' @param x \code{fcFilter} object
#' @param y \code{character} specifying channel name. if \code{missing} then extract priors for all the channels 
#'
#' @inheritParams posteriors
#' @export
#' @aliases priors,fcFilter,ANY-method
setGeneric("priors", function(x, y, ...) standardGeneric("priors"))
setMethod("priors", sig = c("fcFilter", "ANY"),
          definition = function(x, y = "missing") {
  x@priors
})
#' @aliases priors,fcFilter,character-method
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
