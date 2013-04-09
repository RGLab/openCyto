setMethod("show", sig = c("fcFilter"), definition = function(object) {
  cat("fcFilter:")
  callNextMethod(object)
})

setGeneric("posteriors", function(x, y, ...) standardGeneric("posteriors"))
setMethod("posteriors", sig = c("fcFilter", "ANY"), definition = function(x, y = "missing") {
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

setGeneric("priors", function(x, y, ...) standardGeneric("priors"))
setMethod("priors", sig = c("fcFilter", "ANY"), definition = function(x, y = "missing") {
  x@priors
})

setMethod("priors", sig = c("fcFilter", "character"), definition = function(x, y) {
  prior <- priors(x)
  
  priorNames <- names(prior)
  ind <- match(y, priorNames)
  if (is.na(ind)) {
    stop("FcFilter not found for:", y)
  } else {
    prior[[ind]]
  }
}) 
