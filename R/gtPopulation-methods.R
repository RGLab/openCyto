#' get population name
#' @param \code{gtPopulation}
setMethod("names", signature = c("gtPopulation"), definition = function(x) {
  x@name
})

#' get population alias
#' @param \code{gtPopulation}
#' @export 
setMethod("alias", signature = c("gtPopulation"), definition = function(object) {
  object@alias
})

 
