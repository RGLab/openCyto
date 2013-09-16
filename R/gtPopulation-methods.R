#' get population name
#' @param \code{gtPopulation} object
#' @aliases names,gtPopulation-method
setMethod("names", signature = c("gtPopulation"), definition = function(x) {
  x@name
})

#' get population alias
#' @param \code{gtPopulation} object
#' @export
#' @aliases alias,gtPopulation-method 
setMethod("alias", signature = c("gtPopulation"), definition = function(object) {
  object@alias
})

 
