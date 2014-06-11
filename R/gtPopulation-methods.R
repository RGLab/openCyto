#' get population name
#' @param x \code{gtPopulation} object
#' @aliases names,gtPopulation-method
setMethod("names", signature = c("gtPopulation"), definition = function(x) {
  x@name
})

#' get population alias
#' @param object \code{gtPopulation} 
#' @export
#' @aliases alias,gtPopulation-method 
setMethod("alias", signature = c("gtPopulation"), definition = function(object) {
  object@alias
})

 
