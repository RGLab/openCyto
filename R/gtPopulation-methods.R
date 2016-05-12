#' get population name
#' @param x \code{gtPopulation} object
#' @aliases names,gtPopulation-method
setMethod("names", signature = c("gtPopulation"), definition = function(x) {
  x@name
})

# the original S4 method has been suffering from S4 masking issue,thus we keep it as private functioninstead
# get population alias
# @param object \code{gtPopulation} 
#setMethod("alias", signature = c("gtPopulation"), definition = function(object) {
alias <- function(object){
  object@alias
}

 
