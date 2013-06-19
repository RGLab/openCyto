setMethod("show", signature = c("gtMethod"), definition = function(object) {
  # cat('Gating Method: ')
  cat(names(object))
  cat("(")
  chnls <- dims(object)
  cat(paste(paste(names(chnls), chnls, sep = "="), collapse = ","))
  cat(") \n")
#  cat(parameters(object))
})

setMethod("show", signature = c("boolMethod"), definition = function(object) {
  cat(paste(class(object), "(", parameters(object), ")", sep = ""))
  cat("\n")
})

setMethod("names", signature = c("gtMethod"), definition = function(x) {
  x@name
})

setMethod("dims", signature = c("gtMethod"), definition = function(object) {
  dims <- strsplit(object@dims, ",")[[1]]
  if (length(dims) == 1) 
    dims <- c(NA, dims) else if (length(dims) != 2) 
    stop("invalid dimensions!")
  names(dims) <- c("xChannel", "yChannel")
  dims
})

setMethod("parameters", signature = c("gtMethod"), definition = function(object) {
  object@args
})

setGeneric("isCollapse",function(object,...)standardGeneric("isCollapse"))
setMethod("isCollapse", signature = c("gtMethod"), definition = function(object) {
      object@collapse
    })

setGeneric("groupBy",function(object,...)standardGeneric("groupBy"))
setMethod("groupBy", signature = c("gtMethod"), definition = function(object) {
      object@groupBy
    })
