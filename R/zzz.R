#' @name openCyto.options
#' @description  Get/set some global options for openCyto
#' @rdname openCyto.options
#' @aliases openCyto.options
#' @title Some global options for openCyto
#' See examples for the meaning of these options and how to get/set them.
#' @examples 
#' opt <- getOption("openCyto")
#' #the threshold of minimum cell events required for the gating algorithm to proceed
#' opt[["gating"]][["minEvents"]]
#' #to change the threshold
#' opt[["gating"]][["minEvents"]] <- 100
#' options(openCyto = opt)
#' 
#' #switch off the validity check flags(Not recommended)
#' opt[["check.pop"]] <- FALSE
#' options(openCyto = opt)
#' 
NULL

.onLoad <- function(libname, pkgname) 
{
  options("openCyto" = list(gating = list(minEvents = 0 # minimum events needed in order for gating to proceed
                                          )
                            , check.pop = TRUE
                            )
          )
}

