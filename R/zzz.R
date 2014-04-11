#' global options for the package
.onLoad <- function(libname, pkgname) 
{
  options("openCyto" = list(gating = list(minEvents = 0 # minimum events needed in order for gating to proceed
                                          )
                            )
          )
}

