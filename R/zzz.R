#' global options for the package
.onLoad <- function(libname, pkgname) 
{
  options("openCyto" = list(gating = list(minEvents = 0 # minimum events needed in order for gating to proceed
                                          )
                            , default.methods = list(gating1D = NULL #default is mindensity,tailgate,flowClust.1d
                                                     , gating2D = NULL# default is flowClust.2d
                                                     , singlet = NULL) # default is singletGate                                 
                            )
          )
}

