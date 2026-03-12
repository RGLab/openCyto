#' @templateVar old listgtMethods
#' @templateVar new gt_list_methods
#' @template template-depr_pkg
NULL
#'The environment holding the names of registered methods
#' @noRd 
.openCyto_plugin_method_lookup <- new.env()
.openCyto_plugin_method_lookup[["registered_methods"]] <- list(gating = character(0), preprocessing = character(0))

#' Environment holding dynamically registered plugin functions
#' @noRd
.openCyto_plugin_functions <- new.env(parent = emptyenv())

#' Resolve a plugin function by name
#'
#' First checks the plugin functions environment for dynamically registered
#' functions, then falls back to the package namespace for built-in functions.
#' @param name character function name (with leading dot)
#' @return the function object
#' @noRd
.resolve_plugin_function <- function(name) {
  if (exists(name, envir = .openCyto_plugin_functions, inherits = FALSE)) {
    return(get(name, envir = .openCyto_plugin_functions, inherits = FALSE))
  }
  ns <- getNamespace("openCyto")
  if (exists(name, envir = ns, inherits = FALSE)) {
    return(get(name, envir = ns, inherits = FALSE))
  }
  stop(sprintf("Plugin function '%s' not found", name))
}
.DEFAULT_GT <- c("quantileGate", "gate_quantile", "rangeGate","flowClust.2d", "gate_flowclust_2d", "mindensity", "gate_mindensity"
                 , "mindensity2", "gate_mindensity2", "flowClust.1d", "gate_flowclust_1d", "boundary","singletGate"
                 , "quadGate.tmix", "gate_quad_tmix", "quadGate.seq", "gate_quad_sequential", "gate_template", "gate_custom"
                 )
.DEFAULT_PP <- c("prior_flowClust", "warpSet", "standardize_flowset", "pp_gate_template", "pp_gate_custom")
#'Print a list of the registered gating methods
#'@name gt_list_methods
#'@aliases listgtMethods
#'@return Does not return anything. Prints a list of the available gating methods.
#'@export
gt_list_methods <- function(){
  res <- .getPluginMethods()
  cat("Gating Functions:\n")
  for(i in res[["gating"]])
      cat("=== ", i,  "\n")
  cat("Preprocessing Functions:\n")
  for(i in res[["preprocessing"]])
    cat("=== ", i,  "\n")
}

#'@export
listgtMethods <- function(){
  .Deprecated("gt_list_methods")
}

#'Is the method registered
#'
#'This will strip the preceding dot.
#' @noRd 
.isRegistered <- function(name){
  return(gsub("^\\.","", name)%in%unlist(.getPluginMethods()))
}

#'return a list of registered and default gating methods
#' @noRd 
.getPluginMethods <- function(x){
#  ns <- getNamespace("openCyto")
#  plugins <- getFromNamespace(".openCyto_plugin_method_lookup", ns = ns)
  
  plugin_gt <- .openCyto_plugin_method_lookup[["registered_methods"]][["gating"]]
  plugin_pp <- .openCyto_plugin_method_lookup[["registered_methods"]][["preprocessing"]]
  list(gating = c(plugin_gt, .DEFAULT_GT) 
      , preprocessing = c(plugin_pp, .DEFAULT_PP)
      )
  
}
#'Check the formal arguments of a gating function.
#'
#'The formal arguments need to match a certain template
#'We check that they do or do not.
#' @noRd 
.checkFormals <- function(frmls = NA, type = c("gating", "preprocessing")){
  
  type <- match.arg(type, c("gating", "preprocessing"))
  # we don't need to check ... since it is up to either
  # wrapper function or the actual gating function to handle its 
  # own formals
  if(type == "gating"){
    expected <- c("fr","pp_res","channels")
  }else{
    
    expected <- c("fs","gs", "gm", "channels", "groupBy", "isCollapse")
  }
  isMatched <- all(expected %in% names(frmls))
  if(!isMatched)
  {
    message("Formals of function don't match expected template: ")
    message("function(", paste(expected, collapse = ", "), ", ...)")
    return(FALSE)
  }
  return(TRUE)
}

registerGatingFunction <- function(fun=NA,methodName, dep=NA){
  .Defunct("registerPlugins")
}
#' @templateVar old registerPlugins
#' @templateVar new register_plugins
#' @template template-depr_pkg
NULL
#' Register a gating or preprocessing function with OpenCyto
#'
#' Function registers a new gating or preprocessing method with openCyto so that it may be used in the 
#' csv template.
#'
#' @name register_plugins
#' @aliases registerGatingFunction registerPlugins
#' @param fun \code{function} to be registered
#' @param methodName \code{character} name of the gating or preprocessing method
#' @param dep \code{character} name of the library dependency required for the plugin method to work.
#' @param ... other arguments
#'         type \code{character} specifying the type of registering method. Should be either "gating" or "preprocessing".
#' 
#' @return \code{logical} TRUE if successful and prints a message. FALSE otherwise.
#' @useDynLib openCyto,.registration = TRUE
#' @details The \code{fun} argument should be a wrapper function definition for the gating or preprocessing method. 
#'                          Gating method must have formal arguments:
#' 
#'                           fr a \code{flowFrame}
#' 
#'                           pp_res a pre-processing result
#' 
#'                           xChannel \code{character} (optional)
#' 
#'                           yChannel \code{character} (required)
#' 
#'                           filterId \code{character}
#' 
#'                           ... ellipses for the additional parameters.
#' 
#'                          Preprocessing method must have formal arguments:
#' 
#'                          fs a \code{flowSet} that stores the flow data (could be subgrouped data if \code{groupBy} column is defined in the csv template
#' 
#'                          gs a \code{GatingSet}
#'  
#'                          gm a \code{gtMethod} object that stores the information from gating method
#' 
#'                          xChannel \code{character} (required)
#' 
#'                          yChannel \code{character} (required)
#' 
#'                           ... ellipses for the additional parameters.
#' 
#' The gating function must return a filter (i.e. polygonGate or other instance) from flowCore.
#' The preprocessing can return anything and it will be passed on to the gating function. So it is up to gating function to use and interpret the results of preprocessing.
#' Not all formal parameters need to be used. Additional arguments are passed via the ... and can be processed in the wrapper
#' 
#' @import utils
#' @export
register_plugins <- function(fun = NA, methodName, dep = NA, ...){
  
  if(!is.na(dep)){
    if(is.character(dep)){
      if(system.file(package = dep) == ""){
        message(sprintf("Can't register %s with dependency on %s, because dependency is not installed.",methodName,dep))
        return(FALSE)
      }
    }else{
      warning("If provided, dep must be a character naming the package dependency.")
      return(FALSE)
    }
  }
  if(!is.function(fun)){
    warning("You need to put the fun in function! (argument fun is not a function)")
    return(FALSE)
  }else{
    ###Check the formal arguments
    frmls <- formals(fun)
    
    if(.checkFormals(frmls = frmls, ...)){
      if(.register(fun = fun, methodName = methodName, ...)){
        message(sprintf("Registered %s",methodName))
      }
      
    }else{
      warning("Can't register function")
      return(FALSE)
    }
  }
  return(TRUE)
}

#' @export
registerPlugins <- function(fun = NA, methodName, dep = NA, ...){
  .Deprecated("register_plugins")
  register_plugins(fun, methodName, dep, ...)
}


#'Register a gating function for OpenCyto
#' @noRd 
.register <- function(fun = NA,methodName, type = c("gating", "preprocessing")){
  
  type <- match.arg(type, c("gating", "preprocessing"))
  
  methodName <- paste0(".",methodName)
  
  # Store function in the plugin functions environment
  assign(methodName, fun, envir = .openCyto_plugin_functions)
  
  #add to the plugin method list
  
  toAdd <- gsub("^\\.","",methodName)
  current <- .openCyto_plugin_method_lookup[["registered_methods"]][[type]]
  
  if(length(current) == 0)
    found <- FALSE
  else if(any(grepl(toAdd, current)))
    found <- TRUE
  else
    found <- FALSE
  
  if(!found)
    .openCyto_plugin_method_lookup[["registered_methods"]][[type]] <- c(current, toAdd)  
  
  return(TRUE)
}

#' only for internal usage (debug)
#' @noRd 
.unregister <- function(methodName, type = c("gating", "preprocessing")){

  methodName <- paste0(".",methodName)
  type <- match.arg(type, c("gating", "preprocessing"))
  
  
  toRm <- gsub("^\\.","",methodName)
  
  current <- .openCyto_plugin_method_lookup[["registered_methods"]][[type]]
  
  if(length(current) == 0)
    found <- FALSE
  else{
    ind <- grepl(toRm, current)
    if(any(ind))
      found <- TRUE
    else
      found <- FALSE
  } 
  
  if(found)
    .openCyto_plugin_method_lookup[["registered_methods"]][[type]] <- current[!ind]
  
  # Remove from plugin functions environment
  if (exists(methodName, envir = .openCyto_plugin_functions, inherits = FALSE)) {
    rm(list = methodName, envir = .openCyto_plugin_functions)
  }
  
  return(TRUE)
}

