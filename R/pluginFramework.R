#'The environment holding the names of registered methods
.openCyto_gtmethod_lookup<-new.env()
.DEFAULTS <- c("quadrantGate", "quantileGate","rangeGate","flowClust.2d","mindensity","cytokine","flowClust.1d","boundary","singletGate", "tailgate")

#'Print a list of the registered gating methods
#'@return Does not return anything. Prints a list of the available gating methods.
#'@export listgtMethods
listgtMethods<-function(){
  ns<-getNamespace("openCyto")
  tbl<-getFromNamespace(".openCyto_gtmethod_lookup",ns=ns)
  nms<-c(.DEFAULTS,(names(as.list(tbl))))
  if(length(nms)>0){
    cat(nms,sep="\n")
  }else{
    cat("no additional gating methods registered")
  }
}

#'Is the method registered
#'
#'This will strip the preceding dot.
.isRegistered<-function(gtmethod){
  return(gsub("^\\.","",gtmethod)%in%.getgtMethods())
}

#'return a list of registered and default gating methods
.getgtMethods<-function(x){
  ns<-getNamespace("openCyto")
  tbl<-getFromNamespace(".openCyto_gtmethod_lookup",ns=ns)
  nms<-c(names(as.list(tbl)),.DEFAULTS)
  if(length(nms)>0){
    return(nms)
  }else{
    return(NULL)
  }
}
#'Check the formal arguments of a gating function.
#'
#'The formal arguments need to match a certain template
#'We check that they do or do not.
.checkFormals <- function(frmls=NA){
  expected<-c("fr","pp_res","yChannel","filter_id","...")
  posn<-sapply(expected,function(x)which(names(frmls)%in%x))
  if(!(all.equal(posn,c(fr=1,pp_res=2,yChannel=4,filter_id=5,"..."=6))|isTRUE(all.equal(posn,c(fr=1,pp_res=2,yChannel=3,filter_id=4,"..."=5))))){
    message("Formals of function don't match expected template.")
    return(FALSE)
  }else{
    if(any(names(frmls)%in%"xChannel")){
      if(which(names(frmls)%in%"xChannel")!=3){
        message("Formals of function don't match expected template")
        return(FALSE)
      }
    }
  }
  return(TRUE)
}

#'Register a gating function with OpenCyto
#'
#'Function registers a new gating method with openCyto so that it may be used in the 
#'csv template.
#'@param fun \code{function} to be registered
#'@param methodName \code{character} name of the gating method
#'@param dep \code{character} name of the library dependency required for the plugin gating method to work.
#'@return \code{logical} TRUE if successful and prints a message. FALSE otherwise.
#'@export
#'@useDynLib openCyto
#'@details The \code{fun} argument should be a wrapper function definition for the gating method, with formal arguments
#'fr a flowFrame, pp_res a pre-processing result, xChannel (optional character vector), yChannel (required character vector), and a filter_id (character), and an ellipsis ... for additional parameters.
#'The function must return a filter (i.e. polygonGate or other instance) from flowCore.
#'Not all formal parameters need to be used. Additional arguments are passed via the ... and can be processed in the wrapper
#'@import utils
#'@importFrom R.utils isPackageInstalled
registerGatingFunction<-function(fun=NA,methodName="myGatingMethod",dep=NA){
  if(!is.na(dep)){
    if(is.character(dep)){
      if(!isPackageInstalled(dep)){
        message(sprintf("Can't register %s with dependency on %s, because dependency is not installed.",methodName,dep))
        invisible(FALSE)
      }
    }else{
      warning("If provided, dep must be a character naming the package dependency.")
      invisible(FALSE)
    }
  }
  if(!is.function(fun)){
    warning("You need to put the fun in function! (argument fun is not a function)")
    invisible(FALSE)
  }else{
    ###Check the formal arguments
    frmls<-formals(fun)
    if(.checkFormals(frmls=frmls)){
      if(.register(fun=fun,methodName=methodName)){
        message(sprintf("Registered %s",methodName))
      }
      
    }else{
      warning("Can't register function")
      invisible(FALSE)
    }
  }
}

#'Register a gating function for OpenCyto
.register<-function(fun=NA,methodName="myGatingMethod"){
  methodName <- paste0(".",methodName)
  ENV <- getNamespace("openCyto")
  openCyto:::unlockNamespace(ENV)  
  try(unlockBinding(methodName,ENV),silent=TRUE)
  assign(methodName,fun,ENV)
  e<-getFromNamespace(".openCyto_gtmethod_lookup",ns=getNamespace("openCyto"))
  assign(gsub("^\\.","",methodName)," ",envir=e)
  lockBinding(methodName,env=ENV)
  lockEnvironment(ENV)
  return(TRUE)
}
