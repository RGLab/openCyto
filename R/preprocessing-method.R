#' apply a \link{ppMethod} to the \code{GatingSet}
#' 
#' @param x \code{ppMethod}
#' @param x \code{GatingSet}
#' ... other arguments
#' 
#' @inheritParams .preprocessing
#' 
#' @aliases preprocessing,ppMethod,GatingSet-method
setMethod("preprocessing", signature = c("ppMethod", "GatingSet"),
    definition = function(x, y, ...) {
      .preprocessing(x,y,...)
    })
setMethod("preprocessing", signature = c("ppMethod", "GatingSetList"),
    definition = function(x, y, ...) {
      .preprocessing(x,y,...)
    })
#' internal function (preprocessing)
#' 
#' @inheritParams .gating_gtMethod
#' @param gm: \code{gtMethod} object
.preprocessing <- function(x, y, gtPop, parent, gm
                            , mc.cores = 1, parallel_type = c("none", "multicore", "cluster"), cl = NULL
                            , ...) {
#  require("parallel")
  
  args <- parameters(x)
  # overwrite the golbal args with the local one
  args <- lattice:::updateList(args,list(...))
  
  
  ppm <- paste0(".", names(x))
  groupBy <- groupBy(x)
  dims <- dims(x)
  xChannel <- unname(dims["xChannel"])
  yChannel <- unname(dims["yChannel"])
  is_1d_gate <- any(is.na(dims))
  
  popAlias <- alias(gtPop)
  popName <- names(gtPop)
  popId <- gtPop@id
  
  gs_nodes <- basename(getChildren(y[[1]], parent, isPath = TRUE))
  if (length(gs_nodes) == 0 || !popAlias %in% gs_nodes) {
    message("Preprocessing for '", popAlias, "'")
    
    parent_data <- getData(y, parent)
#    parallel_type <- match.arg(parallel_type)
    ## get the accurate channel name by matching to the fr
    if (!is.na(xChannel)) {
      xParam <- getChannelMarker(parent_data[[1]], xChannel)
      xChannel <- as.character(xParam$name)
    }
    yParam <- getChannelMarker(parent_data[[1]], yChannel)
    yChannel <- as.character(yParam$name)
    
#    browser()
    if (nchar(groupBy) > 0) {
      
      split_by <- as.character(groupBy)
      split_by_num <- as.numeric(groupBy)
      #split by every N samples
      if(!is.na(split_by_num)){
        nSamples <- length(parent_data)
        if(nSamples==1){
          split_by <- 1
        }else{
          split_by <-  sample(rep(1:nSamples, each = split_by_num, length.out= nSamples))  
        }
        
      }else{
        #split by study variables
        pd <- pData(parent_data)
        split_by <- strsplit(split_by, ":")[[1]]
        split_by <- apply(pd[, split_by, drop = FALSE], 1, function(i)paste(i, collapse = ":"))
        split_by <- as.character(split_by)
      }
      fslist <- split(parent_data, split_by)
    }else 
    {
      fslist <- list(parent_data)  
    } 
    
    
    # construct method call
    thisCall <- substitute(f1())
    thisCall[["X"]] <- quote(fslist)  #set data
    thisCall[["FUN"]] <- as.symbol(ppm)  #set gating method
    thisCall[["xChannel"]] <- xChannel  #set x,y channel
    thisCall[["yChannel"]] <- yChannel
    thisCall[["gs"]] <- y
    thisCall[["gm"]] <- gm
    
    if(!(all(is.na(args)))){
      for (arg in names(args)) {
        thisCall[[arg]] <- args[[arg]]
      } 
    }
    
    thisCall[[1]] <- quote(lapply)  #select loop mode
    res <- eval(thisCall)
#    browser()
#    res <- unlist(res, recursive = FALSE)
  } else {
    message("Skip preprocessing! Population '", paste(popAlias, collapse = ","), "' already exists.")
    res <- NULL
  }
  
  
  
  res
  
}
