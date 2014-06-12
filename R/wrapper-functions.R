# This file contains all wrapper methods for dispatching data and arguments to
# gating/preprocessing algorithms.


#' An adapter to connect the gating wrapper function with the \link{gating} method
#' 
#' It coerce the input (\code{flowSet}) to a single \code{flowFrame} and apply the gating wrapper function
#' then replicate the gates across samples from the \code{flowSet}.
#' 
#' @inheritParams .prior_flowClust
#' @param pp_res preprocessing result produced by the \code{preprocessing} method
#' @param gFunc \code{character} function name of the wrapper function to be invoked
#' @param popAlias \code{character} the population names that are used to determine how many gates to be expected from the gating function 
#' @param ... other arguments to be passed to wrapper function
#' 
#' @return a \code{list} of \code{filter}s
.gating_wrapper <- function(fs, pp_res, gFunc, popAlias, channels, ...){
    require(openCyto)  #since it is going to be invoked by MPI, better load it
    #coercing
    sn <- sampleNames(fs)
    fr <- as(fs,"flowFrame")
    openCyto.options <- getOption("openCyto")
    minEvents <- openCyto.options[["gating"]][["minEvents"]]
    
    if(is.null(minEvents))
      minEvents <- 0
    if(nrow(fr) <= minEvents){
      warning(paste(sn, collapse =","), ": Not enough events to proceed the data-driven gating!Returning a dummy gate instead.")
      
      #create dummy rectangleGate
      #TODO: move channels to ... to deprecate x,y channel
#      dots <- list(...)
#      channels <- dots$channels
      channels <- as.vector(na.omit(channels))
      nDim <- length(channels) 
      
      if(nDim ==  1)
        gate_coordinates <- list(c(-Inf, -Inf))
      else if(nDim ==  2)
        gate_coordinates <- list(c(-Inf, -Inf), c(-Inf, -Inf))
      else
        stop(nDim, " dimensional gating is not supported yet!")

      names(gate_coordinates) <- channels
      filterRes <- rectangleGate(gate_coordinates)
    
      #this is flowClust-specific operation, which
      # be abstracted out of this framework
      
      if(grepl("flowClust\\.[12]d", gFunc))
        filterRes <- fcRectangleGate(filterRes, priors = list(), posts = list())
      
      nPop <- length(popAlias)
      filterResType <- ifelse(nPop == 1, "filter", "filters")
      if(filterResType == "filters"){
        filterRes <- filters(lapply(1:nPop, function(i)filterRes))
      }
      
    }else{
      if(!.isRegistered(gFunc)){
        stop(sprintf("Can't gate using unregistered method %s",gFunc))
      }
      thisCall <- substitute(f(fr = fr, pp_res = pp_res, ...),list(f=as.symbol(gFunc)))
      filterRes <- try(eval(thisCall), silent = TRUE)  
    }
      
    
        
    if(extends(class(filterRes), "filter")||extends(class(filterRes), "filters")){

      #replicate the filter across samples
      list(sapply(sampleNames(fs),function(i)filterRes, simplify = FALSE))      
    }else{
      stop("failed at ",paste0(sn), "\n", filterRes)
    }
#

}
#' wrapper for \link[flowStats:singletGate]{singletGate}
#' 
#' @param pp_res not used
#' @param fr \code{flowFrame} object as a data input 
#' @param ... arguments to be passed to \link[flowStats:singletGate]{singletGate}
#' @inheritParams .prior_flowClust
#' @return a \code{filter} object
#' @importFrom flowStats singletGate
.singletGate <- function(fr, pp_res = NULL, xChannel = "FSC-A", yChannel = "FSC-H", ...) {
  
  fr <- fr[, c(xChannel,yChannel)]
  # Creates a list of polygon gates based on the prediction bands at the minimum
  # and maximum x_channel observation using a robust linear model trained by
  # flowStats.

  singletGate(fr, area = xChannel, height = yChannel, ...)
}

#' boundary gating function
#' 
#' It essentially constructs a rectangle gate from input range (min, max)
#'
#' @description It is useful for filtering out very large signal (e.g. FSC,SSC, Time)
#'  
#' @inheritParams .singletGate
#' @param min,max the range input for constructing the \code{rectangleGate}
#' @param ... other arguments (not used.)
#' @return a \code{filter} object
.boundary <- function(fr, pp_res = NULL, xChannel = NULL, yChannel, min = NULL, max = NULL, ...) {
  
  if (is.na(xChannel)) {
    xChannel <- NULL
  }
  channels <- c(xChannel, yChannel)
  num_channels <- length(channels)
  
  if (is.null(min)) {
    min <- rep(-Inf, num_channels)
  }
  if (is.null(max)) {
    max <- rep(Inf, num_channels)
  }
  
  if (!(num_channels == length(min) && num_channels == length(max))) {
    stop("The lengths of 'min' and 'max' must match the number of 'channels' given.")
  }
  
  gate_coordinates <- lapply(seq_len(num_channels), function(i) {
        c(min[i], max[i])
      })
  names(gate_coordinates) <- channels
  
  rectangleGate(gate_coordinates)
}

#' wrapper for flowClust.1d
#' 
#' It does some parameter preprocessing before calling the flowClust.1d
#' 
#' @param fr \code{flowFrame} object as a data input 
#' @param ... arguments to be passed to \link{flowClust.1d}
#' @param xChannel must be \code{NA}.
#' @param yChannel the dimension used for gating
#'  
#' @inheritParams .gating_wrapper 
#' @inheritParams .prior_flowClust
#' 
#' @return a \code{filter} object
.flowClust.1d <- function(fr, pp_res, xChannel = NA, yChannel, ...) {
  
  
  prior <- pp_res
  
#  sname <- sampleNames(fs)
#  fr <- fs[[sname]]
  priorList <- list()
  if(!is.null(prior))
    prior <- prior[[yChannel]]
  
  args <- list(...)
  # parse arguments for flowClust
  
  # If any of 'K, 'neg' or 'pos' are given, we extract them from the
  # 'args' and coerce them as integers. Otherwise, they are defaulted to
  # NULL.
  # NOTE: the named elements within 'args' are coerced to lowercase.
  K <- neg_cluster <- pos_cluster <- NULL
  
  if ("K" %in% names(args)) {
    K <- as.integer(args["K"])
  }
  if ("neg" %in% names(args)) {
    neg_cluster <- as.integer(args["neg"])
    args[["neg"]] <- NULL
  }
  if ("pos" %in% names(args)) {
    pos_cluster <- as.integer(args["pos"])
    args[["pos"]] <- NULL
  }
  
  # Examines the various cases of pos, neg and K.
  # In the case that K is not given, we automatically select it.
  # Exception: If both 'pos' and 'neg' are given, we set: K <- pos + neg
  if (is.null(K)) {
    if (!is.null(neg_cluster) && !is.null(pos_cluster)) {
      K <- neg_cluster + pos_cluster
    }
  } else {
    # If pos and neg are given, throw a warning and set: K <- pos + neg
    # In the case that one of pos and neg are given and either exceeds K,
    # we throw an error.         
    if (!is.null(neg_cluster) && !is.null(pos_cluster)) {
      warning("Values given for 'K', 'neg', and 'pos'. Setting K = neg + pos")
      K <- neg_cluster + pos_cluster
    } else if (!is.null(pos_cluster) && K < pos_cluster) {
      stop("The number of positive clusters exceeds 'K'.")
    } else if (!is.null(neg_cluster) && K < neg_cluster) {
      stop("The number of negative clusters exceeds 'K'.")
    }
  }
  
  args[["neg_cluster"]] <- neg_cluster
#        args[["pos_cluster"]] <- pos_cluster
  
#  args[["k"]] <- NULL
  args[["K"]] <- K
  
  # If 'min' and/or 'max' are given, we pass this value along to the
  # prior-elicitation method as well as flowClust. Otherwise, these values
  # are set to NULL.
  min_values <- -Inf
  max_values <- Inf
  if ("min" %in% names(args)) {
    min_values <- as.numeric(args["min"])
  }
  if ("max" %in% names(args)) {
    max_values <- as.numeric(args["max"])
  }                          
  
  # If the number of positive clusters is 0 and no cutpoint method has been
  # specified, we use the quantile gate by default.
  if (!is.null(pos_cluster) && pos_cluster == 0 && is.null(args[["cutpoint_method"]])) {
    args[["cutpoint_method"]] <- "quantile"
  }
  
  
  if (is.na(xChannel)) {
    # 1d gate
    do.call("flowClust.1d"
            ,args = c(list(fr = fr
                        ,params = yChannel
                        ,prior = prior
                        )
                      ,args
                    )
           )
  } else {
    stop("flowClust1d does not support 2d gate!")
  }
}
#' wrapper for cytokine
#' 
#' It does some parameter preprocessing before calling the cytokine
#' 
 
#' @param ... arguments to be passed to \link{cytokine}
#' @inheritParams .flowClust.1d 
#' 
#' @return a \code{filter} object
.cytokine <- function(fr, pp_res, xChannel = NA, yChannel = "FSC-A", ...) {
  .Deprecated("tailgate")
  #TODO:standardize data with pp_res
  .tailgate(fr, pp_res = pp_res, yChannel = yChannel, ...)
}

#' @param ... arguments to be passed to \link{tailgate}
#' @inheritParams .flowClust.1d 
#' 
#' @return a \code{filter} object
.tailgate <- function(fr, pp_res, xChannel = NA, yChannel = "FSC-A", ...) {
#  if(keyword(fr)[["$FIL"]] == "DC,2f,MONO,2f,NK 22013_12828_003.fcs")
#    browser()  
  
  
  #pps_res may contains the standardized and collapsed data and transformation
  if(isTRUE(attr(pp_res, "openCyto_preprocessing") == "standardize")){
      transformedData <- pp_res[["flow_frame"]]
     
      #if no flow data passed in, need to tranform original fr first
     if(is.null(transformedData)){
      data <- exprs(fr)[, yChannel]
      exprs(fr)[, yChannel] <- with(pp_res, (data - center)/scale) 
      transformedData <- fr
     }
     
    g <- tailgate(transformedData, channel = yChannel, ...)
    gate_coordinates <- c(g@min, g@max)

    # Backtransforms the gate 
    gate_coordinates <- with(pp_res, center + scale * gate_coordinates)
    gate_coordinates <- list(gate_coordinates)
    names(gate_coordinates) <- parameters(g)
    rectangleGate(gate_coordinates, filterId = g@filterId)
          
  }else
    tailgate(fr, channel = yChannel, ...)

  
  # If a sample has no more than 1 observation when the 'cytokine' gate is
  # attempted, the 'center' and/or 'scale' will result be NA, in which case we
  # replace the resulting NA cutpoints with the average of the remaining
  # cutpoints. If all of the cutpoints are NA, we set the mean to 0, so that
  # all of the cutpoints are 0.
#  cutpoints_unlisted <- unlist(cutpoints)
#  if (sum(!is.na(cutpoints_unlisted)) > 0) {
#    mean_cutpoints <- mean(cutpoints_unlisted, na.rm = TRUE)
#  } else {
#    mean_cutpoints <- 0
#  }
#  cutpoints <- as.list(replace(cutpoints_unlisted, is.na(cutpoints_unlisted),
#          mean_cutpoints))
  
}

#' wrapper for mindensity
#' 
#' It does some parameter preprocessing before calling the mindensity
#' 

#' @param ... arguments to be passed to \link{mindensity}
#' @inheritParams .flowClust.1d 
#' 
#' @return a \code{filter} object
.mindensity <- function(fr, pp_res, yChannel = "FSC-A", ...) {
  
 
  mindensity(flow_frame = fr, channel = yChannel, ...)
}
#' wrapper for flowClust.2d
#' 
#' It does some parameter preprocessing before calling the flowClust.2d
#' 
#' @param ... arguments to be passed to \link{flowClust.2d}
#' @inheritParams .flowClust.1d 
#' 
#' @return a \code{filter} object
.flowClust.2d <- function(fr, pp_res, xChannel, yChannel, ...) {
  
  
  args <- list(...)
  # get the value of neg and pos
  if (is.element(c("K"), names(args))) {
    K <- as.integer(args["K"])
  } else {
    message("'K' argument is missing! Using default setting: K = 2")
    K <- 2
  
  }
  
  args[["K"]] <- K

  if (is.null(pp_res)) {
    usePrior <- "no"
    pp_res <- list(NA)
  } else {
    usePrior <- "yes"
  }
  
  do.call("flowClust.2d"
      ,args = c(list(fr = fr
                    , xChannel = xChannel
                    , yChannel = yChannel
                    , usePrior = usePrior
                    ,prior = pp_res
                    )
                ,args
              )
        )
  
}
#' wrapper for rangeGate (deprecated)
#' 
#' It does some parameter preprocessing before calling the rangeGate
#' 
#' @param pp_res not used
#' @param ... arguments to be passed to \link[flowStats:rangeGate]{rangeGate}
#' @inheritParams .flowClust.1d 
#' 
#' @return a \code{filter} object
#' @importFrom flowStats rangeGate
.rangeGate <- function(fr, pp_res, xChannel = NA, yChannel,  ...) {
  
  rangeGate(x = fr, stain = yChannel,  ...)
}
#' wrapper for quantileGate
#' 
#' It does some parameter preprocessing before calling the quantileGate
#' 
#' @param ... arguments to be passed to \link{quantileGate}
#' @inheritParams .flowClust.1d 
#' 
#' @return a \code{filter} object
.quantileGate <- function(fr, pp_res, xChannel = NA, yChannel, ...) {
  
  quantileGate(fr = fr, stain = yChannel, ...)
}

.quadGate.tmix <- function(fr, pp_res, xChannel = NA, yChannel, ...) {
  quadGate.tmix(fr, c(xChannel, yChannel), ...)
}
#' deprecated
#' @importFrom flowStats quadrantGate
.quadrantGate <- function(fr, pp_res, xChannel = NA, yChannel, ...) {
  
  qfilter <- quadrantGate(fr, stains = c(xChannel, yChannel), absolute = FALSE, 
                 inBetween = TRUE, ...)
  
  ###############################################################     
  #construct rectangleGates based on the cuts and popNames,clock-wise
  ###############################################################
  cut.x <- qfilter@boundary[xChannel]
  cut.y <- qfilter@boundary[yChannel]
  gateList <- new("filters")
  
  chnls <- c(xChannel, yChannel)
  markers <- chnls
  
  coord <- list(c(-Inf, cut.x), c(cut.y, Inf))
  names(coord) <- as.character(chnls)
  gateList[[paste(paste0(markers, c("-", "+")), collapse = "")]] <- rectangleGate(coord)
  
  coord <- list(c(cut.x, Inf), c(cut.y, Inf))
  names(coord) <- as.character(chnls)
  gateList[[paste(paste0(markers, c("+", "+")), collapse = "")]] <- rectangleGate(coord)
  
  coord <- list(c(cut.x, Inf), c(-Inf, cut.y))
  names(coord) <- as.character(chnls)
  gateList[[paste(paste0(markers, c("+", "-")), collapse = "")]] <- rectangleGate(coord)
  
  coord <- list(c(-Inf, cut.x), c(-Inf, cut.y))
  names(coord) <- as.character(chnls)
  gateList[[paste(paste0(markers, c("-", "-")), collapse = "")]] <- rectangleGate(coord)
  
  gateList
}

############################
# preprocessing wrappers
#########################
#'  wrapper for \link[flowStats:warpSet]{warpSet}
#' 
#' @param stains \code{character} passed to \link[flowStats:warpSet]{warpSet} 
#' @inheritParams .prior_flowClust 
#' 
#' @return \code{NULL}
#' @importFrom flowStats warpSet
.warpSet <- function(fs, gs, gm, xChannel, yChannel, groupBy, isCollapse, stains, ...){
  fs <- fs[, stains]
  if(class(fs) == "ncdfFlowSet")
    flowStats:::warpSetNCDF(fs, stains = stains, ...)
  else
    warpSet(fs, stains = stains, ...)
  return (NULL)
 }
#'  wrapper for prior_flowClust
#' 
#'  This wrapper does some parameter preprocessing before calls \link{prior_flowClust}
#' 
#' @param fs \code{flowSet} or \code{ncdfFlowSet} object
#' @param gs \code{GatingSet}
#' @param gm \code{gtMethod}
#' @param xChannel,yChannel \code{character} specifying the dimensions of flow data used by \code{prior_flowClust}
#' @param prior_source \code{character} specifying the ancester node from where the prior is elicited.
#' @param neg,pos \code{numeric} specifying how many peaks are expected on positive and negative sides from 1d density profile
#' @inheritParams .prior_flowClust1d 
#' 
#' @return a \code{list} of priors, see \link{prior_flowClust} for more details
.prior_flowClust <- function(fs, gs, gm, xChannel, yChannel, groupBy, isCollapse
     , prior_source = NULL
     , K = NULL
     , neg, pos
     , min, max
     , ...){
   prior_list <- list()
   
   # prior estimation is done separately from flowClust routine because
   # .prior_flowClust1d requires the entire parent flowSet yet flowClust only
   # takes one flowFrame
   
   if (is.null(prior_source)) {
     prior_data <- fs
   } else {
     prior_data <- getData(gs, prior_source)
   }
#    browser()
   
   if (names(gm) == "flowClust.1d") {
     # If any of 'K, 'neg' or 'pos' are given, we extract them from the
     # 'args' and coerce them as integers. Otherwise, they are defaulted to
     # NULL.
     
     neg_cluster <- pos_cluster <- NULL
     
     if(!is.null(K)) {
       K <- as.integer(K)
     }
     if(!missing(neg)) {
       neg_cluster <- as.integer(neg)
     }
     if(!missing(pos)) {
       pos_cluster <- as.integer(pos)
     }
     
     # Examines the various cases of pos, neg and K.
     # In the case that K is not given, we automatically select it.
     # Exception: If both 'pos' and 'neg' are given, we set: K <- pos + neg
     if (is.null(K)) {
       if (!is.null(neg_cluster) && !is.null(pos_cluster)) {
         K <- neg_cluster + pos_cluster
       }
     } else {
       # If pos and neg are given, throw a warning and set: K <- pos + neg
       # In the case that one of pos and neg are given and either exceeds K,
       # we throw an error.         
       if (!is.null(neg_cluster) && !is.null(pos_cluster)) {
         warning("Values given for 'K', 'neg', and 'pos'. Setting K = neg + pos")
         K <- neg_cluster + pos_cluster
       } else if (!is.null(pos_cluster) && K < pos_cluster) {
         stop("The number of positive clusters exceeds 'K'.")
       } else if (!is.null(neg_cluster) && K < neg_cluster) {
         stop("The number of negative clusters exceeds 'K'.")
       }
     }
     
     # If 'min' and/or 'max' are given, we pass this value along to the
     # prior-elicitation method as well as flowClust. Otherwise, these values
     # are set to NULL.
     min_values <- -Inf
     max_values <- Inf
     if (!missing(min)) {
       min_values <- as.numeric(min)
     }
     if (!missing(max)) {
       max_values <- as.numeric(max)
     }                          
     
     if (!is.na(xChannel)) {
       prior_list[[xChannel]] <- prior_flowClust(flow_set = prior_data,
           channels = xChannel, K = K,
           min = min_values, max = max_values, ...)
     }
     
     prior_list[[yChannel]] <- prior_flowClust(flow_set = prior_data,
         channels = yChannel, K = K,
         min = min_values, max = max_values, ...)
     
   } else {
     # get the value of neg and pos
     if (!is.null(K)) {
       K <- as.integer(K)
     } else {
       message("'K' argument is missing in prior_flowClust! Using default setting: K = 2.\nYou should set this to the same value as 'K' in the call to flowClust.")
       K <- 2
     }
     
     prior_list <- prior_flowClust(flow_set = prior_data, channels = c(xChannel, 
             yChannel), K = K, ...)
   }
   if(isCollapse)
     prior_list
   else# when collapse is FALSE, gating is done at sample level, thus we need to replicate priors across samples so that it matches with the gating input
     sapply(sampleNames(prior_data), function(i)prior_list, simplify = FALSE)
   
}

 #' preprocessing wrapper for .standardize_flowFrame
#' 
#' @param fs a \code{flowSet} object
#' @param xChannel not used
#' @param yChannel the channel to standardize
#' @param groupBy \code{character} tells whether 'groupBy' column in csv template has been set
#' @param isCollapse \code{logical} indicates whether the gating is done on collapsed data. If so, the prior does not need to be replicated across samples.
 .standardize_flowset <- function(fs, gs, gm, xChannel = NA, yChannel, groupBy, isCollapse, ...) {
   
   if(isCollapse)
     stop("'collapse = TRUE' is not applicable to 'standardize_flowset'!")
   
   isGroup <- nchar(groupBy) > 0
   
   if(isGroup){
#     if(pData(fs)$Center == "UCLA" && pData(fs)$Sample == "12828")
#       browser()
     #standardize multiple flowFrames within the sub-group
     transform_out <- fsApply(fs, .standardize_flowFrame, channel = yChannel, data = TRUE)
     
     # Creates a flowSet object from the transformed flowFrame objects
     fs <- flowSet(lapply(transform_out, function(x) x$flow_frame))
     collapsedFr <- as(fs, "flowFrame")
     
     # update each flow_frame with collapsed flowFrame
     transform_out <- lapply(transform_out, function(x) {
           x$flow_frame <- collapsedFr
           x
         })
     
   }else{
     #standarize within the single flowFrame
     # fs is now the entire data set
     # and we only pass the transformation back to gating function without the flow data (to save memory)
     transform_out <- fsApply(fs, .standardize_flowFrame, channel = yChannel, data = FALSE)
   }
   
   
   transform_out
 }
 