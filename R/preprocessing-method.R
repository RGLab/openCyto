setGeneric("preprocessing", function(x, y, ...) standardGeneric("preprocessing"))
setMethod("preprocessing", signature = c("ppMethod", "GatingSet"),
    definition = function(x, y, ...) {
      .preprocessing(x,y,...)
    })

.preprocessing <- function(x, y, gtPop, parent, mc.cores = 1,
    parallel_type = c("none", "multicore", "cluster"), cl = NULL, plot = FALSE,
    xbin = 128, prior_group = NULL, ...) {
  require("parallel")
#  browser()
  args <- parameters(x)
  ppm <- paste0(".", names(x))
  
  dims <- dims(x)
  xChannel <- unname(dims["xChannel"])
  yChannel <- unname(dims["yChannel"])
  is_1d_gate <- any(is.na(dims))
  
  popAlias <- alias(gtPop)
  popName <- names(gtPop)
  popId <- gtPop@id
  gs_nodes <- getChildren(y[[1]], getNodes(y[[1]], parent))
  
  if (length(gs_nodes) == 0 || !popAlias %in% gs_nodes) {
    message("Preprocessing for '", popAlias, "'")
    
    parent_data <- getData(y, parent)
    parallel_type <- match.arg(parallel_type)
    ## get the accurate channel name by matching to the fr
    if (!is.na(xChannel)) {
      xParam <- getChannelMarker(parent_data[[1]], xChannel)
      xChannel <- as.character(xParam$name)
    }
    yParam <- getChannelMarker(parent_data[[1]], yChannel)
    yChannel <- as.character(yParam$name)
    
    # Splits the flow set into a list.
    # By default, each element in the list is a flowSet containg one flow frame,
    # corresponding to the invidual sample names.
    # If 'split' is given, we split using the unique combinations within pData.
    # In this case 'split' is specified with column names of the pData.
    # For example, "PTID:VISITNO"
    # when split is numeric, do the grouping by every N samples
    if ("split" %in% names(args)) {
      
      split_by <- as.character(args["split"])
      split_by_num <- as.numeric(split_by)
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
        split_by <- strsplit(split_by, ":")[[1]]
        split_by <- apply(pData(parent_data)[, split_by], 1, paste, collapse = ":")
        split_by <- as.character(split_by)
      }
      args[["split"]] <- NULL
    } else {
      split_by <- sampleNames(parent_data)
    }
    
    fslist <- split(parent_data, split_by)
    
    # construct method call
    thisCall <- substitute(f1())
    thisCall[["X"]] <- quote(fslist)  #set data
    thisCall[["FUN"]] <- as.symbol(ppm)  #set gating method
    thisCall[["xChannel"]] <- xChannel  #set x,y channel
    thisCall[["yChannel"]] <- yChannel
    
    if (is_1d_gate) {
      if (grepl("-$", popName)) {
        positive <- FALSE
      } else if (grepl("+$", popName)) {
        positive <- TRUE
      } else {
        stop("Invalid population name! Name should end with '+' or '-' symbol.")
      }
      
      thisCall[["positive"]] <- positive
    }
    
    prior_list <- list()
    
    # prior estimation is done separately from flowClust routine because
    # prior_flowClust1d requires the entire parent flowSet yet flowClust only
    # takes one flowFrame
    if (grepl("^\\.flowClust\\.[12]d$", gm)) {
      local_prior_group <- args[["prior_group"]]
      args[["prior_group"]] <- NULL
      
      # overwrite the golbal one if the local is specified
      if (!is.null(local_prior_group)) 
        prior_group <- local_prior_group
      
      prior_source <- args[["prior_source"]]
      args[["prior_source"]] <- NULL
      
      if (is.null(prior_source)) {
        prior_data <- parent_data
      } else {
        prior_data <- getData(y, prior_source)
      }
      
      if (gm == ".flowClust.1d") {
        # If any of 'K, 'neg' or 'pos' are given, we extract them from the
        # 'args' and coerce them as integers. Otherwise, they are defaulted to
        # NULL.
        # NOTE: the named elements within 'args' are coerced to lowercase.
        K <- neg_cluster <- pos_cluster <- NULL
        
        if ("k" %in% names(args)) {
          K <- as.integer(args["k"])
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
        
        args[["k"]] <- NULL
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
        
        if (!is.na(xChannel)) {
          prior_list[[xChannel]] <- .prior_flowClust(flow_set = prior_data,
              channels = xChannel, K = K, prior_group = prior_group,
              min = min_values, max = max_values, ...)
        }
        
        prior_list[[yChannel]] <- .prior_flowClust(flow_set = prior_data,
            channels = yChannel, K = K, prior_group = prior_group,
            min = min_values, max = max_values, ...)
        
        res <- prior_list
        
       
      } else {
        # get the value of neg and pos
        if (is.element(c("k"), names(args))) {
          K <- as.integer(args["k"])
        } else {
          message("'K' argument is missing! Using default setting: K = 2")
          K <- 2
        }
        
        args[["k"]] <- K
        names(args)[match("k", names(args))] <- "K"  #restore K to capital letter
        names(args)[match("useprior",names(args))]<-"usePrior"
        
        prior_list <- .prior_flowClust(flow_set = prior_data, channels = c(xChannel, 
                yChannel), K = K, prior_group = prior_group, ...)
        res <- prior_list
        # thisCall[['positive']]<-NULL #remove positive arg since 2D gate doesn't
        # understand it
      }
      
    }else{
      if(!(all(is.na(args)))){
        for (arg in names(args)) {
          thisCall[[arg]] <- args[[arg]]
        }  
      }
      
      ## choose serial or parallel mode
      if (parallel_type == "multicore") {
        message("Running in parallel mode with ", mc.cores, " cores.")
        thisCall[[1]] <- quote(mclapply)
        thisCall[["mc.cores"]] <- mc.cores
        flist <- eval(thisCall)
        
      }else if(parallel_type == "cluster"){
        if(is.null(cl))
          stop("cluster object 'cl' is empty!")
        thisCall[[1]] <- quote(parLapply)
        thisCall[["cl"]] <- cl
        # replace FUN with fun for parLapply
        thisCall["fun"] <- thisCall["FUN"]
        thisCall["FUN"] <- NULL
        flist <- eval(thisCall)
      }else {
        thisCall[[1]] <- quote(lapply)  #select loop mode
        flist <- eval(thisCall)
      }
      # Handles the case that 'flist' is a list of lists.
      #   The outer lists correspond to the split by pData factors.
      #   The inner lists contain the actual gates.
      #     The order do not necessarily match up with sampleNames()
      # Unforunately, we cannot simply use 'unlist' because the list element names
      # are mangled.
      if (all(sapply(flist, is.list))) {
        flist <- do.call(c, unname(flist))
      }
      
      if (extends(class(flist[[1]]), "fcFilter")) {
        flist <- fcFilterList(flist)
      } else {
        flist <- filterList(flist)
      }
      
      filterObj <- flist
      gs_node_id <- add(y, flist, parent = parent, name = popAlias)
      recompute(y, gs_node_id)
      
      message("done.")
    }

    
    
  } else {
    message("Skip preprocessing! Population '", paste(popAlias, collapse = ","), "' already exists.")
    res <- NULL
  }
  
  
  
  res
  
}