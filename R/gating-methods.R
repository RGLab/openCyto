setGeneric("gating", function(x, y, ...) standardGeneric("gating"))

.gating_gatingTemplate <- function(x, y, env_fct = NULL, ...) {
  gt <- x
  if (!is.null(env_fct)) {
    # use the fcTree if already exists
    if (exists("fct", env_fct)) {
      fct <- get("fct", env_fct)
    } else {
      # create one from gt if not
      fct <- fcTree(gt)
      assign("fct", fct, env_fct)
    }
    
  }
  
  # gate each node by the topological order
  gt_node_ids <- tsort(gt)
  
  # maintain the mapping between template node ID and gating set node ID in order
  # to refer gating set node ID back to the template ID and find the parent gs
  # node ID
  node_ids <- cbind(gt = gt_node_ids, gs = NA)
  node_ids[1, "gs"] <- 1  #fill out default gsid for root node
  
  for (i in 2:nrow(node_ids)) {
    
    # get parent node to gate
    gt_node_id <- node_ids[i, "gt"]
    gt_node_pop <- getNodes(gt, gt_node_id)
    # parent node in graph is used as reference node
    gt_ref_ids <- getParent(gt, gt_node_id)
    # the parent to be used in gs is from parent slot of pop object
    gt_parent_id <- as.character(gt_node_pop@parentID)
    
    # extract gate method from one edge(since multiple edge to the same node is
    # redudant)
    this_gate <- getGate(gt, gt_ref_ids[1], gt_node_id)
    
    parentInd <- match(gt_parent_id, node_ids[, "gt"])
    if (is.na(parentInd)) 
      stop("parent node '", names(getNodes(gt, gt_parent_id)), "' not gated yet!")
    gs_parent_id <- node_ids[parentInd, "gs"]
    # pass the pops and gate to gating routine
    res <- gating(x = this_gate, y, parent = as.integer(gs_parent_id), gtPop = gt_node_pop, 
      ...)
    gs_node_id <- res[["gs_node_id"]]
    filterObj <- res[["filterObj"]]
    
    # upodate gs node ids
    if (!is.null(gs_node_id)) 
      node_ids[i, "gs"] <- gs_node_id
    # update fct
    if (!is.null(env_fct) && !is.null(filterObj)) {
      nodeData(env_fct$fct, gt_node_id, "fList")[[1]] <- filterObj
    }
  }
  message("finished.")
}

setMethod("gating", signature = c("gatingTemplate", "GatingSetInternal"),
    definition = function(x, y, env_fct = NULL, ...) {
      .gating_gatingTemplate(x,y,env_fct,...)
    })

setMethod("gating", signature = c("gtMethod", "GatingSet"),
    definition = function(x, y, ...) {
      .gating_gtMethod(x,y,...)
    })

.gating_gtMethod <- function(x, y, gtPop, parent, num_cores = 1,
            parallel_type = c("multicore", "SOCK", "MPI"), plot = FALSE,
            xbin = 128, prior_group = NULL, ...) {
  
  require("parallel")
  
  args <- parameters(x)
  gm <- paste0(".", names(x))
  
  dims <- dims(x)
  xChannel <- unname(dims["xChannel"])
  yChannel <- unname(dims["yChannel"])
  is_1d_gate <- any(is.na(dims))
  
  popAlias <- alias(gtPop)
  popName <- names(gtPop)
  popId <- gtPop@id
  gs_nodes <- getChildren(y[[1]], getNodes(y[[1]], parent))
  
  if (length(gs_nodes) == 0 || !popAlias %in% gs_nodes) {
    message("Population '", popAlias, "'")
    
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
    thisCall[["FUN"]] <- as.symbol(gm)  #set gating method
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
        min_values <- max_values <- NULL
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

        args[["prior"]] <- prior_list
        
        # If the number of positive clusters is 0 and no cutpoint method has been
        # specified, we use the quantile gate by default.
        if (!is.null(pos_cluster) && pos_cluster == 0 && is.null(args[["cutpoint_method"]])) {
          args[["cutpoint_method"]] <- "quantile"
        }
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
        args[["prior"]] <- prior_list
        # thisCall[['positive']]<-NULL #remove positive arg since 2D gate doesn't
        # understand it
      }
      
    }
    # update arg_names
    for (arg in names(args)) {
      thisCall[[arg]] <- args[[arg]]
    }
    
    ## choose serial or parallel mode
    if (num_cores > 1) {
      message("Running in parallel mode with ", num_cores, " cores.")
      if (parallel_type == "multicore") {
        thisCall[[1]] <- quote(mclapply)
        thisCall[["mc.cores"]] <- num_cores
        flist <- eval(thisCall)
        
      } else {
        cl <- makeCluster(num_cores, type = parallel_type)
        thisCall[[1]] <- quote(parLapply)
        thisCall[["cl"]] <- cl
        # replace FUN with fun for parLapply
        thisCall["fun"] <- thisCall["FUN"]
        thisCall["FUN"] <- NULL
        flist <- eval(thisCall)
        stopCluster(cl)
      }
    } else {
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
    
  } else {
    message("Skip gating! Population '", paste(popAlias, collapse = ","), "' already exists.")
    gs_node_id <- getChildren(y[[1]], parent)
    # select the corresponding gs node id by matching the node names
    gs_node_name <- getNodes(y[[1]])[gs_node_id]
    gs_node_id <- gs_node_id[match(popAlias, gs_node_name)]
    filterObj <- NULL
  }
  
  if (plot) {
    print(plotGate(y, gs_node_id, xbin = xbin, pos = c(0.5, 0.8)))
  }
  
  list(gs_node_id = gs_node_id, filterObj = filterObj)
}

setMethod("gating", signature = c("boolMethod", "GatingSet"),
          definition = function(x, y, gtPop, parent, ...) {
  
  args <- parameters(x)[[1]]
  gm <- paste0(".", names(x))
  popAlias <- alias(gtPop)
  popName <- names(gtPop)
  popId <- gtPop@id
  
  gs_nodes <- getChildren(y[[1]], getNodes(y[[1]], parent))
  
  tNodes <- deparse(args)
  if (!(tNodes %in% gs_nodes)) {
    message(tNodes, " gating...")
    bf <- eval(substitute(booleanFilter(x), list(x = args)))
    bf@filterId <- tNodes
    invisible(gs_node_id <- add(y, bf, parent = parent))
    invisible(recompute(y, gs_node_id))
    message("done.")
  } else {
    message("Skip gating! Population '", tNodes, "' already exists.")
    gs_node_id <- getChildren(y[[1]], parent)

    # select the corresponding gs node id by matching the node names
    gs_node_name <- getNodes(y[[1]])[gs_node_id]
    gs_node_id <- gs_node_id[match(tNodes, gs_node_name)]
  }
  
  # gs_node_id
  list(gs_node_id = gs_node_id)
})

setMethod("gating", signature = c("polyFunctions", "GatingSet"),
          definition = function(x, y, gtPop, parent, ...) {
  
  refNodes <- x@refNodes
  gm <- paste0(".", names(x))
  popAlias <- alias(gtPop)
  popName <- names(gtPop)
  
  gs_nodes <- getChildren(y[[1]], getNodes(y[[1]], parent))
  
  message("Population '", paste(popAlias, collapse = ","), "'")
  
  nMarkers <- length(refNodes)

  ## all the comibnations of A & B & C
  opList <- permutations(n = 1, r = nMarkers - 1, c("&"), repeats = TRUE)
  isNotList <- permutations(n = 2, r = nMarkers, c("!", ""), repeats = TRUE)
  polyExprsList <- apply(opList, 1, function(curOps) {
    apply(isNotList, 1, function(curIsNot) {
      polyExprs <- curIsNot
      polyExprs[-1] <- paste0(curOps, curIsNot[-1])
      
      paste(paste0(polyExprs, refNodes), collapse = "")
    })
  })
  polyExprsList <- as.vector(polyExprsList)
  
  # actual gating
  lapply(polyExprsList, function(polyExpr) {
    bgt <- new("boolMethod", name = polyExpr, args = list(as.symbol(polyExpr)))
    gating(bgt, y, parent = parent, gtPop = gtPop, ...)
  })
  
  message("done.")
  
  list()
})
setMethod("gating", signature = c("refGate", "GatingSet"),
    definition = function(x, y, ...) {
      .gating_refGate(x, y, ...)
    })
.gating_refGate <- function(x, y, gtPop, parent, plot = FALSE, xbin = 128,
            ...) {
  
  refNodes <- x@refNodes
  gm <- paste0(".", names(x))
  popAlias <- alias(gtPop)
  popName <- names(gtPop)
  dims <- dims(x)
  xChannel <- dims[["xChannel"]]
  yChannel <- dims[["yChannel"]]
  
  gs_nodes <- getChildren(y[[1]], getNodes(y[[1]], parent))
  
  if (length(gs_nodes) == 0 || !popAlias %in% gs_nodes) {
    
    message("Population '", paste(popAlias, collapse = ","), "'")
    if (length(refNodes) > 2) {
      stop("Not sure how to construct gate from more than 2 reference nodes!")
    }
    
    fr <- getData(y[[1]])
    
    flist <- flowWorkspace::lapply(y, function(gh) {
          
       glist <- lapply(refNodes, function(refNode) {
          node_names <- getNodes(gh)
          node_ind <- match(refNode, node_names)
         if (is.na(node_ind)) {
            # match to path
            node_paths <- getNodes(gh, isPath = T)
            toMatch <- gsub("\\+", "\\\\+", refNode)
            toMatch <- paste(toMatch, "$", sep = "")
            node_ind <- grep(toMatch, node_paths)
            if (length(node_ind) == 0) {
            stop(refNode, " not found in gating set!")
            } else if (length(node_ind) > 1) {
            stop("Multiple ", refNode, " found in gating set!")
            }
          }
          getGate(gh, node_ind)
        })

      # standardize the names for the gate parameters and dims
      gate_params <- unlist(lapply(glist, function(g) {
        cur_param <- parameters(g)
        getChannelMarker(fr, cur_param)["name"]
      }))

      if(length(glist)==1){
        #1d ref gate
        dims <- dims[!is.na(dims)]
        nDims <- length(dims)
        if(nDims!=1){
          stop("Can't do 1d gating on ",nDims, " dimensions!")
        }
        dim_params <-  getChannelMarker(fr, dims)["name"]
        y_g <- glist[[1]]               
        y_coord <- c(y_g@min, y_g@max)
        cut.y <- y_coord[!is.infinite(y_coord)]
        
        pos_token <- "[\\+]"
        neg_token <- "[\\-]"
        pop_name_pat <- "[^\\+-]+"
        pos_pop_pat <- paste(pop_name_pat, pos_token, sep = "")
        neg_pop_pat <- paste(pop_name_pat, neg_token, sep = "")
        if (grepl(pos_pop_pat,popName)) {
          gate_coordinates <- list(c(cut.y, Inf))
        } else if(grepl(neg_pop_pat,popName)){
          gate_coordinates <- list(c(-Inf, cut.y))
        }else{
          stop("unknown population pattern, ",popName)
        }
        names(gate_coordinates) <- as.character(dim_params)
        
        rectangleGate(gate_coordinates)
      }else{
          #2d quad gate
          dim_params <- unlist(lapply(dims, function(dim) {
                    getChannelMarker(fr, dim)["name"]
                  }))
          
          # match the gate param to dims to find gates for x, y dimensions
          x_g <- glist[[match(dim_params[1], gate_params)]]
          y_g <- glist[[match(dim_params[2], gate_params)]]
          
          # pick the non-infinite coordinate as the cut points
          x_coord <- c(x_g@min, x_g@max)
          y_coord <- c(y_g@min, y_g@max)
          cut.x <- x_coord[!is.infinite(x_coord)]
          cut.y <- y_coord[!is.infinite(y_coord)]
          
          # In order, the following vector has the patterns:
          # 1. top left (-+)
          # 2. top right (++)
          # 3. bottom right (+-)
          # 4. bottom left (--)
          quadPatterns <- c(".+-.+\\+$", ".+\\+.+\\+$", ".+\\+.+-$", ".+-.+-$")
          
          # check if popname is give as Y[*]X[*]
          YX_pattern <- paste0(dims["yChannel"], ".+", dims["xChannel"], ".+")
          XY_pattern <- paste0(dims["xChannel"], ".+", dims["yChannel"], ".+")
          
          # do the flipping if YX
          if (grepl(YX_pattern, popName)) {
            pos <- regexpr(dims["xChannel"], popName)
            xterm <- substring(popName, pos, nchar(popName))
            yterm <- substring(popName, 1, pos - 1)
            toMatch <- paste(xterm, yterm, sep = "")
          } else if (grepl(XY_pattern, popName)) {
            toMatch <- popName  #keep as it is if XY
          }
          else {
            stop("X,Y axis do not match between 'dims'(", paste(dims, collapse = ","), 
                ") and 'pop'(", popName, ")")
          }
          quadInd <- which(unlist(lapply(quadPatterns, grepl, toMatch)))
          
          # construct rectangleGate from reference cuts
          if (quadInd == 1) {
            coord <- list(c(-Inf, cut.x), c(cut.y, Inf))
          } else if (quadInd == 2) {
            coord <- list(c(cut.x, Inf), c(cut.y, Inf))
          } else if (quadInd == 3) {
            coord <- list(c(cut.x, Inf), c(-Inf, cut.y))
          } else if (quadInd == 4) {
            coord <- list(c(-Inf, cut.x), c(-Inf, cut.y))
          } else stop("Pop names does not match to any quadrant pattern!")
          
          names(coord) <- as.character(dim_params)
          rectangleGate(coord)
      }
      
    })
    
    flist <- filterList(flist)
    gs_node_id <- add(y, flist, parent = parent, name = popAlias)
    recompute(y, gs_node_id)
    
    if (plot) {
      print(plotGate(y, gs_node_id, xbin = xbin, pos = c(0.5, 0.8)))
    }
  } else {
    message("Skip gating! Population '", popAlias, "' already exists.")
    gs_node_id <- getChildren(y[[1]], parent)

    # select the corresponding gs node id by matching the node names
    gs_node_name <- getNodes(y[[1]])[gs_node_id]
    gs_node_id <- gs_node_id[match(popAlias, gs_node_name)]
    flist <- NULL
  }

  message("done.")
  list(gs_node_id = gs_node_id, filterObj = flist)
}

# wrapper for prior_flowClust to return sample-specific priors
.prior_flowClust <- function(flow_set, prior_group = NULL, ...) {
  if (is.null(prior_group)) {
    prior_list <- prior_flowClust(flow_set, ...)
    # replicate the prior for each sample
    sapply(sampleNames(flow_set), function(i) prior_list, simplify = FALSE)
  } else {
    splitBy1 <- factor(pData(flow_set)[, prior_group])
    fs_list <- split(flow_set, splitBy1)
    names(fs_list) <- NULL  #strip group name
    plist <- lapply(fs_list, function(fs) {
      prior_list <- prior_flowClust(fs, ...)

      # replicate the prior for each sample within this prior group
      sapply(sampleNames(fs), function(i) prior_list, simplify = FALSE)
    })
    unlist(plist, recursive = FALSE)
  }
}

## wrappers for the different gating routines
.singletGate <- function(fs, xChannel = "FSC-A", yChannel = "FSC-H",
                         prediction_level = 0.99, ...) {
  require(openCyto)
  # Creates a list of polygon gates based on the prediction bands at the minimum
  # and maximum x_channel observation using a robust linear model trained by
  # flowStats.
  singletGate(fs[[1]], area = xChannel, height = yChannel,
              prediction_level = prediction_level)
}

.flowClust.1d <- function(fs, xChannel = NA, yChannel, tol = 1e-5, prior = NULL,
                          filterId = "", split = TRUE, ...) {
  require(openCyto)

  sname <- sampleNames(fs)
  fr <- fs[[sname]]
  priorList <- list()
  priorList[[yChannel]] <- prior[[yChannel]][[sname]]
  
  if (is.na(xChannel)) {
    # 1d gate
    flowClust.1d(fr = fr, params = yChannel, tol = tol, filterId = filterId, 
      prior = priorList[[yChannel]], ...)
  } else {
    stop("flowClust1d does not support 2d gate!")
  }
}

.cytokine <- function(fs, xChannel = NA, yChannel = "FSC-A", filterId = "",
                      ...) {
  require(openCyto)
  cytokine(flow_set = fs, channel = yChannel, filter_id = filterId, ...)
}

.mindensity <- function(fs, yChannel = "FSC-A", filterId = "", ...) {
  require(openCyto)
  # TODO: Iterate through the flowFrames within 'fs', given that 'fs' may
  # contain more than one flowFrame if 'split' is specified in the CSV file.
  mindensity(flow_frame = fs[[1]], channel = yChannel, filter_id = filterId, ...)
}

.flowClust.2d <- function(fs, xChannel, yChannel, usePrior = "yes", prior = NULL,
                          ...) {
  require(openCyto)
  # TODO: Iterate through the flowFrames within 'fs', given that 'fs' may
  # contain more than one flowFrame if 'split' is specified in the CSV file.
  sname <- sampleNames(fs)
  fr <- fs[[sname]]
  flowClust.2d(fr = fr, xChannel = xChannel, yChannel = yChannel, usePrior = usePrior,
               prior = prior[[sname]], ...)
}

.rangeGate <- function(fs, xChannel = NA, yChannel, absolute = FALSE, filterId = "", 
                       ...) {
  require(openCyto)
  # TODO: Iterate through the flowFrames within 'fs', given that 'fs' may
  # contain more than one flowFrame if 'split' is specified in the CSV file.
  fr <- fs[[1]]
  rangeGate(x = fr, stain = yChannel, inBetween = TRUE, absolute = absolute,
            filterId = filterId, ...)
}

.quantileGate <- function(fs, xChannel = NA, yChannel, probs = 0.999, filterId = "",
                          ...) {
  require(openCyto)
  # TODO: Iterate through the flowFrames within 'fs', given that 'fs' may
  # contain more than one flowFrame if 'split' is specified in the CSV file.
  fr <- fs[[1]]
  quantileGate(fr = fr, probs = probs, stain = yChannel, filterId = filterId, ...)
}

.quadrantGate <- function(fs, xChannel = NA, yChannel, ...) {
  require(openCyto)
  # TODO: Iterate through the flowFrames within 'fs', given that 'fs' may
  # contain more than one flowFrame if 'split' is specified in the CSV file.
  fr <- fs[[1]]
  
  qfilter <- quadrantGate(fr, stain = c(xChannel, yChannel), absolute = FALSE, 
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
 
