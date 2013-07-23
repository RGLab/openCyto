setGeneric("gating", function(x, y, ...) standardGeneric("gating"))
#' Applies gatingTemplate to 1 GatingSetInternal.
#'
#'
#'
#' @param x a \code{gatingTemplate} object
#' @param y a \code{GatingSetInternal} object
#' @param env_fct a \code{environment} that contains \code{fcTree} object named as 'fct'
#' @param stop.at a \code{character} that specifies the population (correspoding to 'alias' column in csv template) where the gating prcoess will stop at.

setMethod("gating", signature = c("gatingTemplate", "GatingSetInternal"),
    definition = function(x, y, env_fct = NULL, ...) {
      .gating_gatingTemplate(x,y,env_fct,...)
    })

.gating_gatingTemplate <- function(x, y, env_fct = NULL, stop.at = NULL, ...) {
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
  #validity check for stop.at argument
  if(!is.null(stop.at)){
    if(is.na(match(stop.at,sapply(getNodes(gt),alias))))
      stop("Can't find stop point: ", stop.at)
  }
  # gate each node 
#  gt_node_ids <- tsort(gt)#by the topological order
  gt_node_ids <- bfs(gt)#by the bfs order
  # maintain the mapping between template node ID and gating set node ID in order
  # to refer gating set node ID back to the template ID and find the parent gs
  # node ID
  node_ids <- cbind(gt = gt_node_ids, gs = NA)
  node_ids[1, "gs"] <- 1  #fill out default gsid for root node
  
  for (i in 2:nrow(node_ids)) {
    
    # get parent node to gate
    gt_node_id <- node_ids[i, "gt"]
    gt_node_pop <- getNodes(gt, gt_node_id)
    if(!is.null(stop.at)){
    
      if(alias(gt_node_pop) == stop.at)
      {
        message("stop at: ",stop.at)
        break
      }
        
    }
#    browser()
    # parent node in graph is used as reference node
    gt_ref_ids <- getParent(gt, gt_node_id)
    # the parent to be used in gs is from parent slot of pop object
    gt_parent_id <- as.character(gt_node_pop@parentID)
    
    # extract gate method from one edge(since multiple edge to the same node is
    # redudant)
    this_gate <- getGate(gt, gt_ref_ids[1], gt_node_id)
    
    #get preprocessing method
    this_ppm <- ppMethod(gt, gt_ref_ids[1], gt_node_id)
    
    parentInd <- match(gt_parent_id, node_ids[, "gt"])
    if (is.na(parentInd)) 
      stop("parent node '", names(getNodes(gt, gt_parent_id)), "' not gated yet!")
    gs_parent_id <- node_ids[parentInd, "gs"]
    if (is.na(gs_parent_id)) 
      stop("parent node '", names(getNodes(gt, gt_parent_id)), "' not gated yet!")
    
    #preprocessing
    pp_res <- NULL
#    browser()
    if(class(this_ppm) == "ppMethod")
      pp_res <- preprocessing(x = this_ppm, y, parent = as.integer(gs_parent_id), gtPop = gt_node_pop, gm = this_gate, ...)
    
    # pass the pops and gate to gating routine
    res <- gating(x = this_gate, y, parent = as.integer(gs_parent_id), gtPop = gt_node_pop, pp_res = pp_res, ...)
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



setMethod("gating", signature = c("gtMethod", "GatingSet"),
    definition = function(x, y, ...) {
      .gating_gtMethod(x,y,...)
    })

.gating_gtMethod <- function(x, y, gtPop, parent, pp_res 
            , mc.cores = 1, parallel_type = c("none", "multicore", "cluster"), cl = NULL
            , plot = FALSE, xbin = 128,  ...) {
  
  require("parallel")
#  browser()
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
#  browser()
  if (length(gs_nodes) == 0 || !popAlias %in% gs_nodes) {
    message("Gating for '", popAlias, "'")
    
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
    groupBy <- groupBy(x)
    if (groupBy != "" && x@collapse) { #when x@collapse == FALSE,then ignore groupBy argument since grouping is only used for collapsed gating
      
      split_by <- as.character(groupBy)
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
        split_by <- apply(pData(parent_data)[, split_by, drop = FALSE], 1, paste, collapse = ":")
        split_by <- as.character(split_by)
      }
    } else {
      split_by <- sampleNames(parent_data)
    }

    fslist <- split(parent_data, split_by)
    
    if(is.null(pp_res))
      pp_res <- sapply(names(fslist),function(i)pp_res)
    # construct method call
    thisCall <- substitute(f1(fslist,pp_res))
    thisCall[["FUN"]] <- as.symbol(".gating_wrapper")
    args[["gFunc"]] <- as.symbol(gm)  #set gating method
    
    args[["xChannel"]] <- xChannel  #set x,y channel
    args[["yChannel"]] <- yChannel
    
    if (is_1d_gate) {
      if (grepl("-$", popName)) {
        positive <- FALSE
      } else if (grepl("+$", popName)) {
        positive <- TRUE
      } else {
        stop("Invalid population name! Name should end with '+' or '-' symbol.")
      }
      
      args[["positive"]] <- positive
    }

    
    thisCall[["MoreArgs"]] <- args
    
    ## choose serial or parallel mode
    
      
    if (parallel_type == "multicore") {
      message("Running in parallel mode with ", mc.cores, " cores.")
      thisCall[[1]] <- quote(mcmapply)
      thisCall[["mc.cores"]] <- mc.cores
    }else if(parallel_type == "cluster"){
      if(is.null(cl))
          stop("cluster object 'cl' is empty!")
        thisCall[[1]] <- quote(clusterMap)
        thisCall[["cl"]] <- cl
     }else {
      thisCall[[1]] <- quote(mapply)  #select loop mode
  
     }
#browser()         
    flist <- eval(thisCall)

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
    definition = function(x, y, ...) {
      .gating_boolMethod(x,y,...)      
    })

.gating_boolMethod <- function(x, y, gtPop, parent, ...) {
  
  args <- parameters(x)[[1]]
  gm <- paste0(".", names(x))
  popAlias <- alias(gtPop)
  popName <- names(gtPop)
  popId <- gtPop@id
  
  gs_nodes <- getChildren(y[[1]], getNodes(y[[1]], parent))
  
  tNodes <- deparse(args)
  if (!(popAlias %in% gs_nodes)) {
    message(tNodes, " gating...")
    bf <- eval(substitute(booleanFilter(x), list(x = args)))
    bf@filterId <- tNodes
    invisible(gs_node_id <- add(y, bf, parent = parent, name = popAlias))
    invisible(recompute(y, gs_node_id))
    message("done.")
  } else {
    message("Skip gating! Population '", popAlias, "' already exists.")
    gs_node_id <- getChildren(y[[1]], parent)

    # select the corresponding gs node id by matching the node names
    gs_node_name <- getNodes(y[[1]])[gs_node_id]
    gs_node_id <- gs_node_id[match(popAlias, gs_node_name)]
  }
  
  # gs_node_id
  list(gs_node_id = gs_node_id)
}
setMethod("gating", signature = c("polyFunctions", "GatingSet"),
    definition = function(x, y, ...) {
      .gating_polyFunctions(x,y,...)      
    })

.gating_polyFunctions <- function(x, y, gtPop, parent, ...) {
  
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
    alias(gtPop) <- polyExpr
    gating(bgt, y, parent = parent, gtPop = gtPop, ...)
  })
  
  message("done.")
  
  list()
}
setMethod("gating", signature = c("refGate", "GatingSet"),
    definition = function(x, y, ...) {
      .gating_refGate(x, y, ...)
    })
.gating_refGate <- function(x, y, gtPop, parent, plot = FALSE, xbin = 128,
            ...) {
#  negated <- FALSE
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
        
        pos_token <- "[\\+]"
        neg_token <- "[\\-]"
        pop_name_pat <- "[^\\+-]+"
        pos_pop_pat <- paste(pop_name_pat, pos_token, sep = "")
        neg_pop_pat <- paste(pop_name_pat, neg_token, sep = "")
        
        if(nDims==2){
          #set negated flag of outer scope when necessary
#          if (grepl(neg_pop_pat,popName)&&!negated) {
#            negated <<- TRUE
#          }
          #before flowCore and flowViz support the negated filter
          #we use refGate+boolGate in csv template as the workaround
          if (grepl(neg_pop_pat,popName)) {
            stop("negated 2d gate is not supported yet!")
          }
          #pass the gate as it is 
          glist[[1]]
        }else{
          dim_params <-  getChannelMarker(fr, dims)["name"]
          y_g <- glist[[1]]               
          y_coord <- c(y_g@min, y_g@max)
          cut.y <- y_coord[!is.infinite(y_coord)]
          
          
          if (grepl(pos_pop_pat,popName)) {
            gate_coordinates <- list(c(cut.y, Inf))
          } else if(grepl(neg_pop_pat,popName)){
            gate_coordinates <- list(c(-Inf, cut.y))
          }else{
            stop("unknown population pattern, ",popName)
          }
          names(gate_coordinates) <- as.character(dim_params)
          
          rectangleGate(gate_coordinates)  
        }
        
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

