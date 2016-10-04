#' @include gating-methods.R
NULL

#' automated gating 
#' 
#' The idea is pretty simple: start with the marker with the best separation between -/+ peaks, 
#' then iterate over all other markers until you can't find a good separation between peaks. 
#' Different methods can be used for determining good/bad separation, such as the ICL or the
#' heuristic method based on peak distance.
#'  
#' @param x a GatingSet
#' @param y either missing or a character specifying the starting node
#' @export
#' @examples 
#' \dontrun{
#' dataDir <- system.file("extdata",package="flowWorkspaceData")
#' gs <- load_gs(list.files(dataDir, pattern = "gs_manual",full = TRUE))
#' 
#' #clone the gating tree strucuture of the existing manually gated gs(through flowJo)
#' gs2 <- clone(gs, isNew = FALSE, isEmpty = FALSE)
#' 
#' #erase all the acestor nodes of CD3+
#' Rm("CD3+", gs2)
#' 
#' # run the exaustive gating from "singlets" node
#' gating("singlets", gs2, min.count = 2000, min.percent = 0.3)
#' 
#' # or proceed from any leaf nodes of the existing gating tree
#' gating(gs2, debug.mode = TRUE) #turn on the debug mode to generate plots and messages (particularly in knitr report) for troubleshooting
#' 
#' }
#' @rdname gating
setMethod("gating", signature = c("GatingSet","missing"),
          definition = function(x, y, ...) {
            
            gating.leafnodes(x,y,...)
          })

init.openCyto.exaustive <- function(gsid){
  
  ops <- getOption("openCyto")
  if(!"exaustive"%in%names(ops))
    ops[["exaustive"]] <- list()
  if(!gsid%in%names(ops[["exaustive"]]))
    ops[["exaustive"]][[gsid]] <- new.env(parent = emptyenv())
  options("openCyto" = ops)
  
}
set.openCyto.exaustive <- function(gsid, node, winner, plotEnv, metrics){
  if(is.null(node)||is.null(gsid))
    stop("Can't record the marker selection process because gs id and node are not set!")
  
  #init global option if needed
  init.openCyto.exaustive(gsid)
  
  #annotate the winner plot
  if(!is.null(winner)) {
    plotEnv[[winner]] <- plotEnv[[winner]] + theme_dark()  
  }
  
  ops <- getOption("openCyto")
  ops[["exaustive"]][[gsid]][[node]] <- list(winner = winner, plotEnv = plotEnv, metrics = metrics)
}
gating.leafnodes <- function(gs, ...){
  
  
  nodes <- getLeafNode(gs)
  for(node in nodes){
    gating(node, gs, ...)
  }
}



setMethod("gating", signature = c("character", "GatingSet"),
          definition = function(x, y, ...) {
            gating.subnode(x,y,...)
          })

#' @rdname gating
#' @param marker.selection a function that selectes the best marker that has the best bi-module separation based on the various peak statistics (e.g. peak area ratio, peak distance, peak/valley ratio, etc...)
#' @param gating.function the 1d gating function used for exaustive gating
#' @param min.count the minimum number of cells that allows the gating proceed further
#' @param max.depth the maximum depths of gating path. Default is -1, which is no limits.
#' @param ... other arguments passed to \link{density} function.
gating.subnode <- function(parent, gs
                           , gating.function = mindensity
                           , marker.selection = best.separation
                           , min.count = 1000, min.percent = 0.2
                           , debug.mode = FALSE, max.depth = -1
                           , mc.cores = getOption("mc.cores", 2L)
                           , parallel_type = c("none", "multicore", "cluster"), cl = NULL
                           , ...){
  parallel_type <- match.arg(parallel_type)
  parent <- file.path(getParent(gs, parent), basename(parent))#ensure to get full path
  depths <- length(strsplit(parent, "/")[[1]])
  if(max.depth >0 && depths >= max.depth)
    message("stop gating at ", parent, ". Reaching the maximum gating depths: ", max.depth)
  else
  {
    fs <- getData(gs, parent)
    #for now we use the first sample
    fr <- fs[[1, use.exprs = FALSE]]
    #exclude the non-stained channels
    pd <- pData(parameters(fr))
    pd <- pd[!is.na(pd[["desc"]]),]
    channels <- pd[["name"]]
    fr <- fs[[1,channels]]
    
    nCell <- nrow(fr)
    if(nCell > min.count){#TODO: use options("openCyto")[[1]][["gating"]][["minEvents"]]
      # apply the gating on each channel
      
      # construct method call
      
      thisFunc <- function(channel){
                            marker <- getChannelMarker(fr, channel)[, "desc"]
                            #avoid gating on the same marker repeately on the same path
                            gated.markers <- strsplit(parent, split = "/")[[1]]
                            gated.markers <- gated.markers[-1] #rm the first empty string
                            is.gated <- sapply(gated.markers, function(i){
                              i <- sub("[\\+\\-]$", "", i)
                              grepl(i, marker)
                              })
                            if(any(is.gated)){
                              g <- NULL
                            }else{
                              g <- gating.function(fr, channel, ...)
                            }
                              
                            g
                          }
      
      thisCall <- substitute(sapply(channels))
      thisCall[["FUN"]] <- as.symbol("thisFunc")
      
      
      ## choose serial or parallel mode

      if (parallel_type == "multicore") {
        message("Running in parallel mode with ", mc.cores, " cores.")
        thisCall[[1]] <- quote(mcmapply)
        thisCall[["mc.cores"]] <- mc.cores
        # thisCall[["mc.preschedule"]] <- FALSE #must set to FALSE to retain the names of the result list
        
      }else if(parallel_type == "cluster"){
        if(is.null(cl))
          stop("cluster object 'cl' is empty!")
        thisCall[[1]] <- quote(parSapply)
        thisCall[["cl"]] <- cl
        # thisCall[["fun"]] <- thisCall[["FUN"]] 
        # thisCall[["FUN"]] <- NULL
        thisCall[["SIMPLIFY"]] <- FALSE
      }
      
      flist <- eval(thisCall)
      flist <- flowWorkspace:::compact(flist)
      if(length(flist) > 0){
        message("parent: ", parent)
        #get measurements for the cutpoints
        plotEnv <- new.env(parent = emptyenv())
        metrics <- marker.selection(flist, fr, debug.mode = debug.mode, plotEnv = plotEnv)
        metrics <- metrics[area.ratio >= min.percent, ] 
          
        
        if(nrow(metrics)==0){
          chnl.selected <- NULL
        }else{
          ind <- which.max(metrics[, score])
          chnl.selected <- metrics[ind, channel]  
        }
        if(debug.mode){
          set.openCyto.exaustive(gs@guid, parent, winner = chnl.selected, plotEnv = plotEnv, metrics = metrics)
          }
        
        if(!is.null(chnl.selected)){
          gate <- flist[[chnl.selected]]  
          # channel <- as.vector(parameters(gate))
          marker <- getChannelMarker(fr, chnl.selected)[["desc"]]
          #clean the marker name(sometime it is in the form of 'antigen isotypecontrol',e.g 'CD38 APC')
          marker <- strsplit(marker, " ")[[1]][[1]]
          message("selected marker: ", marker)
          #add the gates and move on the children of the new node
          for(sub in c("+", "-")){
            negated <- sub == "-"
            child <- paste0(marker, sub)
            suppressMessages(add(gs, gate, name = child, parent = parent, negated = negated)) 
            node.path <- file.path(parent,child)
            suppressMessages(recompute(gs, node.path))
            #resursively gate the sub nodes
            gating(node.path, gs
                   , gating.function = gating.function
                   , marker.selection = marker.selection
                   , min.count = min.count
                   , min.percent = min.percent
                   , debug.mode = debug.mode
                   , max.depth = max.depth
                   , ...)  
          }
        }  
      }
     
    }else{
      if(debug.mode)
        message("skip node '", parent, "' due to the low cell count : ", nCell)
    }
  }
  
}

#' A function that performs the marker selection.
#' 
#' It combines all the following factors to determine which marker has the best the peak separation 
#' 1. distance between two peaks
#' 2. height ratio between two peaks
#' 3. valley to peak(lower of the two peaks) height ratio
#' 4. area ratio between two peaks
#' 
#' @param flist a list of filters/gates generated by mindensity method
#' @param fr flowFrame
#' @param min.percent the minimum ratio of two peak areas. The marker will be disqualified automatically when its measurement is lower than this threshold
#' @param ... passed to density function
#' @return an integer as the index of the winner gate. When it returns 0, it indicates that no further gating should be proceeded
#' @export
#' @examples 
#' data(GvHD)
#' fr <- GvHD[[1]]
#' 
#'  #  transform the raw data
#' channels <- c("FL1-H", "FL2-H", "FL3-H", "FL4-H")
#' trans <- estimateLogicle(fr, channels)
#' fr <- transform(fr, trans)
#' 

#' #run the 1d gating on these channels
#' flist <- sapply(channels, function(channel)mindensity(fr, channel))
#' 
#' #select the best marker
#' ind <- best.separation(flist, fr)
#' channels[ind] #winner channel
#' 
#' objs <- best.separation(flist, fr, debug.mode = TRUE) #turn on the debug mode to display plots
#' objs[["metrics"]]
#' plot(objs[["plotObjs"]])
best.separation <- function(flist, fr, debug.mode = FALSE, plotEnv = new.env(parent = emptyenv()), ...){
  mat <- exprs(fr)
  
  #calculate the peak separation measurements
  res <- lapply(flist, function(gate){
    channel <- as.vector(parameters(gate))
    marker <- getChannelMarker(fr, channel)[, "desc"]
    valley.x <- gate@min #assuming it is always positive gate
    vec <- mat[, channel]
    peaks.result <- .find_peaks(vec, plot = FALSE,...)
    peaks.info <- peaks.result[["peaks"]]
    dens <- peaks.result[["dens"]]
    valley.ind <- which.min(abs(dens[, x] - valley.x))
    valley.y <- dens[valley.ind, y]
    #get the max peak from either side of valley
    left.peak <- peaks.info[x <= valley.x, ][1, ]
    right.peak <- peaks.info[x > valley.x, ][1, ]
    
    if(debug.mode){
      require(ggplot2)
      p <- ggplot(dens, aes(x = x, y)) + geom_line()
      p <- p + geom_point(data = left.peak, col = "red")
      p <- p + geom_point(data = right.peak, col = "blue")
      p <- p + geom_point(data = data.frame(x = valley.x, y = valley.y), pch = "V", col = "brown", cex = 3)
      p <- p + xlab(marker) + theme(axis.text.y= element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank())
      plotEnv[[channel]] <- p
    }
    
    if(is.na(left.peak[,x])||is.na(right.peak[, x]))
    {
      peak.dist <- 0
      peak.ratio <- 0
      valley.peak.ratio <- 1
      area.ratio <- 0
    }else{
      
      left.area <- dens[1:valley.ind, sum(y)]
      total.area <- dens[, sum(y)]
      right.area <- total.area - left.area
      area.ratio <- min(left.area, right.area)/total.area
      peak.dist <- (right.peak[, x] - left.peak[, x])/diff(range(vec))
      
      peak.height <- c(left.peak[, y], right.peak[, y])
      peak.height.min <- min(peak.height)
      peak.height.max <- max(peak.height)
      peak.ratio <- peak.height.min/peak.height.max
      
      valley.peak.ratio <- valley.y/peak.height.min
    }
    data.table(peak.dist, peak.ratio, valley.peak.ratio, area.ratio)
  })
  res <- rbindlist(res, idcol = "channel")
  
  #standardize each metrics
  metrics <- sapply(colnames(res), function(cn){
    
    col <- res[[cn]]
    if(cn%in%c("channel", "area.ratio")) #skip these two columns
      val <- col
    else{
      min.col <- min(col)
      delta <- (max(col) - min.col)
      if(delta > 0)
        val <- (col- min.col)/delta
      else
        val <- rep(0, length(col))  
    }
    
    return(val)
  }, simplify = FALSE)
  metrics <- as.data.table(metrics)
  metrics[, valley.peak.ratio:= 1- valley.peak.ratio] # invert v/p score
  # metrics[, valley.peak.ratio:= valley.peak.ratio] # reduce the weight of valley/peak ratio since it can be overweighted by the extreme low valley
  metrics[, score := rowMeans(.SD), .SDcols = c("valley.peak.ratio", "peak.dist", "peak.ratio", "area.ratio")]
  

  return(metrics)  
}

best.ICL <- function(flist, fr, debug.mode = FALSE, plotEnv = new.env(parent = emptyenv()), ...){
  res <- sapply(flist, function(gate){
    res <- posteriors(gate)[[1]]
    score <- res[["ICL"]]
    
    gate_pct <- res[["gate_pct"]]
    area.ratio <- ifelse(gate_pct > 0.5, 1 - gate_pct, gate_pct)
    data.table(score, area.ratio)
  }, simplify = FALSE)
  metrics <- rbindlist(res, idcol = "channel")
    
  
  if(debug.mode){
    #save plots and metrics to global hash table
    for(cn in colnames(fr)){
      p <- autoplot(fr,cn) + labs_cyto("marker") + labs(title = NULL)
      plotEnv[[cn]] <- as.ggplot(p)
    }
      
  }
   
  return(metrics)
  
}
