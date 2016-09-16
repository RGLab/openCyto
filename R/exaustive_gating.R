#' @include gating-methods.R
NULL

#' automated gating 
#' 
#' The idea is pretty simple: start with the marker with the best separation between -/+ peaks, 
#' then iterate over all other markers until you can't find a good separation between peaks. 
#' Different methods can be used for determining good/bad separation, such as the ICL or the
#' heuristic method based on peak distance.
#'  
#' @param gs a GatingSet
#' @param peak.separation 
#' @export
#' @examples 
#' dataDir <- system.file("extdata",package="flowWorkspaceData")
#' gs <- load_gs(list.files(dataDir, pattern = "gs_manual",full = TRUE))
#' 
#' #clone the gating tree strucuture of the existing manually gated gs(through flowJo)
#' gs2 <- clone(gs, isNew = FALSE, isEmpty = FALSE)
#' 
#' #erase all the manual gates and start from root
#' Rm("not debris", gs2)
#' 
#' getNodes(gs2)
#' 
#' gating(gs2)
#' @rdname gating
setMethod("gating", signature = c("GatingSet","missing"),
          definition = function(x, y, ...) {
            gating.leafnode(x,y,...)
          })
gating.leafnode <- function(gs, ...){
  
  #get used channels from existing gates
  #copied from flowWorkspace::dropRedundantChannels
  #TODO: separate this logic to be reused
  # nodes <- getNodes(gs, showHidden = TRUE)[-1]
  # gh <- gs[[1]]
  # chnls.used <- unlist(lapply(nodes, function(node){
  #   g <- getGate(gh, node)
  #   if(class(g) != "booleanFilter"){
  #     as.vector(parameters(g))  
  #   }
  #   
  # }))
  # chnls.used <- unique(chnls.used)
  # 
  # #remove the used chnls
  # chnls <- setdiff(chnls, chnls.used)
  
  nodes <- getNodes(gs)
  ##loop through all the existing terminal nodes
  isTerminal <- sapply(nodes,function(thisNode){
    length(getChildren(gs,thisNode))==0
  })
  nodes <- nodes[isTerminal]
  for(node in nodes){
  
    gating(node, gs, ...)
  }
}


#' @rdname gating
#' @examples 
#' # specify a particular existing node to start with
#' gs2 <- clone(gs, isNew = FALSE, isEmpty = FALSE)
#' Rm("CD3+", gs2)
#' gating("singlets", gs2, min.count = 2000,min.percent = 0.3)
setMethod("gating", signature = c("character", "GatingSet"),
          definition = function(x, y, ...) {
            gating.subnode(x,y,...)
          })
gating.subnode <- function(parent, gs
                           , gating.function = mindensity
                           , marker.selection = best.separation
                           , min.count = 100, min.percent = 0.1
                           , debug.mode = FALSE, max.depth = -1, ...){
  depths <- length(strsplit(getParent(gs, parent), "/")[[1]]) + 1 #ensure to get full path
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
      # if(debug.mode){
      #   message("start gating on node: ", parent)
      #   require(ggcyto)
      #   print(autoplot(fr))
      # } 
      # apply the gating on each channel
      flist <- sapply(channels, function(channel){
        marker <- getChannelMarker(fr, channel)[, "desc"]
        #avoid gating on the same marker repeately on the same path
        gated.markers <- strsplit(parent, split = "/")[[1]]
        is.gated <- sapply(gated.markers, function(i){
          i <- sub("[\\+\\-]", "", i)
          grepl(i, marker)
          })
        if(any(is.gated)){
          g <- NULL
        }else
          g <- gating.function(fr, channel, ...)
        g
      })
      flist <- flowWorkspace:::compact(flist)
      if(length(flist) > 0){
        message("parent: ", parent)
        #get measurements for the cutpoints
        ind <- best.separation(flist, fr, debug.mode = debug.mode, min.percent = min.percent)
        
        if(ind > 0){
          gate <- flist[[ind]]  
          channel <- as.vector(parameters(gate))
          marker <- getChannelMarker(fr, channel)[["desc"]]
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

#' Combine all three factors (peak height, valley location and height) to determine the best the peak separation 
#' 
#' @param flist a list of filters/gates generated by mindensity method
#' @param fr flowFrame
#' @param peak.ratio the minimum ratio of two peak heights in order for the gate to be considered
#' @param ... passed to density function
#' @return an integer as the index of the winner gate. When it returns 0, it indicates that no further gating should be proceeded
best.separation <- function(flist, fr, min.percent = 0.05, debug.mode = FALSE, ...){
  mat <- exprs(fr)
  plotObjs <- new.env(parent = emptyenv())
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
      plotObjs[[marker]] <- p
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
    data.table(marker, peak.dist, peak.ratio, valley.peak.ratio, area.ratio)
  })
  res <- rbindlist(res)
  
  #standardize each metrics
  metrics <- sapply(colnames(res), function(cn){
    
    col <- res[[cn]]
    if(cn%in%c("marker", "area.ratio")) #skip these two columns
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
  metrics[, valley.peak.ratio:= valley.peak.ratio] # reduce the weight of valley/peak ratio since it can be overweighted by the extreme low valley
  metrics[, score := rowMeans(.SD), .SDcols = c("valley.peak.ratio", "peak.dist", "peak.ratio", "area.ratio")]
  
  # metrics[peak.ratio < min.peak.ratio, score:=0] #automatically rule out the one that has the metrics below the theshold
  
  metrics[area.ratio < min.percent, score:=0] #automatically rule out the one that has the metrics below the theshold
  # metrics[, marker:=res[, marker]]
  # if(debug.mode){
  #   message("peak separation metrics:")
  #   print(metrics)  
  # }
  
  scores <- metrics[, score]
  if(all(scores==0))
    return(0)
  else{
    ind <- which.max(scores)
    select <- metrics[ind, marker]
    if(debug.mode){
      require(gridExtra)
      plotObjs <- as.list(plotObjs)
      
      for(this_marker in names(plotObjs)){
        # get metrics for the current marker
        # this_metrics <- metrics[marker == this_marker, ]
        # this_metrics <- round(this_metrics[, -1, with = FALSE], digits = 3)
        # # this_metrics <-  t(this_metrics)
        # # add metrics as table to top of the plot
        # tbl <- tableGrob(this_metrics, theme = ttheme_minimal(base_size = 5)
        #                  , rows = NULL
        #                  )
        
        p <- plotObjs[[this_marker]]
        # x.range <- layer_scales(p)[["x"]][["range"]][["range"]]
        # x.range <- ggplot_build(p)$panel$ranges[[1]]$x.range
        # x.diff <- diff(x.range)
        # p <- p + annotation_custom(tbl
        #                            # , ymin = max(p[["data"]][["y"]])
        #                            # , xmin = x.range[2] - x.diff * 0.2
        #                            # , xmax = x.range[2]
        #                            )
        
        #highlight the winner
        if(this_marker == select)
          p <- p + theme_dark()  
          
        
        plotObjs[[this_marker]] <- p
      }
      
      
      tbl <- round(metrics[, -1, with = FALSE], digits = 3)
      tbl[, marker := metrics[, marker]]
      # add metrics as table to top of the plot
      tbl <- tableGrob(tbl
                       , theme = ttheme_minimal(base_size = 7)
                       , rows = NULL
                       )
      plotObjs[["metrics"]] <- tbl
      plotObjs <- arrangeGrob(grobs = plotObjs)
      
      plot(plotObjs)
      # plot(gridExtra::arrangeGrob(plotObjs, tbl))
      
    }
    
    return(ind)
  }
    
}
