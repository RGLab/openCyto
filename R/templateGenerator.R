#' generate a partially complete csv template from the existing gating hierarchy 
#' 
#' To ease the process of replicating the existing (usually a manual one) gating schemes, 
#' this function populate an empty gating template with the 'alias', 'pop', 'parent' and 'dims' 
#' columns that exacted from an \code{GatingHierarchy}, and leave the other columns (e.g. `gating_method`) blank.
#' So users can make changes to that template instead of writing from scratch.
#' 
#' @param gh a \code{GatingHierarchy} likely parsed from a xml workspace
#' @return a gating template in \code{data.frame} format that requires further edition after output to csv 
#' @export 
templateGen <- function(gs, ...){
  nodes <- getNodes(gs, order = "tsort", path = "auto")[-1]
  gh <- gs[[1]]
  mergedNodes <- flowWorkspace:::.mergeGates(gh, nodes, bool = FALSE, merge = TRUE, projections = list())
  names(mergedNodes) <- NULL
  ldply(mergedNodes, function(nodes){
         if(is.list(nodes))
           getGatingMethod.multiPopulations(nodes, gs, ...)
         else
           singlePopulation.singlePopulation(as.character(nodes), gs, ...)
  })

  
}

#' single population
getGatingMethod.singlePopulation <- function(node, gs, ...){
   
  alias <- basename(node)
  pop <- alias
  gh <- gs[[1]]
  parent <- getParent(gh, node, path = "auto")
  
  thisGate <- getGate(gh, node)
  
  gateType <- class(thisGate)
  if(gateType == "rectangleGate")
    res <- getGatingMethod.rectangleGates(node, gs, ...)
  else if(gateType == "polygonGate")
    res <- getGatingMethod.polygonGates(node, gs, ...)
  else if(gateType == "ellipsoidGate")
    res <- getGatingMethod.ellipsoidGates(node, gs, ...)
  else
    stop("unkown gate type: ", gateType)
#  browser()
  res <- c(alias = alias
          , pop = pop
          , parent = parent
          , res
          , collapseDataForGating = NA
          , groupBy = NA
          , preprocessing_method = NA
          , preprocessing_args = NA
          )
          
  as.data.frame(t(res))  
}

getGatingMethod.rectangleGates<- function(node, gs, ...){
  cut.points <- sapply(sampleNames(gs), function(sn){
                        thisGate <- getGate(gs[[sn]], node)
                        dims <- parameters(thisGate)
                        fr <- getData(gs, parent)[[sn,dims]]
                        getGatingMethod.rectangleGate(thisGate, fr, ...)
                      })
 null.points <- sapply(cut.points, is.null)
 sum(null.points)/length(cut.points) < 0.5 #less than 50% of null points, then use the majority to derive 1d gate
  
}
#' derive gating method from a \code{rectangleGate}
#' trying to downgrade to 1d gate when applicable
getGatingMethod.rectangleGate <- function(gate, fr, overlap_tol = 0.95, ...){
  
  dims <- parameters(gate)
  data <- exprs(fr)[, dims]
  
  gate_range <- rbind(gate@min,  gate@max)
  data_range <- apply(data, 2, range)
  
  isFullDataRange <- sapply(dims, function(dim){
        
        d_r <- data_range[, dim]
        g_r <- gate_range[, dim]
        
        #get overlap range
        overlap_r <- c(max(d_r[1], g_r[1]), min(d_r[2], g_r[2]))
        
        diff(overlap_r)/diff(d_r) >= overlap_tol
      })
  cut.point <- NULL
  if(sum(isFullDataRange) == 2)
  {
    warning("Can't derive 1d gating method because the gate coordinates cover the full range of both dimesions (a dummy gate)!")
#    getGatingMethod.default(gate, fr)
    
  }else if(sum(isFullDataRange) == 0){
#    getGatingMethod.default(gate, fr)
  }else{
    dims <- dims[!isFullDataRange] #keep one dimension
    #check if either side of the gate reaches(or beyond) to the data boundary
    gate_range <- gate_range[, dims]
    data_range <- data_range[, dims]
    if(gate_range[1] <= data_range[1] || gate_range[2] >= data_range[2]){
      
      cut.point <- ifelse(gate_range[1] <= data_range[1], gate_range[2], gate_range[1])
      
#      getGatingMethod.1d(gate, fr[, dims], cut.point,...)      
      
      
    }else{
      warning("Can't derive 1d gating method because neither edge of ", dims,  " from the gate is at data boundary!")
#      getGatingMethod.default(gate, fr)
    }
    
    
  }
  return(cut.point)
}

#' determine 1d gating method
#' It checks the 1d density and uses 'mindensity' when cut point falls between two modes or 'tailgate' otherwise.
getGatingMethod.1d <- function(gate, fr, cut.point, ...){
  dims <- colnames(fr)
  opt <- options("openCyto")
  gatingMethodName <- opt[["openCyto"]][["default.methods"]][["gating1D"]]
  
  if(is.null(gatingMethodName)){#using default openCyto builin toolsets
    
    x <- exprs(fr)[,1]
#    browser()
    peaks <- sort(openCyto:::.find_peaks(x)[, "x"])
    nPeaks <- length(peaks)
    
    if (nPeaks == 1) {
      getGatingMethod.tailgate(dims, x, cut.point, peaks) 
    }else{
      
      if(cut.point < peaks[1])
        getGatingMethod.tailgate(dims, x, cut.point, peaks[1])
      else if(cut.point > peaks[nPeaks])
        getGatingMethod.tailgate(dims, x, cut.point, peaks[nPeaks])
      else
        getGatingMethod.mindensity(dims, fr, cut.point, ...)
      
    }
    
  }else{#user defined method
    func <- paste('getGatingMethod',gatingMethodName, sep = ".")
    func <- as.symbol(func)
    thisCall <- substitute(f(gate,fr, mu, ...), list(f = func))
    eval(thisCall)  
  }
  
}

getGatingMethod.tailgate<- function(dims, x, cut.point, peak, ...){
  nEvents <- length(x)
  if(cut.point < peak){
    tol <- length(which(x < cut.point))/nEvents
    side <- "left"
  }else{
    tol <- length(which(x > cut.point))/nEvents
    side <- "right"
  }
  
  gating_args <- paste0("side='", side, "',tol=", tol, ")")
  c(dims = dims
      , gating_method = "tailgate"
      , gating_args = gating_args)
}

getGatingMethod.mindensity<- function(dims, fr, cut.point, gate_tol = 0.1, ...){
  data_range <- apply(exprs(fr), 2, range)
  lbound <- cut.point - diff(data_range) * gate_tol
  rbound <- cut.point + diff(data_range) * gate_tol
  gate_range <- paste0("gate_range=c(", lbound, ",", rbound, ")")
  c(dims = dims
      , gating_method = "mindensity"
      , gating_args = gate_range)
}

getGatingMethod.polygonGate <- function(gate,fr, ...){
  rg <- as(gate, "rectangleGate")
  if(!is.null(rg))
    getGatingMethod(rg, fr, ...)
  else
  {
    dims <- parameters(gate)
    if(all(grepl("FSC", dims, ignore.case = TRUE)) || all(grepl("SSC", dims, ignore.case = TRUE))){
      #singlet gate
      c(dims = paste(as.vector(dims), collapse = ",")
          , gating_method = "singletGate"
          , gating_args = NA)
    }else{
      browser()
      getGatingMethod.default(gate, fr)
    }
    
  }
  
}

#' dispatch to the default 2D gating method
getGatingMethod.ellipsoidGate <- function(gate,fr, K = 2,...){
  
  opt <- options("openCyto")
  gatingMethodName <- opt[["openCyto"]][["default.methods"]][["gating2D"]]
  if(is.null(gatingMethodName))
    gatingMethodName <- "flowClust"
  func <- paste('getGatingMethod',gatingMethodName, sep = ".")
  func <- as.symbol(func)
  thisCall <- substitute(f(gate,fr, mu, K, ...), list(f = func, mu = gate@mean, K = K))
  eval(thisCall)
}

#' multiple populations that shares the same parent and projections
getGatingMethod.multiPopulations <- function(nodes, gs, ...){
  gh <- gs[[1]]
  parent <- nodes[["parentId"]]
  children <- nodes[["popIds"]]
  # get all gates (coerce to rectangelgate when applicable)
  gates <- lapply(children, function(node){
        gate <- getGate(gh, node)
        if(class(gate) == "polygonGate"){
          rg <- as(gate, "rectangleGate")
          if(!is.null(rg))
            gate <- rg 
        }
        gate
      })
  
  if(all(sapply(gates, class) == "rectangleGate")){
    #check coordinates to see if they can be merged
    coord_orig <- lapply(gates, function(gate){
          gate_range <- rbind(gate@min,  gate@max)
        })
    coord <- apply(ldply(coord_orig), 2, unique)
    if(nrow(coord) == 3){
      dims <- colnames(coord)
      data <- getData(gs, parent)[[1, dims]]
      data_range <- apply(exprs(data), 2, range)
      
      #check each dims to see if gate coordinates can be coerced to single cutpoint 
      cut.points <- sapply(dims, function(dim){
  
                          d_r <- data_range[, dim]
                          g_r <- sort(coord[, dim])
                          l_bound <- min(d_r)
                          u_bound <- max(d_r)
                          
                          is_edgePoint_out <- g_r[1] <= l_bound && g_r[3] >= u_bound  
                          is_middlePoint_in <- g_r[2] > l_bound && g_r[2] < u_bound
                
                          cut_point <- ifelse(is_edgePoint_out && is_middlePoint_in, g_r[2], NA)
                          cut_point  
                        })
      #if both have valid cut points
      #then construct quad gate 
      if(all(!is.na(cut.points))){
        
        #two 1d helper gates on either dim
        helper_gates <- ldply(dims, function(dim){
                    helper_g <- getGatingMethod.1d(gate = "dummy", data[, dim], cut.points[dim], ...)
                    
                    alias <- pop <- paste(dim, "gate", sep = "_")
                    c(alias = alias
                        , pop = pop
                        , parent = parent
                        , helper_g
                        , collapseDataForGating = NA
                        , groupBy = NA
                        , preprocessing_method = NA
                        , preprocessing_args = NA
                    )
                  })
        #four reference gates based on 1d gating
        #determine the quadrant locations of the original pops by the means of each gate coord
        gates.mu <- lapply(coord_orig, colMeans)
        names(gates.mu) <- children  
        #get 4 corners (as order of -+,++,+-,--)
        corners <- list(top.left = c(min(coord[,1]), max(coord[,2]))
                        , top.right = c(max(coord[,1]), max(coord[,2]))
                        , bottom.right = c(max(coord[,1]), min(coord[,2]))
                        , bottom.left = c(min(coord[,1]), min(coord[,2]))
                        )
        
        # compute the distance from mu to corners 
        pop_order <- sapply(gates.mu, function(mu){
              which.min(sapply(corners,function(corner)dist(rbind(corner,mu))))
            })
        #reorder the children nodes by clock-wise (start from -+)
        pops <- basename(children[pop_order])
        refNodes <- file.path(parent,helper_gates[["alias"]])
        #reference gating
        ref_gates <- c(alias = paste(pops, collapse = ",")
                      , pop = "A+/-B+/-"
                      , parent = parent
                      , dims = paste(dims, collapse = ",")
                      , gating_method = "refGate"
                      , gating_args = paste(refNodes, collapse = ":")
                      , collapseDataForGating = NA
                      , groupBy = NA
                      , preprocessing_method = NA
                      , preprocessing_args = NA
                  )
                  
        res <- rbind(helper_gates, ref_gates)          
        return(data.frame(res))
          
      }                  
                  
    }
  }
  
  res <- ldply(children, function(node){
                  S3Class(node) <- "singlePopulation"
                  getGatingMethod(node, gs, K = 4, ...)
              })

  res
  
}
#' default method leaves 'gating_method' column as NA
getGatingMethod.default <- function(gate,gs){
  dims <- parameters(gate)
  dims <- paste(as.vector(dims), collapse = ",")
  c(dims = dims, gating_method = NA, gating_args = NA)
}



setAs("polygonGate", "rectangleGate", function(from){
      dims <- parameters(from)
      boundaries <- from@boundaries
      
      if(nrow(boundaries) == 4){
        #try to coerce to rectangleGate
        coord <- lapply(dims, function(i)unique(boundaries[,i]))
        if(all(sapply(coord, length) == 2))
          res <- rectangleGate(coord) 
        else
          res <- NULL
      }else
        res <- NULL
      return (res)
    })

getGatingMethod.2d <- function(gate, fr, ...){
  
}

#' generate flowClust based 2d gating method definition
#' default K is 2
getGatingMethod.flowClust<- function(gate, fr, mu, K, ...){
  dims <- parameters(gate)
  c(dims = paste(as.vector(dims), collapse = ",")
      , gating_method = "flowClust"
      , gating_args = paste0("K=", K, ",target=c(", paste0(mu, collapse = ","), ")")
  )
}