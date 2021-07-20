#' @templateVar old getNodes
#' @templateVar new gt_get_nodes
#' @template template-depr_pkg
NULL
#' get nodes from \link[openCyto:gatingTemplate-class]{gatingTemplate} object
#' 
#' @name gt_get_nodes
#' @aliases getNodes getNodes,gatingTemplate-method
#' @param x \code{gatingTemplate}
#' @param y \code{character} node index. When \code{missing}, return all the nodes
#' @param order \code{character} specifying the order of nodes. options are "default", "bfs", "dfs", "tsort"
#' @param only.names \code{logical} specifiying whether user wants to get the entire \code{gtPopulation} object or just the name of the population node
#'   
#' @examples 
#' \dontrun{
#' gt <- gatingTemplate(system.file("extdata/gating_template/tcell.csv",package = "openCyto"))
#' gt_get_nodes(gt)[1:2]
#' gt_get_nodes(gt, only.names = TRUE)
#' gt_get_nodes(gt, "/nonDebris")
#' }
#' @export
gt_get_nodes <- function(x, y
                             , order = c("default", "bfs", "dfs", "tsort")
                             , only.names = FALSE) {
  
  if (missing(y)){
    res <- nodeData(x, attr = "pop")
    order <- match.arg(order)
    if(order != "default"){
      nodeIds <- eval(substitute(f1(x),list(f1=as.symbol(order))))
      if(order == "dfs")
        nodeIds <- nodeIds$discovered
      res <- res[nodeIds]
    }
  }else
  {
    res <- nodeData(x, y, "pop")
  }
  if(only.names){
    res <- sapply(res,alias)
  }
  if(length(res) == 1 && class(res) == "list")
    res <- res[[1]]
  res
}

#' @export
setMethod("getNodes", signature = c("gatingTemplate"),
          definition = function(x, y
                                  , order = c("default", "bfs", "dfs", "tsort")
                                  , only.names = FALSE) {
            .Deprecated("gt_get_nodes")
            gt_get_nodes(x, y, order, only.names)
})

#' @templateVar old getChildren
#' @templateVar new gt_get_children
#' @template template-depr_pkg
NULL
#' get children nodes
#'
#' @name gt_get_children
#' @aliases getChildren getChildren,gatingTemplate,character-method
#' @param obj \code{gatingTemplate}
#' @param y \code{character} parent node path  
#' @export 
#' @examples 
#' \dontrun{
#' gt <- gatingTemplate(system.file("extdata/gating_template/tcell.csv",package = "openCyto"))
#' 
#' gt_get_nodes(gt, "/nonDebris")
#' gt_get_children(gt, "/nonDebris") 
#' }
#' @importClassesFrom methods character ANY data.frame environment list logical matrix missing numeric oldClass
gt_get_children <- function(obj, y){
  edges(obj, y)[[1]]
}

#' @export
setMethod("getChildren", signature = c("gatingTemplate", "character"),
          definition = function(obj, y) {
            .Deprecated("gt_get_children")
            gt_get_children(obj, y)
})

#' @templateVar old getParent
#' @templateVar new gt_get_parent
#' @template template-depr_pkg
NULL
#' get parent nodes
#' 
#' @param obj \code{gatingTemplate}
#' @param y \code{character} child node path
#' @param isRef \code{logical} whether show the reference node besides the parent node
#' 
#' @name gt_get_parent
#' @aliases getParent getParent,gatingTemplate,character-method
#' @importFrom graph inEdges
#' @aliases getParent,gatingTemplate,character-method
#' @examples 
#' \dontrun{
#' gt <- gatingTemplate(system.file("extdata/gating_template/tcell.csv",package = "openCyto"))
#' 
#' gt_get_nodes(gt, "/nonDebris")
#' gt_get_parent(gt, "/nonDebris/singlets") 
#' }
#' @export 
gt_get_parent <- function(obj, y, isRef = FALSE) {
  #            browser()
  src <- inEdges(y, obj)[[1]]
  isRefs <- unlist(lapply(src,function(thisSrc){
    thisEdge <- edgeData(obj,thisSrc,y)
    thisEdge[[1]]$isReference 
  }), recursive = F)
  #  browser()
  if(isRef)
    src[isRefs]
  else
    src[!isRefs]
}

#' @export
setMethod("getParent", signature = c("gatingTemplate", "character"),
          definition = function(obj, y, isRef = FALSE) {
            .Deprecated("gt_get_parent")
            gt_get_parent(obj, y, isRef)
})

#' @templateVar old getGate
#' @templateVar new gt_get_gate
#' @template template-depr_pkg
NULL
#' get gating method from the node
#' 
#' @name gt_get_gate
#' @aliases getGate,gatingTemplate,character-method
#' @param obj \code{gatingTemplate}
#' @param y \code{character} parent node path
#' @param z \code{character} child node path
#' @examples 
#' \dontrun{
#' gt <- gatingTemplate(system.file("extdata/gating_template/tcell.csv",package = "openCyto"))
#' gt_get_nodes(gt, only.names = TRUE)
#' gt_get_nodes(gt, "/nonDebris")
#' gt_get_children(gt, "/nonDebris")
#' gt_get_gate(gt, "/nonDebris", "/nonDebris/singlets")
#' }
#' @export
gt_get_gate <- function(obj, y, z) {
  edgeData(obj, from = y, to = z, attr = "gtMethod")[[1]]
}

#' @export
setMethod("getGate", signature = c("gatingTemplate", "character"),
          definition = function(obj, y, z) {
            .Deprecated("gt_get_gate")
            gt_get_gate(obj, y, z)
})

setGeneric("ppMethod", function(obj, y, ...) standardGeneric("ppMethod"))
#' get preprocessing method from the node
#' 
#' get preprocessing method from the node
#' @param obj \code{gatingTemplate}
#' @param y \code{character} parent node path
#' @param z \code{character} child node path
#' @aliases ppMethod,gatingTemplate,character-method
#' @examples 
#' \dontrun{
#' gt <- gatingTemplate(system.file("extdata/gating_template/tcell.csv",package = "openCyto"))
#' ppMethod(gt,  "/nonDebris/singlets",  "/nonDebris/singlets/lymph")
#' }
setMethod("ppMethod", signature = c("gatingTemplate", "character"),
    definition = function(obj, y, z) {
      edgeData(obj, from = y, to = z, attr = "ppMethod")[[1]]
    })

#' show method for gatingTemplate
#' 
#' show method for gatingTemplate
#' 
#' @param object \code{gatingTemplate}
#' @aliases 
#' show,gatingTemplate-method
#' @export 
#' @importFrom methods show
setMethod("show", signature = c("gatingTemplate"),
          definition = function(object) {
  cat("--- Gating Template: ")
  cat(object@name)
  cat("\n")
  cat("\twith ", length(object@nodes), " populations defined\n")
})

#' plot the gating scheme
#' 
#' plot the gating scheme using Rgraphviz
#' 
#' @param x \code{gatingTemplate} object
#' @param y either \code{character} specifying the root node which can be used to visualize only the subgraph 
#'              or \code{missing} which display the entire gating scheme
#' @param ... other arguments 
#'  
#'      graphAttr, nodeAttr:  graph rendering attributes passed to \link[Rgraphviz:renderGraph]{renderGraph}
#'      showRef \code{logical}: whether to display the reference gates. Sometime it maybe helpful to 
#'                              hide all those reference gates which are not the cell population of interest 
#'                              and used primarily for generating other population nodes.
#' @export 
#' @importFrom RColorBrewer brewer.pal
#' @aliases 
#' plot,gatingTemplate,missing-method
#' plot,gatingTemplate,character-method
#' plot,gatingTemplate-method
#' plot,gatingTemplate,ANY-method
#' @examples 
#' \dontrun{
#' gt <- gatingTemplate(system.file("extdata/gating_template/tcell.csv",package = "openCyto"))
#' plot(gt) #plot entire tree
#' plot(gt, "lymph") #only plot the subtree rooted from "lymph"
#' }
setMethod("plot",c("gatingTemplate","missing"),function(x,y,...){
      .plotTree(x,...)
      
    })

#FIXME:somehow the edge attributes are not correctly assigned to each edge respectively
.plotTree<-function(x
                      , graphAttr = list(rankdir = "LR", page = c(8.5, 11)) 
                      , nodeAttr = list(fixedsize = FALSE, shape = "ellipse")
                      , showRef = TRUE
                      , ...){
#          browser()          
  # get gating method name attr from edges
  hasEdge <- length(edgeData(x)) > 0
  if(hasEdge){
    gm.names <- unlist(lapply(edgeData(x, attr = "gtMethod"), names))
    gm.types <- unique(gm.names)
    
    ref.edges <- unlist(edgeData(x, attr = "isReference"))
    ref.edges <- ref.edges[ref.edges]
    
    
    # fix the name attr
    e.colnames <- gsub("\\|", "~", names(gm.names))
    ref.colnames <- gsub("\\|", "~", names(ref.edges))
    
    #only care about the ref edges that are not conflicting with gm edges
    ref.colnames <- ref.colnames[!ref.colnames%in%e.colnames]
    
    
    
    # encode the method name with color
    nMethods <- length(gm.types)
    gm.col <- brewer.pal(nMethods, name = "Dark2")
    names(gm.col) <- gm.types
    eCol.gm <- gm.col[gm.names]
    
    # restore names
    names(eCol.gm) <- e.colnames  
    
    eCols <- eCol.gm
    
    #add ref edges
    eCol.ref <- sapply(ref.colnames,function(i)return("grey"))
    if(length(eCol.ref)>0)
      eCols <- c(eCols, eCol.ref)
    
    # encode edge style
    gm.isPoly <- unlist(lapply(gm.names, function(y) {
              ifelse(y == "polyFunctions", "poly", "regular")
            }))
    edge.styles <- c("solid", "dashed")
    names(edge.styles) <- unique(gm.isPoly)
    eStyles.gm <- edge.styles[gm.isPoly]
    names(eStyles.gm) <- e.colnames
    
    eStyles <- eStyles.gm
    
    eStyles.ref <- sapply(ref.colnames,function(i)ifelse(showRef,"dotted","blank"))
    if(length(eStyles.ref)>0)
      eStyles <- c(eStyles,eStyles.ref)
    
    eAttrs <- list(color = eCols, lty = eStyles)
    
  }else{
    eAttrs <- list()
  }
    
  # encode the node shape and color with isSubsets
  nodeTypes <- unlist(lapply(nodeData(x, attr = "pop"), function(y) {
            ifelse(class(y) == "gtSubsets", "subset", "pop")
          }))
  n.colnames <- names(nodeTypes)
  nodeType <- unique(nodeTypes)
  
  # color
  subset.col <- brewer.pal(9, name = "Set3")[c(9, 7)]
  names(subset.col) <- nodeType
  nFillCol <- subset.col[nodeTypes]
  names(nFillCol) <- n.colnames
  
  # line type (somehow it doesn't work here)
  LineType <- c("solid", "dotted")
  names(LineType) <- nodeType
  nLtys <- LineType[nodeTypes]
  names(nLtys) <- n.colnames
  
  nLabels <- unlist(lapply(nodeData(x, attr = "pop"), function(thisPop)paste(alias(thisPop), collapse = ",")))
  nAttrs <- list(label = nLabels, fillcolor = nFillCol, lty = nLtys)
  
  plot(as(x,"graphNEL"), nodeAttrs = nAttrs, edgeAttrs = eAttrs
                  ,attrs = list(graph = graphAttr,node = nodeAttr)
    )
  
  if(hasEdge){
    plot.space = par()[["usr"]]
    x1 = plot.space[1]
    x2 = plot.space[2]
    y1 = plot.space[3]
    y2 = plot.space[4]
    
    legend.lty <- ifelse(gm.types == "polyFunctions", "poly", "regular")
    legend(x1 + 100, y2 - 100, legend = gm.types, title = "Gating Methods",
        col = gm.col, lty = edge.styles[legend.lty], cex = 0.8)
  }
}

setMethod("plot", signature = c("gatingTemplate"),
          definition = function(x, y = missing,...) {
  
})

setMethod("plot",c("gatingTemplate","character"),function(x,y,...){
      
      #convert alias to nodeID
      allNodes <- gt_get_nodes(x)
      allAlias <- unlist(lapply(allNodes, alias), recursive = F)
      nodeInd <- match(y, allAlias)
      thisNode <- allNodes[nodeInd]
      thisId <- names(thisNode)
#      browser()
      if(length(y)==1){#use it as the root
        
        nodelist <- new.env(parent=emptyenv())
        nodelist$v <- character()
        flowWorkspace:::.getAllDescendants (x, thisId, nodelist)  
        
        
#        browser()
        nodelist$v <- c(nodelist$v,thisId)
        nodelist$v <- unique(nodelist$v)
        #assume the number y is consistent with  R graph node name: N_x
        subNode_Ind <- nodelist$v
        
      }else{
        #when y is a vector, use it to subset the graph
        subNode_Ind <- thisId
      }
      
#      subNodes <- paste("N",subNode_Ind-1,sep="_")
      if(length(subNode_Ind)<=1){
        stop("Rgraphviz doesn't know how to plot leaf node!")
      }
      x <- subGraph(subNode_Ind, x)
      .plotTree(x,...)
      
    })
