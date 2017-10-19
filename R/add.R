
      

#' the class that carries event indices as well

setClass("ocRectangleGate", contains = "rectangleGate", representation(ind = "raw"))

#' bypass the default flowWorkspace:::.addGate 
#' 
#' to support adding gate along with indices without loading flow data and computing
#' 
#' however it is proven that logical indices are too big to be efficiently passed around
#' 
#' @param wf \code{GatingHierarchy} see \link[flowWorkspace]{add} in \code{flowWorkspace} package
#' @param action \code{ocRectangleGate} or \code{logicalFilterResult}
#' @param recompute \code{logical} see \link[flowWorkspace]{add} in \code{flowWorkspace} package
#' @param ... see \link[flowWorkspace]{add} in \code{flowWorkspace} package
#' @export 
#' @rdname add
setMethod("add",
    signature=c("GatingHierarchy", "ocRectangleGate"),
    definition=function(wf, action, recompute, ... )
    {
      #unpack the bit vector
      indices <- ncdfFlow:::toLogical(action@ind)
      #ignore the recompute flag and force it to be skipped
      nodeID <- flowWorkspace:::.addGate(wf, filterObject(action), recompute = FALSE, ...)
      sn <- sampleNames(wf)
      ptr <- wf@pointer
      flowWorkspace:::.cpp_setIndices(ptr, sn, nodeID-1, indices)
      
      
    })


#' special gate type that mix the rectangleGate with boolean gate
setClass("ocRectRefGate", contains = c("rectangleGate", "booleanFilter"))
#' constructor for ocRectRefGate
#' @param rectGate \code{rectangleGate}
#' @param boolExprs \code{character} boolean expression of reference nodes
ocRectRefGate <- function(rectGate, boolExprs){
  
  bf <- eval(substitute(booleanFilter(x), list(x = as.symbol(boolExprs))))
  g <- as(rectGate, "ocRectRefGate")
  g@expr <- bf@expr
  g@deparse <- bf@deparse
  g
}

#' byPass the default .addGate 
#' 
#' to support adding rectangleGate yet gating through boolean operations 
#' without loading flow data
#' 
#' @export 
#' @rdname add
setMethod("add",
    signature=c("GatingHierarchy", "ocRectRefGate"),
    definition=function(wf, action, recompute, ... )
    {
      rectFilterObj <- selectMethod("filterObject", signature = c("rectangleGate"))(action)
      boolFilterObj <- selectMethod("filterObject", signature = c("booleanFilter"))(action)
      #ignore the recompute flag and force it to be skipped
      nodeID <- flowWorkspace:::.addGate(wf, rectFilterObj, recompute = FALSE, ...)
      sn <- sampleNames(wf)
      ptr <- wf@pointer
      flowWorkspace:::.cpp_boolGating(ptr, sn, boolFilterObj, nodeID - 1)
      nodeID
    })




#' fast version of add gates to gatingset (bypassing some R checks)
#' 
#' used by gating_polyFunctions
.addGate_fast <- function(gs, filter, name = NULL, parent = "root", negated = FALSE){
  
  #preprocess filter
  filterObj <- filterObject(filter)
#  browser()
  if(is.null(name))
    name<-filterObj$filterId
  #replace the slash with colon 
  #since forward slash is reserved for gating path
  if(grepl("/",name)){
    old_name <- name
    name <- gsub("/",":",name)
    warning(old_name, " is replaced with ", name)
  }
  
  
  gh<-gs[[1]]
  ##get node ID
  
  filterObj$negated<-negated
  
  
  if(class(gs) == "GatingSetList"){
    nodeIDs <- lapply(gs, function(thisGS){
          samples <- sampleNames(thisGS)
          lapply(samples,function(sample){
                
                nodeID <- flowWorkspace:::.cpp_addGate(thisGS@pointer,sample,filterObj, parent,name)
                nodeID+1
              })
        }, level = 1)
    nodeIDs <- unlist(nodeIDs)
    
  }else{
    samples <- sampleNames(gs)
    nodeIDs<-lapply(samples,function(sample){
          
          nodeID <- flowWorkspace:::.cpp_addGate(gs@pointer,sample,filterObj, parent,name)
          nodeID+1
        })  
  }
  
  
  nodeID<-nodeIDs[[1]]
  
  if(!all(sapply(nodeIDs[-1],function(x)identical(x,nodeID))))
    stop("nodeID are not identical across samples!")
  nodeID
}
