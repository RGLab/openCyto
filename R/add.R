#' @export 
#' @rdname add
setMethod("add",
    signature=c("GatingHierarchy", "logicalFilterResult"),
    definition=function(wf, action,... )
    {
      g <- filterDetails(action)[[1]][["filter"]]
      ind <- action@subSet
      .addGate_no_data(wf,filterObject(g), indices = ind, ...)
    })

#' the class that carries event indices as well

setClass("ocRectangleGate", contains = "rectangleGate", representation(ind = "raw"))

#' @export 
#' @rdname add
setMethod("add",
    signature=c("GatingHierarchy", "ocRectangleGate"),
    definition=function(wf, action,... )
    {
      #unpack the bit vector
      ind <- ncdfFlow:::.getBitStatus(action@ind)
      .addGate_no_data(wf,filterObject(action), indices = ind, ...)
    })


#' special gate type that mix the rectangleGate with boolean gate
setClass("ocRectRefGate", contains = c("rectangleGate", "booleanFilter"))
#' constructor for ocRectRefGate
#' @param rectGate \code{rectangleGate}
#' @param boolExprs \code{character} boolean expression of reference nodes
ocRectRefGate <- function(rectGate, boolExprs){
  
  bf <- eval(substitute(booleanFilter(x), list(x = boolExprs)))
  g <- as(rectGate, "ocRectRefGate")
  g@expr <- bf@expr
  g@deparse <- bf@deparse
  g
}

#' @export 
#' @rdname add
setMethod("add",
    signature=c("GatingHierarchy", "ocRectRefGate"),
    definition=function(wf, action,... )
    {
      #unpack the bit vector
      flowWorkspace:::.addGate(wf,filterObject(action), indices = ind, ...)
    })


#setMethod("filterObject",signature=c("ocRectRefGate"),function(x){
#
#    rg_res <- selectMethod("filterObject", signature = c("rectangleGate"))(x)
#    bg_res <- selectMethod("filterObject", signature = c("booleanFilter"))(x)
#    list(rect = rg_res, bool = bg_res, filterId = rg_res[["filterId"]])
#    })

#' extend flowWorkspace:::.addGate to support adding indices along with gate yet without loading flow data
#' 
#' @param indices \code{logical} vector used to pass the node indices directly
#'                  or \code{expression} 
.addGate_no_data <- function(gh, filterObject, indices, recompute, ...){
  #ignore the recompute flag and force it to be skipped
  
  nodeID <- flowWorkspace:::.addGate(gh, filterObject, recompute = FALSE, ...)
  sn <- sampleNames(gh)
  ptr <- gh@pointer
  if(is.logical(indices))
    .Call("R_setIndices", ptr, sn, nodeID-1, indices, PACKAGE = "flowWorkspace")
  else
  {
    stop("not supported!")
    #for openCyto::ocRectRefGate
#    .Call("R_boolGating", ptr, sn, ,nodeID, indices)
  }
    
    
  nodeID
}
#' fast version of add gates to gatingset (bypassing some R checks)
.addGate_fast <- function(gs, filter, name = NULL, parent = NULL, negated = FALSE){
  
  
  
  
  #preprocess filter
  filterObj <- flowWorkspace:::filterObject(filter)
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
  if(is.null(parent))
    pid<-1
  else
  {
    if(is.numeric(parent))
      pid <- parent
    else
      pid <- flowWorkspace:::.getNodeInd(gh,parent)
  }
  filterObj$negated<-negated
  
  
  if(class(gs) == "GatingSetList"){
    nodeIDs <- lapply(gs, function(thisGS){
          samples <- sampleNames(thisGS)
          lapply(samples,function(sample){
                
                nodeID <- .Call("R_addGate",thisGS@pointer,sample,filterObj,as.integer(pid-1),name)
                nodeID+1
              })
        }, level = 1)
    nodeIDs <- unlist(nodeIDs)
    
  }else{
    samples <- sampleNames(gs)
    nodeIDs<-lapply(samples,function(sample){
          
          nodeID <- .Call("R_addGate",gs@pointer,sample,filterObj,as.integer(pid-1),name)
          nodeID+1
        })  
  }
  
  
  nodeID<-nodeIDs[[1]]
  
  if(!all(sapply(nodeIDs[-1],function(x)identical(x,nodeID))))
    stop("nodeID are not identical across samples!")
  nodeID
}
