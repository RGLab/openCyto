setMethod("getNodes", signature = c("fcTree"), definition = function(x, y) {
  if (missing(y)) {
    nodeData(x)
  } else {
    nodeData(x, y)
  }
})

setMethod("getGate", sig = c("fcTree", "character"),
    definition = function(x, y,  ...) {
      # get filterList
      nodes <- getNodes(x)
      
      allAlias <- lapply(nodes, function(curNode) alias(curNode$pop))
      ind <- which(allAlias %in% y)
      if (length(ind) > 1) {
        stop("Population '", y, "' is ambiguous!")
      } else if (length(ind) == 0) {
        stop("Population '", y, "' is not found!")
      } else {
        matchedNode <- nodes[[ind]]
      }
      
      flist <- matchedNode[["fList"]]
      flist
    })

setMethod("plot", sig = c("fcTree", "character"),
          definition = function(x, y, channel = NULL, ...) {
  # get filterList
    flist <- getGate(x,y,...)
  plot(x = flist, y = channel, main = matchedNode$pop@name, ...)
})

#setMethod("plot", sig = c("fcTree", "numeric"),
#          definition = function(x, y, channel = NULL, ...) {
#  y <- as.character(y)
#  nodes <- getNodes(x)
#  matchedNode <- nodes[[y]]
#  fcObj <- matchedNode$fcObj
#  
#  plot(x = fcObj, y = channel, main = matchedNode$pop@name, ...)
#}) 
