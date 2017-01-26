#' toggle the hidden flag of the helper gates
#' 
#' The helper gates are defined as the referred gates in csv template. And all the chidlren of referred gates are also 
#' referred gates thus they are considered the helper gates and can usually be hidden to simply the final gating tree.
#' 
#' @param gt gatingTemplate object
#' @param gs GatingSet
#' @export
#' @importFrom flowWorkspace setNode getNodes getChildren
#' @examples 
#' \dontrun{
#' gt <- gatingTemplate(gtFile)
#' #run the gating
#' gating(gt, gs)
#' #hide the gates that are not of interest
#' toggle.helperGates(gt, gs) 
#' }
toggle.helperGates <- function(gt, gs){
  nodes <- names(getNodes(gt, only.names = T))[-1]
  refNodes <- sapply(nodes, function(node){
                          getParent(gt, node, isRef = TRUE)
                    }, USE.NAMES = FALSE)
  refNodes <- unique(unlist(refNodes))
  isHelper <- sapply(refNodes, function(node){
          children <- getChildren(gs, node)
          if(length(children) == 0)
            return(TRUE)
          else{
            all(children %in% refNodes)
          }
          
        }, USE.NAMES = FALSE)
    
  helperGates <- refNodes[isHelper]
  nonHiddenNodes <- getNodes(gs, showHidden = FALSE, path = "full")
  for(i in helperGates){
      if(i%in%nonHiddenNodes)
        setNode(gs, i, FALSE)
      else
        setNode(gs, i, TRUE)
  }
}
