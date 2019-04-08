#' toggle/delete the hidden flag of the helper gates
#' 
#' The helper gates are defined as the referred gates in csv template. And all the chidlren of referred gates are also 
#' referred gates thus they are considered the helper gates and can usually be hidden to simply the final gating tree.
#' 
#' Note that delete action is NOT reversible.
#' 
#' @param gt gatingTemplate object
#' @param gs GatingSet
#' @export
#' @importFrom flowWorkspace gh_set_node_visible gs_get_pop_paths gs_get_children
#' @examples 
#' \dontrun{
#' gt <- gatingTemplate(gtFile)
#' #run the gating
#' gating(gt, gs)
#' #hide the gates that are not of interest
#' toggle.helperGates(gt, gs) 
#' #or simply remove them if you are sure they will not be useful in future
#' delete.helperGates(gt, gs) 
#' }
#' @rdname toggle.helperGates
toggle.helperGates <- function(gt, gs){
  helperGates <- get.helperGates(gt, gs)
  nonHiddenNodes <- gs_get_pop_paths(gs, showHidden = FALSE, path = "full")
  for(i in helperGates){
      if(i%in%nonHiddenNodes)
		  gh_set_node_visible(gs, i, FALSE)
      else
		  gh_set_node_visible(gs, i, TRUE)
  }
}

#' @rdname toggle.helperGates
#' @export
get.helperGates <- function(gt, gs){
  gated.nodes <- gs_get_pop_paths(gs, showHidden = TRUE, path = "full")
  referror <- names(getNodes(gt, only.names = T))[-1]
  referror <- intersect(gated.nodes, referror)#restrict to the gated nodes
  
  referree <- sapply(referror, function(node){
    getParent(gt, node, isRef = TRUE)
  }, USE.NAMES = FALSE)
  referree <- unique(unlist(referree))
  referree <- intersect(gated.nodes, referree)#restrict to the gated nodes
  gh <- gs[[1]]
  isHelper <- sapply(referree, function(node){
    children <- gs_get_children(gh, node)
    if(length(children) == 0)
      return(TRUE)
    else{
      all(children %in% referree)
    }
    
  }, USE.NAMES = FALSE)
  if(length(isHelper) > 0)
  {
    referree[isHelper]  
  }else
  character(0)
  
}

#' @rdname toggle.helperGates
#' @export
delete.helperGates <- function(gt, gs){
  helperGates <- get.helperGates(gt, gs)
  
  for(i in helperGates){
      if(i%in%gs_get_pop_paths(gs))
        Rm(i, gs)
  }
}
