#' @templateVar old toggle.helperGates
#' @templateVar new gt_toggle_helpergates
#' @template template-depr_pkg
NULL
#' toggle/delete the hidden flag of the helper gates
#' 
#' The helper gates are defined as the referred gates in csv template. And all the chidlren of referred gates are also 
#' referred gates thus they are considered the helper gates and can usually be hidden to simply the final gating tree.
#' 
#' Note that delete action is NOT reversible.
#' 
#' @name gt_toggle_helpergates
#' @aliases toggle.helperGates
#' @param gt gatingTemplate object
#' @param gs GatingSet
#' @export
#' @importFrom flowWorkspace gh_pop_set_visibility gs_get_pop_paths gs_pop_get_children
#' @examples 
#' \dontrun{
#' gt <- gatingTemplate(gtFile)
#' #run the gating
#' gt_gating(gt, gs)
#' #hide the gates that are not of interest
#' gt_toggle_helpergates(gt, gs) 
#' #or simply remove them if you are sure they will not be useful in future
#' gt_delete_helpergates(gt, gs) 
#' }
gt_toggle_helpergates <- function(gt, gs){
  helperGates <- gt_get_helpergates(gt, gs)
  nonHiddenNodes <- gs_get_pop_paths(gs, showHidden = FALSE, path = "full")
  for(i in helperGates){
      if(i%in%nonHiddenNodes)
		  gs_pop_set_visibility(gs, i, FALSE)
      else
		  gs_pop_set_visibility(gs, i, TRUE)
  }
}

#' @export
toggle.helperGates <- function(gt, gs){
  .Deprecated("gt_toggle_helpergates")
  gt_toggle_helpergates(gt, gs)
}

#' @templateVar old get.helperGates
#' @templateVar new gt_get_helpergates
#' @template template-depr_pkg
NULL
#' @rdname gt_toggle_helpergates
#' @aliases get.helperGates
#' @export
gt_get_helpergates <- function(gt, gs){
  gated.nodes <- gs_get_pop_paths(gs, showHidden = TRUE, path = "full")
  referror <- names(gt_get_nodes(gt, only.names = T))[-1]
  referror <- intersect(gated.nodes, referror)#restrict to the gated nodes
  
  referree <- sapply(referror, function(node){
    gt_get_parent(gt, node, isRef = TRUE)
  }, USE.NAMES = FALSE)
  referree <- unique(unlist(referree))
  referree <- intersect(gated.nodes, referree)#restrict to the gated nodes
  gh <- gs[[1]]
  isHelper <- sapply(referree, function(node){
    children <- gs_pop_get_children(gh, node)
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

#' @export
get.helperGates <- function(gt, gs){
  .Deprecated("gt_get_helpergates")
  gt_get_helpergates(gt, gs)
}

#' @templateVar old delete.helperGates
#' @templateVar new gt_delete_helpergates
#' @template template-depr_pkg
NULL
#' @rdname gt_toggle_helpergates
#' @aliases delete.helperGates
#' @export
gt_delete_helpergates <- function(gt, gs){
  helperGates <- gt_get_helpergates(gt, gs)
  
  for(i in helperGates){
      if(i%in%gs_get_pop_paths(gs))
        gs_pop_remove(i, gs = gs)
  }
}

#' @export
delete.helperGates <- function(gt, gs){
  .Deprecated("gt_delete_helpergates")
  gt_delete_helpergates(gt, gs)
}
