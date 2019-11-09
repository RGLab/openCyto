#' CytoExploreR exports
#' 
#' Exported wrappers of internal functions for use by CytoExploreR
#' 
#' @name CytoExploreR_exports
NULL

#' @name CytoExploreR_.argDeparser
#' @keywords internal
#' @rdname CytoExploreR_exports
#' @export
CytoExploreR_.argDeparser <- function(args, split = TRUE){
  .argDeparser(args, split = split)
}

#' @name CytoExploreR_.preprocess_csv
#' @keywords internal
#' @rdname CytoExploreR_exports
#' @export
CytoExploreR_.preprocess_csv <- function(dt, strict = TRUE){
  .preprocess_csv(dt, strict = strict)
}