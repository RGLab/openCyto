#' OpenCyto plugin for custom gating functions
#'
#' \code{.gate_custom()} is an openCyto plugin that accepts any custom gating
#' function through the \code{FUN} argument to allow the use of any custom
#' gating function without the need to write additional openCyto plugins.
#'
#' @param fr a \code{flowFrame} or \code{cytoframe} object containing the data
#'   to be gated by \code{FUN}.
#' @param pp_res the output of the \code{.pp_gate_custom()} pre-processing
#'   method.
#' @param channels vector of channels names to indicate the parameters to use
#'   for gating. \code{channels} are parsed by OpenCyto using the \code{dims}
#'   specified in the \code{gatingTemplate} and therefore \code{channels} can
#'   only include one or two parameters. The \code{params} argument is reserved
#'   for cases where more than two channels are to be passed to the gating
#'   function.
#' @param FUN a gating function or name of a gating function supplied as a
#'   character string. \code{FUN} is applied to the data following
#'   \code{inverse} transformations and/or scaling along with any additional
#'   arguments specified through \code{...}.
#' @param input indicates how the data extracted from the \code{cytoframe}
#'   should be formatted prior to passing it to the gating function supplied to
#'   \code{FUN}, set to \code{"cytoframe"} by default. Options include
#'   \code{"matrix"}, \code{"data.table"} or \code{"cytoframe"}.
#' @param inverse logical indicating whether inverse data transformations should
#'   be applied to the extracted data prior to passing it to the gating
#'   function, set to FALSE by default.
#' @param scale can be either a logical or the name of a scaling function to
#'   apply to each channel of the \code{cytoframe} prior to applying the gating
#'   function supplied through \code{FUN}, set to FALSE by default. Min/max
#'   scaling can be applied by setting \code{scale = "range'}, which is the
#'   default behavior when \code{scale = TRUE}.
#' @param slot provides a mechanism by which gate objects can be extracted from
#'   the outputs of gating functions that return complex data structures (e.g.
#'   list objects with multiple elements). \code{slot} can either be the name of
#'   a slot or a the name of a function extracts or formats the output of the
#'   gating function to return a valid gate object.
#' @param ... additional arguments passed to \code{FUN} for gating.
#'
#' @return a gate object as returned by \code{FUN}.
#' 
#' @importFrom flowCore exprs exprs<-
#' @importFrom flowWorkspace realize_view cf_unlock cf_lock transform
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@ozette.com}
#'
#' @noRd
.gate_custom <- function(fr,
                         pp_res,
                         channels,
                         FUN,
                         input = "cytoframe",
                         inverse = FALSE,
                         scale = FALSE,
                         slot  = NULL,
                         ...) {
  
  # NOTE: openCyto currently supports 1D or 2D gates only!
  
  # extract transformers from pp_res
  trans <- pp_res$trans
  
  # data requires scaling or inverse transformations
  if(!scale %in% FALSE | (inverse & any(channels %in% names(trans)))) {
    # cytoframe must be copied and unlocked for editing
    if(inherits(fr, "cytoframe")) {
      fr <- realize_view(fr)
      cf_unlock(fr)
    }
    # inverse transformations
    if(inverse & any(channels %in% names(trans))) {
      fr <- transform(
        fr,
        lapply(trans, `[[`, "inverse")
      )
    }
    # scaling
    if(!scale %in% FALSE) {
      # default range scale
      if(scale %in% TRUE) {
        scale <- "range"
      }
      # range (min/max) scaling
      if(scale %in% "range") {
        scale <- function(x) {
          (x - min(x))/diff(range(x))
        }
      # scaling function 
      } else {
        scale <- .fun_match(
          scale
        )
      }
      # apply scaling function
      exprs(fr)[, channels] <- apply(
        exprs(fr)[, channels],
        2,
        scale
      )
    }
    # lock copied cytoframe for editing
    if(inherits(fr, "cytoframe")) {
      cf_lock(fr)
    }
  }
  
  # input - data extraction required
  if(!input %in% "cytoframe") {
    # extract matrix
    fr <- exprs(fr[, channels])
    # coerce to vector | data.table | data.frame | matrix
    fr <- do.call(
      "as",
      list(
        fr,
        input
      )
    )
  }
  
  # extract gating arguments
  gating_args <- list(
    fr,
    ...
  )
  
  # get gating function
  FUN <- .fun_match(
    FUN
  )
  
  # apply gating function
  gate <- do.call(
    FUN,
    gating_args
  )
  
  # extract gate object
  if(!is.null(slot)) {
    # slot name
    if(slot %in% names(gate)) {
      gate <- gate[[slot]]
    # slot function
    } else {
      gate <- do.call(
        .fun_match(slot),
        list(gate)
      )
    }
  }
  
  # return gate object
  return(gate)
  
}

#' OpenCyto pre-processing plugin to transfer transformers to gating function
#'
#' @param fs a cytoset containing samples for the current sample group.
#' @param gs a GatingSet.
#' @param gm a gating method.
#' @param channels the names of the channels to use or gating.
#' @param groupBy grouing variables separated by colon to be used to separate
#'   samples into groups.
#' @param isCollapse logical indicating whether the data is collapsed prior to
#'   gating.
#' @param trans an optional \code{transformerList} containing the transformers
#'   applied to the GatingSet, automatically extracted from the\code{GatingSet}
#'   directly.
#' @param ... not in use.
#'
#' @return returns a list of pre-processing arguments including \code{trans}
#'   which contains the transformers applied to the \code{GatiingSet}.
#'
#' @importFrom flowWorkspace gh_get_transformations transformerList
#'
#' @author Dillon Hammill, \email{Dillon.Hammill@ozette.com}
#'
#' @noRd
.pp_gate_custom <- function(fs,
                            gs,
                            gm,
                            channels,
                            groupBy = NA,
                            isCollapse = NA,
                            trans = list(),
                            ...) {
  
  # extract transformers
  if(!inherits(trans, "transformerList")) {
    trans <- gh_get_transformations(
      gs[[1]],
      only.function = FALSE
    )
    if(length(trans) > 0) {
      trans <- transformerList(
        names(trans),
        trans
      )
    }
  }
  
  # return required arguments
  return(
    list(
      "trans" = trans
    )
  )
  
}

#' Internal function to source a namespaced function provided by name
#'
#' @param FUN character vector specifying the name of the function to source,
#'   may be prefixed with a namespace using either \code{pkg::FUN} for exported
#'   or \code{pkg:::FUN} for internal functions respectively.
#' @param ... additional arguments passed to \code{match.fun()}.
#'
#' @return the sourced function.
#' 
#' @author Dillon Hammill, \email{Dillon.Hammill@ozette.com}
#' 
#' @examples
#' FUN <- .fun_match(
#'   "stats::median"
#' )
#' 
#' @noRd
.fun_match <- function(FUN,
                       ...) {
  
  # name of function -> function
  if(is.character(FUN)) {
    # namespaced function from another package
    if(grepl(":{2,3}", FUN)) {
      # split into namespace and function
      FUN <- unlist(
        strsplit(
          FUN, 
          ":{2,3}"
        )
      )
      # source function from namespace
      FUN <- get(
        FUN[2],
        envir = asNamespace(FUN[1]),
        mode = "function"
      )
    }
    # get function
    FUN <- match.fun(
      FUN, 
      descend = TRUE
    )
  }
  
  # return a function to call
  return(FUN)
  
}
