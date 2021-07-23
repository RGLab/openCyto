#adapted from plyr::list_to_array
id_var<-
function (x, drop = FALSE) 
{
  if (length(x) == 0) 
    return(structure(integer(), n = 0L))
  if (!is.null(attr(x, "n")) && !drop) 
    return(x)
  if (is.factor(x) && !drop) {
    x <- addNA(x, ifany = TRUE)
    id <- as.integer(x)
    n <- length(levels(x))
  }
  else {
    levels <- sort(unique(x), na.last = TRUE)
    id <- match(x, levels)
    n <- max(id)
  }
  structure(id, n = n)
}
id <-
function (.variables, drop = FALSE) 
{
  lengths <- vapply(.variables, length, integer(1))
  .variables <- .variables[lengths != 0]
  if (length(.variables) == 0) {
    n <- nrow(.variables) %||% 0L
    return(structure(seq_len(n), n = n))
  }
  if (length(.variables) == 1) {
    return(id_var(.variables[[1]], drop = drop))
  }
  ids <- rev(lapply(.variables, id_var, drop = drop))
  p <- length(ids)
  ndistinct <- vapply(ids, attr, "n", FUN.VALUE = numeric(1), 
                      USE.NAMES = FALSE)
  n <- prod(ndistinct)
  if (n > 2^31) {
    char_id <- do.call("paste", c(ids, sep = "\r"))
    res <- match(char_id, unique(char_id))
  }
  else {
    combs <- c(1, cumprod(ndistinct[-p]))
    mat <- do.call("cbind", ids)
    res <- c((mat - 1L) %*% combs + 1L)
  }
  attr(res, "n") <- n
  if (drop) {
    id_var(res, drop = TRUE)
  }
  else {
    structure(as.integer(res), n = attr(res, "n"))
  }
}
.dims <- function (x) {length(amv_dim(x))}

amv_dim <- function (x) {
  if (is.vector(x)) length(x) else dim(x)}

amv_dimnames <-function (x) 
{
  d <- if (is.vector(x)) 
    list(names(x))
  else dimnames(x)
  if (is.null(d)) 
    d <- rep(list(NULL), .dims(x))
  null_names <- which(unlist(lapply(d, is.null)))
  d[null_names] <- lapply(null_names, function(i) seq_len(amv_dim(x)[i]))
  d
}
list_to_array <- function (res, labels = NULL) 
  {
    if (length(res) == 0) 
      return(vector())
    n <- length(res)
    # atomic <- sapply(res, is.atomic)
    # if (all(atomic) || all(!atomic)) {
      # dlength <- unique.default(llply(res, .dims))
      # if (length(dlength) != 1)
      #   stop("Results must have the same number of dimensions.")
      # dims <- unique(do.call("rbind", llply(res, amv_dim)))
      # if (is.null(dims))
      #   stop("Results must have one or more dimensions.",
      #        call. = FALSE)
      # if (nrow(dims) != 1)
      #   stop("Results must have the same dimensions.", call. = FALSE)
      res_dim <- amv_dim(res[[1]])
      res_labels <- amv_dimnames(res[[1]])
      if (any(vapply(res_labels, anyDuplicated, integer(1)) != 
              0)) {
        warning("Duplicate names not supported.", call. = FALSE)
      }
      res_index <- expand.grid(res_labels)
      res <- unlist(res, use.names = FALSE, recursive = FALSE)
    # }
    # else {
    #   stop("Results must have compatible types.")
    # }
    if (is.null(labels)) {
      labels <- data.frame(X = seq_len(n))
      in_labels <- list(NULL)
      in_dim <- n
    }
    else {
      in_labels <- lapply(labels, function(x) if (is.factor(x)) 
        levels(x)
        else sort(unique(x), na.last = TRUE))
      in_dim <- sapply(in_labels, length)
    }
    index_old <- rep(id(rev(labels)), each = nrow(res_index))
    index_new <- rep(id(rev(res_index)), nrow(labels))
    index <- (index_new - 1) * prod(in_dim) + index_old
    out_dim <- unname(c(in_dim, res_dim))
    out_labels <- c(in_labels, res_labels)
    n <- prod(out_dim)
    if (length(index) < n) {
      overall <- match(1:n, index, nomatch = NA)
    }
    else {
      overall <- order(index)
    }
    out_array <- res[overall]
    dim(out_array) <- out_dim
    dimnames(out_array) <- out_labels
    out_array
  }


.aaply <- function(x, margin, FUNC){
  res <- unlist(apply(x, margin, function(i) {
    
    list(FUNC(i))
  }), recursive = F) 
  list_to_array(res)
}