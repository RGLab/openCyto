#' Creates a singlet polygon gate using the prediction bands from a robust linear model
#'
#' We construct a singlet gate by applying a robust linear model. By default, we model the forward-scatter height
#' (FSC-H)as a function of forward-scatter area (FSC-A). If \code{sidescatter}
#' is given, forward-scatter height is as a function of \code{area} +
#' \code{sidescatter} + \code{sidescatter / area}.
#'
#' Because \code{rlm} relies on iteratively reweighted least
#' squares (IRLS), the runtime to construct a singlet gate is dependent in part
#' on the number of observations in \code{x}. To improve the runtime, we provide
#' an option to subsample randomly a subset of \code{x}. A percentage of
#' observations to subsample can be given in \code{subsample_pct}. By default, no
#' subsampling is applied.
#' 
#' @param x a \code{\link[flowCore:flowFrame-class]{flowFrame}} object
#' @param area character giving the channel name that records the signal
#' intensity as peak area
#' @param height character giving the channel name that records the signal
#' intensity as peak heightchannel name of height
#' @param sidescatter character giving an optional channel name for the
#' sidescatter signal. By default, ignored.
#' @param prediction_level a numeric value between 0 and 1 specifying the level
#' to use for the prediction bands
#' @param subsample_pct a numeric value between 0 and 1 indicating the percentage
#' of observations that should be randomly selected from \code{x} to construct
#' the gate. By default, no subsampling is performed.
#' @param wider_gate logical value. If \code{TRUE}, the prediction bands used to
#' construct the singlet gate use the robust fitted weights, which increase
#' prediction uncertainty, especially for large FSC-A. This leads to wider gates,
#' which are sometimes desired.
#' @param filterId the name for the filter that is returned
#' @param maxit the limit on the number of IWLS iterations 
#' @param ... additional arguments (not used)
#' @return a \code{\link[flowCore]{polygonGate}} object with the singlet gate
#' @rdname gate_singlet
gate_singlet <- function(x, area = "FSC-A", height = "FSC-H", sidescatter = NULL,
                        prediction_level = 0.99, subsample_pct = NULL,
                        wider_gate = FALSE, filterId = "singlet", maxit = 5, ...) {
  flowCore:::checkClass(x, "flowFrame")
  flowCore:::checkClass(area, "character")
  flowCore:::checkClass(height, "character")
  if (!is.null(sidescatter)) {
    flowCore:::checkClass(sidescatter, "character")
  }
  if (length(area) + length(height) != 2) {
    stop("Each of 'area' and 'height' must be 'character' vectors of length 1.")
  }

  x <- exprs(x[, c(area, height, sidescatter)])
  channel_names <- c(area, height)
  area <- make.names(area)
  height <- make.names(height)
  colnames(x)[1:2]  <- c(area, height)
  if (!is.null(sidescatter)) 
  {
    sidescatter <- make.names(sidescatter)
    colnames(x)[3] <- sidescatter
  }
    
  # If specified, subsample from 'x'
  if (!is.null(subsample_pct)) {
    subsample_pct <- as.numeric(subsample_pct)
    if (subsample_pct <= 0 || subsample_pct > 1) {
      warning("The subsampling percentage must be between 0 and 1. ",
              "Setting percentage to default value...", call. = FALSE)
      subsample_pct <- 1
    }
    n <- nrow(x)
    x <- x[sample(x = seq_len(nrow(x)), size = subsample_pct * n), ]
  }
  # Creates polygon gate based on the prediction bands at the minimum and maximum
  # 'area' observation using the trained robust linear model.
  which_min <- which.min(x[,area])
  which_max <- which.max(x[,area])
  x_extrema <- x[c(which_min, which_max), ]

  y <- x[, height]
  x <- x[, c(area, sidescatter), drop = FALSE]
  x <- cbind(1, x)#add default weight column
  rlm_formula <- paste(make.names(height), make.names(area), sep = " ~ ")
  
  if (!is.null(sidescatter)) {
    sidescatter <- make.names(sidescatter)
    x <- cbind(x, x[, sidescatter] / x[, area])
    ssc_ratio <- paste0("I(", sidescatter, " / ", area, ")")
    colnames(x)[4] <- ssc_ratio
    
    rlm_formula <- paste(rlm_formula, sidescatter, ssc_ratio, sep = " + ")
  }
  rlm_formula <- as.formula(rlm_formula)
  
  rlm_fit <- withCallingHandlers(fast_rlm(x, y, maxit = maxit, ...),
                                 warning = function(w){
                                   if(grepl("failed to converge", conditionMessage(w))){
                                     invokeRestart("muffleWarning")
                                   }
                                 })

  # if (!rlm_fit$converged) {
  #   warning("The IRLS algorithm employed in 'rlm' did not converge.")
  # }

  
  # Prediction weights. By default, these are weighted equally. If 'wider_gate'
  # is set, the weights from 'rlm' are used. The weight for the maximum FSC-A
  # is usually much smaller and yields a more uncertain prediction interval.
  # This leads to a wider gate.
  prediction_weights <- c(1, 1)
  if (wider_gate) {
    prediction_weights <- rlm_fit$w[c(which_min, which_max)]
  }

  #not really needed since default weights are 1s
  #but somehow this recompute qr yield slightly different predicted
  #vals (5e-7 < err < 2e-6), which is actually small enough for gates to be considered as ignorable
  #but for the sake of comparison to MASS version, we keep this step for now
  rlm_fit$weights <- rep(1, nrow(x))
  rlm_fit$qr <- qr(sqrt(rlm_fit$weights) * rlm_fit$x)
  
  #strip rlm class to dispatch stats::predict.lm
  class(rlm_fit) <- "lm"
  rlm_fit$terms <- rlm_formula
   predictions <- predict(rlm_fit, scale = rlm_fit$s, data.frame(x_extrema), interval = "prediction",
                         level = prediction_level, weights = prediction_weights)

  # Create a matrix of the vertices using the prediction bands at the minimum
  # and maximum values of x. The ordering matters. Must go clockwise.
  # Otherwise, the polygon is not convex and makes an X-shape.
  gate_vertices <- rbind(cbind(x_extrema[, area][1], predictions[1, "lwr"]),
                         cbind(x_extrema[, area][1], predictions[1, "upr"]),
                         cbind(x_extrema[, area][2], predictions[2, "upr"]),
                         cbind(x_extrema[, area][2], predictions[2, "lwr"]))
  colnames(gate_vertices) <- channel_names

  polygonGate(gate_vertices, filterId = filterId)
}

singletGate <- gate_singlet