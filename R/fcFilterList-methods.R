#' plot a \code{fcFilterList}
#' 
#' It is usually called by \code{plot} method for \code{fcTree} instead of directly by users.
#' 
#' @param x \code{fcFilterList}
#' @param y \code{character} channel name
#' @param samples \code{character} a vector of sample names to be plotted
#' @param posteriors \code{logical} indicating whether posteriors should be plotted
#' @param xlim,ylim scale settings for x,y axises
#' @param lwd line width
#' @param breaks passed to \link{hist}
#' @param data \code{GatingSet} object
#' @param node \code{character} population name associated with the \code{fcFilterList}
#' @param ... other arguments passed to base \code{plot}
#' @examples
#' \dontrun{ 
#'  env1<-new.env(parent=emptyenv())
#'  #gt is a gatingTemplate, gs is a GatingSet
#'  gt_gating(gt,gs,env1) #the flowClust gating results are stored in env1 
#'  plot(env1$fct,"nonDebris",post=T) #plot the priors as well as posteriors for the "nonDebris" gate
#' }
setMethod("plot", sig = c("fcFilterList", "ANY"),
          definition = function(x, y, samples = NULL, posteriors = FALSE
        , xlim = NULL, ylim = NULL, node = NULL, data = NULL, breaks = 20, lwd = 1, ...) {
#            browser()
  prior1 <- priors(x[[1]], y)
  
  if (is.null(y)) {
    priorNames <- names(prior1)
    if (length(prior1) == 1) {
      y <- priorNames
    } else {
      stop("Need to specify which channel:", paste(priorNames, collapse = "or"))
    }
  }

  # refetch prior according to modified y
  prior1 <- priors(x[[1]], y)  

  # try to set x,y lim from
  minX <- min(unlist(lapply(x, function(curFilter) posteriors(x = curFilter, y = y)$min)))
  maxX <- max(unlist(lapply(x, function(curFilter) posteriors(x = curFilter, y = y)$max)))
  x_dens <- seq(minX, maxX, length = 1000)

  # TODO: right now the priors from the first sample is used to determine y scale
  K <- nrow(prior1$Mu0)
  prior_density <- lapply(seq_len(K), function(k) dnorm(x_dens, mean = prior1$Mu0[k], 
    sd = sqrt(prior1$Omega0[k])))
  minY <- min(unlist(lapply(prior_density, min)))
  maxY <- max(unlist(lapply(prior_density, max)))
  if(is.null(xlim)){
    xlim <- c(minX, maxX)
  }
  if(is.null(ylim)){
    ylim <- c(minY, maxY)
  }
  plot(x = NULL, type = "n", xlim = xlim, ylim = ylim, xlab = y, 
    ylab = "", ...)
  
  # plot post
  if (posteriors) {
    if (is.null(samples)) 
      samples <- names(x)
    for (samp in samples) {
      curFilter <- x[[samp]]
#      browser()
      if(!is.null(data)){
        fr <- gh_pop_get_data(data[[samp]],node)
        thisChnl <- priorNames
        hist(exprs(fr)[,thisChnl],add=T,prob=T,breaks=breaks)  
      }
      curPost <- posteriors(x = curFilter, y = y)
      
      curPrior <- priors(x = curFilter, y = y)
      K <- nrow(curPrior$Mu0)
      prior_density <- lapply(seq_len(K), function(k) dnorm(x_dens, mean = curPrior$Mu0[k], 
        sd = sqrt(curPrior$Omega0[k])))
#      browser()
      for (k in seq_len(K)) {
        lines(x_dens, prior_density[[k]], col = rainbow(K)[k], lty = 2, lwd = lwd)
        
        posterior_density <- curPost$w[k]*flowClust::dmvt(x_dens, mu = curPost$mu[k, ], sigma = curPost$sigma[k, 
          , ], nu = curPost$nu, lambda = curPost$lambda)$value
        lines(x_dens, posterior_density, col = rainbow(K)[k], lwd = lwd)

        # plot gate
        g <- c(curFilter@min[y], curFilter@max[y])
        ind <- which(!is.infinite(g))
        abline(v = g[ind], col = "grey")
      }
    }
  }
  
  
  
})

