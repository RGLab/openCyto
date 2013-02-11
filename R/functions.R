#' @param cutpoint_method How should the cutpoint be chosen from the fitted
#' \code{flowClust} model? See Details.
#' @param neg_cluster integer. The index of the negative cluster. The cutpoint
#' is computed between clusters \code{neg_cluster} and \code{neg_cluster + 1}.
#' @param truncate_min Truncate observations less than this minimum value. By
#' default, this value is \code{NULL} and is ignored.
#' @param truncate_max Truncate observations greater than this maximum value. By
#' default, this value is \code{NULL} and is ignored.
#' @param quantile the quantile for which we will find the cutpoint using
#' the quantile \code{cutpoint_method}. If the \code{cutpoint_method} is not set
#' to \code{quantile}, this argument is ignored.
#' @param quantile_interval a vector of length 2 containing the end-points of the
#' interval of values to find the quantile cutpoint. If the
#' \code{cutpoint_method} is not setto \code{quantile}, this argument is ignored.
#' @param ... additional arguments that are passed to \code{flowClust}
#' @return a \code{rectangleGate} object ranging consisting of all values greater
#' than the cutpoint determined
flowClust.1d <- function(fr, params, filterId = "", K = 2,  adjust = 1, trans = 0,
		positive = TRUE, plot = FALSE, usePrior = 'no', prior = NULL,
		cutpoint_method = c("boundary", "min_density", "quantile", "posterior_mean"),
		neg_cluster = 1, truncate_min = NULL, truncate_max = NULL,
		quantile = 0.99, quantile_interval = c(0, 10), ...) {
	
	cutpoint_method <- match.arg(cutpoint_method)
	
	if (neg_cluster + 1 > K) {
		stop("The value for K specified is larger than the index of the positive cluster.")
	}
	
	# If appropriate, we generate prior parameters for the Bayesian version of flowClust.
	if (usePrior == "yes" && is.null(prior)) {
		prior <- prior_flowClust1d(fr = fr, channel = params[1], K = K, adjust = adjust)
	}
	
	# HACK: Circumvents a bug in flowClust.
	if (missing(prior) || is.null(prior)) {
		prior <- list(NA)
	}
	
	# If a truncation value is specified, we remove all observations less than
	# (greater than) this value for the marker specified to construct the gate.
	# NOTE: These observations are removed from the 'flowFrame' locally and are gated
	# out only for the determining the gate.
	if (!is.null(truncate_min) || !is.null(truncate_max)) {
		fr <- truncate_flowframe(flow_frame = fr, channel = params[1],
				min = truncate_min, max = truncate_max)
	}
	
	# Applies `flowClust` to the feature specified in the `params` argument using
	# the data given in `fr`. We use priors with hyperparameters given by the
	# elements in the list `prior_list`.
	filter1 <- tmixFilter(filterId, params[1], K = K, trans = trans,
			usePrior = usePrior, prior = prior, ...)
	tmixRes1 <- filter(fr, filter1)
	
	# To determine the cutpoint, we first sort the centroids so that we can determine
	# the second largest centroid.
	centroids_sorted <- sort(getEstimates(tmixRes1)$locations)
	
	# Also, because the cluster labels are arbitrary, we determine the cluster
	# the cluster labels, sorted by the ordering of the cluster means.
	labels_sorted <- order(getEstimates(tmixRes1)$locations)
	
	# Grabs the data matrix that is being gated.
	x <- exprs(fr)[, params[1]]
	
	# Determines the cutpoint between clusters 1 and 2.
	if (cutpoint_method == "boundary") {
		# Choose the cutpoint as the boundary between the first two clusters.
		# First, we sort the data.
		order_x <- order(x)
		x_sorted <- x[order_x]
		labels <- Map(tmixRes1, rm.outliers = FALSE)[order_x]
		
		# Determine which observations are between the first two centroids and their
		# corresponding cluster labels.
		which_between <- which(centroids_sorted[neg_cluster] < x_sorted & x_sorted < centroids_sorted[neg_cluster + 1])
		x_between <- x_sorted[which_between]
		labels_between <- labels[which_between]
		
		# For the observations between the first two centroids, we find the last
		# observation that belongs to the first cluster and the first observation
		# that belongs to the second cluster. In the rare occurrence that no
		# observations from one of the clusters is between the labels, we set the
		# max/min index as NA. This results in the cutpoint being set to the max/min
		# observation on the boundary.
		which_cluster1 <- which(labels_between == labels_sorted[neg_cluster])
		which_cluster2 <- which(labels_between == labels_sorted[neg_cluster + 1])
		max_obs_cluster1 <- ifelse(length(which_cluster1) == 0, NA, max(which_cluster1))
		min_obs_cluster2 <- ifelse(length(which_cluster2) == 0, NA, min(which_cluster2))
		
		# We define the cutpoint to be the midpoint between the two clusters.
		cutpoint <- mean(x_between[c(max_obs_cluster1, min_obs_cluster2)], na.rm = TRUE)
		if (is.nan(cutpoint)) {
			cutpoint <- centroids_sorted[neg_cluster]
		}
	} else if (cutpoint_method == "min_density") {
		# Determine the minimum density value of the observations between clusters 1 and 2.
		# Sets the cutpoint at the observation that attained this minimum density.
		x_between <- x[centroids_sorted[neg_cluster] < x & x < centroids_sorted[neg_cluster + 1]]
		x_dens <- dmvtmix(x = as.matrix(x_between), object = tmixRes1)
		cutpoint <- x_between[which.min(x_dens)]
	} else if (cutpoint_method == "quantile") {
		cutpoint <- quantile_flowClust(p = quantile, object = tmixRes1, interval = quantile_interval)
	} else { # cutpoint_method == "posterior_mean"
		cutpoint <- centroids_sorted[neg_cluster]
	}
	
	# After the 1D cutpoint is set, we set the gate coordinates used in the
	# rectangleGate that is returned. If the `positive` argument is set to TRUE,
	# then the gate consists of the entire real line to the right of the cut point.
	# Otherwise, the gate is the entire real line to the left of the cut point.
	if (positive) {
		gate_coordinates <- list(c(cutpoint, Inf))
	} else {
		gate_coordinates <- list(c(-Inf, cutpoint))
	}
	
	names(gate_coordinates) <- params
	
	fres <- rectangleGate(gate_coordinates, filterId = filterId)
	
	if (plot) {
		gate_pct <- round(100 * mean(x > cutpoint), 3)
		plot_title <- paste0(filterId, " (", gate_pct, "%)")
		plot(fr, tmixRes1, main = plot_title, labels = FALSE)
		abline(v = centroids_sorted, col = rainbow(K))
		abline(v = cutpoint, col = "black", lwd = 3, lty = 2)
		
		if (!is.null(prior)) {
			x_dens <- seq(min(x), max(x), length = 1000)
			
			for(k in seq_len(K)) {
				prior_density <- with(prior, dnorm(x_dens, mean = Mu0[k], sd = sqrt(Omega0[k])))
				
				# Grab posterior estimates for the degrees of freedom (nu) and the
				# transformation parameters. Because these can either be of length 1 or
				# of length K, we grab the appropriate posterior estimates for the
				# current mixture component.
				nu <- ifelse(length(tmixRes1@nu) == 1, tmixRes1@nu, tmixRes1@nu[k])
				lambda <- ifelse(length(tmixRes1@lambda) == 1, tmixRes1@lambda, tmixRes1@lambda[k])        
				posterior_density <- dmvt(x_dens, mu = tmixRes1@mu[k,],
						sigma = tmixRes1@sigma[k,,], nu = nu,
						lambda = lambda)$value
				lines(x_dens, prior_density, col = rainbow(K)[k], lty = 2, lwd = 1)
				lines(x_dens, posterior_density, col = rainbow(K)[k], lwd = 1)
			}
		}
	}
	
	fres
}

#' Applies flowClust to two features in a flowFrame to construct an elliptical
#' gate.
#'
#' We cluster the observations in \code{fr} into \code{K} clusters. We set the
#' cutpoint to be the point at which the density between the first and second
#' smallest cluster centroids is minimum.
#'
#' @param fr a \code{flowFrame} object
#' @param xChannel TODO
#' @param yChannel TODO
#' @param filterId TODO
#' @param K the number of clusters to find
#' @param usePrior Should we use the Bayesian version of \code{flowClust}?
#' Answers are "yes", "no", or "vague". The answer is passed along to
#' \code{flowClust}.
#' @param prior_list list of prior parameters for the Bayesian \code{flowClust}.
#' If \code{usePrior} is set to 'no', then the list is unused.
#' @param trans numeric indicating whether the Box-Cox transformation parameter
#' is estimated from the data. May take 0 (no estimation), 1 (estimation) or 2
#' (cluster-speciﬁc estimation). NOTE: For the Bayesian version of
#' \code{flowClust}, this value cannot be 2.
#' @param plot a logical value indicating if the fitted mixture model should be
#' plotted. By default, no.
#' @param truncate_min A vector of length 2. Truncate observations less than this
#' minimum value. The first value truncates the \code{xChannel}, and the second
#' value truncates the \code{yChannel}. By default, this vector is \code{NULL}
#' and is ignored.
#' @param truncate_max A vector of length 2. Truncate observations greater than
#' this maximum value. The first value truncates the \code{xChannel}, and the
#' second value truncates the \code{yChannel}. By default, this vector is
#' \code{NULL} and is ignored.
#' @param ... additional arguments that are passed to \code{flowClust}
#' @return a \code{polygonGate} object containing the contour (ellipse) for 2D
#' gating.
flowClust.2d <- function(fr, xChannel, yChannel, filterId = "", K = 2,
		usePrior = 'no', prior_list = list(NA), trans = 0,
		plot = FALSE, which_gate = c("bottom", "top", "left", "right"),
		gate_type = c("ellipse", "axis"), quantile = 0.995,
		truncate_min = NULL, truncate_max = NULL, ...) {
	
	which_gate <- match.arg(which_gate)
	gate_type <- match.arg(gate_type)
	
	# If a truncation value is specified, we remove all observations less than this
	# value for the marker specified to construct the gate.
	# NOTE: These observations are removed from the 'flowFrame' locally and are gated
	# out only for the determining the gate.
	if (!is.null(truncate_min)) {
		exprs(fr) <- exprs(fr)[exprs(fr)[, xChannel] >= truncate_min[1], ]
		exprs(fr) <- exprs(fr)[exprs(fr)[, yChannel] >= truncate_min[2], ]
	}
	
	if (!is.null(truncate_max)) {
		exprs(fr) <- exprs(fr)[exprs(fr)[, xChannel] <= truncate_max[1], ]
		exprs(fr) <- exprs(fr)[exprs(fr)[, yChannel] <= truncate_max[2], ]    
	}
	
	x <- exprs(fr)[, xChannel]
	y <- exprs(fr)[, yChannel]
	
	# If appropriate, we generate prior parameters for the Bayesian version of flowClust.
	if (usePrior == "yes" && identical(prior_list, list(NA))) {
		prior_list <- prior_flowClust2d(fr = fr, xChannel = xChannel, yChannel = yChannel, K = K)
	}
	
	
	# Applies `flowClust` to the feature specified in the `params` argument using
	# the data given in `fr`. We use priors with hyperparameters given by the
	# elements in the list `prior_list`.
	filter1 <- tmixFilter(filterId, c(xChannel, yChannel), K = K, trans = trans,
			usePrior = usePrior, prior = prior_list, ...)
	tmixRes1 <- filter(fr, filter1)
	
	# Converts the tmixFilterResult object to a polygonGate.
	# We select the cluster with the minimum 'yChannel' to be the subpopulation from
	# which we obtain the contour (ellipse) to generate the polygon gate.
	fitted_means <- getEstimates(tmixRes1)$locations
	cluster_selected <- switch(which_gate,
			bottom = which.min(fitted_means[, 2]),
			top = which.max(fitted_means[, 2]),
			left = which.min(fitted_means[, 1]),
			right = which.max(fitted_means[, 1]))    
	
	if (gate_type == "ellipse") {
		contour_ellipse <- .getEllipse(filter = tmixRes1,
				include = cluster_selected)
		flowClust_gate <- polygonGate(.gate = matrix(contour_ellipse, ncol = 2,
						dimnames = list(NULL, tmixRes1@varNames)),
				filterId = filterId)
	} else if (gate_type == "axis") {
		chisq_quantile <- qchisq(quantile, df = 2)
		
		xbar <- tmixRes1@mu[cluster_selected, ]
		Sigma <- tmixRes1@sigma[cluster_selected, , ]
		
		Sigma_eigen <- eigen(Sigma, symmetric = TRUE)
		u1 <- Sigma_eigen$vectors[, 1]
		u2 <- Sigma_eigen$vectors[, 2]
		lambda1 <- Sigma_eigen$values[1]
		lambda2 <- Sigma_eigen$values[2]
		
		# Determines which eigenvector points towards the first quadrant (has
		# positive slope). Because both eigenvectors can potentially point in the
		# negative direction, we also check to see the negated eigenvectors point
		# towards the first quadrant.
		if (all(u1 >= 0)) {
			axis <- sqrt(lambda1 * chisq_quantile) * u1
			axis_perp <- sqrt(lambda2 * chisq_quantile) * u2
		} else if (all(-u1 >= 0)) {
			axis <- -sqrt(lambda1 * chisq_quantile) * u1
			axis_perp <- sqrt(lambda2 * chisq_quantile) * u2
		} else if (all(u2 >= 0)) {
			axis <- sqrt(lambda2 * chisq_quantile) * u2
			axis_perp <- sqrt(lambda1 * chisq_quantile) * u1
		} else if (all(-u2 >= 0)) {
			axis <- -sqrt(lambda2 * chisq_quantile) * u2
			axis_perp <- sqrt(lambda1 * chisq_quantile) * u1
		}
		
		# The gate location is the frame of reference for the gate. If it is xbar,
		# then the frame of reference is the cross-section along the eigenvector
		# from top-left to bottom-right. We translate this reference as a function
		# of the appropriate chi-squared coefficient.
		gate_location <- xbar + 0.25 * axis
		
		# To construct the gate, we have effectively shifted the eigenvector with the
		# negative slope. We then extend the gate both horizontally and vertically
		# to the maximum observed values in the horizontal and vertical directions.
		# NOTE: We extend the gate one standard deviation beyond the maximum values
		# observed to mimic a gate that extends without limit in the positive
		# directions. However, because flowCore cannot handle such a shape, we force
		# the gate to have the same shape for the observed data.
		x_max <- max(x) + sd(x)
		y_max <- max(y) + sd(y)
		gate_x <- c(x_max, (gate_location - axis_perp)[2])
		gate_y <- c((gate_location + axis_perp)[1], y_max)
		
		# We extend the gate to the min and max values
		polygon_gate <- rbind(gate_location + axis_perp,
				gate_y,
				c(x_max, y_max),
				gate_x,
				gate_location - axis_perp)
		
		colnames(polygon_gate) <- c(xChannel, yChannel)
		flowClust_gate <- polygonGate(filterId = filterId, boundaries = polygon_gate)
	}
	
	if (plot) {
		plot(fr, tmixRes1, main = filterId)
		
		if (gate_type == "axis") {
			# The major and minor axes (eigenvectors) scaled by their respective
			# eigenvalues and the chi-squared quantile.
			lines(rbind(xbar, xbar + sqrt(lambda1 * chisq_quantile) * u1), col = "darkgreen")
			lines(rbind(xbar, xbar - sqrt(lambda1 * chisq_quantile) * u1), col = "darkgreen")
			lines(rbind(xbar, xbar + sqrt(lambda2 * chisq_quantile) * u2), col = "darkgreen")
			lines(rbind(xbar, xbar - sqrt(lambda2 * chisq_quantile) * u2), col = "darkgreen")
			
			# Draws the polygon gate.
			lines(polygon_gate, col = "red")
			lines(rbind(gate_location - axis_perp, gate_location + axis_perp), col = "red")
			
			# Also, draws points at the vertices of the polygon gate.
			points(polygon_gate, col = "red", pch = 16)
		}
	}
	
	flowClust_gate
}

read.FCS.csv<-function(file,stains=NA){
	mat<-as.matrix(read.csv(file,check.names=FALSE))
	
	fr<-new("flowFrame",exprs=mat)
	
	pd<-pData(parameters(fr))
	pd$desc<-as.character(pd$desc)
	pd$name<-as.character(pd$name)
#	browser()
	##update the desc with marker name
	if(!is.na(stains))
	{
		ind<-match(names(stains),pd$name)
		pd[ind,]$desc<-as.character(stains)
		##update SSC and FSC description with NA
		ind<-grepl("[F|S]SC",pd$desc)
		pd[ind,]$desc<-NA
	}else
		pd$desc<-NA
#	{
#		#when statins is NA,indicates that marker name is concatenated with channel name
#		#then we strip marker name from colnames
#		lapply(pd$name[1],function(x){
#					browser()
#					substr(x)
#				})
#	}
	
	
	
	##update minRange with -111 for proper display of the data
#	browser()
	pd$minRange[pd$minRange<(-111)]<--111
	pData(parameters(fr))<-pd
	fr
}

read.flowSet.csv<-function(files,...){
#	browser()
	fs<-flowSet(lapply(files,read.FCS.csv,...))
	sampleNames(fs)<-basename(files)
	fs
}
spillover1<-function(x, unstained=NULL, cols=NULL, fsc="FSC-A",
				ssc="SSC-A", method="median", useNormFilt=FALSE,isOrdered=FALSE)
		{
			
			if(is.null(unstained)) {
				stop("Sorry, we don't yet support unstained cells blended ",
						"with stained cells", call.=FALSE)
			} else {
				## We often only want spillover for a subset of the columns 
				allcols <- colnames(x)
#				cols <- if(is.null(patt)) allcols else grep(patt, allcols,
#									value=TRUE)
				
				## Ignore these guys if they somehow got into cols.
				## cols <- cols[-match(c(fsc,ssc),cols)]
				cols <- cols[!(cols %in% c(fsc,ssc))]
				
				## There has got to be a better way of doing this...
				if(!is.numeric(unstained)) {
					unstained <- match(unstained,sampleNames(x))
					if(is.na(unstained))
						stop("Baseline not in this set.", call.=FALSE)
				}
				## Check to see if the unstained sample is in the list of
				## stains. If not, we need to add it, making it the first
				## row and adjust the unstained index accordingly.
				## If it is there we adjust to the appropriate index.
#				browser()
				## pdh: you shouldn't use the nor2Filter as a default without telling people!
				if(useNormFilt){
					if(is.numeric(fsc)) fsc <- allcols[fsc]
					if(is.numeric(ssc)) ssc <- allcols[ssc]
					
					if(is.na(match(fsc,allcols)))
						stop("Could not find forward scatter parameter. ",
								"Please set the fsc parameter", call.=FALSE)
					if(is.na(match(ssc,allcols)))
						stop("Could not find side scatter parameter. ",
								"Please set the ssc parameter", call.=FALSE)
					n2f <- norm2Filter(fsc, ssc, scale.factor=1.5)
					x <- Subset(x,n2f)
				}
#				browser()
		
				#select positive population on its stained channel
				newX<-fsApply(x[-unstained],function(curFr){

							filName<-basename(keyword(curFr)$FIL)
							##validity check by matching channels with filenames
#							
							
							ind<-which(unlist(lapply(cols,function(y){
#														browser()
														#strip the -X from the end
														y<-substr(y,1,nchar(y)-2)
														y<-paste(y,"")#append space at the end
														grepl(y,filName)
													})
												)
										)
							if(length(ind)==0)
								stop(filName,"does not match any of the channels!")
							if(length(ind)>1)
								stop(filName,"matches more than one channels!")
							curChannel<-cols[ind]
#							browser()
							##rangeGate to select positive pop
							fres<-rangeGate(curFr,stain=curChannel
										,inBetween=T
#										,plot=T
										,borderQuant=0
										,absolute=F
										)
#							browser()										
							newFr<-Subset(curFr,fres)
							
#							densityplot(~.,newFr)
#							densityplot(~.,curFr)
#							median(exprs(newFr)[,curChannel])
							
							newFr
						})
				newX<-rbind2(newX,x[unstained])
#				browser()
#				lgcl_cont<-estimateLogicle(newX[[3]],channels=channels)
#				xyplot(`PerCP-Cy5-5-A`~`FITC-A`,transform(x[[3]],lgcl_cont),smooth=F,xbin=128)
#				xyplot(`SSC-A`~`FSC-A`,x,smooth=F,xbin=128)
#				CairoX11()
				
				
#				grid.arrange(
#						densityplot(~.,transform(x[[unstained]],lgcl_cont))
#						,densityplot(~.,transform(newX[[3]],lgcl_cont))
#				##							,xyplot(`APC-Cy7-A`~`APC-A`,transform(x[[1]],lgcl_1),smooth=F,xbin=128)
##						
#				)
				
				if(method=="mode")
				{
					inten<-fsApply(newX,function(curFr){
	#								browser()
									modes<-sapply(cols,function(curStain){
	#											browser()
												sig<-exprs(curFr)[, curStain]
												
												
	#											fres<-filter(curFr,curv1Filter(curStain,bwFac=1.2))
	#											bnds <- flowStats:::curvPeaks(fres,sig,borderQuant=0)
	#											curMode<-as.numeric(bnds$peaks[which.max(bnds$peaks[,"y"]),"x"])
												
												res<-density(sig)
												curMode<-res$x[which.max(res$y)]
												
	#											print(densityplot(as.formula(paste("~`",curStain,"`",sep="")),curFr,refline=curMode))
	#											hist(sig,breaks=1000)
	#											abline(v=curMode,col="red")
												curMode
	
											},USE.NAMES=T)
	#								filterList<-sapply(names(modes),function(x){
	#											g<-rectangleGate(list(x=c(modes[x],Inf)))
	#											parameters(g)<-x
	#											g
	#										})
									modes
								})
					
						
				}else
				{
					inten <- fsApply(newX, each_col,method)[, cols]	
				}
				
				#background correction
				inten <- pmax(sweep(inten[-unstained,], 2,inten[unstained,]), 0)
				#normalize by max of each row
				inten <- sweep(inten, 1,apply(inten, 1, max), "/")
#				browser()
				#if the files is already ordered by channels
				#which means we know which channel is stained for each control file
				#we don't need to guess it by pmax 
				if(isOrdered)
				{
					row.names(inten) <- channels
				}else
				{
					#guessing row names by picking the colname which has maximun value  
					row.names(inten) <- colnames(inten)[apply(inten ,1,which.max)]	
					
				}
				
				
				inten[colnames(inten),]
			}
		}



#' Creates a matrix of points on an ellipse from a fitted flowClust model.
#'
#' The ellipse is constructed from a contour from the fitted flowClust model.
#' The contour level is specified with \code{level}.
#'
#' @param filter object containing the fitted flowClust model.
#' @param include the mixture component in the fitted flowClust model for which
#' the contour (ellipse) is returned
#' @param ecol TODO
#' @param elty TODO
#' @param level the contour level of the ellipse
#' @param npoints the number of points on the ellipse
#' @param subset the dimensions of the mixture component to return
#' @return matrix containing the points of the ellipse from the flowClust contour
.getEllipse <- function(filter = NULL, include = seq_len(filter@K), ecol = 1,
                        elty = 1, level = NULL, npoints = 501, subset = c(1, 2)) {
  py <- 2
  ecol <- matrix(ecol, length(include))
  elty <- matrix(elty, length(include))

  if (all(filter@nu != Inf)) {
    if (filter@ruleOutliers[1] == 0) { # 0 means quantile
      cc <- py * qf(filter@ruleOutliers[2], py, filter@nu)
    } else { # 1 means u.cutoff
      cc <- ((filter@nu + py) / filter@ruleOutliers[2] - filter@nu)    
    }
  } else {
    cc <- qchisq(filter@ruleOutliers[2], py)
  }
	
	j <- 0
	if (length(filter@lambda) > 0) {
    lambda <- rep(filter@lambda, length.out = filter@K)
  }
	else {
    lambda <- numeric(0)
  }
	cc <- rep(cc, length.out = filter@K)
  for (i in include) {
    eigenPair <- eigen(filter@sigma[i, subset, subset])
    l1 <- sqrt(eigenPair$values[1]) * sqrt(cc)
    l2 <- sqrt(eigenPair$values[2]) * sqrt(cc)
    angle <- atan(eigenPair$vectors[2, 1] / eigenPair$vectors[1, 1]) * 180 / pi

    if (length(lambda) > 0) {
      res <- rbox(flowClust:::.ellipsePoints(a = l1[i], b = l2[i], alpha = angle,
                                             loc = filter@mu[i, subset],
                                             n = npoints),
                  lambda[i])
    } else {
      res <- flowClust:::.ellipsePoints(a = l1[i], b = l2[i], alpha = angle,
                                        loc = filter@mu[i, subset], n = npoints)
    }
  }
  res
}



getChannelMarker<-function(frm,name,fix=FALSE)
{
	#try stain name
	pd<-pData(parameters(frm))
	pname<-paste(name,"([ ]|$)",sep="")
#	browser()
	if(fix)
		ind<-which(toupper(pd$name)%in%toupper(name))
	else
		ind<-which(grepl(pname,pd$name,ignore.case=T))
		
	
	
		
	if(length(ind)==0)
	{
		#try marker name
		ind<-which(unlist(lapply(pd$des,function(x){
								#split by white space and then match each individual string
#									browser()
									if(fix)
										any(unlist(lapply(strsplit(x," "),function(y)toupper(y)%in%toupper(name))))
									else
										grepl(pattern=pname,x,ignore.case=T)
							})
						)
					)
		if(length(ind)==0)
			stop("can't find ",name)
		if(length(ind)>1)
			stop("multiple markers matched: ",name)
	}
	
	pd[ind,c("name","desc")]
}







#' For the given workflow, we look up the given markers and return the
#' corresponding channels.
#'
#' @param flow_frame object of type \code{flowFrame}
#' @param markers the markers from which we obtain the corresponding channel names
#' @return vector of channel names
markers2channels <- function(flow_frame, markers) {
  # First, we build a lookup table for the channels and markers.
  channel_markers <- lapply(colnames(flow_frame), function(channel) {
    marker_name_desc <- getChannelMarker(flow_frame, channel)
    marker <- with(marker_name_desc, ifelse(is.na(desc), name, desc))
    cbind(channel, marker = unname(marker))
  })
  channel_markers <- data.frame(do.call(rbind, channel_markers),
                                stringsAsFactors = FALSE)

  # Now, we query the channels for the specified markers.
  channels <- sapply(markers, function(marker) {
    channel_markers$channel[grepl(marker, channel_markers$marker)]
  })

  as.vector(channels)
}

#' For the given flow frame, we look up the given markers and return the
#' corresponding channels.
#'
#' @param flow_frame object of type \code{flowFrame}
#' @param channels the channels from which we obtain the corresponding markers
#' @return vector of markers
channels2markers <- function(flow_frame, channels) {
  markers <- sapply(channels, function(channel) {
    marker <- getChannelMarker(flow_frame, channel)
    with(marker, ifelse(is.na(desc), name, desc))
  })
  unname(markers)
}

#' Removes any observation from the given flow frame that has values less than
#' (greater than) the minimum (maxim) value.
#'
#' The minimum/maximum values are ignored if \code{NULL}.
#'
#' @param flow_frame an object of type \code{flowFrame}
truncate_flowframe <- function(flow_frame, channel, min = NULL, max = NULL) {
  if (is.null(min) && is.null(max)) {
    warning("No truncation value was provided. Returning the original 'flow_frame'.")
  }
  x_channel <- exprs(flow_frame)[, channel]

  # For comparison purposes, we update the min and max values to -Inf and Inf, respectively,
  # if one is NULL.
  min <- ifelse(is.null(min), -Inf, min)
  max <- ifelse(is.null(max), Inf, max)

  # Removes any observation that has an observation outside of the min and max
  # values specified.
  exprs(flow_frame) <- exprs(flow_frame)[min < x_channel & x_channel < max, ]

  flow_frame
}

#' Computes the quantile from flowClust for a given vector of probabilties
#' 
#' We estimate the quantile from a \code{flowClust} fit with a combination of
#' numerical integration and a root-finding method. We are effectively
#' estimating the cumulative distribution function (CDF) of the mixture density
#' estimated by \code{flowClust}.
#'
#' Because we are using numerical methods, we also need an \code{interval} of
#' values in which we will attempt to find the specified quantile.
#'
#' @param p vector of probabilities
#' @param object an object containing the \code{flowClust} fit
#' @param interval a vector of length 2 containing the end-points of the interval
#' of values to find the quantile
#' @param ... Additional arguments that are passed to \code{uniroot} to find the
#' quantile.
#' @return the quantile corresponding to the specified probabilities
quantile_flowClust <- function(p, object, interval, ...) {
  cdf_target <- function(x, p, object) {
    cdf_values <- sapply(seq_len(object@K), function(k) {
      nu <- ifelse(length(object@nu) == 1, object@nu, object@nu[k])
      lambda <- ifelse(length(object@lambda) == 1, object@lambda, object@lambda[k])        

      # TODO: Incorporate the Box-Cox transformation (i.e., box(qt(...), lambda = lambda)) into quantile
      # The case of 'lambda = 1' may be not be trivial -- this case is largely ignored in flowClust.
      pt((x - object@mu[k]) / sqrt(object@sigma[k]), df = nu)
    })
    weighted.mean(cdf_values, w = object@w) - p
  }

  uniroot(cdf_target, interval = interval, p = p, object = object, ...)$root
}

#' Extracts the quadrants of a quadGate as a list of rectangleGates
#'
#' The quadrants are numbered in a clockwise manner with the top-left quadrant
#' numbered 1, the top-right quadrant numbered 2, and so on.
#' 
#' @param quad_gate a \code{quadGate} object
#' @param markers character vector of the marker names for the x- and y-axes
#' @param channels character vector of the channel names for the x- and y-axes
#' @param quadrants a vector indicating the quadrants to extract
#' @return a \code{filters} object containing a list of the rectangle gates
quadGate2rectangleGates <- function(quad_gate, markers, channels, quadrants = 1:4) {
  x_gate <- quad_gate@boundary[1]
  y_gate <- quad_gate@boundary[2]

  gates_list <- list()
  
  # Top-left quadrant
  gates <- list(c(-Inf, x_gate), c(y_gate, Inf))
  names(gates) <- channels
  gates_list[[paste0(markers[1], "-", markers[2], "+")]] <- rectangleGate(gates)

  # Top-right quadrant
  gates <- list(c(x_gate, Inf), c(y_gate, Inf))
  names(gates) <- channels
  gates_list[[paste0(markers[1], "+", markers[2], "+")]] <- rectangleGate(gates)

  # Lower-right quadrant
  gates <- list(c(x_gate, Inf), c(-Inf, y_gate))
  names(gates) <- channels
  gates_list[[paste0(markers[1], "+", markers[2], "-")]] <- rectangleGate(gates)

  # Lower-left quadrant
  gates <- list(c(-Inf, x_gate), c(-Inf, y_gate))
  names(gates) <- channels
  gates_list[[paste0(markers[1], "-", markers[2], "-")]] <- rectangleGate(gates)

  filters(gates_list[quadrants])
}

#' Constructs a vector of all the combinations of A & B & C
#'
#' The \code{permutations} function is from the \code{gregmisc} package on CRAN.
#' @param markers character vector of marker names
#' @return vector containing all combinations of the markers
#' @examples
#' polyfunction_nodes(c("IFNg", "IL2", "TNFa", "GzB", "CD57"))
polyfunction_nodes <- function(markers) {
  require('gregmisc')
  markers <- paste0(markers, "+")
  num_markers <- length(markers)
  and_list <- as.vector(permutations(n = 1, r = num_markers - 1, c("&"), repeats = TRUE))
  isnot_list <- permutations(n = 2, r = num_markers, c("!", ""),
                             repeats = TRUE)
  apply(isnot_list, 1, function(isnot_row) {
    isnot_row[-1] <- paste0(and_list, isnot_row[-1])
    paste(paste0(isnot_row, markers), collapse = "")
  })
}
