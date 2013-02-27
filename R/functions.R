
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
#'
#' By default, the contour level is extracted from the \code{ruleOutliers} slot
#' in the \code{filter}. The user can override this level by specifying value
#' between 0 and 1 in \code{quantile}.
#'
#' @param filter object containing the fitted \code{flowClust} model.
#' @param include the mixture component in the fitted \code{flowClust} model for
#' which the contour (ellipse) is returned
#' @param ecol TODO
#' @param elty TODO
#' @param quantile the contour level of the ellipse. See details.
#' @param npoints the number of points on the ellipse
#' @param subset the dimensions of the mixture component to return
#' @return matrix containing the points of the ellipse from the flowClust contour
.getEllipse <- function(filter = NULL, include = seq_len(filter@K), ecol = 1,
                        elty = 1, quantile = NULL, npoints = 501, subset = c(1, 2)) {

  # Sets the quantile of the ellipse.
  if (is.null(quantile)) {
    quantile <- filter@ruleOutliers[2]
  } else {
    if (!is.numeric(quantile) || quantile < 0 || quantile > 1) {
      stop("The 'quantile' must be a numeric value between 0 and 1.")
    }
  }

  # py is the degrees of freedom?
  py <- 2
  ecol <- matrix(ecol, length(include))
  elty <- matrix(elty, length(include))

  if (all(filter@nu != Inf)) {
    if (filter@ruleOutliers[1] == 0) { # 0 means quantile
      cc <- py * qf(p = quantile, py, filter@nu)
    } else { # 1 means u.cutoff
      cc <- ((filter@nu + py) / quantile - filter@nu)
    }
  } else {
    cc <- qchisq(p = quantile, py)
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

getChannelMarker <- function(frm, name, fix = FALSE) {
  #try stain name
  pd<-pData(parameters(frm))
  pname<-paste(name,"([ ]|$)",sep="")
  if (fix) {
    ind <- which(toupper(pd$name) %in% toupper(name))
  } else {
    ind <- which(grepl(pname, pd$name, ignore.case = TRUE))
  }

  if(length(ind) == 0) {
    # try marker name
    ind <- which(unlist(lapply(pd$des, function(x) {
      # split by white space and then match each individual string
      if (fix) {
        any(unlist(lapply(strsplit(x, " "), function(y) {
          toupper(y) %in% toupper(name)
        })))
      } else {
        grepl(pattern=pname,x,ignore.case=T)
      }
    })))

    if (length(ind) == 0) {
      stop("can't find ", name)
    } else if(length(ind) > 1) {
      stop("multiple markers matched: ", name)
    }
  }

  pd[ind, c("name", "desc")]
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
