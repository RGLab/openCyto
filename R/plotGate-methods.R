# TODO: Add comment
# 
# Author: wjiang2
###############################################################################



setMethod("plotGate",signature(x="workFlow",y="missing"),function(x,y,...){
						
#		for(curSample in sampleNames(Data(x[["base view"]])))
#		{
#			plotGate(x,curSample,...)
#		}
		plotGate(x,sampleNames(Data(x[["base view"]])),...)
	})
				
setMethod("plotGate",signature(x="workFlow",y="numeric"),function(x,y,...){
			plotGate(x,y=sampleNames(Data(x[["base view"]]))[y],...)
		})
setMethod("plotGate",signature(x="workFlow",y="character"),function(x,y,z,main=NULL,arrange=TRUE,smooth=FALSE,lattice=FALSE,...){

  actNames <- actions(x)

	# Excludes the compensation, transformation, and boundary filter actions.
  ind <- !grepl("compensation", actNames, ignore.case = TRUE)
  ind <- ind & !grepl("transformation", actNames, ignore.case = TRUE)
  ind <- ind & !grepl("boundary", actNames, ignore.case = TRUE)

  # Also, excludes intersection filter actions. The intersection filter uses the
  # word "and" in the action. So, we exclude actions that include "and".
  ind <- ind & !grepl("and", actNames, ignore.case = TRUE)

  actNames <- actNames[ind]
  
  ##filter by the provided index z
  if (!missing(z)) {
    if (is.numeric(z)) {
      actNames<-actNames[z]
    } else {
      actNames <- z
    }
  }
			
			if(lattice)
			{
			
				plotObjs<-sapply(actNames,function(curActName){
							
							
							curAct<-x[[curActName]]
							pView<-parent(curAct)
							curData<-Data(pView)[y]
							curFres<-gate(curAct)
							smooth<-ifelse(nrow(curData[[1]])<100,TRUE,smooth)
							if(is.null(main))
								main<-names(pView)
#						browser()
							plotGate(x=curData
									,y=curFres[y]
									,smooth=smooth
									,main.panel = main
									,...
							)
							
						},simplify=F)
				plotObjs
				
			}else
			{
			
				lapply(y,function(curSample){
							if(is.null(main))
								main<-paste(rev(pData(Data(x[["base view"]]))[curSample,][-2]),collapse=":")
							plotObjs<-sapply(actNames,function(curActName){
										
										
										curAct<-x[[curActName]]
										pView<-parent(curAct)
										curData<-Data(pView)
										curFres<-gate(curAct)
										smooth<-ifelse(nrow(curData[[curSample]])<100,TRUE,smooth)
#						browser()
										plotGate(x=curData[curSample]
												,y=curFres[curSample]
												,smooth=smooth
												,main.panel = names(pView)
												,...
										)
										
									},simplify=F)
							if(arrange)
								do.call(grid.arrange
										,c(plotObjs[actNames],main=main)
										)		
							else
								invisible(plotObjs[actNames])	
						})
					
			}
			
		})

setMethod("plotGate",signature(x="workFlow",y="list"),function(x,y,...){
			
		})

setMethod("plotGate",signature(x="flowSet",y="filter"),function(x,y,...){
#		
			.plotGate(obj=x,filterlist=y,params=parameters(y),...)
		})


setMethod("plotGate",signature(x="flowSet",y="list"),function(x,y,...){
#			browser()
			
			.plotGate(obj=x,filterlist=y,params=parameters(y[[1]]),...)
		})

setMethod("plotGate",signature(x="flowSet",y="filtersList"),function(x,y,...){
#			browser()
			
			.plotGate(obj=x,filterlist=y,params=parameters(y[[1]][[1]]),...)
		})

.plotGate<-function(obj,filterlist,params
		,x.axis=NULL,y.axis=NULL
		,groupBy=NULL
		,main.panel=NULL
		,cond="name"
		,stat=TRUE
		,smooth=FALSE
		,margin=FALSE
		,...){
#	browser()
	if(is.null(x.axis)&&is.null(y.axis))#get channel info from filterlist
	{
		
#		if(is(filterlist,"list"))
#			if(is(filterlist,"filtersList"))
#				params<-parameters(filterlist[[1]][[1]])
#			else
#				params<-parameters(filterlist[[1]])
#		else
#			params<-parameters(filterlist)
		if(length(params)==1)
		{
			x=params
			y=NULL
		}else
		{
			x=params[1]
			y=params[2]
		}
	}else
	{
		x<-x.axis
		y<-y.axis
	}
	#guess the x channel name
	xChnl<-getChannelMarker(obj[[1]],x)
	
	#guess the y channel name
	if(is.null(y))
		y<-"SSC-A"
	yChnl<-getChannelMarker(obj[[1]],y)
	
	
#	browser()
	formula1<-paste("`",yChnl$name,"`~`",xChnl$name,"`",sep="")
	if(is.null(groupBy))
	{

		if(length(obj)==1)
		{
			obj<-obj[[1]]
			rg<-range(obj)
			if(is.list(filterlist))
				filterlist<-filterlist[[1]]
		}else
		{
			formula1<-paste(formula1,cond,sep="|")
			rg<-range(obj[[1]])
		}
#		browser()
	
		
		d1<-xyplot(as.formula(formula1)
					,data=obj
					,filter=filterlist
					,main=main.panel
					,xlab=sub("NA","",paste(unlist(xChnl),collapse=" "))
					,ylab=sub("NA","",paste(unlist(yChnl),collapse=" "))
					,stat=stat
					,margin=margin
					,smooth=smooth
					,ylim=c(min(-1,rg[yChnl$name][1,]),max(5,rg[yChnl$name][2,]))
					
					,...
					)
		d1
	}else
	{
#		browser()
		f1<-eval(substitute(pData(obj)$f,list(f=groupBy)))
#	browser()
		if(class(obj)=="ncdfFlowSet")
		{
			fslist<-split(obj,f1)@datalist
		}else
		{
			fslist<-flowCore:::split(obj,f1)
		}
		
		plotObjs<-lapply(names(fslist),function(pid)
		{
			obj<-fslist[[pid]]
			if(!is.null(x))
			{
			
				if(length(obj)==1)
				{
					obj<-obj[[1]]
					rg<-range(obj)
					if(is.list(filterlist))
						filterlist<-filterlist[[1]]
				}else
				{
					formula1<-paste(formula1,cond,sep="|")
					rg<-range(obj[[1]])
				}
			
				
				d1<-xyplot(as.formula(formula1)
							,data=obj
							,filter=filterlist
							,main=paste(pid,main.panel)
							,xlab=sub("NA","",paste(unlist(xChnl),collapse=" "))
							,ylab=sub("NA","",paste(unlist(yChnl),collapse=" "))
							,stat=stat
							,margin=margin
							,smooth=smooth
							,ylim=c(min(-1,rg[yChnl$name][1,]),max(5,rg[yChnl$name][2,]))
							
							,...
							)
			}else
			{
				d1<-densityplot(`antigen`~`<PE Tx RD-A>`,data=obj[[pid]],filter=filterlist[[pid]]
						,abline=TRUE
						,main=main
				)
			}
			d1
		})
		plotObjs
	}
#	
}
##############################################################################
#copied from flowViz::panel.xyplot.flowset,the change is:
#1.call panel.flowClust.xyplot.flowframe instead of default one panel.xyplot.flowframe
#2.parse include list 
#3.use named argument x=x,y=y instead simply paass x,y
##############################################################################
panel.flowClust.xyplot.flowset <- function(x,
		frames,
		filter=NULL,
		channel.x,
		channel.y
		,xbins=0 #passed to hexbin routine
		,binTrans=sqrt	  
		,include
		,...)
{
	nm <- as.character(x)
	if (length(nm) < 1) return()
	## 'filter' either has to be a single filter, or a list of filters matching
	## the flowSet's sample names, or a filterResultList.
	if(!is.null(filter)){
		if(!is.list(filter)){
			if(is(filter, "filter")){
				filter <- lapply(seq_along(nm), function(x) filter)
				names(filter) <- nm
			}
		}else if(!is(filter, "filterResultList"))
			filter <- as(filter, "filterResultList")
		if(!nm %in% names(filter) || !is(filter[[nm]] ,"filter")){
			warning("'filter' must either be a filterResultList, a single\n",
					"filter object or a named list of filter objects.",
					call.=FALSE)
			filter <- NULL
		}
	}
	x <- flowViz:::evalInFlowFrame(channel.x, frames[[nm]])
	y <- flowViz:::evalInFlowFrame(channel.y, frames[[nm]])
#	browser()
	panel.flowClust.xyplot.flowframe(x=x, y=y, frame=frames[[nm]], filter=filter[[nm]],include=include[[nm]],xbins=xbins,binTrans=binTrans, ...)
}

##############################################################################
#TODO:match the parameter names to make sure the order is right   
#panel function to plot elipse for tmixFilter returned by flowClust   ##
##############################################################################
panel.flowClust.xyplot.flowframe <- function(ellipse=T
		,filter=NULL
		, include=1:(filter@K)
		, ecol=1, elty=1
		,pch=gpar$flow.symbol$pch
		,alpha=gpar$flow.symbol$alpha
		,cex=gpar$flow.symbol$cex
		,col=gpar$flow.symbol$col
		,gp=NULL
		,level=NULL
		, npoints=501
		,subset=c(1,2)
		, ...)
{
#	browser()	
	flowViz:::panel.xyplot.flowframe(gp=gp,...)
	
	gpar <- flowViz.par.get()
	
	if(!is.null(gp))
		gpar <- lattice:::updateList(gpar, gp)
	if(is.null(gpar$gate$cex))
		gpar$gate$cex <- cex
	if(is.null(gpar$gate$pch))
		gpar$gate$pch <- pch
	if(is.null(gpar$gate$plotType))
		gpar$gate$plotType<-"l"
	if(is.null(gpar$density))
		gpar$density<-TRUE
	if(!is.null(level))
		filter@ruleOutliers[2]<-level
	py<-2
	j <- 0
	# plot ellipses
	if (ellipse) {
		ecol <- matrix(ecol, length(include))
		elty <- matrix(elty, length(include))
		
		if (all(filter@nu!=Inf)) {
			if (filter@ruleOutliers[1]==0) {     # 0 means quantile
				cc <- py * qf(filter@ruleOutliers[2], py, filter@nu)
			}  else {     # 1 means u.cutoff
				cc <- ((filter@nu+py)/filter@ruleOutliers[2] - filter@nu)    
			}
		}  else cc <- qchisq(filter@ruleOutliers[2], py)
		
		j <- 0
		if (length(filter@lambda)>0)
			lambda <-rep(filter@lambda, length.out=filter@K)
		else
			lambda <-numeric(0)
		cc <- rep(cc, length.out=filter@K)
		for (i in include) {
			eigenPair <- eigen(filter@sigma[i,subset,subset])
			l1 <- sqrt(eigenPair$values[1]) * sqrt(cc)
			l2 <- sqrt(eigenPair$values[2]) * sqrt(cc)
			angle <- atan(eigenPair$vectors[2,1] / eigenPair$vectors[1,1]) * 180/pi
			
			if (length(lambda)>0) {
#				browser()
				lpoints(rbox(flowClust:::.ellipsePoints(a=l1[i], b=l2[i]
										, alpha=angle
										, loc=filter@mu[i,subset]
										, n=npoints)
								, lambda[i])
						
						,type=gpar$gate$plotType
						,lty=gpar$gate$lty
						,col=gpar$gate$col
						,alpha=gpar$gate$alpha
						,lwd=gpar$gate$lwd
				)
			} else {
				lpoints(flowClust:::.ellipsePoints(a=l1[i], b=l2[i]
								, alpha=angle
								, loc=filter@mu[i,subset]
								, n=npoints)
						,type=gpar$gate$plotType
						,lty=gpar$gate$lty
						,col=gpar$gate$col
						,alpha=gpar$gate$alpha
						,lwd=gpar$gate$lwd
				)
			}
		}  
	}
	
}


###############################################################################
#plot manual gates
###############################################################################
#plotGate_gh_tcell<-function(gh){
#	
#	
#	nodePath<-getNodes(gh,isPath=T)
#	plotObjs<-NULL
#	nLength<-length(nodePath)
##	nLength<-3w
#	for(curNodeInd in 2:nLength)
#	{
#		
#		curPath<-nodePath[curNodeInd]
#		curNodeName<-basename(curPath)
#		
#		
#		pos<-c(1.2,1.2)
#		abs=F
#		
#		
#		
#		if(curNodeName=="Viable")
#		{
#			pos<-0.5
#		}
#		if(curNodeName=="Singlets")
#		{
#			pos<-c(0.3,0.2)
#			abs=T
#		}
#		if(curNodeName=="Lymph")
#		{
#			pos<-c(0.3,0)
#			abs=T
#		}
##		if(curNodeName=="CD3+")
##		{
##			pos<-c(1.2,1.2)
##		}
##		if(curNodeName=="CD4+")
##		{
##			pos<-c(1.2,1.2)
##		}
##		if(curNodeName=="activated CD4")
##		{
##			pos<-c(0.1,0.9)
##		}
#		if(curNodeName=="central memory")
#		{
#			pos<-c(0.9,0.1)
#		}
#		if(curNodeName=="naive")
#		{
#			pos<-c(0.9,0.9)
#		}
#		if(curNodeName=="effector")
#		{
#			pos<-c(0.2,0.9)
#		}
#		
#		plotObjs[[curNodeInd]]<-flowWorkspace:::plotGate(gh,curNodeInd
#				,smooth=FALSE
#				,margin=FALSE
#				#				,colramp=colorRampPalette(rev(brewer.pal(11, "Spectral")))
##								,xbin=64
#				,pos=pos
#				,abs=abs
#				,stat=T
#		)
#		
##			browser()
#		
#	}
#	plotObjs
#	
#}
plotGate_G<-function(G){
	
	plotObjs<-NULL
	
	pd<-pData(G)
	nodePath<-getNodes(G[[1]],isPath=T)
	nLength<-length(nodePath)
#	nLength<-3
	for(curNodeInd in 2:nLength)
	{
		
		pNodeInd<-getParent(G[[1]],curNodeInd)
		parentData<-getData(G,pNodeInd)
		pData(parentData)<-pd
		filters<-getGate(G,curNodeInd)
		curFres<-filterList(filters)
		
		pos<-c(1.2,1.2)
		abs=F
		curPath<-nodePath[curNodeInd]
		curNodeName<-basename(curPath)
		
		if(curNodeName=="Viable")
		{
			pos<-0.5
		}
		if(curNodeName=="Singlets")
		{
			pos<-c(0.3,0.2)
			abs=T
		}
		if(curNodeName=="Lymph")
		{
			pos<-c(0.3,0)
			abs=T
		}
#		if(curNodeName=="CD3+")
#		{
#			pos<-c(1.2,1.2)
#		}
#		if(curNodeName=="CD4+")
#		{
#			pos<-c(1.2,1.2)
#		}
#		if(curNodeName=="activated CD4")
#		{
#			pos<-c(0.1,0.9)
#		}
		if(curNodeName=="central memory")
		{
			pos<-c(0.9,0.1)
		}
		if(curNodeName=="naive")
		{
			pos<-c(0.9,0.9)
		}
		if(curNodeName=="effector")
		{
			pos<-c(0.2,0.9)
		}
#		
#		
#				
#		browser()
		plotObjs[[curNodeInd]]<-plotGate(obj=parentData
				,filterlist=curFres
				,smooth=FALSE
				#						,par.settings=list(par.main.text=list(cex=0.8))
				#						,margin=F
				#				,colramp=colorRampPalette(rev(brewer.pal(11, "Spectral")))
				,dataRange=F
				,pos=pos
				,abs=abs
				,cond="factor(Center):factor(Donors):factor(Replicate)"
				,main=curPath
		)
		
		
		
	}
	plotObjs
	
}

pairPlot<-function(x,...){
	
	chnls<-colnames(x)
	chnlPairs<-combn(chnls,2)
#		browser()
	apply(chnlPairs,2,function(curPair){
#				browser()
				f1<-paste("`",curPair[1],"`~`",curPair[2],"`",sep="")
				xyplot(as.formula(f1),x,...)
				
			})
}
