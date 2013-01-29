##It is a hack,using parent and child names to determine the unique ID for child
##might not work if the path of parent/child itself is not unqiue 
#getNodeID<-function(gh,parent,nodeName,...){
#	browser()
#	if()
#	if(parent=="root")
#		tPath<-paste("/",nodeName)
#	else
#		tPath<-paste("/",nodeName)
#	getNodes(wf[[1]],isPath=T)
#}
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



#################################
##stats extraction 
#################################
extractStatsFromG<-function(G,cytotrol=TRUE,pops=NULL){
	
	pd<-pData(G)
	pd$sample<-getSamples(G)#sometime samplename in G is differnt from pd
	res<-lapply(getSamples(G),function(curSample)
								{
									#select the terminal nodes only								
									
									stats<-getPopStats(G[[curSample]])[,4,drop=F]
									stats$node<-rownames(stats)
#									browser()
									if(!is.null(pops))
									{
										if(is.integer(pops))
											stats<-stats[pops,]
										else if(is.character(pops))
											stats<-stats[grepl(pops,stats$node),]
									}
											
									
									
									stats$sample<-pd$name[match(curSample,pd$sample)]
									rownames(stats)<-NULL

									##rename the pop names
									popNames<-basename(as.character(stats$node))
									parentNames<-basename(dirname(as.character(stats$node)))
#									browser()									
									popNames<-gsub("[a,A]ctivated[ ]*[CD]*[4,8]*","Activated",popNames)
									popNames<-paste(popNames,parentNames)
									stats$"Cell Population"<-popNames
									names(stats)[1]<-"percent"
									stats[,-2]	
								})
	res<-do.call("rbind",res)
	
	res<-merge(res,pd,by.x="sample",by.y="name")
	if(cytotrol)
		selectedCol<-c("percent","Replicate","Panel","Center","Cell Population")
	else
		selectedCol<-c("percent","Donors","Replicate","Panel","Center","Cell Population")
	res<-res[,selectedCol]
	res
}

extractStats<-function(wf,popList){
	
	###collect % from gate action view
	res<-NULL
		browser()
	for(gActName in rev(actions(wf))[1:2])
	{
		
		actView<-wf[[gActName]]
		localRes<-summary(actView)[,c(1:3)]
		localRes$parent<-names(parent(actView))
		res<-rbind(res,localRes)
	}
	
#		head(res)
#		nrow(res)
	
	##remove unwanted pops
	res<-subset(res,population%in%popList)
	
	##rename the pop names
	popNames<-as.character(res$population)
#	browser()	
#	popNames<-gsub("CD38\\+HLA\\-DR\\+","Activated",popNames)
	popNames<-gsub("activated CD[4,8]\\+","Activated",popNames)
	popNames<-gsub("CCR7\\+CD45RA\\+","Naive",popNames)
	popNames<-gsub("CCR7\\+CD45RA\\-","Central Memory",popNames)
	popNames<-gsub("CCR7\\-CD45RA\\-","Effector Memory",popNames)
	popNames<-gsub("CCR7\\-CD45RA\\+","Effector",popNames)
	res$population<-popNames
	
	##rename the parent
	pNames<-as.character(res$parent)
	pNames<-gsub("CD4\\+CD8\\-","CD4+",pNames)
	pNames<-gsub("CD4\\-CD8\\+","CD8+",pNames)
	res$parent<-pNames
	
	
	res$"Cell Population"<-paste(res$population,res$parent)
	pd<-pData(Data(wf[["base view"]]))
#		browser()
	res<-merge(res,pd,by.x="sample",by.y="name")
	
	selectedCol<-c("percent","Donors","Replicate","Panel","Center","Cell Population")
	res<-res[,selectedCol]
	
}#convert tmixFilterResult to polygon gate

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
		ind<-which(grepl(pname,pd$name))
		
	
	
		
	if(length(ind)==0)
	{
		#try marker name
		ind<-which(unlist(lapply(pd$des,function(x){
								#split by white space and then match each individual string
#									browser()
									if(fix)
										any(unlist(lapply(strsplit(x," "),function(y)toupper(y)%in%toupper(name))))
									else
										grepl(pattern=pname,x)
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







#######plot legend for density scale
ff<-function(){
#	cR<-cols
	cols<-c("blue","green","yellow","red")
	cR<-IDPcolorRamp(21,
			t(col2hsv(cols)),
			fr=c(0.7,0))
	colramps <- colorRampPalette(cR,bias=3)
	showColors(cR)
	
	fs1<-Data(vLymph)
#	mat<-exprs(Data(vSinglet)[[1]])[,c("FSC-A","SSC-A")]
	
#	cols <- densCols(mat, colramp=colramps)
	cols<-sort(unique(cols))
	
	ind<-as.integer(unlist(lapply(strsplit(levels(cut(1:length(cols),10)),","),function(x)substr(x[[2]],1,nchar(x[[2]])-1))))
	
	unloadNamespace("flowStats")
	unloadNamespace("QUALIFIER")
	unloadNamespace("flowWorkspace")
	unloadNamespace("ncdfFlow")
	unloadNamespace("flowViz")
	library(flowViz)
	
	time1<-Sys.time()
	png("xyplot_hexbin.png")
	Rprof()
	xyplot(`<FITC-A>`~`<Pacific Blue-A>`,fs1[1]
#						,colramp=colramps
						,smooth=F
						,xbins=512
#					,key=list(cex.title=1
#							,space="right"
#							,title="density scale"
#							,rect=list(col=c(cols[ind][-1],"red"),border=F)
#							,text=list(as.character(c(1:10)))
#							)


				)
	Rprof(NULL)
	summaryRprof()
	
				
	dev.off()
	Sys.time()-time1

	time1<-Sys.time()
	png("xyplot.png")
	xyplot(`SSC-A`~`FSC-A`,fs1,colramp=colramps,smooth=F
#						,nbin = 1000
#			,key=list(cex.title=1
#					,space="right"
#					,title="density scale"
#					,rect=list(col=c(cols[ind][-1],"red"),border=F)
#					,text=list(as.character(c(1:10)))
#					)
	)
	dev.off()
	Sys.time()-time1
	
	
}

####################################
##7.cytokine gating by rangeGate
#####################################

ckGating<-function(ckPop.list,ckChnls
					,groupBy="patient",isNcdf=FALSE
					,sdvalue=4,collapsed=TRUE,isRare=FALSE)
{
	
	time1<-Sys.time()
	cl <- makeCluster(10,type = "SOCK")

#	clusterExport(cl,c("ckChnls","dataCollapse","sdvalue","isRare","groupBy"))


	CkFresList<-lapply(ckPop.list,function(curPop){
		#####################################################################
		#apply rangeGate on each cytokine channels for each patient group
		#####################################################################
#		browser()
		
#		lapply(curPop,function(nc,ckChnls,collapsed,sdvalue,isRare,groupBy){
		parLapply(cl,curPop,function(nc,ckChnls,collapsed,sdvalue,isRare,groupBy){
#					library(ncdfFlow)
					library(flowStats)
					curfresList<-list()
					pd<-pData(parameters(nc[[1]]))
					for(curChannel in ckChnls)
					{
						
						curCk<-as.character(subset(pd,name==curChannel)$desc)
						
						if(groupBy=="sample")
						{
							curfresList[[curCk]]<-fsApply(nc,function(frm){
															g<-rangeGate(frm,stain=curChannel
																			,absolute=F,sd=sdvalue
																			,ref=0,rare=isRare)
															filter(frm,g)
														})
						}
						
						if(groupBy=="patient")
						{
							if(collapsed)
							{###
								g<-rangeGate(nc,stain=curChannel,absolute=F,sd=sdvalue,ref=0,rare=isRare)
							}else
							{	#select postive control to get the gate
								sebctrl<-as.character(subset(pData(nc),antigen=="sebctrl")$name)
								g<-rangeGate(nc[[sebctrl]],stain=curChannel,absolute=F,sd=sdvalue,ref=0,rare=TRUE)
							}
							curfresList[[curCk]]<-filter(nc,g)
						}
						
					}
					curfresList
					
				},ckChnls,collapsed,sdvalue,isRare,groupBy)
		
		})
		
	
	stopCluster(cl)
	print(Sys.time()-time1)
	
	CkFresList
}	

ckGatePlot<-function(cd4cd8Pos.list,CkFresList,plotType="xyplot"
						,groupBy="patient",sdvalue=4,collapsed=TRUE,ckChnls=ckChnls,isRare=FALSE){
	
	#####################################################################
	#plot results of gated population
	#####################################################################
	pdf(paste("RangeGateCk_sd",sdvalue,groupBy
					,ifelse(groupBy=="sample",""
								,ifelse(collapsed
										,"collapse","sebctrl")
							)
					,ifelse(isRare,"Rare","")
					,".pdf",sep="")
#			,width=10
		)
	for(parentPopName in names(CkFresList))
	{
		nc<-cd4cd8Pos.list[[parentPopName]]
		for(pid in names(nc))
		{
			###update meta info
#				curMeta<-pData(nc[[pid]])
#				curMeta$antigen<-as.character(curMeta$antigen)
#				curMeta$antigen[which(is.na(curMeta$antigen))]<-c("unknown1","unknown2")
#				curMeta$antigen[which(curMeta$antigen=="negctrl")]<-c("negctrl1","negctrl2")
#				curMeta$antigen<-as.factor(curMeta$antigen)
#				pData(nc[[pid]])<-curMeta
#				
			pd<-pData(parameters(nc[[pid]][[1]]))
			for(i in 1:length(ckChnls))
			{
				curChl<-ckChnls[i]
				curCk<-as.character(subset(pd,name==curChl)$desc)
				curCkFres<-CkFresList[[parentPopName]][[pid]][[curCk]]
#				browser()
				if(plotType=="xyplot")
				{
					
#					d1<-xyplot(x=as.formula(paste("`",curChl,"`~`SSC-A`|antigen",sep="")),nc[[pid]]
#							,filter=curCkFres
#							,main=paste(parentPopName,pid,curCk,"vs SSC-A")
#							,abline=TRUE
#							)
#					plot(d1)
					
#					d1<-xyplot(x=as.formula(paste("`",curChl,"`~`FSC-A`|antigen",sep="")),nc[[pid]]
#							,filter=curCkFres
#							,main=paste(parentPopName,pid,curCk,"vs FSC-A")
#							,abline=TRUE
#					)
#					plot(d1)
#					
#					d1<-xyplot(x=as.formula(paste("`",curChl,"`~`<Pacific Blue-A>`|antigen",sep="")),nc[[pid]]
#							,filter=curCkFres
#							,main=paste(parentPopName,pid,curCk,"vs live")
#							,abline=TRUE
#					)
#					plot(d1)
					
					if(curCk=="IL2")
					{
						d1<-xyplot(x=as.formula(paste("`",curChl,"`~`<PE Cy7-A>`|antigen",sep="")),nc[[pid]]
							,filter=curCkFres
							,main=paste(parentPopName,pid,curCk,"vs IFNg")
							,abline=TRUE
							)
						plot(d1)
					}
					if(curCk=="TNFa")
					{
						d1<-xyplot(x=as.formula(paste("`",curChl,"`~`<APC-A>`|antigen",sep="")),nc[[pid]]
							,filter=curCkFres
							,main=paste(parentPopName,pid,curCk,"vs IL4")
							,abline=TRUE
						)
						plot(d1)
					}
					
				}else
				{
					d1<-densityplot(x=as.formula(paste("`antigen`~`",curChl,"`",sep="")),nc[[pid]]
							,filter=curCkFres
							,main=paste(ifelse(i==1,parentPopName,""),
									ifelse(i==1,pid,""),
									curCk)
							,abline=TRUE
					)
					plot(d1)
				}
				

				
			}
			
			
		}
	}
	
	dev.off()
	
	
}



plot2dCk<-function(parentPopName,nc)
{
	for(pid in names(nc))
	{
		
			d1<-xyplot(x=`<PE Green laser-A>`~`<PE Cy7-A>`|antigen,nc[[pid]]
					,xlab="IFNg",ylab="IL2"
	#				,colramp=cols
	#				,smooth=FALSE
					,main=paste(parentPopName,pid,"IL2 vs IFNg")
	#				,panel=panel.xyplot.flowsetEx
					
					)
			plot(d1)
		
			d1<-xyplot(x=`<Alexa 680-A>`~`<APC-A>`|antigen,nc[[pid]]
					,xlab="TNFa",ylab="IL4"
					,main=paste(parentPopName,pid,"TNFa vs IL4")
					
			)
			plot(d1)
		
		
	}
}






featureExtraction<-function(CkFresList,ckPop.list,meta,groupBy
								,sdvalue=4,collapsed=TRUE
								,ckLogic,isRare=FALSE,negSelFun
								,MFI=FALSE
								,isTraining=TRUE)
{
	time1<-Sys.time()
	
	library(gtools)
#	browser()
	
	ckfeatures<-NULL
	##extract proportion
	for(pid in names(CkFresList[[1]]))
	{
		curPidfeatures<-NULL
		for(parentPopName in names(CkFresList))
		{
			
			
			curCkFresList<-CkFresList[[parentPopName]]
			
			parentCount<-as.matrix(sapply(curCkFresList[[pid]][[1]],function(curCkFres){
								ind<-curCkFres@subSet
								length(ind)
							}))
			
			###add  proportions of 4 cks and their polyfunctional combinations
			curfeatures<-do.call("rbind",lapply(1:2,function(k)
			{
				pfs<-combinations(n=4,r=k,v=names(curCkFresList[[pid]]))
					
				do.call("rbind",apply(pfs,1,function(curPf){
								
							curPfvalues<-sapply(names(curCkFresList[[pid]][[1]]),function(curFcs,curPf){
					
											ind<-rep(TRUE,length(curCkFresList[[pid]][[1]][[curFcs]]@subSet))
											for(n in 1:length(curPf))
											{
												curpfCk<-curPf[n]
												curInd<-curCkFresList[[pid]][[curpfCk]][[curFcs]]@subSet
												if(n==1)
												{
													ind<-ind&curInd
													if(MFI)
													{
#														browser()
														curFrm<-ckPop.list[[parentPopName]][[pid]][[curFcs]]						
														curChnl<-as.character(subset(pData(parameters(curFrm)),desc==curpfCk)$name)
														curMFI<-mean(exprs(curFrm)[ind,curChnl])
														#normailze it by mu and s
														est<-huber(exprs(curFrm)[,curChnl])
														curMFI<-(curMFI-est$mu)/est$s
													}
												}else
												{
													if(ckLogic=="and")
													{
														ind<-ind&curInd
													}else
													{
														ind<-ind|curInd	
													}		
													
												}
											}
											
											c(curPfCount=length(which(ind)),curMFI=curMFI)
										}
									,curPf)
																
							curPfvalues<-t(curPfvalues)

									
							ret<-data.frame(pos=curPfvalues[,"curPfCount"]
											,neg=parentCount-curPfvalues[,"curPfCount"]
											,MFI=curPfvalues[,"curMFI"]
											,fname=paste(curPf,collapse="."))
							ret$fcsfile<-rownames(ret)
							ret
							})
					)
				
			}))
			
			curfeatures$parent<-parentPopName
			
			curPidfeatures<-rbind(curPidfeatures,curfeatures)
		}
		ckfeatures<-rbind(ckfeatures,curPidfeatures)
		
	}
#	browser()
	
	features1<-merge(ckfeatures,meta,by.x="fcsfile",by.y="fcsfile")
#	browser()
	write.csv(features1[,c(colnames(ckfeatures),"antigen","PTID","pub_id1")]
				,file=paste("counts_sd",sdvalue,negSelFun,"_",ckLogic,".csv",sep="")
				,row.names = F)
	print(Sys.time()-time1)
		

	
	
}






featurecd4cd8<-function(meta,negSelFun="min",MFI=FALSE,isTraining=TRUE)
{
	library(gtools)
#	browser()

	ckfeatures<-do.call("rbind",lapply(names(cd4cd8.norm.list),function(curPid){
					
					parentCount<-fsApply(cd4cd8.norm.list[[curPid]],nrow)
					curfeatures<-do.call("rbind",sapply(names(cd4cd8Pos.list),function(parentPopName){
											curPfvalues<-fsApply(cd4cd8Pos.list[[parentPopName]][[curPid]],nrow)
											ret<-data.frame(pos=curPfvalues
													,neg=parentCount-curPfvalues
													,fname=parentPopName)
											ret$fcsfile<-rownames(ret)
											ret
										},simplify=F))
				
				
					curfeatures
					})
			)
	
#	browser()
	
	features1<-merge(ckfeatures,meta,by.x="fcsfile",by.y="fcsfile")
#	browser()
	write.csv(features1[,c(colnames(ckfeatures),"antigen","pub_id")]
			,file=paste("counts_cd4cd8.csv",sep="")
			,row.names = F)
	
	
	training<-ckfeatures
	training$proportion<-training$pos/(training$pos+training$neg)
#	training$proportion<-training$proportion#*training$MFI
	training$fname<-paste(training$parent,training$fname,sep=".")
#	training<-training[,colnames(cdFeature)]
#	training<-rbind(training,cdFeature)
	training<-merge(training,meta,by.x="fcsfile",by.y="fcsfile")
	training<-training[,c("proportion","fname","antigen","pub_id")]
	
	
	
	#construct training table
	
	#get negative control
	T_neg<-subset(training,antigen%in%c("negctrl1","negctrl2"))
	T_neg<-T_neg[order(T_neg$fname),]
	T_neglist<-split(T_neg,T_neg$fname)
	T_neglist_sel<-lapply(T_neglist,function(curT_neg){
				curT_neg$pub_id<-as.factor(as.character(curT_neg$pub_id))
				curT_neg<-split(curT_neg,curT_neg$pub_id)
#		browser()
				curT_neg<-lapply(curT_neg,function(x)eval(parse(text=negSelFun))(x$proportion))
				curT_neg<-do.call("rbind",curT_neg)
				data.frame(proportion=curT_neg
				)
			})
	T_neg<-do.call("cbind",T_neglist_sel)
	colnames(T_neg)<-names(T_neglist_sel)		
	
	
	# subtract negative control value from training data
	#ENV
	if(isTraining)
	{
		T_env<-subset(training,antigen=="ENV-1-PTEG")
	}else
	{
		T_env<-subset(training,antigen=="unknown1")
	}
	
	T_env<-T_env[order(T_env$fname),]
	T_envlist<-split(T_env,T_env$fname)
	T_envlist<-lapply(T_envlist,function(x){
				x$pub_id<-as.character(x$pub_id)
				x<-x[order(x$pub_id),]
				rownames(x)<-x$pub_id
				x[,"proportion",drop=FALSE];
			})
	T_env<-do.call("cbind",T_envlist)
	colnames(T_env)<-names(T_envlist)		
	T_env<-T_env-T_neg
	T_env$antigen<-"ENV-1-PTEG"
	T_env$pub_id<-rownames(T_env)
	
	#GAG
	if(isTraining)
	{
		T_gag<-subset(training,antigen=="GAG-1-PTEG")
	}else
	{
		T_gag<-subset(training,antigen=="unknown2")
	}
	
	T_gag<-T_gag[order(T_gag$fname),]
	T_gaglist<-split(T_gag,T_gag$fname)
	T_gaglist<-lapply(T_gaglist,function(x){
				x$pub_id<-as.character(x$pub_id)
				x<-x[order(x$pub_id),]
				rownames(x)<-x$pub_id
				x[,"proportion",drop=FALSE];
			})
	T_gag<-do.call("cbind",T_gaglist)
	colnames(T_gag)<-names(T_gaglist)		
	T_gag<-T_gag-T_neg
	T_gag$antigen<-"GAG-1-PTEG"
	T_gag$pub_id<-rownames(T_gag)
	
	##put two subset together and run classifier
	tt<-rbind(T_env,T_gag)
	tt$antigen<-factor(tt$antigen)
	
	##classification
	
	library(RWeka)
	library(RColorBrewer)
	mypalette<-brewer.pal(9,"Greens")
	tt<-tt[order(tt$antigen),]
	
	data<-tt[,c(-ncol(tt),-ncol(tt)+1)]
	rownames(data)<-paste(tt$pub_id,substr(tt$antigen,1,1))
	data[data<0]<-0
	data[is.na(data)]<-0
#	pdf(paste("heatmap_sd",sdvalue,groupBy
#					,ifelse(groupBy=="sample",""
#							,ifelse(collapsed
#									,"collapse","sebctrl")
#					)
#					,ifelse(isRare,"Rare","")
#					,negSelFun
#					,".pdf",sep=""))
	
#	colnamesel<-colnames(data)[grep("cd4cd8Pos",colnames(data))]
#	data<-data[,colnamesel]
	data1<-log(as.matrix(data.frame(scale(data,center=FALSE))))
	data1[is.infinite(data1)]<-NA
	heatmap(x=data1,Colv=NA,Rowv=NA,col=mypalette,scale="none")
#	dev.off()
#	
	
	data$antigen<-tt$antigen
	write.csv(data,paste("trainging_",negSelFun,".csv",sep=""))
#	
	##	browser()
#	
	m1 <- J48(antigen ~ ., data = data)
	e1 <- evaluate_Weka_classifier(m1,
			cost = matrix(c(0,2,1,0), ncol = 2),
			numFolds = 10, complexity = TRUE,
			seed = 123, class = TRUE)
#	
	e1
#	
#	print(e1)
	
#	write.csv(cbind(data,tt[c("fcsfile","pub_id","response.cd4.","response.cd8.")])
#				,file="trainingRatio.csv")
	
	
	
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
