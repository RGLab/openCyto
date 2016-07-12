#' lightweight version of flowClust::flowClust
#' Speed up by removing some layers to R methods
#' and validity checks as well as changing some default parameters
.flowClustFast<-function(x,varNames=NULL, K, B=500, tol=1e-5, nu=4, lambda=1
                     , nu.est=0, trans=1, min.count=10, max.count=10
                     , min=NULL, max=NULL, level=0.9, u.cutoff=NULL
                     , z.cutoff=0, randomStart=0, B.init=B, tol.init=1e-2
                     , seed=1, criterion="BIC", control=NULL, prior=NULL,usePrior="no"
                     , nstart = 1
                     , parallel = TRUE)
{
#	if (is(x, "flowFrame")) {
#		if (length(varNames)==0) {
#			y <- exprs(x)
#			varNames <- colnames(y)
#		}
#		else {
			y <- as.matrix(exprs(x)[, varNames,drop=FALSE])
#		}
#	}
#	else if (is(x, "matrix")) {
#		if (length(varNames)==0) {
#			y <- x
#			if (length(colnames(x))==0)
#				varNames <- "Not Available"
#			else varNames <- colnames(x)
#		}
#		else {
#			y <- as.matrix(x[, varNames,drop=FALSE])
#		}
#	}
#	else if (is(x, "data.frame")) {
#		if (length(varNames)==0) {
#			y <- as.matrix(x)
#			varNames <- colnames(x)
#		}
#		else {
#			y <- as.matrix(x[, varNames,drop=FALSE])
#		}
#	}
#	else if (is(x, "vector")) {
#		y <- matrix(x)
#		if (length(varNames)==0)
#			varNames <- "Not Available"
#	}
#	else {
#		stop(paste("Object ", as.character(x), " is not of class flowFrame / matrix / data frame!"))
#	}
#
#	# finding filtered observations
#	rm.max <- rm.min <- rep(FALSE, nrow(y))
#	if (max.count > -1) {
#		if (is.null(max)[1])
#			max <- apply(y, 2, max)
#		for (k in 1:ncol(y))  if (sum(y[,k]>=max[k]) >= max.count)
#				rm.max <- rm.max | (y[,k] >= max[k])
#	}
#	if (min.count > -1) {
#		if (is.null(min)[1])
#			min <- apply(y, 2, min)
#		for (k in 1:ncol(y))  if (sum(y[,k]<=min[k]) >= min.count)
#				rm.min <- rm.min | (y[,k] <= min[k])
#	}
#	include <- !rm.max & !rm.min

	usePrior=match.arg(as.character(usePrior),c("yes","no"))
	if(usePrior=="yes"){
		if(is.null(prior)){
			stop("You must specify a prior with usePrior=\"yes\"");
		}
		if (length(K)>1)
		{
			stop("You can only use a prior if the number of cluster is fixed!")
		}
		if(randomStart){
			message("randomStart>0 has no effect when using a prior. Labels are initialized from the prior.");
		}
		if(trans==2)
		{
			stop("You are using a prior with cluster specific transformations.\n  This is not recommended.")
		}
                if(any(is.infinite(c(nu,prior$nu0)))){
                    stop("If usePrior='yes', nu or nu0 may not be Inf");
                }
                if(lambda!=1&trans==1){
                    warning("Use of a prior with transformation estimation and lambda!=1 requires the prior means to be on the transformed scale.")
                }
		# TODO Add tests for validity of w0. Set a default. Same for oorder.
		if(!is.null(prior)&length(K)==1){
			#Check that the prior dimensions match the model being fit.
			mp<-ncol(prior$Mu0);
			mk<-nrow(prior$Mu0);
			ld<-dim(prior$Lambda0)
			lomega<-dim(prior$Omega0);
			#message(lomega);
			lnu0<-length(prior$nu0);
			py<-ncol(y);
			warn.o<-options("warn")
			options(warn=-1);
			if(any(is.null(prior$w0)))stop("w0 prior should be defined")
			if(length(prior$w0)!=K|(length(ld)!=3)|((lnu0!=1)&(lnu0!=mk))|(mp!=py)|(mk!=K)|(ld[1]!=K)|(ld[2]!=py)|(ld[3]!=py)|((any(lomega!=c(K,py,py)))&(any(lomega!=c(py,py))))){
				stop("Prior dimensions must match the number of clusters and the number of dimensions being analyzed.")
			}
			if(lnu0==1){
				#Extend nu0 to be a vector of length k. We allow cluster specific weight for covariance matrices.
				prior$nu0<-rep(prior$nu0,mk);
			}
			if(length(prior$w0)==1){
					#Extend w0 to be a vector of length k.
					prior$w0<-rep(prior$w0,mk);
			}
			options(warn=warn.o$warn)
		}
	}
	if (usePrior=="no")
	{
		if(!is.null(prior)){
			message("The prior specification has no effect when usePrior=",usePrior);
		}
		prior<-list(NA);
	}
	if(length(grep("parallel",loadedNamespaces()))==0||!parallel)
    {
		message("Using the serial version of flowClust")
		# C version
		result<-lapply(as.list(1:length(K)),.flowClustKfast, y
		               , varNames=varNames, K=K, B=B, tol=tol, nu=nu, lambda=lambda
		               , nu.est=nu.est, trans=trans, min.count=min.count, max.count=max.count
		               , min=min, max=max, level=level, u.cutoff=u.cutoff, z.cutoff=z.cutoff
		               , randomStart=randomStart, B.init=B.init, tol.init=tol.init, seed=seed
		               , criterion=criterion, control=control, nstart = nstart, prior,usePrior)
	}else 
	{
        require(parallel)
		cores<-getOption("cores")
		if(is.null(cores))
		{
          
			nClust<- detectCores()
		}
		else
		{
			nClust<-cores
		}
		message("Using the parallel (multicore) version of flowClust with ",nClust," cores")
		# Split into nClust segReadsList
		result<-mclapply(as.list(1:length(K)),.flowClustKfast, y
		                 , varNames=varNames, K=K, B=B, tol=tol, nu=nu, lambda=lambda
		                 , nu.est=nu.est, trans=trans, min.count=min.count, max.count=max.count
		                 , min=min, max=max, level=level, u.cutoff=u.cutoff, z.cutoff=z.cutoff
		                 , randomStart=randomStart, B.init=B.init, tol.init=tol.init, seed=seed
		                 , criterion=criterion, control=control, nstart = nstart
		                 ,  prior,usePrior, mc.preschedule=FALSE)
	}
	

	# Simply return a flowClust object
	if (length(K)==1)
	{
		result[[1]]
	}
	# Create a list flowClustList
	else
	{
		result <- new("flowClustList", result, criterion=criterion)
		result@index <- which.max(criterion(result, criterion))
		result
	}
}

.flowClustKfast<-function(i, y
                      , varNames=NULL, K, B=500, tol=1e-5, nu=4, lambda=1, nu.est=0
                      , trans=1, min.count=10, max.count=10, min=NULL, max=NULL, level=0.9
                      , u.cutoff=NULL, z.cutoff=0, randomStart=10, B.init=B, tol.init=1e-2, seed=1
                      , criterion="BIC", control=NULL, prior,usePrior, nstart = 1)
{
  
	oorder<-1:K[i]
	.model<-1; #Tells the C code whether to run ECM with non-conjugate priors, or classic flowClust.'
	match.arg(usePrior,c("yes","no"));
	switch(usePrior,
			yes=priorFlag<-1,
			no=priorFlag<-0)
	
	ly <- nrow(y)
	py <- ncol(y)
	if (min(y)<=0 && lambda<=0)
		stop("lambda must be positive when data contain zero / negative values!")
	else if(usePrior=="yes")
	{
    #If we have a prior lambda.. use it to override the specified lambda.
    if(exists("lambda",prior)){
      if(!is.null(prior$lambda)){
        if(prior$lambda!=0){
          lambda<-prior$lambda
        }
      }
    }
		Mu0<-prior$Mu0
		Lambda0<-prior$Lambda0
		if(length(dim(prior$Omega0))==2){
			Omega0<-array(prior$Omega0,c(py,py,K[i]))
		}else{
			Omega0<-aperm(prior$Omega0,c(2:3,1));
		}

		for(j in 1:K[i]){
			if(!all(as.vector(Omega0[,,j])==0)){
				#message(j)
				#message(Omega0[,,j])
				Omega0[,,j]<-try(solve((Omega0[,,j])))
			}else{
				#message(j);
				#message(Omega0[,,j])
				Omega0[,,j]<-Omega0[,,j];
			}
		}
		nu0<-prior$nu0
		w0<-prior$w0;
		kappa0<-0.0
		Lambda0<-aperm(Lambda0,c(2:3,1))
		initprec<-array(0,c(py,py,K[i]));
		#cat("initprec dim: ",dim(initprec),"\n")
		#cat("Lambda0 dim:",dim(Lambda0),"\n");
		if(!all(as.vector(Lambda0==0))){
			for(j in 1:dim(Lambda0)[3]){
				#cat("Lambda0 value:",(Lambda0[,,j]),"\n")
				#cat("initprec value:",initprec[,,j],"\n")
				initprec[,,j]<-try(chol(solve(Lambda0[,,j])),silent=TRUE);
				#cat("initprec class:",class(initprec[,,]))
				#cat("solve initprec value:",initprec[,,j],"\n")
				if(inherits(initprec[,,j],"try-error")){
					##Should do something here to catch the error.
				}
			}
		}else{
			initprec<-array(var(box(y,lambda))/(K[i]^(2/py)),c(py,py,K[i]));
		}
		priorFlag<-1
		.model<-2;
	}
	else
	{
		# Non informative prior
		Mu0<-matrix(rep(colMeans(y),K[i]),K[i],py,byrow=TRUE)
		nu0<- rep(-py-1,K[i]);
		#Don't need Omega0, we'll use the conjugate prior code.
		Omega0<-rep(0,py*py);
		kappa0<-0;
		initprec<-rep(0,K[i]*py*py);
		#Lambda0, the prior covariance
		Lambda0<-rep(0,K[i]*py*py)
		w0<-rep(0,K[i]);
		.model<-1;
	}

# to determine the rule of calling outliers
	if (nu != Inf) {
		if (is.null(u.cutoff)) {
			if (!nu.est) {
				cc <- py * qf(level, py, nu)
				u.cutoff <- (nu + py) / (nu + cc)
			}
			else {
				u.cutoff <- level     # pass level to .C call to compute u.cutoff
			}
			ruleOutliers <- c(0, level, z.cutoff)     # 0 means quantile
		}
		else {
			ruleOutliers <- c(1, u.cutoff, z.cutoff)     # 1 means cutoff
		}
	}
	else {
		if (level != 1)
			q.cutoff <- qchisq(level, py)
		else
			q.cutoff <- -1     # -1 means no outlier identification
		ruleOutliers <- c(0, level, z.cutoff)
	}


	if (is.null(control$B.lambda)) control$B.lambda <- B    # BSolve=100
	if (is.null(control$B.brent)) control$B.brent <- 10000    # iterSolveMax=50
	if (is.null(control$tol.brent)) control$tol.brent <- 1e-5    # DiffSolve=1e-3
	if (is.null(control$xLow)) control$xLow <- 0.1    # xLow=.1
	if (is.null(control$xUp)) control$xUp <- 10    # xUp=1
	if (is.null(control$nuLow)) control$nuLow <- 2    # nuLow=2
	if (is.null(control$nuUp)) control$nuUp <- 100    # nuUp=30
	if (is.null(control$seed)) control$seed <- TRUE

	ind <- 0
	if (K[i]==1)
	{
		label <- rep(1, ly)
	}
	else if (!randomStart) #kmeans initialization
	{
		if (py==1)
		{
			q <- quantile(y, seq(from=0, to=1, by=1/K[i]))
			label <- rep(0, ly)
			q[1] <- q[1]-1
			for (k in 1:K[i]) label[y>q[k] & y<=q[k+1]] <- k
		}else{
			#label<-try(kmeans(y,Mu0)$cluster,silent=TRUE)
			#if(inherits(label,"try-error"))
			label<-try(kmeans(scale(y),K[i],nstart=nstart,iter.max=100)$cluster,silent=TRUE)
		}
	}
	if(priorFlag==0)
	{ # Initialization based on short EMs with random partitions if randomStart=TRUE
		if (control$seed) set.seed(seed)
		if (randomStart==1)
		{
			label <- sample(1:K[i], ly, replace=T)
		}
		else if(randomStart==0)
		{
			M<-0
			if(inherits(label,"try-error")){
				label <- sample(1:K[i], ly, replace=T)
			}
			ind<-0;
		}
		else
		{
			maxLabel <- vector("list",randomStart)
			maxLogLike <- rep(NA,randomStart)
			for (j in 1:randomStart)
			{
				label <- sample(1:K[i], ly, replace=TRUE)
				if (nu != Inf)
				{
					#cat(initprec,"\n");
					#ordering of the priors.. used to reorder the population names;
					obj <- try(.C("flowClust", as.double(t(y)), as.integer(ly),
									as.integer(py), as.integer(K[i]),
									w=rep(0,K[i]), mu=rep(0,K[i]*py),
									precision=initprec,
									lambda=as.double(rep(lambda, length.out=(if (trans>1) K[i] else 1))),
									nu=as.double(rep(nu,K[i])),
									z=rep(0,ly*K[i]), u=rep(0,ly*K[i]),
									as.integer(label), uncertainty=double(ly),
									as.double(rep(u.cutoff,K[i])), as.double(z.cutoff),
									flagOutliers=integer(ly), as.integer(B.init),
									as.double(tol.init), as.integer(trans),
									as.integer(nu.est), logLike=as.double(0),
									as.integer(control$B.lambda), as.integer(control$B.brent),
									as.double(control$tol.brent), as.double(control$xLow),
									as.double(control$xUp), as.double(control$nuLow),
									as.double(control$nuUp),
									mu0=as.double(t(Mu0)),
									as.double(kappa0),
									nu0=as.double(nu0),
									lambda0=as.double(Lambda0),
									omega0=as.double(Omega0),
									w0=as.double(w0),
									as.integer(.model),
									oorder=as.integer(oorder),
									PACKAGE="flowClust"))
									if (class(obj)=="try-error")
									{
										message("flowClust failed")
									}
				}
				else
				{
					obj <- try(.C("flowClustGaussian", as.double(t(y)), as.integer(ly),
									as.integer(py), as.integer(K[i]),
									w=rep(0,K[i]), mu=rep(0,K[i]*py),
									precision=initprec,
									lambda=as.double(rep(lambda, length.out=(if (trans>1) K[i] else 1))),
									z=rep(0, ly*K[i]), u=rep(0,ly*K[i]),
									as.integer(label), uncertainty=double(ly),
									as.double(q.cutoff), as.double(z.cutoff),
									flagOutliers=integer(ly), as.integer(B.init),
									as.double(tol.init), as.integer(trans),
									logLike=as.double(0),
									as.integer(control$B.lambda), as.integer(control$B.brent),
									as.double(control$tol.brent), as.double(control$xLow),
									as.double(control$xUp),
									as.double(t(Mu0)),
									as.double(kappa0),
									as.double(nu0),
									as.double(Lambda0),
									as.double(Omega0),
									as.integer(.model),
									PACKAGE="flowClust"))
				}
				if (class(obj)!="try-error")
				{
					maxLabel[[j]] <- label
					maxLogLike[j] <- obj$logLike
				}
			}
			ind <- order(maxLogLike, decreasing=T, na.last=NA)
		}
	}
	else
	{
		if(usePrior=="yes")
		{
			M<-0
                        #FIXME priors with transformation is currently borked
			#Assuming Mu0 is on transformed scale
			if(lambda!=1){
				ytmp<-box(y,lambda);
			}else{
				ytmp<-y;
			}
      # We use the prior densities to initialize the cluster labels to the
      # mixture component (cluster) that maximizes the prior probability.
      prob <- sapply(seq_len(K), function(k) {
        with(prior, w0[k] * dmvt(x = ytmp, mu = Mu0[k, ], sigma = Lambda0[k,,],
                                 nu = nu0[k], lambda = lambda)$value)
      })
      label <- apply(prob, 1, which.max)
		}
	}


# long EMs
	for (M in ind)
	{
		if (nu != Inf)
		{
			obj <- try(.C("flowClust", as.double(t(y)), as.integer(ly),
							as.integer(py), as.integer(K[i]),
							w=rep(0,K[i]), mu=rep(0,K[i]*py),
							precision=initprec,
							lambda=as.double(rep(lambda, length.out=(if (trans>1) K[i] else 1))),
							nu=as.double(rep(nu,K[i])),
							z=rep(0,ly*K[i]), u=rep(0,ly*K[i]),
							as.integer(if (M==0) label else maxLabel[[M]]), uncertainty=double(ly),
							as.double(rep(u.cutoff,K[i])), as.double(z.cutoff),
							flagOutliers=integer(ly), as.integer(B),
							as.double(tol), as.integer(trans),
							as.integer(nu.est), logLike=as.double(0),
							as.integer(control$B.lambda), as.integer(control$B.brent),
							as.double(control$tol.brent), as.double(control$xLow),
							as.double(control$xUp), as.double(control$nuLow),
							as.double(control$nuUp),
							mu0=as.double(t(Mu0)),
							as.double(kappa0),
							nu0=as.double(nu0),
							lambda0=as.double(Lambda0),
							omega0=as.double(Omega0),
							w0=as.double(w0),
							as.integer(.model),
							oorder=as.integer(oorder),
							PACKAGE="flowClust"))
							if (class(obj)=="try-error"){
								message("flowClust failed")
							}

		}
		else
		{
			obj <- try(.C("flowClustGaussian", as.double(t(y)), as.integer(ly),
							as.integer(py), as.integer(K[i]),
							w=rep(0,K[i]), mu=rep(0,K[i]*py),
							precision=initprec,
							lambda=as.double(rep(lambda, length.out=(if (trans>1) K[i] else 1))),
							z=rep(0,ly*K[i]), u=rep(0,ly*K[i]),
							as.integer(if (M==0) label else maxLabel[[M]]), uncertainty=double(ly),
							as.double(q.cutoff), as.double(z.cutoff),
							flagOutliers=integer(ly), as.integer(B),
							as.double(tol), as.integer(trans),
							logLike=as.double(0),
							as.integer(control$B.lambda), as.integer(control$B.brent),
							as.double(control$tol.brent), as.double(control$xLow),
							as.double(control$xUp),
							as.double(t(Mu0)),
							as.double(kappa0),
							as.double(nu0),
							as.double(Lambda0),
							as.double(Omega0),
							as.integer(.model),
							PACKAGE="flowClust"))
			if (class(obj)!="try-error"){

				obj$nu <- Inf
			}
		}
		if (class(obj)!="try-error")
			break
	}
	if (class(obj)=="try-error")
		stop(geterrmessage())

	if(usePrior=="vague"){
		trans<-1;
	}
# output obj$precision to sigma
	sigma <- array(0, c(K[i], py, py))
	precision <- matrix(obj$precision, K[i], py * py, byrow=TRUE)
	for (k in 1:K[i])
		sigma[k,,] <- matrix(precision[k,], py, py, byrow = TRUE)

# output BIC & ICL
	BIC <- 2*obj$logLike - log(ly) * (K[i]*(py+1)*py/2 + K[i]*py + K[i]-1 + (if (trans>1) K[i] else trans) + (if (nu.est>1) K[i] else abs(nu.est)))
	z <- matrix(obj$z, ly, K[i], byrow = TRUE)
	ICL <- BIC + 2 * sum(z*log(z), na.rm = TRUE)

# output z, u, label, uncertainty, flagOutliers
	z <- u <- matrix(NA, ly, K[i])
	z <- matrix(obj$z, ly, K[i], byrow=TRUE)
	u <- matrix(obj$u, ly, K[i], byrow=TRUE)
#cat(M);
	tempLabel <- if (M==0) label else maxLabel[[M]]
	label <- uncertainty <- flagOutliers <- rep(NA, ly)
	label <- tempLabel
	uncertainty <- obj$uncertainty
	flagOutliers <- as.logical(obj$flagOutliers)

# output reordered prior
	prior$Mu0<-matrix({if(all(!is.null(obj$mu0))){obj$mu0}else{NA}},K[i],py,byrow=TRUE);
	prior$Lambda0<-aperm(array({if(all(!is.null(obj$lambda0))){obj$lambda0}else{NA}},c(py,py,K[i])),c(3,1:2))
	prior$Omega0<-aperm(array({if(all(!is.null(obj$omega0))){obj$omega0}else{NA}},c(py,py,K[i])),c(3,1:2))
	prior$nu0<-{if(all(!is.null(obj$nu0))){obj$nu0}else{NA}}
	prior$w0<-{if(all(!is.null(obj$w0))){obj$nu0}else{NA}}
#omit
#result<- new("flowClust", expName=expName, varNames=varNames, K=K[i],
#		w=obj$w, mu=matrix(obj$mu, K[i], py, byrow=TRUE), sigma=sigma,
#		lambda=(if (trans>0) obj$lambda else numeric(0)), nu=(if (nu.est>1) obj$nu else obj$nu[1]), z=z,
#		u=u, label=label, uncertainty=uncertainty,
#		ruleOutliers=ruleOutliers, flagOutliers=flagOutliers, rm.min=sum(rm.min),
#		rm.max=sum(rm.max), logLike=obj$logLike, BIC=BIC, ICL=ICL);
class(prior)<-"list";
prior$order<-obj$oorder;
	if(trans==1&obj$lambda==1){
		obj$mu<-rbox(obj$mu,obj$lambda)
	}
	#do nothing in particular if trans>1
        #Not sure if the above does the right thing when trans=1, and the returned lambda=1, and usePrior="no"
	result<- new("flowClust", expName="expName", varNames=varNames, K=K[i],
			w=obj$w, mu=matrix(obj$mu, K[i], py, byrow=TRUE), sigma=sigma,
			lambda= obj$lambda, nu=(if (nu.est>1) obj$nu else obj$nu[1]), z=z,
			u=u, label=label, uncertainty=uncertainty,
			ruleOutliers=ruleOutliers, flagOutliers=flagOutliers, rm.min=0,
			rm.max=0, logLike=obj$logLike, BIC=BIC, ICL=ICL,prior=prior);
	# if(!any(is.na(prior))&usePrior=="yes"&ruleOutliers[1]==0){
		# label<-.fcbMap(result,ruleOutliers[2])
		# result@flagOutliers<-label==0
		# label[result@flagOutliers]<-NA;
		# result@label<-label;
	# }
	result
}
