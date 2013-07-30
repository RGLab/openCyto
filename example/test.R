# 
# ICS
# 
###############################################################################
unloadNamespace("flowIncubator")
unloadNamespace("openCyto")
unloadNamespace("QUALIFIER")
unloadNamespace("flowStats")
unloadNamespace("flowWorkspace")


library(openCyto)
#library(flowIncubator)
library(flowWorkspace)
library(flowClust)
library(flowStats)
library(MASS)
library(plyr)
library(clue)
library(gtools)

source("/home/wjiang2/rglab/workspace/openCyto/R/AllClasses.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gatingTemplate-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gating-cytokines.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gating-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gating-functions.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gtMethod-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/gtPopulation-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/fcTree-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/fcFilterList-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/fcFilter-methods.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/bayes-flowClust.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/functions.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/bayes-flowClust.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/wrapper-functions.R")
source("/home/wjiang2/rglab/workspace/openCyto/R/preprocessing-method.R")

path<-"/home/wjiang2/rglab/workspace/openCyto"

###############
##ICS
###############
gt<-gatingTemplate(file.path(path,"data/ICS.csv"),"ICS")
gt
plot(gt)

##transform the ICS data
load(file.path(path,"data/fs_080.rda"))
#paramters<-colnames(fs[[1]])
#trans <- estimateLogicle(fs[[1]], channels = paramters[!grepl("[F|S]SC|[T|t]ime",paramters)])
#fs_trans<-transform(fs,trans)

gs<-GatingSet(fs)
env1<-new.env(parent=emptyenv())
gating(gt,gs[1],env1)
plot(gs[[1]])
plot(gs[[1]],bool=T)
getNodes(gs[[1]])
plotGate(gs[[1]],xbin=64,margin=T)
xyplot(`SSC-A`~`<FITC-A>`,fr,smooth=F)

getPopStats(gs[[1]])[22,]

#plot priors
plotGate(x=gs,3)
X11()
plot(env1$fct,"v",posteriors=T)
plotGate(x=gs,4)
plot(env1$fct,"nonDebris",posteriors=T)
plot(env1$fct,"cd4",posteriors=T,channel="PE Cy55-A")
plot(env1$fct,"cd4",posteriors=T,channel="FITC-A")
plot(env1$fct,"TNFa",posteriors=T)

###############
#Tcell is already transformed
###############
gt1<-gatingTemplate(file.path(path,"data/Cytotrol_Tcell.csv"),"Tcell")
plot(gt1)
getNodes(gs1[[1]])
Rm("cd3",gs1)
load(file.path(path,"data/fs_tcell.rda"))
gs1<-GatingSet(fs_tcell)
env1<-new.env(parent=emptyenv())
gating(gt1,gs1,env1)
plotGate(gs1[[1]],xbin=64)
plot(env1$fct,"nonDebris",post=T)
plot(env1$fct,"cd3",post=T)
plot(env1$fct,"cd4+",post=T,channel="<B710-A>")
plot(env1$fct,"activated cd4",post=T,channel="<R660-A>")
plot(env1$fct,2,post=T)
plot(env1$fct,"lymph",post=T)

#xyplot(`<B710-A>`~`<R660-A>`,fs_tcell)
#densityplot(~.,fs_tcell[[1]])

###############
#Bcell is already transformed
###############
gt2<-gatingTemplate(file.path(path,"data/Cytotrol_Bcell.csv"),"Bcell")
plot(gt2)

getNodes(gs2[[1]])
load(file.path(path,"data/fs_bcell.rda"))
gs2<-GatingSet(fs_bcell)
env1<-new.env(parent=emptyenv())
gating(gt2,gs2,env1)
getGate(gs2,4)
plot(gs2[[1]],bool=T)
plotGate(gs2[[1]],bool=T,xbin=64)
getData(gs2[[1]],"cd19&cd20")
getNodes(gs2[[1]])
densityplot(~`<G780-A>`,getData(gs2,"cd19&!cd20"))

getGate(env1$fct,"nonDebris")
plot(env1$fct,"nonDebris",post=T)
plot(env1$fct,"cd19",post=T)
plot(env1$fct,"IgD-cd27+",channel="<G780-A>",post=T)
plot(env1$fct,"IgD-cd27+",post=T,channel="<V545-A>")

###debug
library(openCyto)
gs_HVTN065 <- load_gs("/loc/no-backup/ramey/HVTN/065/gating-results")
gating_template <- gatingTemplate("/home/jramey/rglab/papers/paper-opencyto/gt-HVTN065.csv", "HVTN065")
Rm("cd4", gs_HVTN065)
Rm("cd8", gs_HVTN065)
gs <- clone(gs_HVTN065[1:12])
getData(gs)
Rm("cd4", gs)
Rm("cd8", gs)

gating(gating_template, gs, prior_group = 'Stim'
    , num_cores = 6, parallel_type = "MPI"
)

set.seed(42)

n <- 50

x1 <- rnorm(n, mean = 5)
x2 <- rnorm(n, mean = 0)
x3 <- rnorm(n, mean = 20)

plot(density(x1), xlim = c(-5, 25))
lines(density(x2), col = "red")
lines(density(x3), col = "blue")

ecdf1 <- ecdf(x1)
ecdf2 <- ecdf(x2)
ecdf3 <- ecdf(x3)

plot(ecdf1, xlim = c(-5, 25))
lines(ecdf2, col = "red")
lines(ecdf3, col = "blue")

# The statistics are equal. So with this criterion, samples 1 and 2 have the same distance as samples 1 and 3
ks.test(x1, x2)$statistic
ks.test(x1, x3)$statistic



library(Rcpp)
library(inline)

cstring <- '
double * _s = REAL(s);
double * _t1 = REAL(t1);
double * _t2 = REAL(t2);
*_t1 =*_s;
*_t2 =  (float)(*_s);


return R_NilValue;
'

funx <- cfunction(signature(s="numeric",t1="numeric",t2 = "numeric")
    ,cstring
    
    ,Rcpp=TRUE
#    ,cppargs="-I/usr/include"
#    ,libargs="-lgsl -lgslcblas"
)
s=.Machine$double.xmax/100;
t1=0;
t2=0
funx(s,t1,t2)
s
t1
t2

p11 <- histogram( ~ height | voice.part, data = singer, xlab="Height")
p12 <- densityplot( ~ height | voice.part, data = singer, xlab = "Height")
p2 <- histogram( ~ height, data = singer, xlab = "Height")


## simple positioning by split
print(p11, split=c(1,1,1,2), more=TRUE)
print(p2, split=c(1,2,1,2))

## Combining split and position:
print(p11, position = c(0,0,.75,.75), split=c(1,1,1,2), more=TRUE)
print(p12, position = c(0,0,.75,.75), split=c(1,2,1,2), more=TRUE)
print(p2, position = c(.5,.75,1,1), more=FALSE)

## Using seekViewport

## repeat same plot, with different polynomial fits in each panel
xyplot(Armed.Forces ~ Year, longley, index.cond = list(rep(1, 6)),
    layout = c(3, 2),
    panel = function(x, y, ...)
    {
      panel.xyplot(x, y, ...)
      fm <- lm(y ~ poly(x, panel.number()))
      llines(x, predict(fm))
    })

## Not run: 
grid::seekViewport(trellis.vpname("panel", 1, 1))
cat("Click somewhere inside the first panel:\n")
ltext(grid::grid.locator(), lab = "linear")

## End(Not run)

grid::seekViewport(trellis.vpname("panel", 1, 1))
grid::grid.text("linear")

grid::seekViewport(trellis.vpname("panel", 2, 1))
grid::grid.text("quadratic")

grid::seekViewport(trellis.vpname("panel", 3, 1))
grid::grid.text("cubic")

grid::seekViewport(trellis.vpname("panel", 1, 2))
grid::grid.text("degree 4")

grid::seekViewport(trellis.vpname("panel", 2, 2))
grid::grid.text("degree 5")

grid::seekViewport(trellis.vpname("panel", 3, 2))
grid::grid.text("degree 6")
