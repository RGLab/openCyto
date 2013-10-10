# 
# ICS
# 
###############################################################################
#unloadNamespace("flowIncubator")
unloadNamespace("openCyto")
#unloadNamespace("QUALIFIER")
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

#modify functions within package namespace
environment(.gatingTemplate) <- getNamespace("openCyto")
assignInNamespace(".gatingTemplate", .gatingTemplate, ns = "openCyto")

#lapply(list.files("/home/wjiang2/rglab/workspace/openCyto/R/",full=T),source)
#source("/home/wjiang2/rglab/workspace/openCyto/R/AllClasses.R")
#source("/home/wjiang2/rglab/workspace/openCyto/R/gatingTemplate-methods.R")
#source("/home/wjiang2/rglab/workspace/openCyto/R/gating-cytokines.R")
#source("/home/wjiang2/rglab/workspace/openCyto/R/gating-methods.R")
#source("/home/wjiang2/rglab/workspace/openCyto/R/gating-functions.R")
#source("/home/wjiang2/rglab/workspace/openCyto/R/gtMethod-methods.R")
#source("/home/wjiang2/rglab/workspace/openCyto/R/gtPopulation-methods.R")
#source("/home/wjiang2/rglab/workspace/openCyto/R/fcTree-methods.R")
#source("/home/wjiang2/rglab/workspace/openCyto/R/fcFilterList-methods.R")
#source("/home/wjiang2/rglab/workspace/openCyto/R/fcFilter-methods.R")
#source("/home/wjiang2/rglab/workspace/openCyto/R/bayes-flowClust.R")
#source("/home/wjiang2/rglab/workspace/openCyto/R/functions.R")
#source("/home/wjiang2/rglab/workspace/openCyto/R/bayes-flowClust.R")
#source("/home/wjiang2/rglab/workspace/openCyto/R/wrapper-functions.R")
#source("/home/wjiang2/rglab/workspace/openCyto/R/preprocessing-method.R")

path<-"/home/wjiang2/rglab/workspace/openCyto"
library(Cairo)
CairoX11()
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

gs<-GatingSet(fs[1:2])
env1<-new.env(parent=emptyenv())
gating(gt,gs,env1)
plot(gs[[1]])
plot(gs[[1]],bool=F)
getNodes(gs[[1]])
plotGate(gs[[1]],xbin=32,margin=T)
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
gating(gt1,gs1,env1,mc.cores=4,parallel_type = "multicore")
plotGate(gs1[[1]],xbin=32)
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
plotGate(gs2[[1]],bool=T,xbin=32)
getData(gs2[[1]],"cd19&cd20")
getNodes(gs2[[1]])
densityplot(~`<G780-A>`,getData(gs2,"cd19&!cd20"))

getGate(env1$fct,"nonDebris")
plot(env1$fct,"nonDebris",post=T)
plot(env1$fct,"cd19",post=T)
plot(env1$fct,"IgD-cd27+",channel="<G780-A>",post=T)
plot(env1$fct,"IgD-cd27+",post=T,channel="<V545-A>")


##debug John's gt

#bcell

gating_template <- gatingTemplate(file.path(path,"data/gt-bcell.csv"))
gs <- load_gs(path = file.path(path,"data/gs-bcell"))
plot(gating_template)
getNodes(gs[[1]])
Rm("boundary",gs)
gating(gating_template, gs
    , mc.cores = 3, parallel_type = "multicore"
)

plotGate(gs[[1]],xbin=32,margin=T)
#t-reg

#gs <- load_gs("/loc/no-backup/ramey/Lyoplate/gating-sets/gs-treg")
gs <- load_gs(file.path(path,"data/gs-treg"))
gating_template <- gatingTemplate(file.path(path,"data/gt-treg.csv"))
#gs_sub <- gs[1:3]
#gs_sub <- clone(gs_sub)
#save_gs(gs_sub,path = file.path(path,"data/gs-treg"))

getNodes(gs[[1]])
#Rm("boundary",gs)
Rm("CD4",gs)
#Rm("memory",gs)
plot(gs[[1]])
gating(gating_template, gs
    , mc.cores = 3, parallel_type = "multicore"
)

plotGate(gs[[1]],xbin=32,margin=T,bool =T)

plot(gating_template
       , graph =list(rankdir ="TB")
#    ,y=c(   "CD25+CD127-"
#            ,"Memory"
#            ,"CCR4+CD45RO+"
#            ,"CCR4+"
#            ,"CD45RO+")
    ,y = "CD4"
#    , showRef = F
    )
dev.off()

#TODO: HVTN065 (plotGate doesn't work because range slot is invalid( character) 
#gs <- load_gs("/loc/no-backup/ramey/HVTN/065/gating-set/")
#gs_sub <- clone(gs[1])
#getNodes(gs_sub[[1]])
#getData(gs_sub[[1]],"cd3")
#fr <- getData(gs_sub[[1]])
#hist(exprs(fr)[,"PE Tx RD-A"])
#
#densityplot(~`SSC-A`,fr)
#
#xyplot(`SSC-A`~`PE Tx RD-A`
#    , getData(gs_sub[[1]]
##            ,"cd3"
#              )
##    ,filter =getGate(gs_sub[[1]],"cd4")
#)
library(Rcpp)
library(inline)
getNodes(gs[[1]])
plotGate(gs[[1]],"cd4+",xbin=32)
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

 