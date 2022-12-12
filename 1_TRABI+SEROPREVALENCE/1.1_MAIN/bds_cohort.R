require(data.table)
####################################################################
require(data.table)
require(ggplot2)
system("mkdir Figs/kim_cohort/")
system("mkdir Results/kim_cohort/")
system("mkdir interimData/kim_cohort/")

source("Code/prepareDataHelper.R")
source("Code/helperFunctions.R")


aa=load("interimData/kim_cohort/summary_fit.RData")

uscTbl=fread("Results/usz_cohort/fullUscTblblood.txt")

uscTbl[,Pseudo_ID:=Patient_or_Control_ID]
uscTbl[neglog< -3,neglog:= -3]


tbl1=uscTbl[mydate<200200,]
tbl1[,myy:=(mydate==-1)+0]


##tbl1[sample_cat==categ,myy:=1]
myl=list()
myl[[1]]=tbl1[,list(target="averaged",neglog=mean(neglog),myy=myy[1]),by=Pseudo_ID]
myl[[2]]=tbl1[target!="NC_IgG",list(target="averagedSpikeRBD",neglog=mean(neglog),myy=myy[1]),by=Pseudo_ID]
myl[[3]]=tbl1[,list(target="min",neglog=min(neglog),myy=myy[1]),by=Pseudo_ID]
myl[[4]]=tbl1[target!="NC_IgG",list(target="minSpikeRBD",neglog=min(neglog),myy=myy[1]),by=Pseudo_ID]
myl[[5]]=getPredict(tbl1,kimfit)
tbl2=tbl1[,list(Pseudo_ID,target,neglog,myy)]
myr=do.call("rbind",myl)
myr=rbind(tbl2,myr)
setkey(myr,neglog)

tt=lapply(unique(myr[,target]),function(mytarget){
	xx=getRocValsInternal(myr[target==mytarget,myy], myr[target==mytarget,neglog])
	xx[,target:=mytarget]
	return(xx)
})

tbl=do.call("rbind",tt)
write.table(tbl,file="interimData/bds_cohort/bds_cohort_RocData.R",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
ggplot(tbl, aes(x=fpr,y=tpr, group=target,colour=target)) +geom_line(colour="black",stat="identity",position="dodge") + theme_bw()
ggplot(tbl, aes(x=fpr,y=tpr, group=target,colour=target)) +geom_line() +theme_bw() + xlim(0,0.05)





####################################################################
uscTbl=fread("Results/usz_cohort/fullUscTblblood.txt")
uscTbl[neglog< -3,neglog:=-Inf]
uscTbl=uscTbl[target=="Spike_IgG",]


#uscTbl=uscTbl[target!="NC_IgG",list(target="averagedSpikeRBD",mydate=mydate[1],neglog=mean(neglog)),by=Patient_or_Control_ID]
valTbl=uscTbl[mydate==210000,]
valTblUsed=valTbl[neglog> -3,]
valTblNotused=valTbl[neglog<= -3,]
tbl1=uscTbl[mydate<200200,]
#tbl1=uscTbl[mydate<=180000,]
tbl1[,myy:=(mydate==-1)+0]



tbl1[myy==1,neglog]
sd2=sd(tbl1[myy==1,neglog])
mean2=mean(tbl1[myy==1,neglog])

sd1=sd(tbl1[myy==0 & neglog> -3,neglog])
mean1=mean(tbl1[myy==0 & neglog> -3,neglog])

myvals=valTblUsed[,neglog]
ee=estimateFrac(myvals,mean1,sd1,mean2,sd2)

ss=data.table(prob=ee[[2]],myy=valTblUsed[,neglog])
ggplot(ss,aes(x=myy, group=prob>0.5,fill=prob>0.5))+geom_histogram(alpha=0.5) +facet_grid((prob>0.5)~.)
ggsave("Figs/estimatedProb.pdf")
dim(valTblUsed)[1]

valTblNotused=valTbl[neglog<= -3,]

potential=nrow(valTblUsed)
impossible=nrow(valTblNotused)
estimated_prob=(potential*(1-ee[[1]]))/(potential+impossible)

setkey(ss,myy)
mycutoff=tail(ss[prob>0.5,],1)[,myy]
uscTbl[neglog>mycutoff & mydate!=-1 & mydate!=210000,]




####################################################################

####################################################################
ONLYHISTORICAL=TRUE
uscTbl=fread("Results/usz_cohort/fullUscTblblood.txt")
uscTbl[neglog< -3,neglog:=-Inf]
##uscTbl=uscTbl[target=="Spike_IgG",]
###########
uscTbl=uscTbl[target!="NC_IgG",list(mydate=mydate[1], target="averagedSpikeRBD",neglog=mean(neglog)*sqrt(2)),by=Patient_or_Control_ID]
if(ONLYHISTORICAL){
	posNegTbl=uscTbl[mydate<=180000,]
	valTbl=uscTbl[mydate>180000,]
	}else{
	posNegTbl=uscTbl[mydate<200200,]
	valTbl=uscTbl[mydate>=200200,]
}
posNegTbl[,myy:=(mydate==-1)+0]
posNegTbl=posNegTbl
#valTbl=uscTbl[mydate==210000,]

print("number of positives.")
print(posNegTbl[,sum(myy)])
#bloodPrevalenceFits=getPrevalence(posNegTbl, valTbl, sureNegLevel= -3)
bloodPrevalenceFits=getPrevalenceFull(posNegTbl, valTbl, sureNegLevel= -3)

myvalTbl=bloodPrevalenceFits$resTbl
myvalTbl[,mymonth:=floor(mydate/100)]
prevalenceTbls=myvalTbl[,list(prevalenceEstimate=mean(posteriorProb),len=length(posteriorProb)),by=mymonth]
prevalenceTbls[,method:="gaussianMethod"]
save(bloodPrevalenceFits,file="Results/bloodPrevalenceFitsGaussian.RData")
####################################################################
####################################################################
if(ONLYHISTORICAL){
	fprVsCutoffTbl=fread("interimData/cutoffsblood_onlyHistorical.txt")
	}else{
	fprVsCutoffTbl=fread("interimData/cutoffsBlood.txt")

}

uscTbl=fread("Results/usz_cohort/fullUscTblblood.txt")

uscTbl=uscTbl[target!="NC_IgG",list(mydate=mydate[1], target="averagedSpikeRBD",neglog=mean(neglog)*sqrt(2)),by=Patient_or_Control_ID]
mycutoff=fprVsCutoffTbl[method=="neglog" & fpr==0.001,cutoff]
uscTbl[,mymonth:=round(mydate/100)]
monthTbl=uscTbl[mydate!=-1,list(prevalenceEstimate=mean(neglog>mycutoff),len=length(neglog)),by=mymonth]
monthTbl[,method:="fdrMethod"]
outTbl=rbind(monthTbl,prevalenceTbls)
print(outTbl)
write.table(outTbl,file=paste("Results/bloodPrevalenceFitsBoth",("OnlyHistorical")[ONLYHISTORICAL],".txt",sep=""),sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)



outTbl[,Month:=factor(mymonth)]
ggplot(outTbl[mymonth!=1800,],aes(x=Month,y=prevalenceEstimate,group=method,fill=method)) + geom_bar(stat="identity",position="dodge") + ylim(0,0.017) + theme_bw()
ggsave(paste("Figs/bdsPrevalences",("OnlyHistorical")[ONLYHISTORICAL],".pdf",sep=""),width=4,height=4)





