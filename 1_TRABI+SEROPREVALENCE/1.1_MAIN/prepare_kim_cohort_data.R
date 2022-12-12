require(data.table)
####################################################################
require(data.table)
require(ggplot2)
system("mkdir Figs/kim_cohort/")
system("mkdir Results/kim_cohort/")
system("mkdir interimData/kim_cohort/")

source("Code/prepareDataHelper.R")
source("Code/helperFunctions.R")

USZNORMED=TRUE


filepaths=do.call("c",lapply(c("RBD","NC","Spike"), getInHouseFitPaths,useImunology=FALSE,useKim=TRUE))

inhouseTbl=do.call("rbind",lapply(filepaths,fread))
inhouseTbl[,neglog:=-log10(IC50)]

inhouseTbl[,idWithoutMonth:=sub("^000","", sub("\\d\\d$","", Patient_or_Control_ID))]
inhouseTbl[,neglog:=-log10(IC50)]
inhouseTbl[,Pseudo_ID:=Patient_or_Control_ID]

kimTbl=fread("Data/SARS_CoV_2_share/20200506_KIM_reference_dataset/COVID19_samples_Dropbox.csv")
kimTbl[,timepoint:=`time point (days since first symptoms)`]
hasAgeButNoTimepoint=kimTbl[,!is.na(Age) & is.na(timepoint)]
kimTbl=kimTbl[!hasAgeButNoTimepoint,]
setkey(inhouseTbl,Pseudo_ID)
setkey(kimTbl,Pseudo_ID)


mergedTbl=merge(inhouseTbl,kimTbl)
mergedTbl[neglog< -3,neglog:= -3 ]
mytbl2=mergedTbl[!is.na(timepoint),]

mytbl2[,neglogAbove0:=neglog]
mytbl2[neglogAbove0<0,neglogAbove0:=0]
ggplot(mytbl2[target=="Spike_IgG",],aes(y=timepoint,x=neglog)) + ylim(0,4)+ theme_bw() + geom_smooth(method="lm") + xlab("-logEC50 spike")
ggsave("Figs/kim_cohort/spike_vs_time.pdf", width=4, height=4) 


ggplot(mytbl2[target=="Spike_IgG",],aes(x=timepoint,y=neglogAbove0,group=Pseudo_Pat_ID,colour=Pseudo_Pat_ID)) + geom_point() + geom_line(linetype=1) + theme_bw()  + ylab("-logEC50 spike")
ggsave("Figs/tmp.pdf", width=8, height=4) 


ggplot(mytbl2[target=="Spike_IgG",],aes(x=timepoint,y=neglogAbove0,group=Pseudo_Pat_ID,colour=Pseudo_Pat_ID)) + geom_point() + geom_line(linetype=1) + theme_bw() + theme(legend.position = "none") + ylab("-logEC50 spike")
ggsave("Figs/kim_cohort/spike_vs_time_perPat.pdf", width=4, height=4) 
ggsave("Figs/Preprint/Fig1D.pdf")

ggplot(mytbl2[target=="RBD_IgG",],aes(y=timepoint,x=neglog)) + geom_point() + theme_bw() + geom_smooth(method="lm") + xlab("-logEC50 RBD")
ggsave("Figs/kim_cohort/RBD_vs_time.pdf", width=4, height=4)

ggplot(mytbl2[target=="NC_IgG",],aes(y=timepoint,x=neglog)) + geom_point() + theme_bw() + geom_smooth(method="lm") + xlab("-logEC50 NC")
ggsave("Figs/kim_cohort/NC_vs_time.pdf", width=4, height=4)


castTbl=data.table(dcast(mytbl2,timepoint+Pseudo_ID~target,value.var="neglog"))
castTbl[,mypred:=RBD_IgG+Spike_IgG]

ggplot(castTbl,aes(y=timepoint,x=mypred)) + geom_point() + theme_bw() + geom_smooth(method="lm") + xlab("-logEC50 Combined")
ggsave("Figs/kim_cohort/combined_vs_time.pdf", width=4, height=4)



#castTbl2=data.table(dcast(mergedTbl,timepoint+Pseudo_ID~target,value.var="neglog"))
mergedTblSub=mergedTbl[,list(Pseudo_ID,timepoint,target,neglog)]
subt=mergedTblSub[target!="NC_IgG",list(timepoint=timepoint[1],target="averagedSpikeRBD",neglog=mean(neglog)),by=Pseudo_ID]
subt2=mergedTblSub[,list(timepoint=timepoint[1],target="averaged",neglog=mean(neglog)),by=Pseudo_ID]
subt3=mergedTblSub[,list(timepoint=timepoint[1],target="min",neglog=min(neglog)),by=Pseudo_ID]
subt4=mergedTblSub[target!="NC_IgG",list(timepoint=timepoint[1],target="minSpikeRBD",neglog=min(neglog)),by=Pseudo_ID]
mergedTblSub=rbind(mergedTblSub,subt)
mergedTblSub=rbind(mergedTblSub,subt2)
mergedTblSub=rbind(mergedTblSub,subt3)
mergedTblSub=rbind(mergedTblSub,subt4)
#subt=mergedTblSub[,list(timepoint,target="averaged_sqrt",neglog=mean(sqrt(abs(neglog))*sign(neglog))),by=Pseudo_ID]
#mergedTblSub=rbind(mergedTblSub,subt)

if(USZNORMED){
	myg=fread("Results/usz_cohort/uscTbl.txt")
	negTbl=myg[myy==0,list(Pseudo_ID=editedId,timepoint=NA, target=target,neglog=neglog,myy)]
	tmpTbl=negTbl[,list(timepoint=timepoint[1],target="averaged",neglog=mean(neglog),myy=0),by=Pseudo_ID]
	tmpTbl2=negTbl[target!="NC_IgG",list(timepoint=timepoint[1],target="averagedSpikeRBD",neglog=mean(neglog),myy=0),by=Pseudo_ID]
	tmpTbl3=negTbl[,list(timepoint=timepoint[1],target="min",neglog=min(neglog),myy=0),by=Pseudo_ID]
	tmpTbl4=negTbl[target!="NC_IgG",list(timepoint=timepoint[1],target="minSpikeRBD",neglog=min(neglog),myy=0),by=Pseudo_ID]
	negTbl=rbind(negTbl,tmpTbl)
	negTbl=rbind(negTbl,tmpTbl2)
	negTbl=rbind(negTbl,tmpTbl3)
	negTbl=rbind(negTbl,tmpTbl4)
	negTbl=negTbl[!duplicated(paste(Pseudo_ID,target)),]
	#tmpTbl=negTbl[,list(timepoint,target="averaged_sqrt ",neglog=mean(sqrt(abs(neglog))*sign(neglog)),myy=0),by=Pseudo_ID]
	#negTbl=rbind(negTbl,tmpTbl)	
}else{
	negTbl=mergedTblSub[is.na(timepoint),]
	negTbl[,myy:=FALSE]
}



pos1Tbl=mergedTblSub[timepoint<7,]
pos1Tbl[,myy:=TRUE]
pos2Tbl=mergedTblSub[timepoint>6 & timepoint<14,]
pos2Tbl[,myy:=TRUE]
pos3Tbl=mergedTblSub[timepoint>13 ,]
pos3Tbl[,myy:=TRUE]
pos4Tbl=mergedTblSub[timepoint>6 ,]
pos4Tbl[,myy:=TRUE]
pos5Tbl=mergedTblSub[timepoint>13 ,]
pos5Tbl[,myy:=TRUE]

tbl1=rbind(negTbl,pos1Tbl)
setkey(tbl1,neglog)
tbl2=rbind(negTbl,pos2Tbl)
setkey(tbl2,neglog)
tbl3=rbind(negTbl,pos3Tbl)
setkey(tbl3,neglog)
tbl4=rbind(negTbl,pos4Tbl)
setkey(tbl4,neglog)
tbl5=rbind(negTbl,pos5Tbl)
setkey(tbl5,neglog)

##negTbl[,category:="negative"]
tbl1[,category:="<7 days"]
tbl2[,category:="7-13  days"]
tbl3[,category:=">13 days"]
tbl4[,category:=">6  days"]
################################################
################################################
fullTbl=rbind(rbind(rbind(tbl1,tbl2),tbl3),tbl4)
write.table(fullTbl,file="interimData/kim_cohort/kim_cohort_inputData.txt",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")


kimfit=getLogModel(fullTbl[category==">6  days",])
save(kimfit, file = "interimData/kim_cohort/summary_fit.RData")
####################################################################
####################################################################

uscTbl=fread("Results/usz_cohort/fullUscTbl.txt")
targetL=c("averagedSpikeRBD","Spike_IgG")
uszPrevalenceFits=lapply(c("averagedSpikeRBD","Spike_IgG"), function(mytarget){
#	posNegTbl=tbl3[target==mytarget,]
	posNegTbl=tbl5[target==mytarget,]
	uscTbl[neglog< -3,neglog:=-Inf]
	if(mytarget=="averagedSpikeRBD"){
		uscTbl=uscTbl[target!="NC_IgG",list(mydate=mydate[1], target="averagedSpikeRBD",neglog=mean(neglog)),by=Patient_or_Control_ID]
	}
	if(mytarget=="Spike_IgG"){
		uscTbl=uscTbl[target==mytarget,list(mydate=mydate[1], target=target[1],neglog=neglog[1]),by=Patient_or_Control_ID]
	}
	valTbl=uscTbl[mydate>180000,]
	print("number of positives.")
	print(posNegTbl[,sum(myy)])
	tt2=getPrevalenceFull(posNegTbl, valTbl, sureNegLevel= -3)
	#tt3=getPrevalence(posNegTbl, valTbl, sureNegLevel= -3)	
	outTbl=tt2[[1]]
	outTbl[,mymonth:=round(mydate/100)]
	min(outTbl[posteriorProb>0.5,neglog])
	print(tt2[[2]])
	monthTbl=outTbl[,list(prevalenceEstimate=mean(posteriorProb),len=length(posteriorProb)),by=mymonth]
	outl=list()
	outl[["monthTbl"]]=monthTbl
	outl[["prevalenceRes"]]=tt2
	return(outl)
})

names(uszPrevalenceFits)=targetL
save(uszPrevalenceFits,file="Results/uszPrevalenceFitsGaussian.RData")
####################################################################
fprVsCutoffTbl=fread("interimData/cutoffs.txt")


uscTbl=fread("Results/usz_cohort/fullUscTbl.txt")
uscTbl=uscTbl[target!="NC_IgG",list(mydate=mydate[1], target="averagedSpikeRBD",neglog=mean(neglog)*sqrt(2)),by=Patient_or_Control_ID]
mycutoff=fprVsCutoffTbl[method=="neglog" & fpr==0.001,cutoff]
uscTbl[,mymonth:=round(mydate/100)]
monthTbl=uscTbl[,list(prevalenceEstimate=mean(neglog>mycutoff),len=length(neglog)),by=mymonth]
monthTbl[,method:="fdrMethod"]
prevalenceTbls=uszPrevalenceFits[[1]]$monthTbl
prevalenceTbls[,method:="gaussianMethod"]
outTbl=rbind(monthTbl,prevalenceTbls)
print(outTbl)
write.table(outTbl,file="Results/uszPrevalenceFitsBoth.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


outTbl[,Month:=factor(mymonth)]
ggplot(outTbl[mymonth!=1800,],aes(x=Month,y=prevalenceEstimate,group=method,fill=method)) + geom_bar(stat="identity",position="dodge") + theme_bw()
ggsave("Figs/uszPrevalences.pdf")



