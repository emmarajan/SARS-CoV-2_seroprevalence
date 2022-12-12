require(data.table)
####################################################################
require(data.table)
require(ggplot2)
system("mkdir Figs/kim_cohort/")
system("mkdir Results/kim_cohort/")
system("mkdir interimData/kim_cohort/")

source("Code/prepareDataHelper.R")
source("Code/helperFunctions.R")

filepaths=do.call("c",lapply(c("RBD","NC","Spike"), getInHouseFitPaths,useImunology=FALSE,useKim=TRUE))

inhouseTbl=do.call("rbind",lapply(filepaths,fread))
inhouseTbl[,neglog:=-log10(IC50)]
inhouseTbl[,Pseudo_ID:=Patient_or_Control_ID]

kimTbl=fread("Data/SARS_CoV_2_share/20200506_KIM_reference_dataset/COVID19_samples_Dropbox.csv")
kimTbl[,timepoint:=`time point (days since first symptoms)`]
setkey(inhouseTbl,Pseudo_ID)
setkey(kimTbl,Pseudo_ID)
mergedKimTbl=merge(inhouseTbl,kimTbl)
kimPosForMerge=mergedKimTbl[!is.na(Age),list(Pseudo_ID,dataset="kim",timepoint,target,valueUsed=neglog)]
myg=fread("Results/usz_cohort/uscTbl.txt")
uscNegForMerge=myg[mydate<190000,list(Pseudo_ID=Patient_or_Control_ID,dataset="uscNeg",timepoint=NA,target=target,valueUsed=neglog)]

mergedTbl=rbind(kimPosForMerge,uscNegForMerge)

myl = list()
myl[[1]]=mergedTbl[,list(target="averaged",valueUsed=mean(valueUsed)),by=list(dataset,Pseudo_ID,timepoint)]
myl[[2]]=mergedTbl[target!="NC_IgG",list(target="averagedSpikeRBD",valueUsed=mean(valueUsed)),by=list(dataset,Pseudo_ID,timepoint)]
myl[[3]]=mergedTbl[,list(target="min",valueUsed=min(valueUsed)),by=list(dataset,Pseudo_ID,timepoint)]
myl[[4]]=mergedTbl[target!="NC_IgG",list(target="min",valueUsed=min(valueUsed)),by=list(dataset,Pseudo_ID,timepoint)]

completeTbl=rbind(mergedTbl,do.call("rbind",myl))
completeTbl[,hasCovid:=dataset=="kim"]

write.table(completeTbl,file="interimData/kim_cohort/DLdataForKimROC.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

########################################################
########################################################
########################################################

myf=function(subtbl,categorystr){
	setkey(subtbl,valueUsed)
	rr=lapply(unique(completeTbl[,target]),function(mytarget){
		xx1=getRocValsInternal(subtbl[target==mytarget,hasCovid+0], subtbl[target==mytarget,valueUsed])
		xx1[,target:=mytarget]
		return(xx1)
	})
	rr2=do.call("rbind",rr)
	rr2[,category:=categorystr]
	return(rr2)
}

myl2=list()
myl2[[1]]=myf(completeTbl[timepoint<7 | hasCovid==FALSE,],"timepoint < 7")
myl2[[2]]=myf(completeTbl[(timepoint>6 & timepoint<14)| hasCovid==FALSE,],"timepoint > 6 & timepoint < 14")
myl2[[3]]=myf(completeTbl[timepoint>13 | hasCovid==FALSE,],"timepoint > 13")
plotTbl=do.call("rbind",myl2)


write.table(plotTbl,file="interimData/kim_cohort/DLForKimROC_res.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


ggplot(plotTbl,aes(y=tpr,x=fpr,group=target,  colour=target)) + geom_line(size=1) + xlim(0,1) + ylim(0,1) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + geom_abline() +facet_grid(category~.)
ggsave(paste("Figs/kim_cohort/DL_overviewRocCombined_uszNeg.pdf",sep=""), width=4, height=4)

ggplot(plotTbl,aes(y=tpr,x=fpr,group=target,  colour=target)) + geom_line(size=1) + xlim(0,1) + ylim(0,1) + xlab("False Positive Rate") + ylab("True Positive Rate") + ylim(0,0.5) + xlim(0,0.05) + theme_bw() + geom_abline() +facet_grid(category~.)
ggsave(paste("Figs/kim_cohort/DL_overviewRocCombined_uszNeg_cropped.pdf",sep=""), width=4, height=4)

##########################
##########################
## compare rj to dl results.
##########################
##########################
require(data.table)
####################################################################
require(data.table)
require(ggplot2)
system("mkdir Figs/kim_cohort/")
system("mkdir Results/kim_cohort/")
system("mkdir interimData/kim_cohort/")

source("Code/prepareDataHelper.R")
source("Code/helperFunctions.R")


myf=function(subtbl,categorystr){
	setkey(subtbl,valueUsed)
	rr=lapply(unique(completeTbl[,target]),function(mytarget){
		xx1=getRocValsInternal(subtbl[target==mytarget,hasCovid+0], subtbl[target==mytarget,valueUsed])
		xx1[,target:=mytarget]
		return(xx1)
	})
	rr2=do.call("rbind",rr)
	rr2[,category:=categorystr]
	return(rr2)
}


xx=myf(completeTbl[(timepoint>6 & timepoint<14)| hasCovid==FALSE,],"timepoint>6 & timepoint<14")
vv=fread("Results/RJ/table_min_mean.csv")
vv[is.na(hasCovid),hasCovid:=0]
vv2=vv[!(hasCovid==0 & cohort=="KIM"),]
vv2[,valueUsed:=data]
xx2=myf(vv2[(timepoint>6 & timepoint<14)| hasCovid==0,],"timepoint>6 & timepoint<14")

ggplot(xx,aes(y=tpr,x=fpr,group=target,  colour=target)) + geom_line(size=1) + xlim(0,1) + ylim(0,1) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + geom_abline() +facet_grid(category~.)
ggsave("blub.pdf")
ggplot(xx2,aes(y=tpr,x=fpr,group=target,  colour=target)) + geom_line(size=1) + xlim(0,1) + ylim(0,1) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + geom_abline() +facet_grid(category~.)
ggsave("blub2.pdf")

