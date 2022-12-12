require(data.table)
####################################################################
require(reshape2)
require(ggplot2)
system("mkdir Figs/kim_cohort/")
system("mkdir Results/kim_cohort/")
system("mkdir interimData/kim_cohort/")

source("Code/prepareDataHelper.R")
source("Code/helperFunctions.R")


#commercialTbl=fread("Data/SARS_CoV_2_share/20200506_KIM_reference_dataset/Results_commercial_assays/COVID19_samples_comparison_20200517.csv")


commercialTbl=fread("Data/SARS_CoV_2_share/20200506_KIM_reference_dataset/Results_commercial_assays/Comparison_20200530_complete.csv")
#commercialTbl["COVID19_24"==Pseudo_ID,]


commercialTbl=commercialTbl[,c(1:30),with=FALSE]
commercialTbl[,diasorin:=as.numeric(sub("<","",`DIASORIN Extinktion`))]
commercialTbl[,roche:=as.numeric(`Roche Extinktion`)]
commercialTbl[,oxford:=as.numeric(`Oxford 1:100`)]

commercialTbl[,abbott:=as.numeric(`Abbott [S/Co]`)]
commercialTbl[,euroimmunIgA:=as.numeric(`Euroimmun IgA [S/Co]`)]
commercialTbl[,euroimmunIgG:=as.numeric(`Euroimmun IgG [S/Co]`)]

commercialTbl[grepl(">",`Euroimmun IgA [S/Co]`),euroimmunIgA:=11.7]


commercialTbl=commercialTbl[,list(Pseudo_ID,roche,rochePos=`ROche Bewertung`=="reactive",diasorin,diasorinPos=`DIASORIN Bewertung`=="pos",oxford,abbott,euroimmunIgA,euroimmunIgG)]
#commercialTbl=commercialTbl[!is.na(diasorin) & !is.na(roche),]
#commercialTbl2=commercialTbl[!is.na(diasorin) & !is.na(roche),]
#commercialTbl=commercialTbl[!is.na(diasorin) & !is.na(roche) & !is.na(abbott) & !is.na(euroimmunIgA) & !is.na(euroimmunIgG),]
commercialTbl=commercialTbl[!is.na(diasorin) & !is.na(roche) & !is.na(abbott) & !is.na(euroimmunIgG),]

filepaths=do.call("c",lapply(c("RBD","NC","Spike"), getInHouseFitPaths,useImunology=FALSE,useKim=TRUE))

inhouseTbl=do.call("rbind",lapply(filepaths,fread))
inhouseTbl[,neglog:=-log10(IC50)]

inhouseTbl[,idWithoutMonth:=sub("^000","", sub("\\d\\d$","", Patient_or_Control_ID))]
inhouseTbl[,neglog:=-log10(IC50)]
inhouseTbl[,Pseudo_ID:=Patient_or_Control_ID]
inhouseTbl=inhouseTbl[,list(Pseudo_ID,neglog,target)]
inhouseTbl[neglog< -3,neglog:= -Inf ]
inhouseTblCast=data.table(dcast(inhouseTbl , Pseudo_ID~target,value.var="neglog"))

kimTbl=fread("Data/SARS_CoV_2_share/20200506_KIM_reference_dataset/COVID19_samples_Dropbox.csv")
kimTbl[,timepoint:=`time point (days since first symptoms)`]
hasAgeButNoTimepoint=kimTbl[,!is.na(Age) & is.na(timepoint)]
kimTbl=kimTbl[!hasAgeButNoTimepoint,]
setkey(inhouseTblCast,Pseudo_ID)
setkey(kimTbl,Pseudo_ID)

mergedTbl=merge(inhouseTblCast,kimTbl)
setkey(commercialTbl,Pseudo_ID)
mergedTbl=merge(mergedTbl,commercialTbl)

mergedTbl[is.na(timepoint),timepoint:=-1]
removed=mergedTbl[,NC_IgG== -Inf | RBD_IgG== -Inf  | Spike_IgG== -Inf]
castTblUnused=mergedTbl[removed,]
castTblUsed=mergedTbl[!removed,]

#fit=estimateMVNFrac(as.matrix(castTblUsed[,list(NC_IgG,RBD_IgG,`Spike_IgG`)]),castTblUsed[,timepoint>13],castTblUsed[,timepoint==-1], nround=50)
knownNeg=castTblUsed[,timepoint==-1]
knownNegW =which(knownNeg)
knownPos=castTblUsed[,timepoint>13]
knownPosW=which(knownPos)
#########
posIndexVec=sample(rep(rep(c(1:10)),ceiling(length(knownPosW)/10))[1:length(knownPosW)])
negIndexVec=sample(rep(rep(c(1:10)),ceiling(length(knownNegW)/10))[1:length(knownNegW)])
mymat=as.matrix(castTblUsed[,list(NC_IgG,RBD_IgG,Spike_IgG)])
allVals=rep(NA,nrow(castTblUsed))
for(myi in c(1:10)){
	removedPos=knownPosW[posIndexVec==myi]
	removedNeg=knownNegW[negIndexVec==myi]
	torm=c(removedPos,removedNeg)
	mymat2=mymat[-torm,]
	fit=QDA(mymat2,knownPos[-torm],knownNeg[-torm])
	myposPosterior=getQDAFromFit(valueMat=mymat[removedPos,],fit)
	mynegPosterior=getQDAFromFit(valueMat=mymat[removedNeg,],fit)
	allVals[removedPos]=myposPosterior
	allVals[removedNeg]=mynegPosterior
}
fullfit=QDA(mymat,knownPos,knownNeg)
unknown=!(knownPos | knownNeg)
allVals[unknown]=fullfit[[1]][unknown]


castTblUsed[,posteriorProb:=allVals]
castTblUsed[,averagedSpikeRBD:=(Spike_IgG+RBD_IgG)/sqrt(2)]
castTblUnused[,posteriorProb:=-Inf]
castTblUnused[,averagedSpikeRBD:=(Spike_IgG+RBD_IgG)/sqrt(2)]

mergedTblSub=rbind(castTblUsed,castTblUnused)[,list(Pseudo_ID,timepoint,NC_IgG,RBD_IgG,Spike_IgG,roche,diasorin,oxford,posteriorProb,averagedSpikeRBD,abbott,euroimmunIgG)]


melted=data.table(melt(mergedTblSub, id.vars=c("Pseudo_ID","timepoint")))
setnames(melted, c("Pseudo_ID","timepoint","target","myvalue"))
negTbl=melted[timepoint==-1,]
negTbl[,myy:=FALSE]

pos1Tbl=melted[timepoint<7,]
pos1Tbl[,myy:=TRUE]
pos2Tbl=melted[timepoint>6 & timepoint<14,]
pos2Tbl[,myy:=TRUE]
pos3Tbl=melted[timepoint>13 ,]
pos3Tbl[,myy:=TRUE]
pos4Tbl=melted[timepoint>6 ,]
pos4Tbl[,myy:=TRUE]




tbl1=rbind(negTbl,pos1Tbl)
setkey(tbl1,myvalue)
tbl2=rbind(negTbl,pos2Tbl)
setkey(tbl2,myvalue)
tbl3=rbind(negTbl,pos3Tbl)
setkey(tbl3,myvalue)
tbl4=rbind(negTbl,pos4Tbl)
setkey(tbl4,myvalue)


tbl1[,category:="<7 days"]
tbl2[,category:="7-13 days"]
tbl3[,category:=">13 days"]
tbl4[,category:=">6 days"]

#fullTbl=rbind(rbind(rbind(rbind(negTbl,pos1Tbl),pos2Tbl),pos3Tbl),pos4Tbl)





fullTbl=rbind(rbind(rbind(tbl1,tbl2),tbl3),tbl4)
fullTbl[,timepoint:=NULL]
write.table(fullTbl,file="interimData/kim_cohort_WCommercial_inputData.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
#mytbl= melt(fullTbl, id.vars = c("Pseudo_ID", "myy","category"))
#setnames(mytbl,c("Pseudo_ID","myy","category","target","myvalue"))
setkey(fullTbl,myvalue)
tt=lapply(unique(fullTbl[,target]),function(mytarget){
	# subTbl=fullTbl[category=="<7 days",]
	# xx1=getRocValsInternal(subTbl[target==mytarget,myy], subTbl[target==mytarget,myvalue])
	# xx1[,category:="<7 days"]

	subTbl=fullTbl[category=="7-13 days",]
	xx2=getRocValsInternal(subTbl[target==mytarget,myy], subTbl[target==mytarget,myvalue])
	xx2[,category:="7-13 days"]

	subTbl=fullTbl[category==">13 days",]
	xx3=getRocValsInternal(subTbl[target==mytarget,myy], subTbl[target==mytarget,myvalue])
	xx3[,category:=">13 days"]

	# subTbl=fullTbl[category==">6 days",]
	# xx4=getRocValsInternal(subTbl[target==mytarget,myy], subTbl[target==mytarget,myvalue])
	# xx4[,category:=">6 days"]
	# xx=rbind(rbind(rbind(xx1,xx2),xx3),xx4)
	xx=rbind(xx2,xx3)
	xx[,target:=mytarget]
	return(xx)
})
tt3=do.call("rbind",tt)
write.table(tt3, file="interimData/kimComercialROCData.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

ggplot(tt3,aes(y=tpr,x=fpr,group=target,  colour=target)) + geom_line(size=1) + xlim(0,1) + ylim(0,1) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + theme(legend.position = c(0.7 , 0.2))+geom_abline() + facet_grid(category~.)
ggsave("Figs/commercialComparison.pdf")

ggplot(tt3,aes(y=tpr,x=fpr,group=target,  colour=target)) + geom_line(size=1) + ylim(0,1) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + theme(legend.position = c(0.7 , 0.2))+geom_abline() + facet_grid(category~.)

ggplot(tt3[is.element(target,c("NC_IgG","Spike_IgG","RBD_IgG","posteriorProb")),],aes(y=tpr,x=fpr,group=target,  colour=target)) + geom_line(size=1) + ylim(0,1) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + theme(legend.position = c(0.7 , 0.2))+geom_abline() + facet_grid(category~.)
ggsave("Figs/KIMinternalComparison.pdf")
ggsave("Figs/Preprint/Fig1B.pdf")

ggplot(tt3[is.element(target,c("oxford","diasorin","roche","abbott","euroimmunIgG","posteriorProb","Spike_IgG")),],aes(y=tpr,x=fpr,group=target,  colour=target)) + geom_line(size=1) + ylim(0,1) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + theme(legend.position = c(0.7 , 0.2))+geom_abline() + facet_grid(category~.)
ggsave("Figs/KIMcommercialComparison.pdf")
ggsave("Figs/Preprint/Fig1C.pdf")

