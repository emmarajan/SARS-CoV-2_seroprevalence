require(data.table)
####################################################################
require(data.table)
require(reshape2)
require(ggplot2)

system("mkdir Figs/kim_cohort/")
system("mkdir Results/kim_cohort/")
system("mkdir interimData/kim_cohort/")
system("mkdir interimData/bds_cohort/")
source("Code/helperFunctions.R")


USEBLOOD=FALSE
USENORM=FALSE
USESTRINGENT=TRUE

USET=TRUE
USEMIX=FALSE

bootstrapSample=function(castTbl){
	mysamples=castTbl[,sample(Patient_or_Control_ID,replace=TRUE)]	
	matcher=match(mysamples,castTbl[,Patient_or_Control_ID])
	castTblOut=castTbl[matcher,]
	return(castTblOut)
}


runFunction=function(USET, USEMIX, USEBLOOD, USELDA=FALSE, RUNBOOTSTRAP=FALSE,GETFULLFIT=FALSE){
	if(GETFULLFIT + RUNBOOTSTRAP > 1){
		stop("only one allowed")
	}
	if(USET + USEMIX + USELDA> 1){
		stop("only one allowed(2)")
	}
	if(RUNBOOTSTRAP==FALSE){
		set.seed(11)
	}
	if(USEBLOOD){
		dataTbl=fread("Results/usz_cohort/fulluscTblblood.txt")
		}else{
		dataTbl=fread("Results/usz_cohort/fulluscTbl.txt")
	}

	dataTbl[,Pseudo_ID:=Patient_or_Control_ID]
	dataTbl[neglog< -3,neglog:= -Inf]

	dataTbl[,neglogNormed:=neglog-med]
	if(USENORM){
		valStr="neglogNormed"
	}else{
		valStr="neglog"
	}
	castTbl=as.data.table(dcast(dataTbl , mydate+Patient_or_Control_ID~target,value.var=valStr))
	castTbl=castTbl[!is.na(NC_IgG) & !is.na(RBD_IgG) & !is.na(Spike_IgG),]


	if(RUNBOOTSTRAP){
		castTbl=bootstrapSample(castTbl)
	}

	removed=castTbl[,NC_IgG== -Inf | RBD_IgG== -Inf  | Spike_IgG== -Inf]
	castTblUnused=castTbl[removed,]
	castTblUsed=castTbl[!removed,]
	if(USEBLOOD){
		knownNegNames=castTbl[mydate>-1 & mydate<200200,Patient_or_Control_ID]
		knownPosNames=castTbl[mydate==-1,Patient_or_Control_ID]
	}else{
		 meta=fread("Data/SARS_CoV_2_share/Additional_information/20200520_COVID_Spike_reduced.csv")
		if(USESTRINGENT){
			knownPosNames=meta[`CoV2 time since onset`==">14" & grepl("positive",`CoV2-PCR`),PatientIDList]
			}else{
			knownPosNames=meta[`CoV2-PCR`=="positive",PatientIDList]
		}
		knownNegNames=castTbl[mydate==180000,Patient_or_Control_ID]
	}
	####
	WUnused=castTblUnused[,list(Pseudo_ID=Patient_or_Control_ID,NC_IgG,RBD_IgG,Spike_IgG)]
	WUnused[,avgSpikeRBD:= (RBD_IgG + Spike_IgG)/sqrt(2)]
	WUnused[,knownPos:=is.element(Pseudo_ID,knownPosNames)]
	WUnused[,knownNeg:=is.element(Pseudo_ID,knownNegNames)]
	WUnused[,posteriorProb:=0]
	####
	knownNeg=castTblUsed[,is.element(Patient_or_Control_ID,knownNegNames)]
	knownPos=castTblUsed[,is.element(Patient_or_Control_ID,knownPosNames)]
	####
	knownPosW=which(knownPos)
	knownNegW=which(knownNeg)
	posIndexVec=sample(rep(rep(c(1:10)),ceiling(length(knownPosW)/10))[1:length(knownPosW)])
	negIndexVec=sample(rep(rep(c(1:10)),ceiling(length(knownNegW)/10))[1:length(knownNegW)])
	####
	####
	mymat=as.matrix(castTblUsed[,list(NC_IgG,RBD_IgG,Spike_IgG)])
	allVals=rep(NA,nrow(castTblUsed))
	for(myi in c(1:10)){
		removedPos=knownPosW[posIndexVec==myi]
		removedNeg=knownNegW[negIndexVec==myi]
		torm=c(removedPos,removedNeg)
		mymat2=mymat[-torm,]
	#####
		if(!USEMIX & !USET){
			fit=estimateMVNFrac(mymat2,knownPos[-torm],knownNeg[-torm], nround=50)
			myposPosterior=getPosteriorFromMVNFull(valueMat=mymat[removedPos,],fit)
			mynegPosterior=getPosteriorFromMVNFull(valueMat=mymat[removedNeg,],fit)
		}
		if(USET){
			fit=estimateMVTFrac(mymat2,knownPos[-torm],knownNeg[-torm], nround=50)
			myposPosterior=getPosteriorFromMVTFull(valueMat=mymat[removedPos,],fit)
			mynegPosterior=getPosteriorFromMVTFull(valueMat=mymat[removedNeg,],fit)
		}
		if(USEMIX){
			fit=estimateMVMixedFrac(mymat2,knownPos[-torm],knownNeg[-torm], nround=50)
			myposPosterior=getPosteriorFromVMixedFull(valueMat=mymat[removedPos,],fit)
			mynegPosterior=getPosteriorFromVMixedFull(valueMat=mymat[removedNeg,],fit)
		}
		if(USELDA){
			ldaCorr=TRUE
			fit=estimateLDAFrac(mymat2,knownPos[-torm],knownNeg[-torm], nround=50)
			myposPosterior=getPosteriorFromLDAFull(valueMat=mymat[removedPos,],fit,use1dimCor=ldaCorr)
			mynegPosterior=getPosteriorFromLDAFull(valueMat=mymat[removedNeg,],fit,use1dimCor=ldaCorr)
			if(ldaCorr){
				browser()
				myposPosterior=myposPosterior[removedPos]
				mynegPosterior=mynegPosterior[removedNeg]
			}
		}
		allVals[removedPos]=myposPosterior
		allVals[removedNeg]=mynegPosterior
	}
	###to estimate proportion fit all data: add values to predictor
	if(!USEMIX & !USET & !USELDA){
		fullfit=estimateMVNFrac(mymat,knownPos,knownNeg, nround=50)
	}
	if(USET){
		fullfit=estimateMVTFrac(mymat,knownPos,knownNeg, nround=50)
	}
	if(USEMIX){
		fullfit=estimateMVMixedFrac(mymat,knownPos,knownNeg, nround=50)
		print("degrees of freedom:")
		print(fullfit$paraml$df)
	
	}
	if(USELDA){
		fullfit=estimateLDAFrac(mymat,knownPos,knownNeg, nround=50)
	}
	if(GETFULLFIT){
		return(fullfit)
	}
	unknown=!(knownPos | knownNeg)
	allVals[unknown]=fullfit[[2]][unknown]
	allValsTotFit=fullfit[[2]]
	###
	castTblUsed[,avgSpikeRBD:= (RBD_IgG + Spike_IgG)/sqrt(2)]
	WW=data.table(Pseudo_ID=castTblUsed[,Patient_or_Control_ID],mymat,avgSpikeRBD=castTblUsed[,avgSpikeRBD],knownPos=knownPos,knownNeg=knownNeg,posteriorProb=allVals)
	WTotal=rbind(WW,WUnused)
	mymatcher=match(WTotal[,Pseudo_ID],castTbl[,Patient_or_Control_ID])
	WTotal[,mydate:=castTbl[mymatcher,mydate]]
	WTotal[,mymonth:=floor(mydate/100)]
	prevTbl=getPrevTbl(WTotal,USEBLOOD=USEBLOOD)
	if(RUNBOOTSTRAP){
		return(prevTbl)
	}
	write.table(WTotal,file=paste("interimData/bds_cohort/cohort_combined_outTbl",("_USELDA_")[USELDA],("_USEDMIX_")[USEMIX],("_USEDT_")[USET],("_Blood")[USEBLOOD],("_Normed")[USENORM],("_Above14")[USESTRINGENT & !USEBLOOD],".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
	ggplot(prevTbl,aes(x=month,y=prevalence,fill=knownPos))+geom_bar(stat="identity") + theme_bw()
	ggsave(paste("Figs/cohort_combined_Prevalence",("_USELDA_")[USELDA],("_USEDMIX_")[USEMIX],("_USEDT_")[USET],("_Blood")[USEBLOOD],("_Normed")[USENORM],("_Above14")[USESTRINGENT & !USEBLOOD],".pdf",sep=""))
	WTotalPosNeg=WTotal[knownPos | knownNeg,]
	WTotalPosNeg[,myy:=knownPos+0]
	WTotalPosNeg[,knownPos:=NULL]
	WTotalPosNeg[,knownNeg:=NULL]
	WTotalPosNeg[,mydate:=NULL]
	WTotalPosNeg[,mymonth:=NULL]
	mytbl= melt(WTotalPosNeg, id.vars = c("Pseudo_ID", "myy"))
	setkey(mytbl,value)
	tt=lapply(unique(mytbl[,variable]),function(mytarget){
		xx=getRocValsInternal(mytbl[variable==mytarget,myy], mytbl[variable==mytarget,value])
		xx[,target:=mytarget]
		return(xx)
	})
	tbl=do.call("rbind",tt)
	write.table(tbl,file=paste("interimData/bds_cohort/cohort_combined_RocData",("_USELDA_")[USELDA],("_USEDMIX_")[USEMIX],("_USEDT_")[USET],("_Blood")[USEBLOOD],("_Normed")[USENORM],("_Above14")[USESTRINGENT & !USEBLOOD],".txt",sep=""),quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

	ggplot(tbl, aes(x=fpr,y=tpr, group=target,colour=target)) +geom_line() +theme_bw() + xlim(0,0.01) + theme_bw()
	ggsave(paste("Figs/cohort_combined_RocData",("_USELDA_")[USELDA],("_USEDMIX_")[USEMIX],("_USEDT_")[USET],("_Blood")[USEBLOOD],("_Normed")[USENORM],("_Above14")[USESTRINGENT & !USEBLOOD],".pdf",sep=""))
	return(WTotal)
}


getPrevTbl=function(WTotal,USEBLOOD){
	if(USEBLOOD){
		prevTbl=WTotal[mydate!=-1 & knownNeg==FALSE,list(mylen=length(posteriorProb),nrs=sum(posteriorProb)),by=c("mymonth","knownPos")]
	}else{
		prevTbl=WTotal[knownNeg==FALSE,list(mylen=length(posteriorProb),nrs=sum(posteriorProb)),by=c("mymonth","knownPos")]

	}
	fullLength=prevTbl[,list(fulllen=sum(mylen)),by=mymonth]
	setkey(fullLength,mymonth)
	setkey(prevTbl,mymonth)
	prevTbl=merge(prevTbl,fullLength)
	prevTbl[,prevalence:=(nrs/fulllen)]
	prevTbl[,month:=as.factor(mymonth)]
	return(prevTbl)
}


for(USELDA in c(TRUE,FALSE)){
	for(USET in c(TRUE,FALSE)){
		for(USEMIX in c(TRUE,FALSE)){
			for(USEBLOOD in c(TRUE,FALSE)){	

			set.seed(11)

			if(USET + USEMIX + USELDA > 1){
				next
			}
				runFunction(USET, USEMIX, USEBLOOD)
			}
		}
	}
}



fullfitUsz=runFunction(USET=FALSE, USEMIX=FALSE, USEBLOOD=FALSE,GETFULLFIT=TRUE)
fullfitBlood=runFunction(USET=FALSE, USEMIX=FALSE, USEBLOOD=TRUE,GETFULLFIT=TRUE)

save(fullfitUsz,file="interimData/fullfitUsz.RData")
save(fullfitBlood,file="interimData/fullfitBlood.RData")



aa1=runFunction(USET=FALSE, USEMIX=FALSE, USEBLOOD=FALSE)
aa2=runFunction(USET=FALSE, USEMIX=TRUE, USEBLOOD=FALSE)
aa3=runFunction(USET=FALSE, USEMIX=FALSE, USEBLOOD=FALSE,USELDA=TRUE)

bb1=runFunction(USET=FALSE, USEMIX=FALSE, USEBLOOD=TRUE)
bb2=runFunction(USET=FALSE, USEMIX=TRUE, USEBLOOD=TRUE)
bb3=runFunction(USET=FALSE, USEMIX=FALSE, USEBLOOD=TRUE,USELDA=TRUE)


set.seed(11)
myl=lapply(c(1:1000),function(i){
	runFunction(USET=FALSE, USEMIX=FALSE, USEBLOOD=FALSE, RUNBOOTSTRAP=TRUE)
})
save(myl, file="interimData/bootstrappingUSZ.RData")


set.seed(11)
myl=lapply(c(1:1000),function(i){
	runFunction(USET=FALSE, USEMIX=TRUE, USEBLOOD=FALSE, RUNBOOTSTRAP=TRUE)
})
save(myl, file="interimData/bootstrappingUSZ_MIX.RData")



set.seed(11)
myl=lapply(c(1:1000),function(i){
	runFunction(USET=FALSE, USEMIX=FALSE, USEBLOOD=TRUE, RUNBOOTSTRAP=TRUE)
})
save(myl, file="interimData/bootstrappingBlood.RData")


set.seed(11)
myl=lapply(c(1:1000),function(i){
	runFunction(USET=FALSE, USEMIX=TRUE, USEBLOOD=TRUE, RUNBOOTSTRAP=TRUE)
})
save(myl, file="interimData/bootstrappingBlood_MIX.RData")



ggl=list()

load("interimData/bootstrappingUSZ.RData")
gg=lapply(myl,function(x){
	resl=x[,sum(prevalence),by=month]
})
gg=do.call("rbind",gg)
print("USZ gauss")


ggl[[1]]=gg[,list(cohort="USZ gaussian",lower025=quantile(V1,0.025),upper975=quantile(V1,0.975)),by=month]

load("interimData/bootstrappingUSZ_MIX.RData")
gg=lapply(myl,function(x){
	resl=x[,sum(prevalence),by=month]
})
gg=do.call("rbind",gg)
print("USZ mix")
ggl[[2]]=gg[,list(cohort="USZ mix",lower025=quantile(V1,0.025),upper975=quantile(V1,0.975)),by=month]

load("interimData/bootstrappingBlood.RData")
gg=do.call("rbind",myl)
print("BDS gauss")
quantile(gg[mymonth==2100,prevalence],c(0.025,0.975))
ggl[[3]]=gg[,list(cohort="BDS gauss",lower025=quantile(prevalence,0.025),upper975=quantile(prevalence,0.975)),by=month]

load("interimData/bootstrappingBlood_MIX.RData")
gg=do.call("rbind",myl)
print("BDS mix")
quantile(gg[mymonth==2100,prevalence],c(0.025,0.975))
ggl[[4]]=gg[,list(cohort="BDS mix",lower025=quantile(prevalence,0.025),upper975=quantile(prevalence,0.975)),by=month]


ggM=do.call("rbind",ggl)
write.table(ggM,file="interimData/confidenceIntervals.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

tbl1=getPrevTbl(aa1,USEBLOOD=FALSE)
tbl2=getPrevTbl(aa2,USEBLOOD=FALSE)
tbl3=getPrevTbl(aa3,USEBLOOD=FALSE)
tbl1[,method:="gaussian neg"]
tbl2[,method:="t neg"]
tbl3[,method:="lda"]
tbl4=rbind(rbind(tbl1,tbl2),tbl3)
write.table(tbl4, file="interimData/prevalenceTbl_USZ.txt",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")

tbl32=tbl4[,list(prevalence=sum(prevalence)),by=list(month,method)]
ggplot(tbl32,aes(x=month,y=prevalence,fill=method,group=method))+geom_bar(stat="identity",position="dodge") + theme_bw()
ggsave("Figs/PrevalencesBothMethodsUSZ.pdf")

tbl1=getPrevTbl(bb1,USEBLOOD=TRUE)
tbl2=getPrevTbl(bb2,USEBLOOD=TRUE)
tbl3=getPrevTbl(bb3,USEBLOOD=TRUE)
tbl1[,method:="gaussian neg"]
tbl2[,method:="t neg"]
tbl3[,method:="lda"]
tbl5=rbind(rbind(tbl1,tbl2),tbl3)
write.table(tbl5, file="interimData/prevalenceTbl_BDS.txt",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")
tbl42=tbl5[,list(prevalence=sum(prevalence)),by=list(month,method)]
ggplot(tbl42[month!=-1,],aes(x=month,y=prevalence,fill=method,group=method))+geom_bar(stat="identity",position="dodge") + theme_bw()
ggsave("Figs/PrevalencesBothMethodsBDS.pdf")



ggplot(tbl2,aes(x=month,y=prevalence,fill=knownPos))+geom_bar(stat="identity") + theme_bw()


###########
###########
tbl=fread("interimData/bds_cohort/cohort_combined_RocData_Blood.txt")
ggplot(tbl[target!="avgSpikeRBD",], aes(x=fpr,y=tpr, group=target,colour=target)) +geom_line() +theme_bw() + xlim(0,0.01) + theme_bw()
ggsave("Figs/cohort_combined_RocData_Blood_edit.pdf")
ggsave("Figs/Preprint/S1FigF.pdf")


tbl2=fread("interimData/bds_cohort/cohort_combined_RocData_Above14.txt")
ggplot(tbl2[target!="avgSpikeRBD",], aes(x=fpr,y=tpr, group=target,colour=target)) +geom_line() +theme_bw() + xlim(0,0.01) + theme_bw()
ggsave("Figs/cohort_combined_RocData_Above14_edit.pdf")
ggsave("Figs/Preprint/S1FigE.pdf")

tbl=fread("interimData/bds_cohort/cohort_combined_RocData_USELDA__Above14.txt")
ggplot(tbl[target!="avgSpikeRBD",], aes(x=fpr,y=tpr, group=target,colour=target)) +geom_line() +theme_bw() + xlim(0,0.01) + theme_bw()
#ggsave("Figs/cohort_combined_RocData_Blood_edit.pdf")
#ggsave("Figs/Preprint/S1FigF.pdf")
tbl=fread("interimData/bds_cohort/cohort_combined_RocData_USELDA__Blood.txt")
ggplot(tbl[target!="avgSpikeRBD",], aes(x=fpr,y=tpr, group=target,colour=target)) +geom_line() +theme_bw() + xlim(0,0.01) + theme_bw()



#####
#df: BLOOD:[1] 5.670253

#df: USZ:[1] 6.31327
#####
# USEBLOOD=FALSE
# USENORM=FALSE
# USESTRINGENT=TRUE

# USET=TRUE
# USEMIX=FALSE

aa=runFunction(USET=FALSE, USEMIX=FALSE, USEBLOOD=FALSE)
bb=runFunction(USET=TRUE, USEMIX=FALSE, USEBLOOD=FALSE)
cc=runFunction(USET=FALSE, USEMIX=TRUE, USEBLOOD=FALSE)

# setkey(aa,Pseudo_ID)
# setkey(bb,Pseudo_ID)
# setkey(cc,Pseudo_ID)
addgroup=function(bb){
	bb[,mygroup:="unannotated"]
	bb[knownNeg==TRUE,mygroup:="historical samples"]
	bb[knownPos==TRUE,mygroup:="known positives"]
	return(bb)
}

addgroup(cc)
addgroup(bb)
addgroup(aa)
ggplot(bb,aes(x=RBD_IgG,y=Spike_IgG,shape=mygroup,colour=posteriorProb)) + geom_point(alpha=0.5) + theme_bw()

ggplot(aa,aes(x=RBD_IgG,y=Spike_IgG,shape=mygroup,colour=posteriorProb)) + geom_point(alpha=0.5) + theme_bw()

ggplot(cc,aes(x=RBD_IgG,y=Spike_IgG,shape=mygroup,colour=posteriorProb)) + geom_point(alpha=0.5) + theme_bw()


aa2=runFunction(USET=FALSE, USEMIX=FALSE, USEBLOOD=TRUE)
addgroup(aa2)
aa2[,`Posterior Probability`:=posteriorProb]
aa2[,`Category`:=mygroup]
ggplot(aa2,aes(x=RBD_IgG,y=Spike_IgG,shape=Category,colour=`Posterior Probability`)) + geom_point(alpha=0.5,size=2) + geom_point(data=aa2[mygroup=="known positives",], colour="black",shape=3,size=1) + theme_bw() + scale_colour_gradient(low = "lightblue", high = "black")
ggsave("Figs/cohort_corr_Prevalence_Blood.pdf",width=6,height=4)
ggsave("Figs/Preprint/Fig2B.pdf")



aa2=runFunction(USET=FALSE, USEMIX=FALSE, USEBLOOD=FALSE)
addgroup(aa2)
aa2[,`Posterior Probability`:=posteriorProb]
aa2[,`Category`:=mygroup]
ggplot(aa2,aes(x=RBD_IgG,y=Spike_IgG,shape=Category,colour=`Posterior Probability`)) + geom_point(alpha=0.5,size=2) + geom_point(data=aa2[mygroup=="known positives",], colour="black",shape=3,size=1) + theme_bw() + scale_colour_gradient(low = "lightblue", high = "black")
ggsave("Figs/cohort_corr_Prevalence_USZ.pdf",width=6,height=4)
ggsave("Figs/Preprint/Fig2A.pdf")

aa2=runFunction(USET=FALSE, USEMIX=FALSE, USEBLOOD=FALSE,USELDA=TRUE)
addgroup(aa2)
aa2[,`Posterior Probability`:=posteriorProb]
aa2[,`Category`:=mygroup]
ggplot(aa2,aes(x=RBD_IgG,y=Spike_IgG,shape=Category,colour=`Posterior Probability`)) + geom_point(alpha=0.5,size=2) + geom_point(data=aa2[mygroup=="known positives",], colour="black",shape=3,size=1) + theme_bw() + scale_colour_gradient(low = "lightblue", high = "black")
ggsave("Figs/cohort_corr_Prevalence_LDA_USZ.pdf",width=6,height=4)

aa2=runFunction(USET=FALSE, USEMIX=FALSE, USEBLOOD=TRUE,USELDA=TRUE)
addgroup(aa2)
aa2[,`Posterior Probability`:=posteriorProb]
aa2[,`Category`:=mygroup]
ggplot(aa2,aes(x=RBD_IgG,y=Spike_IgG,shape=Category,colour=`Posterior Probability`)) + geom_point(alpha=0.5,size=2) + geom_point(data=aa2[mygroup=="known positives",], colour="black",shape=3,size=1) + theme_bw() + scale_colour_gradient(low = "lightblue", high = "black")
ggsave("Figs/cohort_corr_Prevalence_LDA_Blood.pdf",width=6,height=4)


aa2=runFunction(USET=FALSE, USEMIX=TRUE, USEBLOOD=FALSE)
addgroup(aa2)
aa2[,`Posterior Probability`:=posteriorProb]
aa2[,`Category`:=mygroup]
ggplot(aa2,aes(x=RBD_IgG,y=Spike_IgG,shape=Category,colour=`Posterior Probability`)) + geom_point(alpha=0.5,size=2) + geom_point(data=aa2[mygroup=="known positives",], colour="black",shape=3,size=1) + theme_bw() + scale_colour_gradient(low = "lightblue", high = "black")
ggsave("Figs/cohort_corr_Prevalence_Mix_USZ.pdf",width=6,height=4)
ggsave("Figs/Preprint/S1FigC.pdf")

aa2=runFunction(USET=FALSE, USEMIX=TRUE, USEBLOOD=TRUE)
addgroup(aa2)
aa2[,`Posterior Probability`:=posteriorProb]
aa2[,`Category`:=mygroup]
ggplot(aa2,aes(x=RBD_IgG,y=Spike_IgG,shape=Category,colour=`Posterior Probability`)) + geom_point(alpha=0.5,size=2) + geom_point(data=aa2[mygroup=="known positives",], colour="black",shape=3,size=1) + theme_bw() + scale_colour_gradient(low = "lightblue", high = "black")
ggsave("Figs/cohort_corr_Prevalence_Mix_Blood.pdf",width=6,height=4)
ggsave("Figs/Preprint/S1FigD.pdf")


aa=fread("interimData/prevalenceTbl_USZ.txt")
aa2=aa[mymonth>2002,list(sum(nrs),mean(fulllen)), by=list(month,method)]
aa3=aa2[,sum(V1)/sum(V2),by=method]
aa3[method=="t neg",V1]/aa3[method=="gaussian neg",V1]




aa=fread("interimData/prevalenceTbl_BDS.txt")
aa2=aa[mymonth>2002,list(sum(nrs),mean(fulllen)), by=list(month,method)]
aa3=aa2[,sum(V1)/sum(V2),by=method]
aa3[method=="t neg",V1]/aa3[method=="gaussian neg",V1]

aa3[method=="t neg",V1]/aa3[method=="gaussian neg",V1]
###ggsave("Figs/cohort_corr_Prevalence",("_USEDMIX_")[USEMIX],("_USEDT_")[USET],("_Blood")[USEBLOOD],("_Normed")[USENORM],("_Above14")[USESTRINGENT & !USEBLOOD],".pdf"



fitUszLDA=runFunction(USET=FALSE, USEMIX=FALSE, USEBLOOD=FALSE,USELDA=TRUE,GETFULLFIT=TRUE)
fitBloodLDA=runFunction(USET=FALSE, USEMIX=FALSE, USEBLOOD=TRUE,USELDA=TRUE,GETFULLFIT=TRUE)




myfit=fitUszLDA
#myfit=fitBloodLDA
tstat = myfit[[4]]
annot = myfit[[5]]
post = myfit[[2]]
tbl=data.table(tstat,annot,post)
setkey(tbl,V1)
tbl2=tbl[annot==-1,]
tbl2[,emp:=scale(V1)]
#tbl2[,theor:=qt(c(1:nrow(tbl2))/(nrow(tbl2)+1),100)]
tbl2[,theor:=qnorm(c(1:nrow(tbl2))/(nrow(tbl2)+1))]

ggplot(tbl2,aes(x=theor,y=emp))+geom_point()+geom_abline()


plot(tbl[annot==-1,V1]),sum(tbl[,annot==-1]
########################################
# rank(ww)




# ww=scale(sort(castTbl[NC_IgG!=-Inf & mydate<200100,NC_IgG]))

# ww2=rank(ww)/c(length(ww)+1)

# plot(ww,qnorm(ww2))
# abline(0,1)

# plot(ww,qt(ww2, 11))
# abline(0,1)


# ww=scale(sort(castTbl[Spike_IgG!=-Inf & mydate<200100,Spike_IgG]))
# ww2=rank(ww)/c(length(ww)+1)
# #plot(ww,qt(ww2, 15))
# plot(ww,qnorm(ww2))
# abline(0,1)


# ww=scale(sort(castTbl[Spike_IgG!=-Inf & mydate<200100,Spike_IgG]))
# ww2=rank(ww)/c(length(ww)+1)
# plot(ww,qt(ww2,))
# abline(0,1)

