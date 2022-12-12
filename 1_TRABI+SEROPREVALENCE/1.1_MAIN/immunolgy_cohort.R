require(data.table)
####################################################################
require(data.table)
require(ggplot2)
system("mkdir Figs/immunology_cohort/")
system("mkdir Results/immunology_cohort/")
system("mkdir interimData/immunology_cohort/")

source("Code/prepareDataHelper.R")
source("Code/helperFunctions.R")

filepaths=do.call("c",lapply(c("RBD","NC","Spike"), getInHouseFitPaths,useImunology=TRUE))
filepaths=filepaths[!grepl("repli",filepaths)]

inhouseTbl=do.call("rbind",lapply(filepaths,fread))
inhouseTbl[,idWithoutMonth:=sub("^000","", sub("\\d\\d$","", Patient_or_Control_ID))]
inhouseTbl[,neglog:=-log10(IC50)]

eurImmuneTbl=fread("Data/SARS_CoV_2_share/20200503_Immunology_dataset/Identifiers/Immunology_cohort_clean.csv")
eurImmuneTbl=eurImmuneTbl[!is.na(ID),]


eurImmuneTbl[,idWithoutMonth:=as.character(ID)]
eurImmuneTbl[,EurImmunCall:=as.numeric(Result)>=1.1]
setkey(inhouseTbl,idWithoutMonth)
setkey(eurImmuneTbl,idWithoutMonth)
mergedTbl=merge(inhouseTbl,eurImmuneTbl)
##########################################################
## have 
##########################################################
mergedTbl[,length(unique(Patient_or_Control_ID)),by=target]
xx=mergedTbl[target=="RBD_IgG",]
print(length(setdiff(unique(inhouseTbl[,Patient_or_Control_ID]),xx[,Patient_or_Control_ID])))
## missing: "114784522003"
print(length(setdiff(unique(eurImmuneTbl[,idWithoutMonth]),xx[,idWithoutMonth])))
##########################################################
##########################################################

ee3=mergedTbl[,mean(-log10(IC50)>2),by=list(EurImmunCall,target)]

ggplot(ee3, aes(y=V1,fill=EurImmunCall,x = target))+geom_bar(stat="identity",position="dodge") + ylab("fraction -log10(EC50)>2") +theme_bw()
ggsave("Figs/immunology_cohort/replicationInEurImmune.pdf")


setkey(mergedTbl,neglog)

mytarget="Spike_IgG"
getRocVals = function(mytarget,tbl){	
	mytbl = tbl[target==mytarget,]
	myrev=rev(mytbl[,EurImmunCall])
	myrevVal = rev(mytbl[,neglog])
	tp = cumsum(myrev)
	myp = sum(myrev)
	tpr = tp/myp
	fp = cumsum(!myrev)
	myn = sum(!myrev)
	fpr = fp/myn
	tn = rev(cumsum(rev(!myrev)))
	tnr = tn/myn
	ppv = tp/ c(1:length(myrev))
	data.table(fpr,tpr, tnr,ppv,val=myrevVal,target=mytarget)
}

getRocVals2=function(mytarget,tbl){
	mytbl = tbl[target==mytarget,]
	outtbl = getRocValsInternal(mytbl[,EurImmunCall],mytbl[,neglog])
	outtbl[,val]
	outtbl[,target:=mytarget]
	return(outtbl)
}


## check manually if it gives the same
resl=lapply(c("RBD_IgG","NC_IgG","Spike_IgG"),getRocVals,tbl=mergedTbl)
resl2=lapply(c("RBD_IgG","NC_IgG","Spike_IgG"),getRocVals2,tbl=mergedTbl)
res=do.call("rbind",resl)
res2=do.call("rbind",resl2)
## should give the same

write.table(res, file ="Results/immunology_cohort/immunologyCohortROCData.txt", quote=FALSE,row.names=FALSE,col.names= TRUE,sep="\t")

ggplot(res,aes(y=tpr,x=fpr,colour=target, group=target))+geom_line(size=1) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + theme(legend.position = c(0.8, 0.2))
ggsave("Figs/immunology_cohort/EurImmunRoc.pdf", width=4, height=4)

ggplot(res,aes(y=tpr,x=fpr,colour=target, group=target))+geom_line(size=1) + xlim(0, 0.05)+ xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + theme(legend.position = c(0.8, 0.2))
ggsave("Figs/immunology_cohort/EurImmunRoc_crop.pdf", width=4, height=4)

#ggplot(res,aes(y=ppv,x=tpr,colour=target, group=target))+geom_line(size=1)
#ggsave("Figs/immunology_cohort/EurImmunPrecisionRecall.pdf")
ggplot(res[val>0,],aes(y=tpr,x=val,colour=target, group=target))+geom_line(size=1)+xlab("cutoff") + ylab("False Positive Rate") + theme_bw() + theme(legend.position = c(0.8, 0.5))
ggsave("Figs/immunology_cohort/tprVsCutoff.pdf", width=4, height=4)

ggplot(res[val>0,],aes(y=fpr,x=val,colour=target, group=target))+geom_line(size=1)+xlab("cutoff") + ylab("False Positive Rate") + theme_bw() + theme(legend.position = c(0.8, 0.5))
ggsave("Figs/immunology_cohort/fprVsCutoff.pdf", width=4, height=4)
ggplot(res[val>0,],aes(y=fpr,x=val,colour=target, group=target))+geom_line(size=1) + xlim(2, 3) +  ylim(0, 0.25) +xlab("cutoff") + ylab("False Positive Rate") + theme_bw() + theme(legend.position = c(0.8, 0.7))
ggsave("Figs/immunology_cohort/fprVsCutoff_crop.pdf", width=4, height=4)

################################################################
################################################################

write.table(mergedTbl, file="interimData/immunology_cohort/annotatedData.txt", quote=FALSE,row.names=FALSE,col.names= TRUE,sep="\t")

