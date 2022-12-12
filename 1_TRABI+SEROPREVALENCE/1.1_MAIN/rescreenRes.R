##################################################################
##################################################################
## investigate assay correlation among negatives..

require(data.table)
require(reshape2)
require(ggplot2)
source("Code/prepareDataHelper.R")
source("Code/helperFunctions.R")

system("mkdir Figs/usz_cohort/")
system("mkdir Figs/Preprint/")
system("mkdir Results/usz_cohort/")
system("mkdir interimData/usz_cohort/")

set.seed(11)

bloodMeasured=fread("interimData/bds_cohort/cohort_combined_outTbl_Blood.txt")
uszMeasured=fread("interimData/bds_cohort/cohort_combined_outTbl_Above14.txt")

bloodMeasured2=bloodMeasured[,list(Patient_or_Control_ID=Pseudo_ID,posteriorProb=posteriorProb)]
uszMeasured2=uszMeasured[,list(Patient_or_Control_ID=Pseudo_ID,posteriorProb=posteriorProb)]

setkey(bloodMeasured2,Patient_or_Control_ID)
setkey(uszMeasured2,Patient_or_Control_ID)


filepaths2=getReplicationFitPaths(nr1=FALSE)
tbl2=do.call("rbind",lapply(filepaths2,fread))
tbl2[,replication:=2]
#GETBLOOD=TRUE
filepaths1=getReplicationFitPaths(nr1=TRUE)
tbl1=do.call("rbind",lapply(filepaths1,fread))
tbl1[,replication:=1]
tbl=rbind(tbl1,tbl2)
tbl[,neglog:=-log10(IC50)]
tbl[neglog< -3,neglog:=-Inf]
#####
tblReplication=tbl
qcTbl=fread("Data/SARS_CoV_2_share/QC/Overview_data_transfer_CoV2.csv",header=TRUE)
dim(tbl[!is.element(Plate_ID,qcTbl[,`Destination plate`])])## should be empty
tblReplication=tbl[is.element(Plate_ID,qcTbl[QCUser==1,`Destination plate`])]
#####
tblReplication=tblReplication[grep("^(J|BDS|H)",Patient_or_Control_ID)]
tblReplication=tblReplication[,list(Patient_or_Control_ID,target,neglog,replication)]
#####
tblReplication[,targetRepl:=paste(target,"_",replication,sep="")]

castTbl=data.table(dcast(tblReplication,  Patient_or_Control_ID + target ~ replication,value.var="neglog"))
castTbl[`1`==-Inf,`1`:=0]
castTbl[`2`==-Inf,`2`:=0]
ggplot(castTbl,aes(x=`1`,y=`2`))+geom_point() + geom_abline() + facet_grid(target~.) + ylab("-log(EC50) replicate 2") + xlab("-log(EC50) replicate 1") + theme_bw() 
ggsave("Figs/replication/replicationReplication.pdf",width=6,height=6)
ggsave("Figs/Preprint/S1FigG.pdf")
castTbl[,`On Baseline`:="no"]
castTbl[`1`==0 | `2`==0,`On Baseline`:="yes"]

ggplot(castTbl,aes(x=`1`,y=`2`,colour=`On Baseline`))+geom_point() + geom_abline() + facet_grid(target~.) + ylab("-log(EC50) replicate 2") + xlab("-log(EC50) replicate 1")  + theme_bw() + theme(legend.position = c(0.8, 0.2))
ggsave("Figs/replication/replicationReplication.pdf",width=5,height=12)
ggsave("Figs/Preprint/Fig2D_v2.pdf",width=4,height=12)

ggplot(castTbl,aes(x=`1`,y=`2`))+geom_point() + geom_abline() + facet_grid(target~.) + ylab("-log(EC50) replicate 2") + xlab("-log(EC50) replicate 1")  + theme_bw() + theme(legend.position = c(0.8, 0.2))
ggsave("Figs/Preprint/Fig2D_v1.pdf",width=4,height=12)

WW=data.table(castTbl)
WW2=WW[target=="Spike_IgG" & `1`> -Inf & `2`> -Inf,summary(lm(`1`~`2`))]
print(WW2$r.squared)




castTbl2=data.table(dcast(tblReplication,  Patient_or_Control_ID ~ replication + target,value.var="neglog"))


bloodTbl=fread("Results/usz_cohort/fullUscTblblood.txt")
bloodTbl=bloodTbl[,list(Patient_or_Control_ID,target,neglog,mydate,replication=0,targetRepl=paste(target,0,sep="_"))]
bloodTbl[neglog< -3,neglog:=-Inf]
bloodTbl[,firstValue:=neglog]
setkeyv(bloodTbl,c("Patient_or_Control_ID","target"))


uszTbl=fread("Results/usz_cohort/fullUscTbl.txt")
uszTbl=uszTbl[,list(Patient_or_Control_ID,target,neglog,mydate,replication=0,targetRepl=paste(target,0,sep="_"))]
uszTbl[neglog< -3,neglog:=-Inf]
uszTbl[,firstValue:=neglog]
setkeyv(uszTbl,c("Patient_or_Control_ID","target"))


ggplot(castTbl2,aes(x=`1_RBD_SARS_IgG`,y=`1_RBD_IgG`))+geom_point() + geom_abline() + theme_bw()
ggsave("Figs/replication/repl_SARSRBDvsRBD.pdf")
ggplot(castTbl2,aes(x=`1_RBD_SARS_IgG`,y=`1_Spike_IgG`))+geom_point() + geom_abline() + theme_bw()
ggsave("Figs/replication/repl_SARSRBDvsSpike.pdf")
ggplot(castTbl2,aes(x=`1_Spike_IgG`,y=`1_RBD_IgG`))+geom_point() + geom_abline() + theme_bw()
ggsave("Figs/replication/repl_RBDvsSpike.pdf")



tblReplication[,replicationValue:=neglog]
repl1=tblReplication[replication==1,]
setkeyv(repl1,c("Patient_or_Control_ID","target"))



aa=fread("Results/selectedPositivesAnnot.txt")
bb=fread("Results/selectedNegativesAnnot.txt")
repl1[,inPos:=is.element(Patient_or_Control_ID,aa[,`Patient_or_Control_ID`])]
repl1[,inNeg:=is.element(Patient_or_Control_ID,bb[,`Patient_or_Control_ID`])]
repl1[,condition:="neg"]

myindex=match(repl1[inPos==TRUE,Patient_or_Control_ID],aa[,`Patient_or_Control_ID`])
repl1[inPos==TRUE,condition:=aa[myindex,target]]

bloodMerged=merge(bloodTbl,repl1)


uszMerged=merge(uszTbl,repl1)
uszMerged2=merge(uszMerged,uszMeasured2)
uszMerged2[,cohort:="USZ"]

bloodMerged=merge(bloodTbl,repl1)
bloodMerged2=merge(bloodMerged,bloodMeasured2)
bloodMerged2[,cohort:="BDS"]

bothMerged2=rbind(uszMerged2,bloodMerged2)


ggplot(bothMerged2[target=="Spike_IgG" & is.element(condition,c("neg","averagedSpikeRBD")),], aes(y=firstValue,x=replicationValue-firstValue,colour=posteriorProb,shape=condition))+geom_point() + geom_vline(xintercept=0) + scale_colour_gradient(low = "lightblue", high = "black") + theme_bw() 
ggsave("Figs/SpikeReplication.pdf",width=4,height=4)
ggplot(bothMerged2[target=="NC_IgG" & is.element(condition,c("neg","averagedSpikeRBD")),], aes(y=firstValue,x=replicationValue-firstValue,colour=posteriorProb,shape=condition))+geom_point() + geom_vline(xintercept=0) + scale_colour_gradient(low = "lightblue", high = "black") + theme_bw() 
ggsave("Figs/NCReplication.pdf",width=4,height=4)
ggplot(bothMerged2[target=="RBD_IgG" & is.element(condition,c("neg","averagedSpikeRBD")),], aes(y=firstValue,x=replicationValue-firstValue,colour=posteriorProb,shape=condition))+geom_point() + geom_vline(xintercept=0) + scale_colour_gradient(low = "lightblue", high = "black") + theme_bw() 
ggsave("Figs/RBDReplication.pdf",width=4,height=4)
#ggsave("Figs/Preprint/Fig2D.pdf",width=5.5,height=4) 



AA=bothMerged2[target=="Spike_IgG" & is.element(condition,c("neg","averagedSpikeRBD")),]




ggplot(bothMerged2[target=="NC_IgG" & is.element(condition,c("neg","averagedSpikeRBD")),], aes(y=firstValue,x=replicationValue-firstValue,colour=posteriorProb,shape=condition))+geom_point() + geom_vline(xintercept=0) + scale_colour_gradient(low = "lightblue", high = "black") + theme_bw() 

ggsave("Figs/NCReplication.pdf",width=4,height=4)  + scale_colour_gradient(low = "lightblue", high = "black") + theme_bw() 

ggplot(bothMerged2[target=="RBD_IgG" & is.element(condition,c("neg","averagedSpikeRBD")),], aes(y=firstValue,x=replicationValue-firstValue,colour=posteriorProb,shape=condition))+geom_point() + geom_vline(xintercept=0) + scale_colour_gradient(low = "lightblue", high = "black") + theme_bw()
ggsave("Figs/RBDReplication.pdf",width=4,height=4) + scale_colour_gradient(low = "lightblue", high = "black") + theme_bw() 




setkey(bloodMeasured2,Patient_or_Control_ID)
setkey(uszMeasured2,Patient_or_Control_ID)
bothMeasured2=rbind(bloodMeasured2,uszMeasured2)


setkey(castTbl2,Patient_or_Control_ID)
setkey(bothMeasured2,Patient_or_Control_ID)
castTbl3=merge(castTbl2,bothMeasured2)
### needs condition
conditionTbl=unique(repl1[,list(Patient_or_Control_ID,condition)])
setkey(conditionTbl,Patient_or_Control_ID)
setkey(castTbl3,Patient_or_Control_ID)
castTbl4=merge(castTbl3,conditionTbl)
castTbl5=castTbl4[is.element(condition,c("averagedSpikeRBD","neg")),]

#castTbl6=castTbl5[,list(rbd=(`1_RBD_IgG`+`2_RBD_IgG`)/2,sars=(`1_RBD_SARS_IgG`+`2_RBD_SARS_IgG`)/2)]
castTbl6=castTbl5[,list(rbd=`1_RBD_IgG`,sars=`1_RBD_SARS_IgG`,posteriorProb=posteriorProb)]
ggplot(castTbl6,aes(y=rbd,x=sars))+geom_point()
castTbl6[rbd< -3,rbd:=-3]
castTbl6[sars< -3,sars:=-3]

castTbl6[,mygroup:="<1.5"]
castTbl6[rbd>=1.5,mygroup:=">1.5 & <2.5"]
castTbl6[rbd>=2.5,mygroup:=">2.5"]
castTbl6[,mygroupF:=factor(mygroup,levels=c("<1.5",">1.5 & <2.5",">2.5"))]


ggplot(castTbl6,aes(x=mygroupF, y=sars))+geom_violin(fill="lightgrey") + geom_jitter(aes(colour=posteriorProb),width=0.1,height=0) + ylab("SARS-RBD")+xlab("RBD level") + scale_colour_gradient(low = "lightblue", high = "black") + theme_bw() 
ggsave("Figs/rbd_vs_sars.pdf",width=4,height=4)
ggsave("Figs/Preprint/Fig2E.pdf")

#################################
#################################
bloodMeasured=fread("interimData/bds_cohort/cohort_combined_outTbl_Blood.txt")

uszMeasured=fread("interimData/bds_cohort/cohort_combined_outTbl_Above14.txt")


uszMeasured[,isNeg:=mydate==180000]
uszMeasured[,mymonth:=floor(mydate/100)]
uszMeasured2=uszMeasured[mymonth==1800 | mymonth==2005,]

uszMeasured3=uszMeasured2[,list(Pseudo_ID,NC_IgG,RBD_IgG,Spike_IgG,isNeg,posteriorProb)]
uszMeasured4=melt(uszMeasured3, id.vars = c("Pseudo_ID","isNeg","posteriorProb"))
uszMeasured4[,category:="historical"]
uszMeasured4[isNeg==FALSE,category:="May/20"]
uszMeasured4[,category:=factor(category,levels=c("historical","May/20"))]

ggplot(uszMeasured4[value>2,],aes(y=value,x=category,colour=posteriorProb)) + geom_violin(fill="white") + geom_jitter(aes(colour=posteriorProb),size=1,width=0.1,height=0) + facet_grid(.~variable) + scale_colour_gradient(low = "lightblue", high = "black") +  theme_bw()
ggsave("Figs/compareTopVals.pdf",width=5.5,height=4) 
ggsave("Figs/Preprint/Fig4A.pdf",width=5.5,height=4) 
browser()




ggplot(bloodMeasured[value>2 & (knownNeg==TRUE | mydate==210000),],aes(y=value,x=category,colour=posteriorProb)) + geom_violin(fill="white") + geom_jitter(aes(colour=posteriorProb),size=1,width=0.1,height=0) + facet_grid(.~variable) + scale_colour_gradient(low = "lightblue", high = "black") +  theme_bw()
############
uszMeasured4[value>2,list(fractionPosteriorAboveDot5=mean(posteriorProb>0.5)),by=list(category,variable)]

############


bloodMeasured[,isNeg:=knownNeg]
bloodMeasured[,mymonth:=floor(mydate/100)]
bloodMeasured2=bloodMeasured[knownNeg==TRUE | mymonth==2100,]

bloodMeasured3=bloodMeasured2[,list(Pseudo_ID,NC_IgG,RBD_IgG,Spike_IgG,isNeg,posteriorProb)]
bloodMeasured4=melt(bloodMeasured3, id.vars = c("Pseudo_ID","isNeg","posteriorProb"))
bloodMeasured4[,category:="historical"]
bloodMeasured4[isNeg==FALSE,category:="May/20"]
bloodMeasured4[,category:=factor(category,levels=c("historical","May/20"))]

ggplot(bloodMeasured4[value>1.6 & value<2.3,],aes(y=value,x=category,colour=posteriorProb)) + geom_violin(fill="white") + geom_jitter(aes(colour=posteriorProb),size=1,width=0.1,height=0) + facet_grid(.~variable) + scale_colour_gradient(low = "lightblue", high = "black") +  theme_bw()
####ggsave("Figs/compareTopVals.pdf",width=5.5,height=4) 


