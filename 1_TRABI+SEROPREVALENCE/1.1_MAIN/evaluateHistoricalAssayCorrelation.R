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


#GETBLOOD=TRUE

for(GETBLOOD in c(FALSE,TRUE)){

print("loading all plates, both bsd and usz")
filepaths=do.call("c",lapply(c("NC","Spike","RBD"), getInHouseFitPaths))
tbl=do.call("rbind",lapply(filepaths,fread))
tbl[,neglog:=-log10(IC50)]
xx=tbl[,length(KD),by=Plate_ID]
low_numbers_in_plate=xx[V1< 10,Plate_ID]
print("nr of plates with low numbers")
print(length(low_numbers_in_plate))
tbl=tbl[!is.element(Plate_ID,low_numbers_in_plate),]
print("nr of plates left")
print(length(unique(tbl[,Plate_ID])))

if(GETBLOOD){
	sampleTbl=tbl[grepl("^(H|BDS)",Patient_or_Control_ID),]
	sampleTbl[,editedId:=Patient_or_Control_ID]
	sampleTbl=extractDateBDS(sampleTbl)
	xx=fread("Data/REKO.csv")[,`Patient_or_Control_ID`]
	xx2=sub("(.{13}).+","\\1",xx)
	sampleTbl[is.element(sub("(.{13}).+","\\1",Patient_or_Control_ID),xx2),mydate:=-1]
	}else{
	print("extracting sample dates (only J samples left)")
	sampleTbl=extractDate(tbl)	
	print("nr of plates left")
	print(length(unique(sampleTbl[,Plate_ID])))
	print("nr of Spike samples left")
	print(length(unique(sampleTbl[target=="Spike_IgG",editedId])))
}
blubTbl=sampleTbl[neglog> -5,]
medsTbl=blubTbl[,list(med=median(neglog),iqr=IQR(neglog),len=length(neglog)),by=Plate_ID]
##
##
### remove qc tables
qcTbl=fread("Data/SARS_CoV_2_share/QC/Overview_data_transfer_CoV2.csv",header=TRUE)
dim(sampleTbl[!is.element(Plate_ID,qcTbl[,`Destination plate`])])## should be empty

sampleTbl=sampleTbl[is.element(Plate_ID,qcTbl[QCUser==1,`Destination plate`])]
print("removing plates not passing QC.")
print("nr of plates left")
print(length(unique(sampleTbl[,Plate_ID])))
print("nr of Spike samples left")
print(length(unique(sampleTbl[target=="Spike_IgG",editedId])))

sampleTbl[,toflag:=FALSE]
sampleTbl[neglog <= -5,toflag:=TRUE]
setkey(medsTbl,Plate_ID)
setkey(sampleTbl,Plate_ID)
sampleTbl=merge(sampleTbl,medsTbl)
sampleTbl=sampleTbl[!duplicated(paste(editedId,target)),]

##
##
sampleTbl[,neglogNormed:=(neglog-med)/iqr]
#sampleTbl[,neglogNormed:=(neglog-med)]
sampleTbl[toflag==TRUE,neglogNormed:=0]

hist(sampleTbl[,median(neglogNormed),by=Plate_ID][,V1])
# sampleTbl[,neglogNormed:=(neglog-med)]
# sampleTbl[toflag==TRUE,neglogNormed:=0]
#tbl[grepl("^J",Patient_or_Control_ID),length(IC50),by=Patient_or_Control_ID][,V1]
#sampleTbl[grepl("^J",Patient_or_Control_ID),length(IC50),by=Patient_or_Control_ID][,V1]

if(GETBLOOD){
	negativeTbl = sampleTbl[grepl("BDS",Patient_or_Control_ID),]
	positiveTbl = sampleTbl[grepl("H",Patient_or_Control_ID),]
	newTbl=rbind(negativeTbl,positiveTbl)
	newTbl[,myy:=grepl("H",Patient_or_Control_ID)]
}else{
	negativeTbl = sampleTbl[mydate<190000,]
	print("use before 19 samples for negativeTbl")
	print("nr of plates left")
	print(length(unique(negativeTbl[,Plate_ID])))
	print("nr of Spike samples left")
	print(length(unique(negativeTbl[target=="Spike_IgG",editedId])))
	positiveTbl = sampleTbl[mydate>200400,]
	newTbl=rbind(negativeTbl,positiveTbl)
	newTbl[,myy:=(mydate>200400)+0]
}
##################################################################
##################################################################



numAssaysPerSample=newTbl[,length(med),by=editedId]
allAssays=numAssaysPerSample[V1==3,editedId]
newTbl[is.element(Patient_or_Control_ID,allAssays),]

castTblNormed=data.table(dcast(newTbl,editedId + myy ~target,value.var="neglogNormed"))
castTblUnNormed=data.table(dcast(newTbl,editedId + myy ~target,value.var="neglog"))

castTblNormed[, mypred:=NC_IgG + RBD_IgG + Spike_IgG]
castTblNormed[, mypred2:= RBD_IgG + Spike_IgG]
castTblNormed[, mypred3:=sqrt(abs(RBD_IgG))*sign(RBD_IgG) + sqrt(abs(Spike_IgG))*sign(Spike_IgG)]
castTblNormed[, mypred4:=sqrt(abs(RBD_IgG))*sign(RBD_IgG) + sqrt(abs(Spike_IgG))*sign(Spike_IgG) + sqrt(abs(NC_IgG))*sign(NC_IgG)]


isn=castTblNormed[,!is.na(mypred)]
castTblNormed=castTblNormed[isn,]
castTblUnNormed=castTblUnNormed[isn,]
castTblNormed[,`during Epidemic`:=myy==1]
castTblUnNormed[,`during Epidemic`:=myy==1]



setkey(castTblUnNormed,NC_IgG)
ee1=getRocValsInternal(castTblUnNormed[,myy], castTblUnNormed[,NC_IgG])
ee1[,target:="nc"]
setkey(castTblUnNormed,Spike_IgG)
ee=getRocValsInternal(castTblUnNormed[,myy], castTblUnNormed[,Spike_IgG])
ee[,target:="spike"]
ee1=rbind(ee1,ee)
setkey(castTblUnNormed,RBD_IgG)
ee=getRocValsInternal(castTblUnNormed[,myy], castTblUnNormed[,RBD_IgG])
ee[,target:="RBD"]
ee1=rbind(ee1,ee)
setkey(castTblNormed,mypred)
ee=getRocValsInternal(castTblNormed[,myy], castTblNormed[,mypred])
ee[,target:="normedMeanAll"]
ee1=rbind(ee1,ee)
setkey(castTblNormed,mypred2)
ee=getRocValsInternal(castTblNormed[,myy], castTblNormed[,mypred2])
ee[,target:="normedMeanRBD_Spike"]
ee1=rbind(ee1,ee)
setkey(castTblNormed,mypred3)
ee=getRocValsInternal(castTblNormed[,myy], castTblNormed[,mypred3])
ee[,target:="normedMeanRBD_Spike_sqrt"]
ee1=rbind(ee1,ee)
setkey(castTblNormed,mypred4)
ee=getRocValsInternal(castTblNormed[,myy], castTblNormed[,mypred4])
ee[,target:="normedMeanAll_sqrt"]
ee1=rbind(ee1,ee)
setkey(castTblNormed,Spike_IgG)
ee=getRocValsInternal(castTblNormed[,myy], castTblNormed[,Spike_IgG])
ee[,target:="spike"]



if(GETBLOOD){
	ggplot(ee1,aes(y=tpr,x=fpr,colour=target, group=target)) + geom_line(size=1) + xlim(0,1) + ylim(0,1) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + theme(legend.position = c(0.7 , 0.2))+geom_abline()
	ggsave("Figs/usz_cohort/simpleMeanMetricBlood.pdf", width=4, height=4)
	ggplot(castTblNormed, aes(x=Spike_IgG,y=RBD_IgG,`during Epidemic`,colour=`during Epidemic`)) + geom_point(alpha=0.5,size=1.5) + theme_bw() + theme(legend.position = c(0.7 , 0.2)) + scale_colour_manual(values = c("black", "blue"))
	ggsave("Figs/usz_cohort/corrPlotBlood.pdf", width=4, height=4)
	write.table(newTbl,file="Results/usz_cohort/uscTbllood.t xt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
	write.table(sampleTbl,file="Results/usz_cohort/fullUscTblblood.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
	}else{
	ggplot(ee1,aes(y=tpr,x=fpr,colour=target, group=target)) + geom_line(size=1) + xlim(0,0.04) + ylim(0,0.04) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + theme(legend.position = c(0.7 , 0.2))+geom_abline()
	ggsave("Figs/usz_cohort/simpleMeanMetric.pdf", width=4, height=4)
	ggplot(castTblNormed, aes(x=Spike_IgG,y=RBD_IgG,`during Epidemic`,colour=`during Epidemic`)) + geom_point(alpha=0.5,size=1.5) + theme_bw() + theme(legend.position = c(0.7 , 0.2)) + scale_colour_manual(values = c("black", "blue"))
	ggsave("Figs/usz_cohort/corrPlot.pdf", width=4, height=4)
	write.table(newTbl,file="Results/usz_cohort/uscTbl.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
	write.table(sampleTbl,file="Results/usz_cohort/fullUscTbl.txt",sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
}

}

