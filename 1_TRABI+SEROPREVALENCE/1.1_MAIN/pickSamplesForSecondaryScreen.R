require(data.table)
source("Code/prepareDataHelper.R")


uszBloodTbl0=fread("Results/usz_cohort/fullUscTblblood.txt")
uszBloodTbl0=uszBloodTbl0[,list(Plate_ID,Patient_or_Control_ID,mydate,target,neglog)]
uszBloodTbl0=uszBloodTbl0[mydate!=-1,]
uszBloodTbl0[,cohort:="BDS"]
uszBloodTbl=uszBloodTbl0[target!="NC_IgG",list(Plate_ID=Plate_ID[1],Patient_or_Control_ID=Patient_or_Control_ID[1],mydate=mydate[1],target="averagedSpikeRBD",neglog=(mean(neglog)*sqrt(2)),cohort=cohort[1]),by=Patient_or_Control_ID]

uszTbl0=fread("Results/usz_cohort/fullUscTbl.txt")
uszTbl0=uszTbl0[,list(Plate_ID,Patient_or_Control_ID,mydate,target,neglog)]
uszTbl0[,cohort:="USZ"]
uszTbl=uszTbl0
uszTbl=uszTbl[target!="NC_IgG",list(Plate_ID=Plate_ID[1],Patient_or_Control_ID=Patient_or_Control_ID[1],mydate=mydate[1],target="averagedSpikeRBD",neglog=(mean(neglog)*sqrt(2)),cohort=cohort[1]),by=Patient_or_Control_ID]



bloodTbl=uszBloodTbl[,list(Plate_ID,Patient_or_Control_ID,mydate,target,neglog,cohort)]
uszTbl=uszTbl[,list(Plate_ID,Patient_or_Control_ID,mydate,target,neglog,cohort)]


setkey(uszTbl,neglog)
selectedUszTbl=uszTbl[rev(c(1:nrow(uszTbl))),][c(1:150),]
setkey(bloodTbl,neglog)
selectedBloodTbl=bloodTbl[rev(c(1:nrow(bloodTbl))),][c(1:60),]

selTbl=rbind(selectedUszTbl,selectedBloodTbl)

bothTbl=rbind(uszBloodTbl0,uszTbl0)
uszBothTblDiff=bothTbl[!is.element(Patient_or_Control_ID,selTbl[,Patient_or_Control_ID]),]

setkey(uszBothTblDiff,neglog)

ee=lapply(c("RBD_IgG","NC_IgG","Spike_IgG"),function(mytarget){
	mytbl=uszBothTblDiff[target==mytarget,]
	mytbl[,mypheno:=-neglog]
	setkey(mytbl,mypheno)
	mytbl2=mytbl[1:10,]
	mytbl2[,mypheno:=NULL]
	mytbl2[,cond:=mytarget]
	return(mytbl2)
})
ee2=do.call("rbind",ee)
selTbl[,cond:="averagedSpikeRBD"]

fullSelTbl=rbind(selTbl,ee2)

write.table(fullSelTbl,"Results/selectedPositives.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


########################
########################
filepaths=do.call("c",lapply(c("NC","Spike","RBD"), getInHouseFitPaths, getLayout=TRUE))
tableLayout=do.call("rbind",lapply(filepaths,fread))
tableLayout1=unique(tableLayout[,list(`Source Plate Barcode`,`Source Well`,`patient_or_control_id`)])
setnames(tableLayout1, c("Source Plate Barcode","Source Well","Patient_or_Control_ID"))
setkey(fullSelTbl,Patient_or_Control_ID)
setkey(tableLayout1,Patient_or_Control_ID)
fullSelTblAnnot=merge(fullSelTbl,tableLayout1)
sum(fullSelTblAnnot[,Patient_or_Control_ID]!=fullSelTbl[,Patient_or_Control_ID])
write.table(fullSelTblAnnot,"Results/selectedPositivesAnnot.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

########################
########################
set.seed(13)
uszTblNeg=uszTbl[180000==mydate,]
uszTblNeg=uszTblNeg[!is.element(Patient_or_Control_ID,fullSelTbl[,Patient_or_Control_ID])]
myplates=sample(unique(uszTblNeg[,Plate_ID]),7)
uszTblNegSelPlates=uszTblNeg[is.element(Plate_ID,myplates),]
indices=sample(c(1:nrow(uszTblNegSelPlates)),70)
negUszTbl=uszTblNegSelPlates[indices,]
hist(negUszTbl[neglog>-2,neglog])
#######

bloodTblNeg=bloodTbl[180000==mydate,]
bloodTblNeg=bloodTblNeg[!is.element(Patient_or_Control_ID,fullSelTbl[,Patient_or_Control_ID])]
myplates=sample(unique(bloodTblNeg[,Plate_ID]),5)
bloodTblNegSelPlates=bloodTblNeg[is.element(Plate_ID,myplates),]
indices=sample(c(1:nrow(bloodTblNegSelPlates)),52)
negBloodTbl=bloodTblNegSelPlates[indices,]
hist(negBloodTbl[neglog>-2,neglog])

negBothTbl=rbind(negUszTbl,negBloodTbl)
write.table(negBothTbl,"Results/selectedNegatives.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

########################
########################
setkey(negBothTbl,Patient_or_Control_ID)
negBothTblAnnot=merge(negBothTbl,tableLayout1)
sum(negBothTblAnnot[,Patient_or_Control_ID]!=negBothTbl[,Patient_or_Control_ID])
write.table(negBothTblAnnot,file="Results/selectedNegativesAnnot.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

########################
########################

nc=fread("/Users/dlamparter/Downloads/table/10_highest_NC.csv")
rbd=fread("/Users/dlamparter/Downloads/table/10_highest_RBD.csv")
spike=fread("/Users/dlamparter/Downloads/table/10_highest_Spike.csv")
bds=fread("/Users/dlamparter/Downloads/table/60_highest_BDS_muSpikeRBD.csv")
usz=fread("/Users/dlamparter/Downloads/table/150_highest_usz_muSpikeRBD.csv")
setkey(usz,Patient_or_Control_ID)
setkey(bds,Patient_or_Control_ID)
setkey(nc,Patient_or_Control_ID)
setkey(spike,Patient_or_Control_ID)
setkey(rbd,Patient_or_Control_ID)

setkey(fullSelTblAnnot,Patient_or_Control_ID)

m3=merge(usz,fullSelTblAnnot)
dim(m3)
sum(m3[,`Source Well.y`!=`Source Well.x`])
sum(m3[,`Source Plate Barcode.x`!=`Source Plate Barcode.y`])

m3=merge(bds,fullSelTblAnnot)
dim(m3)
sum(m3[,`Source Well.y`!=`Source Well.x`])
sum(m3[,`Source Plate Barcode.x`!=`Source Plate Barcode.y`])

m3=merge(spike,fullSelTblAnnot)
dim(m3)
sum(m3[,`Source Well.y`!=`Source Well.x`])
sum(m3[,`Source Plate Barcode.x`!=`Source Plate Barcode.y`])

m3=merge(rbd,fullSelTblAnnot)
dim(m3)
sum(m3[,`Source Well.y`!=`Source Well.x`])
sum(m3[,`Source Plate Barcode.x`!=`Source Plate Barcode.y`])

m3=merge(nc,fullSelTblAnnot)
dim(m3)
sum(m3[,`Source Well.y`!=`Source Well.x`])
sum(m3[,`Source Plate Barcode.x`!=`Source Plate Barcode.y`])


