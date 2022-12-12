####################################################################
## Request by Marc as of 4th of May: evaluate plate calibration:
####################################################################
require(data.table)
require(ggplot2)
system("mkdir Figs/evaluateByPlate/")

source("Code/prepareDataHelper.R")

filepaths=do.call("c",lapply(c("RBD","NC","Spike"), getInHouseFitPaths))
tbl=do.call("rbind",lapply(filepaths,fread))
tbl[,neglog:=-log10(IC50)]
xx=tbl[,length(KD),by=Plate_ID]
low_numbers_in_plate=xx[V1< 10,Plate_ID]
print("plates with low numbers")
print(low_numbers_in_plate)
tbl=tbl[!is.element(Plate_ID,low_numbers_in_plate),]

negTbl=tbl[Patient_or_Control_ID=="neg",list(nc_l50=neglog,Plate_ID, target)]
posTbl=tbl[Patient_or_Control_ID=="Positive_serum",list(pc_l50=neglog,Plate_ID, target)]

sampleTbl=extractDate(tbl)
ggplot(sampleTbl, aes(x=neglog)) + geom_histogram()
ggsave("Figs/evaluateByPlate/neglogHist.pdf")
print("only analyse values above -5")

neglogCutoff = -5



lowNrRatioTbl=sampleTbl[,list(lowVals=mean(neglog < neglogCutoff)),by=list(Plate_ID,target)]
collapsedTbl=sampleTbl[neglog> neglogCutoff,list(ml50=mean(neglog),sdl50=sd(neglog),medl50=median(neglog)), by=list(Plate_ID,target)]

setkeyv(collapsedTbl, c("Plate_ID", "target"))
setkeyv(negTbl, c("Plate_ID", "target"))
setkeyv(posTbl, c("Plate_ID", "target"))
setkeyv(lowNrRatioTbl, c("Plate_ID", "target"))

mTbl=merge(negTbl, posTbl)
mTbl=merge(collapsedTbl, mTbl)
mTbl=merge(lowNrRatioTbl, mTbl)

mTbl[,summary(lm(medl50~pc_l50))]
mTbl[,predval:=predict(lm(medl50~pc_l50))]

ggplot(mTbl,aes(y=pc_l50,colour = target, x= Plate_ID)) + geom_point() 
ggsave("Figs/evaluateByPlate/posOverTime.pdf")
ggplot(mTbl,aes(y=nc_l50,colour = target, x= Plate_ID)) + geom_point()
ggsave("Figs/evaluateByPlate/negOverTime.pdf")


ggplot(mTbl, aes(x=medl50,y=lowVals,colour = target)) + geom_point() + ylab("fraction < -5") + xlab("-logEC50 plate median(only above -5)") + geom_smooth(method="lm")
ggsave("Figs/evaluateByPlate/factionLowValsvsMedIC50.pdf")


ggplot(mTbl, aes(x=medl50,y=pc_l50)) + geom_point() + ylab("-logEC50 Positive C") + xlab("-logEC50 plate median(only above -5)") + geom_smooth(method="lm")
ggsave("Figs/evaluateByPlate/posCvsMedIC50.pdf")

ggplot(mTbl, aes(x=medl50,y=pc_l50, colour=target)) + geom_point() + ylab("-logEC50 Positive C") + xlab("-logEC50 plate median(only above -5)") + geom_smooth(method="lm")
ggsave("Figs/evaluateByPlate/posCvsMedIC50_perTarget.pdf")


ggplot(mTbl[lowVals<0.1,], aes(x=medl50,y=pc_l50, colour=target)) + geom_point() + ylab("-logEC50 Positive C") + xlab("-logEC50 plate median(only above -5)") + geom_smooth(method="lm")
ggsave("Figs/evaluateByPlate/posCvsMedIC50_perTarget_fewLowVals.pdf")
####################################################################
####################################################################


####################################################################
####################################################################

