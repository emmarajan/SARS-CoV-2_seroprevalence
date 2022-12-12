
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


fullTbl=fread("interimData/kim_cohort/kim_cohort_inputData.txt")

tbl1=fullTbl[category=="<7 days",]
tbl2=fullTbl[category=="7-13  days",]
tbl3=fullTbl[category==">13 days",]
tbl4=fullTbl[category==">6  days",]

fullTbl[,uniqid:=paste(Pseudo_ID,"_",timepoint,sep="")]

#tt=lapply(unique(tbl3[,target]),function(mytarget){
tt=lapply(unique(tbl3[,target]),function(mytarget){
	myl=list()
	xx1=getRocValsInternal(tbl1[target==mytarget,myy], tbl1[target==mytarget,neglog])
	xx1[,category:="<7 days"]
	xx2=getRocValsInternal(tbl2[target==mytarget,myy], tbl2[target==mytarget,neglog])
	xx2[,category:="7-13  days"]
	xx3=getRocValsInternal(tbl3[target==mytarget,myy], tbl3[target==mytarget,neglog])
	xx3[,category:=">13  days"]
	xx4=getRocValsInternal(tbl4[target==mytarget,myy], tbl4[target==mytarget,neglog])
	xx4[,category:=">6  days"]
	xx=rbind(rbind(rbind(xx1,xx2),xx3),xx4)
	xx[,target:=mytarget]
	return(xx)
})
tt2=do.call("rbind",tt)
write.table(tt2,file="interimData/kim_cohort/kim_cohort_RocData.R",quote=FALSE,row.names=FALSE,col.names=TRUE,sep="\t")



lapply(tt2[,unique(target)],function(mytarget){
	ggplot(tt2[target==mytarget],aes(y=tpr,x=fpr,group=category,  colour=category)) + geom_line(size=1) + xlim(0,1) + ylim(0,1) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + theme(legend.position = c(0.7 , 0.2))+geom_abline()
	ggsave(paste("Figs/kim_cohort/roc_stratifiedpertime_",mytarget,("_uszNeg")[USZNORMED],".pdf",sep=""), width=4, height=4)
	ggplot(tt2[target==mytarget],aes(y=fpr,x=val,group=category,  colour=category)) + geom_line(size=1) + ylim(0,0.05) + xlim(0.5,3) + theme_bw() + theme(legend.position = c(0.7 , 0.2))
	ggsave(paste("Figs/kim_cohort/fprVsVal_stratifiedpertime_",mytarget,("_uszNeg")[USZNORMED],".pdf",sep=""), width=4, height=4)
	ggplot(tt2[target==mytarget],aes(y=tpr,x=fpr,group=category,  colour=category)) + geom_line(size=1) + xlim(0,0.03) + ylim(0,0.5) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + theme(legend.position = c(0.7 , 0.2))+geom_abline()
	ggsave(paste("Figs/kim_cohort/roc_stratifiedpertime_",mytarget,("_uszNeg")[USZNORMED],"_crop.pdf",sep=""), width=4, height=4)

})

ggplot(tt2[!grepl("averaged",target),],aes(y=tpr,x=fpr,group=category,  colour=category)) + geom_line(size=1) + xlim(0,1) + ylim(0,1) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + geom_abline() +facet_grid(target~.)
ggsave(paste("Figs/kim_cohort/overviewRocSingleTargets",("_uszNeg")[USZNORMED],"_crop.pdf",sep=""), width=4, height=4)


ggplot(tt2,aes(y=tpr,x=fpr,group=target,  colour=target)) + geom_line(size=1) + xlim(0,1) + ylim(0,1) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + geom_abline() +facet_grid(category~.)
ggsave(paste("Figs/kim_cohort/overviewRocCombined",("_uszNeg")[USZNORMED],".pdf",sep=""), width=5, height=5)

ggplot(tt2,aes(y=tpr,x=fpr,group=target,  colour=target)) + geom_line(size=1) + xlim(0,1) + ylim(0,1) + xlab("False Positive Rate") + ylab("True Positive Rate") + ylim(0,0.9) + xlim(0,0.02) + theme_bw() + geom_abline() +facet_grid(category~.)
ggsave(paste("Figs/kim_cohort/overviewRocCombined",("_uszNeg")[USZNORMED],"_cropped.pdf",sep=""), width=5, height=5)

fullTbl[,category:=NULL]
castTbl2=data.table(dcast(unique(fullTbl),uniqid + timepoint ~ target,value.var="neglog"))
castTbl2[,myy:=NULL]
castTbl2[,myy:=0]
castTbl2[!is.na(timepoint),myy:=(timepoint>0) + (timepoint>6) + (timepoint>13)]

castTbl2[Spike_IgG < -3,Spike_IgG:=0]
castTbl2[RBD_IgG < -3,RBD_IgG:=0]

ggplot(castTbl2, aes(x=Spike_IgG,y=RBD_IgG,colour=`myy`)) + geom_point(alpha=0.5,size=1.5) + theme_bw() + theme(legend.position = c(0.15 , 0.25)) +geom_abline()
ggsave(paste("Figs/kim_cohort/corrPlot",("_uszNeg")[USZNORMED],".pdf",sep=""), width=4, height=4)

ggplot(castTbl2, aes(x=Spike_IgG,y=RBD_IgG,colour=`myy`)) + geom_point(alpha=0.5,size=1.5) + theme_bw() + theme(legend.position = c(0.15 , 0.25)) + geom_vline(xintercept=1.8,linetype=2,size=0.5) 
ggsave(paste("Figs/kim_cohort/corrPlot",("_uszNeg")[USZNORMED],"single_line.pdf",sep=""), width=4, height=4)

ggplot(castTbl2, aes(x=Spike_IgG,y=RBD_IgG,colour=`myy`)) + geom_point(alpha=0.5,size=1.5) + theme_bw() + theme(legend.position = c(0.15 , 0.25)) + geom_abline(slope=-1, intercept=3.3,linetype=2)
ggsave(paste("Figs/kim_cohort/corrPlot",("_uszNeg")[USZNORMED],"average_line.pdf",sep=""), width=4, height=4)

ggplot(castTbl2, aes(x=Spike_IgG,y=RBD_IgG,colour=`myy`)) + geom_point(alpha=0.5,size=1.5) + theme_bw() + theme(legend.position = c(0.15 , 0.25)) + geom_segment(x=1.5, y=1.5, yend=4,xend=1.5,linetype=2,size=0.5) + geom_segment(x=1.5, y=1.5, yend=1.5,xend=4,linetype=2,size=0.5)
ggsave(paste("Figs/kim_cohort/corrPlot",("_uszNeg")[USZNORMED],"min_line.pdf",sep=""), width=4, height=4)


##########################################################
## have 
##########################################################
print("number of neg samples: Spike")
print(dim(tbl1[target=="Spike_IgG" & myy==0,]))
print("number of pos samples cutoff1: Spike: ")
print(dim(tbl1[target=="Spike_IgG" & myy==1,]))
print("number of pos samples cutoff2: Spike: ")
print(dim(tbl2[target=="Spike_IgG" & myy==1,]))
print("number of pos samples cutoff3: Spike: ")
print(dim(tbl3[target=="Spike_IgG" & myy==1,]))


