require(reshape2)
require(maptree)
require(mgcv)
require(data.table)
require(rpart)
require(earth)

require(ggplot2)

source("Code/helperFunctions.R")

mergedTbl = fread("interimData/immunology_cohort/annotatedData.txt")
mergedTbl[,myy:=EurImmunCall+0]
subtbl=mergedTbl[,list(neglog, idWithoutMonth, myy, EurImmunCall, target)]
castTbl=dcast(subtbl,idWithoutMonth + myy + EurImmunCall~target,value.var="neglog")
castTbl[NC_IgG< -5,NC_IgG:=-0.5]
castTbl[Spike_IgG< -5,Spike_IgG:=-0.5]
castTbl[RBD_IgG< -5,RBD_IgG:=-0.5]
ggplot(castTbl[NC_IgG>-5 & Spike_IgG>-5,],aes(x=RBD_IgG, y=Spike_IgG,colour=myy))+geom_point()



getROCForComposite=function(castTbl, myfunc){
	set.seed(11)
	indexVec=sample(rep(rep(c(1:10)),ceiling(nrow(castTbl)/10))[1:nrow(castTbl)])
	resVec = rep(NA, length(indexVec))
	for (i in c(1:10)){
		sum(castTbl[indexVec==i,myy])
		testDat = castTbl[indexVec==i,]
		trainDat = castTbl[indexVec!=i,]
		tt=myfunc(trainDat, testDat)
		resVec[indexVec==i]=tt
	}
	castTbl[,gamPred:=resVec]
	tableForRoc=castTbl[,list(gamPred,EurImmunCall)]
	setkey(tableForRoc,gamPred)
	myroc=getRocValsInternal(tableForRoc[,EurImmunCall],tableForRoc[,gamPred])
	return(myroc)
}

gamTrainFunc=function(trainDat, testDat){
	gamres = gam(myy ~ s(NC_IgG)+s(RBD_IgG)+s(Spike_IgG),family=binomial(),
     data = trainDat)
	summary(gamres)
	tt=predict.gam(gamres, testDat) 
}

marsTrainFunc=function(trainDat, testDat){
	marsres = earth(myy ~ NC_IgG + RBD_IgG + Spike_IgG,
     data = trainDat)
	summary(marsres)
	tt=predict(marsres, testDat) 
}

treeTrainFunc=function(trainDat, testDat){
	fit = rpart(myy ~ NC_IgG + RBD_IgG + Spike_IgG,
     data = trainDat)
	mincp = fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"]## use se rule
	fitpruned=prune(fit, cp=mincp)
	tt=predict(fitpruned, testDat) 
}

logregTrainFunc=function(trainDat, testDat){

	fit = glm(myy ~ NC_IgG + RBD_IgG + Spike_IgG,
     data = trainDat,family="binomial")
	tt=predict(fit, testDat) 
}

bagTrainFunc=function(trainDat, testDat){
	modelsl = lapply(c(1:30),function(dummy){
		tmp = c(1:nrow(trainDat))
		newSample = sample(tmp,replace=TRUE) 
			fit = rpart(myy ~ NC_IgG + RBD_IgG + Spike_IgG,
	     data = trainDat[newSample,])
			mincp = fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"]## use se rule
		fitpruned=prune(fit, cp=mincp)
	})
	resl= lapply(modelsl, function(x){predict(x, testDat)})
	resMat = do.call("cbind",resl)
	rowSums(resMat)
}


rocList=list()

# gamRoc=getROCForComposite(castTbl, gamTrainFunc)
# gamRoc[,target:="gamComposite"]

# rocList[["gam"]]=gamRoc

# marsRoc=getROCForComposite(castTbl, marsTrainFunc)
# marsRoc[,target:="marsComposite"]
# rocList[["mars"]]=marsRoc

# treeRoc=getROCForComposite(castTbl, treeTrainFunc)
# treeRoc[,target:="treeComposite"]
# rocList[["tree"]]=treeRoc


bagRoc=getROCForComposite(castTbl, bagTrainFunc)
bagRoc[,target:="bagComposite"]
rocList[["bag"]]=bagRoc


logregRoc=getROCForComposite(castTbl, logregTrainFunc)
logregRoc[,target:="logregComposite"]
rocList[["logreg"]]=logregRoc

value=svd(castTbl[, list(NC_IgG,RBD_IgG,Spike_IgG)])$u[,1]
castTbl[,mysvd:=-value]
tableForRoc=castTbl[,list(mysvd,EurImmunCall)]
setkey(tableForRoc,mysvd)
myroc=getRocValsInternal(tableForRoc[,EurImmunCall],tableForRoc[,mysvd])
myroc[,target:="svd"]
rocList[["svd"]]=myroc



simplreRocData =fread("Results/immunology_cohort/immunologyCohortROCData.txt")

allRocs=do.call("rbind",rocList)
allRocs=rbind(allRocs,simplreRocData)


ggplot(allRocs,aes(y=tpr,x=fpr,colour=target, group=target))+geom_line(size=1) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + theme(legend.position = c(0.8, 0.2))
ggsave("Figs/immunology_cohort/EurImmunRocCompositeTryout.pdf", width=4, height=4)

ggplot(allRocs,aes(y=tpr,x=fpr,colour=target, group=target))+geom_line(size=1) + xlim(0, 0.05)+ xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + theme(legend.position = c(0.8, 0.2))
ggsave("Figs/immunology_cohort/EurImmunRocCompositeTryout_crop.pdf", width=4, height=4)



allRocs=do.call("rbind",rocList)
allRocs=rbind(allRocs,simplreRocData)


ggplot(myroc,aes(y=tpr,x=fpr,colour=target, group=target))+geom_line(size=1) + xlab("False Positive Rate") + ylab("True Positive Rate")  + theme_bw() + theme(legend.position = c(0.8, 0.2))



subtbl=mergedTbl[,list(neglog, idWithoutMonth, EurImmunCall, target)]
castTbl=dcast(subtbl,idWithoutMonth + EurImmunCall ~target,value.var="neglog")
require(rpart)
fit=rpart(EurImmunCall ~ NC_IgG + RBD_IgG + Spike_IgG , data= castTbl, method = "anova")
print(plotcp(fit))
mincp = fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"]## use se rule
fitpruned=prune(fit, cp=mincp)


draw.tree(fit)

ggplot(castTbl[NC_IgG>-5 & Spike_IgG>-5,],aes(x=NC_IgG, y=Spike_IgG,colour=myy))+geom_point()


