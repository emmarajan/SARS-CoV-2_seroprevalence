require(data.table)
require(mvtnorm)
require(QRM)
require(MASS)

getPrevalenceFull=function(posNegTbl, valTbl, sureNegLevel= -3){
	valTblOut=valTbl[1:nrow(valTbl),]
	posNegTbl2=posNegTbl[,list(neglog,myy)]
	valTblUsed=valTbl[neglog > sureNegLevel,]
	valTblNotused=valTbl[neglog <= sureNegLevel,]
	valTbl2Used=valTblUsed[,list(neglog)]
	valTbl2Used[,myy:=2]
	combTbl=rbind(valTbl2Used,posNegTbl2[neglog> sureNegLevel,])
	knownPos=combTbl[,myy==1]
	knownNeg=combTbl[,myy==0]
	ee=estimateFracFull(combTbl[,neglog],knownPos,knownNeg)
	gg=ee[[2]][!knownNeg & !knownPos]
	valTblOut[,posteriorProb:=0]
	valTblOut[neglog > sureNegLevel,posteriorProb:=gg]
	potential=nrow(valTblUsed)
	impossible=nrow(valTblNotused)
	estimated_prob=(potential*ee[[1]])/(potential+impossible)
	out=list()
	out$resTbl=valTblOut
	out$estimated_prob=estimated_prob
	out$paraml=ee$paraml
	return(out)
}

# posNegTbl needs myy=(0|1) and neglog for both.
getPrevalence=function(posNegTbl, valTbl, sureNegLevel= -3){
	valTblUsed=valTbl[neglog> sureNegLevel,]
	valTblNotused=valTbl[neglog<= sureNegLevel,]
	posNegTbl[myy==1,neglog]
	sd2=sd(posNegTbl[myy==1,neglog])
	mean2=mean(posNegTbl[myy==1,neglog])
	sd1=sd(posNegTbl[myy==0 & neglog> sureNegLevel,neglog])
	mean1=mean(posNegTbl[myy==0 & neglog> sureNegLevel,neglog])
	myvals=valTblUsed[,neglog]
	ee=estimateFrac(myvals,mean1,sd1,mean2,sd2)
	valTblUsed[,posteriorProb:=1-ee[[2]]]
	meanProb=1-ee[[1]]
	valTblNotused[,posteriorProb:=0]
	potential=nrow(valTblUsed)
	impossible=nrow(valTblNotused)
	estimated_prob=(potential*meanProb)/(potential+impossible)
	valTblOut=rbind(valTblUsed,valTblNotused)
	outl=list()
	outl$resTbl=valTblOut
	outl$estimated_prob=estimated_prob
	outl$paraml=list(mean1,sd1,mean2,sd2)
	return(outl)
}

estimateMVNFracFull=function(valueMat,knownPos,knownNeg,nround=50){
	mean1=colMeans(valueMat[knownPos,,drop=FALSE])
	mean2=colMeans(valueMat[knownNeg,,drop=FALSE])
	cov1=cov(valueMat[knownPos,,drop=FALSE])
	cov2=cov(valueMat[knownNeg,,drop=FALSE])
	vec2Mat=function(vec,mylen){
		kronecker(vec,t(rep(1,mylen)))
	}
	myp=0.02
	for(i in c(1:nround)){
		density1=dmvnorm(valueMat,mean=mean1,sigma=cov1,log=FALSE)
		density2=dmvnorm(valueMat,mean=mean2,sigma=cov2,log=FALSE)
		myt=(myp*density1)/((myp*density1)+((1-myp)*density2))
		myt[knownPos]=1
		myt[knownNeg]=0
		myp=mean(myt[knownPos==FALSE & knownNeg==FALSE])
##		myp=mean(myt)
		myTMat=myt%*%t(rep(1,ncol(valueMat)))
		mean1=colSums(myTMat * valueMat)/sum(myt)
		mean2=colSums((1-myTMat) * valueMat)/sum(1-myt)
		mean1Mat=t(vec2Mat(mean1, length(myt)))
		mean2Mat=t(vec2Mat(mean2, length(myt)))
		mat1=valueMat-mean1Mat
		mat2=valueMat-mean2Mat
		cov1=t(mat1)%*%(myTMat*mat1)/sum(myt)
		cov2=t(mat2)%*%((1-myTMat)*mat2)/sum(1-myt)
		print(myp)
	}
	print("done")
	paraml=list(mean1=mean1,cov1=cov1,mean2=mean2,cov2=cov2)
	outl=list(myp=myp,myt=myt,paraml=paraml)
	return(outl)
}


estimateMVMixedFrac=function(valueMat,knownPos,knownNeg,nround=50){
	gg=fit.mst(valueMat[knownNeg,,drop=FALSE], method = "BFGS")
	mydf=gg$df
	mean1=colMeans(valueMat[knownPos,,drop=FALSE])
	cov1=cov(valueMat[knownPos,,drop=FALSE])
	sigma2=as.matrix(gg$Sigma)
	mean1=colMeans(valueMat[knownPos,,drop=FALSE])
	mean2=colMeans(valueMat[knownNeg,,drop=FALSE])
	myp=0.02
	for(i in c(1:nround)){
		density1=dmvnorm(valueMat,mean=mean1,sigma=cov1,log=FALSE)
		density2=dmvt(valueMat,delta=mean2,sigma=sigma2,df=mydf,log=FALSE)
		myt=(myp*density1)/((myp*density1)+((1-myp)*density2))
		myt[knownPos]=1
		myt[knownNeg]=0
		myp=mean(myt[knownPos==FALSE & knownNeg==FALSE])
		print(myp)
	}
	print("done")
	paraml=list(mean1=mean1,cov1=cov1,mean2=mean2,sigma2=sigma2,df=mydf)
	outl=list(myp=myp,myt=myt,paraml=paraml)
	return(outl)
}
getPosteriorFromVMixedFull=function(valueMat,fittedModel){
	myp=fittedModel[[1]]
	paraml=fittedModel[[3]]
	density1=dmvnorm(valueMat,paraml$mean1,sigma=paraml$cov1,log=FALSE)
	density2=dmvt(valueMat,delta=paraml$mean2,sigma=paraml$sigma2,df=paraml$df,log=FALSE)
	myt=(myp*density1)/((myp*density1)+((1-myp)*density2))
	return(myt)
}

estimateMVTFrac=function(valueMat,knownPos,knownNeg,nround=50){
	gg=fit.mst(valueMat[knownNeg,,drop=FALSE], method = "BFGS")
	mydf=gg$df
	mean1=colMeans(valueMat[knownPos,,drop=FALSE])
	cov1fit=cov.trob(valueMat[knownPos,,drop=FALSE],nu=mydf)
	cov1=cov1fit$cov
	sigma1=(mydf-2)/mydf*cov1
	sigma2=as.matrix(gg$Sigma)
	mean1=cov1fit$center #colMeans(valueMat[knownPos,,drop=FALSE])
	mean2=colMeans(valueMat[knownNeg,,drop=FALSE])	
	myp=0.05
	for(i in c(1:nround)){
		density1=dmvt(valueMat,delta=mean1,sigma=sigma1,df=mydf,log=FALSE)
		density2=dmvt(valueMat,delta=mean2,sigma=sigma2,df=mydf,log=FALSE)
		myt=(myp*density1)/((myp*density1)+((1-myp)*density2))
		myt[knownPos]=1
		myt[knownNeg]=0
		myp=mean(myt[knownPos==FALSE & knownNeg==FALSE])	
		print(myp)
	}
	print("done")
	paraml=list(mean1=mean1,sigma1=sigma1,mean2=mean2,sigma2=sigma2,df=mydf)
	outl=list(myp=myp,myt=myt,paraml=paraml)
	return(outl)
}

getPosteriorFromMVTFull=function(valueMat,fittedModel){
	myp=fittedModel[[1]]
	paraml=fittedModel[[3]]
	density1=dmvt(valueMat,delta=paraml$mean1,sigma=paraml$sigma1,df=paraml$df,log=FALSE)
	density2=dmvt(valueMat,delta=paraml$mean2,sigma=paraml$sigma2,df=paraml$df,log=FALSE)
	myt=(myp*density1)/((myp*density1)+((1-myp)*density2))
	return(myt)
}


QDA=function(valueMat,knownPos,knownNeg){
	mean1=colMeans(valueMat[knownPos,,drop=FALSE])
	mean2=colMeans(valueMat[knownNeg,,drop=FALSE])
	cov1=cov(valueMat[knownPos,,drop=FALSE])
	cov2=cov(valueMat[knownNeg,,drop=FALSE])	
	ldensity1=dmvnorm(valueMat,mean=mean1,sigma=cov1,log=TRUE)
	ldensity2=dmvnorm(valueMat,mean=mean2,sigma=cov2,log=TRUE)
	myt=ldensity1-ldensity2
	paraml=list(mean1=mean1,cov1=cov1,mean2=mean2,cov2=cov2)
	outl=list(myt=myt,paraml=paraml)
	return(outl)
}

getQDAFromFit=function(valueMat,fittedModel){
	paraml=fittedModel[[2]]
	density1=dmvnorm(valueMat,mean=paraml$mean1,sigma=paraml$cov1,log=TRUE)
	density2=dmvnorm(valueMat,mean=paraml$mean2,sigma=paraml$cov2,log=TRUE)
	myt=density1-density2
	return(myt)
}

estimateMVNFrac=function(valueMat,knownPos,knownNeg,nround=50){
	mean1=colMeans(valueMat[knownPos,,drop=FALSE])
	mean2=colMeans(valueMat[knownNeg,,drop=FALSE])
	cov1=cov(valueMat[knownPos,,drop=FALSE])
	cov2=cov(valueMat[knownNeg,,drop=FALSE])
	myp=0.02
	for(i in c(1:nround)){
		density1=dmvnorm(valueMat,mean=mean1,sigma=cov1,log=FALSE)
		density2=dmvnorm(valueMat,mean=mean2,sigma=cov2,log=FALSE)
		myt=(myp*density1)/((myp*density1)+((1-myp)*density2))
		myt[knownPos]=1
		myt[knownNeg]=0
		myp=mean(myt[knownPos==FALSE & knownNeg==FALSE])
		print(myp)
	}
	print("done")
	paraml=list(mean1=mean1,cov1=cov1,mean2=mean2,cov2=cov2)
	outl=list(myp=myp,myt=myt,paraml=paraml)
	return(outl)
}

estimateLDAFrac=function(valueMat,knownPos,knownNeg,nround=50,use1dim=TRUE,use1dimCor=TRUE){
	mean1=colMeans(valueMat[knownPos,,drop=FALSE])
	mean2=colMeans(valueMat[knownNeg,,drop=FALSE])
	##cov1=cov(valueMat[knownPos,,drop=FALSE])
	cov2=cov(valueMat[knownNeg,,drop=FALSE])
	myp=0.02
	annot=rep(0,nrow(valueMat))
	annot[knownPos]=1
	annot[knownNeg]=-1

	linearComb=solve(cov2)%*%(mean1-mean2)
	tstat=valueMat%*%linearComb
	linVar=(t(mean1-mean2)%*%solve(cov2)%*%(mean1-mean2))[1]
	linMean1=(t(mean1-mean2)%*%solve(cov2)%*%(mean1))[1]
	linMean2=(t(mean1-mean2)%*%solve(cov2)%*%(mean2))[1]
	sd1=sqrt(var(tstat[knownPos]))
	sd2=sqrt(var(tstat[knownNeg]))
	for(i in c(1:nround)){
		if(use1dim){
			density1_lin=dnorm(tstat,mean=linMean1,sqrt(linVar))
			density2_lin=dnorm(tstat,mean=linMean2,sqrt(linVar))
			if(use1dimCor){
				density1=dnorm(tstat,mean=mean(tstat[knownPos]),sd=sd1)
				density2=dnorm(tstat,mean=mean(tstat[knownNeg]),sd=sd2)
			}
			density1=density1_lin
			density2=density2_lin
		}else{
			density1=dmvnorm(valueMat,mean=mean1,sigma=cov2,log=FALSE)
			density2=dmvnorm(valueMat,mean=mean2,sigma=cov2,log=FALSE)
		}
		myt=(myp*density1)/((myp*density1)+((1-myp)*density2))
		myt[knownPos]=1
		myt[knownNeg]=0
		myp=mean(myt[knownPos==FALSE & knownNeg==FALSE])
		print(myp)
	}
	print("done")
	paraml=list(mean1=mean1,cov1=NA,mean2=mean2,cov2=cov2,linearComb=linearComb,linMean1=linMean1,linMean2=linMean2,sd1=sd1,sd2=sd2)
	outl=list(myp=myp,myt=myt,paraml=paraml,mytstat=tstat,annot=annot)
	return(outl)
}


getPosteriorFromMVNFull=function(valueMat,fittedModel){
	myp=fittedModel[[1]]
	paraml=fittedModel[[3]]
	density1=dmvnorm(valueMat,mean=paraml$mean1,sigma=paraml$cov1,log=FALSE)
	density2=dmvnorm(valueMat,mean=paraml$mean2,sigma=paraml$cov2,log=FALSE)
	myt=(myp*density1)/((myp*density1)+((1-myp)*density2))
	return(myt)
}

getPosteriorFromLDAFull=function(valueMat,fittedModel,use1dimCor=TRUE){
	myp=fittedModel[[1]]
	paraml=fittedModel[[3]]
	if(use1dimCor){
		mytstat=valueMat%*%paraml$linearComb
		density1=dnorm(mytstat,paraml$linMean1,paraml$sd1)
		density2=dnorm(mytstat,paraml$linMean2,paraml$sd2)
	}else{
		density1=dmvnorm(valueMat,mean=paraml$mean1,sigma=paraml$cov2,log=FALSE)
		density2=dmvnorm(valueMat,mean=paraml$mean2,sigma=paraml$cov2,log=FALSE)
	}
	myt=(myp*density1)/((myp*density1)+((1-myp)*density2))
	return(myt)
}


getLikelihoodRatio=function(posteriorProb,prevalence){
	myt=posteriorProb
	myp=prevalence
	(myt/(1-myt))*((1-myp)/myp)
}

getPosteriorFromRatio=function(lRatio,prevalence){
	postInv=1+((1-prevalence)/prevalence)/lRatio
	return(postInv)
}


estimateFracMVNFullTest=function(){
	mean1 = c(0.7,0.2,-0.5)
	mean2 = c(0.2,-0.1,0.8)
	cov1 = toeplitz(c(1,0.8,0.64))
	cov2 = toeplitz(c(1,0.5,0.25))	
	xx1=rmvnorm(300,mean1,cov1)
	xx2=rmvnorm(300,mean2,cov2)
	valueVec=rbind(xx1,xx2)
	knownPos=rep(FALSE,nrow(valueVec))
	knownPos[1:150]=TRUE
	knownNeg=rep(FALSE,nrow(valueVec))
	knownNeg[400:600]=TRUE
	out=estimateMVNFracFull(valueVec,knownPos,knownNeg)
	return(out)
}

estimateFracFull=function(valueVec,knownPos,knownNeg){
	mean1=mean(valueVec[knownPos])
	mean2=mean(valueVec[knownNeg])
	sd1=sd(valueVec[knownPos])
	sd2=sd(valueVec[knownNeg])
	myp=0.02
	for(i in c(1:50)){
		density1=dnorm(valueVec,mean=mean1,sd=sd1,log=FALSE)
		density2=dnorm(valueVec,mean=mean2,sd=sd2,log=FALSE)
		myt=(myp*density1)/((myp*density1)+((1-myp)*density2))
		myt[knownPos]=1
		myt[knownNeg]=0
		myp=mean(myt[knownPos==FALSE & knownNeg==FALSE])
##		myp=mean(myt)
		mean1=sum(myt*valueVec)/sum(myt)
		mean2=sum((1-myt)*valueVec)/sum(1-myt)
		var1=sum(myt*((valueVec-mean1)^2))/sum(myt)
		var2=sum((1-myt)*((valueVec-mean2)^2))/sum(1-myt)
		sd1=sqrt(var1)
		sd2=sqrt(var2)
		print(myp)
	}
	paraml=list(mean1=mean1,sd1=sd1,mean2=mean2,sd2=sd2)
	outl=list(myp=myp,myt=myt,paraml=paraml)
	return(outl)
}


estimateFracFullTest=function(){
	mean1 = 0.7
	mean2 = -0.1
	sd1 = 0.2
	sd2 = 0.2
	xx1=rnorm(1000,mean1,sd1)
	xx2=rnorm(1000,mean2,sd2)
	valueVec=c(xx1,xx2)
	knownPos=rep(FALSE,length(valueVec))
	knownPos[1:400]=TRUE
	knownNeg=rep(FALSE,length(valueVec))
	knownNeg[1900:2000]=TRUE
	out=estimateFracFull(valueVec,knownPos,knownNeg)
	return(out)
}

estimateFrac=function(valueVec,mean1,sd1,mean2,sd2){
	pInit=0.02
	density1=dnorm(valueVec,mean=mean1,sd=sd1,log=FALSE)
	density2=dnorm(valueVec,mean=mean2,sd=sd2,log=FALSE)
	myp=pInit
	for(i in c(1:300)){
		myt=(myp*density1)/((myp*density1)+((1-myp)*density2))
		myp=mean(myt)
	}
	outl=list(myp,myt)
	return(outl)
}

estimateFracTest=function(){
	mean1 = 0.7
	mean2 = -0.1
	sd1 = 0.5 
	sd2 = 0.3 
	xx1=rnorm(20000,mean1,sd1)
	xx2=rnorm(80000,mean2,sd2)
	valueVec=c(xx1,xx2)
	out=estimateFrac(valueVec,mean1,sd1,mean2,sd2)
	return(out)
}


collapseTbl=function(mytbl){
	mytbl[target!="NC_IgG",list(target="averagedSpikeRBD",pheno=(mean(pheno)*sqrt(2)),myy=0),by=id]
}

getRocValsInternal = function(condSortedByVal, sortedVal){
	myrev=rev(condSortedByVal)
	myrevVal = rev(sortedVal)
	tp = cumsum(myrev)
	myp = sum(myrev)
	tpr = tp/myp
	fp = cumsum(!myrev)
	myn = sum(!myrev)
	fpr = fp/myn
	tn = rev(cumsum(rev(!myrev)))
	tnr = tn/myn
	ppv = tp/ c(1:length(myrev))
	out=data.table(fpr,tpr, tnr,ppv,val=myrevVal)
	out[,torem:=rev(duplicated(rev(val)))]
	out=out[torem==FALSE,]
	out[,torem:=FALSE]
	return(out)
}


getLogModel=function(dataTbl){
	tbl4=dataTbl[is.element(target,c("Spike_IgG","RBD_IgG","NC_IgG")),]
	tbl4[neglog< -3,neglog:=-3]
	castTbl=dcast(tbl4,Pseudo_ID+ myy ~ target ,value.var="neglog")
	fit=glm(myy ~ Spike_IgG + RBD_IgG + NC_IgG, family="binomial",data=castTbl)
	return(fit)
}


getPredictForComposite=function(castTbl, myfunc){
	mycastTbl=data.table(castTbl)
	set.seed(11)
	indexVec=sample(rep(rep(c(1:10)),ceiling(nrow(mycastTbl)/10))[1:nrow(mycastTbl)])
	resVec = rep(NA, length(indexVec))
	for (i in c(1:10)){
		sum(mycastTbl[indexVec==i,myy])
		testDat = mycastTbl[indexVec==i,]
		trainDat = mycastTbl[indexVec!=i,]
		tt=myfunc(trainDat, testDat)
		resVec[indexVec==i]=tt
	}
	outTbl=mycastTbl[,list(Pseudo_id,target="logPredCrossValid")]
	outTbl[,neglog:=resVec]
	outTbl[,myy]=mycastTbl[,myy]
	return(outTbl)
}


getPredictFromExternal=function(castTbl, modelTbl){
	mycastTbl=data.table(castTbl)
	mymatcher=match(modelTbl[,varname],colnames(castTbl))
	mypred=(as.matrix(mycastTbl[,mymatcher,with=FALSE])%*%modelTbl[,weight])[,1]
	outTbl=mycastTbl[,list(Pseudo_ID,target="externalPred")]
	outTbl[,neglog:=mypred]
	newy=mycastTbl[,myy]
	outTbl[,myy:=newy]
	return(outTbl)
}

getPredict=function(dataTbl, glmFit){
	coefs=glmFit$coef
	len=length(coefs)
	modelTbl=data.table(varname=names(coefs)[2:len],weight=coefs[2:len])
	subTbl=dataTbl[is.element(target,modelTbl[,varname]),list(Pseudo_ID,target,neglog,myy)]
	castTbl=dcast(Pseudo_ID+myy~target,data=subTbl,value.var="neglog")
	getPredictFromExternal(castTbl, modelTbl)
}




###########################
## fitting data
###########################


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
	tableForRoc=castTbl[,list(gamPred,myy)]
	setkey(tableForRoc,gamPred)
	myroc=getRocValsInternal(tableForRoc[,myy],tableForRoc[,gamPred])
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


spikeTrainFunc=function(trainDat, testDat){
	fit = glm(myy ~ Spike_IgG,
     data = trainDat,family="binomial")
	tt=predict(fit, testDat) 
}

meanFunc=function(trainDat, testDat){
	tt=testDat$NC_IgG + testDat$RBD_IgG + testDat$Spike_IgG
}