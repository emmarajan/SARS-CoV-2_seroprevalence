require(data.table)
####################################################################
require(data.table)
require(ggplot2)
require(QRM)
system("mkdir Figs/kim_cohort/")
system("mkdir Results/kim_cohort/")
system("mkdir interimData/kim_cohort/")

source("Code/prepareDataHelper.R")
source("Code/helperFunctions.R")

myvar="neglog"
vars=c("neglog","neglogNormed")

#USEUSZBACKGROUND=FALSE
USEBLOOD=TRUE
USEUSZBACKGROUND=FALSE

ONLYHISTORICAL=TRUE



ggl=lapply(c(FALSE,TRUE),function(USEUSZBACKGROUND){

resl=lapply(c(TRUE,FALSE),function(USEBLOOD){

uscTbl=fread(paste("Results/usz_cohort/fullUscTbl",("blood")[USEBLOOD],".txt",sep=""))
uscTbl=uscTbl[mydate> -1,]
fprVsCutoffTbl=lapply(vars, function(myvar){	
	if(USEBLOOD){
		if(ONLYHISTORICAL){
			mytbl=uscTbl[180000==mydate,c(c("mydate","Patient_or_Control_ID","target"),myvar),with=FALSE]
			}else{
			mytbl=uscTbl[200100>=mydate,c(c("mydate","Patient_or_Control_ID","target"),myvar),with=FALSE]
		}
	}else{
		mytbl=uscTbl[year<19,c(c("mydate","Patient_or_Control_ID","target"),myvar),with=FALSE]
	}
	setnames(mytbl,c("mydate","id","target","pheno"))
	tmpTbl2=mytbl[target!="NC_IgG",list(mydate=mydate[1],target="averagedSpikeRBD",pheno=(mean(pheno)*sqrt(2)),myy=0),by=id]
	##tmpTbl2=mytbl[,list(target="averaged",pheno=(mean(pheno)*sqrt(3)),myy=0),by=id]
	setkey(tmpTbl2,pheno)
	fprVsCutoffTbl=getRocValsInternal(tmpTbl2[,myy], tmpTbl2[,pheno])
	fprs=c(0.001, 0.005,0.01, 0.02,0.05)
	res = sapply(fprs, function(x){
		tail(fprVsCutoffTbl[fpr<x,val],1)
	})
	tbl = data.table(method=myvar,fpr=fprs, cutoff=res)
	return(tbl)
})
fprVsCutoffTbl=do.call("rbind",fprVsCutoffTbl)
write.table(fprVsCutoffTbl, file = paste("interimData/cutoffs",("blood")[USEBLOOD],("_onlyHistorical")[ONLYHISTORICAL],".txt",sep=""), sep= "\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
if(USEUSZBACKGROUND){
	fprVsCutoffTbl=fread("interimData/cutoffs.txt")
}
####################################################################
####################################################################

#myvar="neglog"
myvars=c("neglog","neglogNormed")


getList=function(uscTbl, fprVsCutoffTbl,myvar, myfpr){
	mytbl=uscTbl[,c(c("Patient_or_Control_ID","target","mydate"),myvar),with=FALSE]
	setnames(mytbl,c("id","target","mydate","pheno"))
	mycutoff=fprVsCutoffTbl[fpr==myfpr & method==myvar,cutoff]
	dataTbl=mytbl[target!="NC_IgG",list(mydate=mydate[1],target="averagedSpikeRBD",pheno=(mean(pheno)*sqrt(2))),by=id]
	#dataTbl=mytbl[,list(mydate=mydate[1],target="averaged",pheno=(mean(pheno)*sqrt(3))),by=id]
	dataTbl[pheno>mycutoff,]
}

tblList=lapply(myvars, function(myvar){
	rr=fprVsCutoffTbl[method==myvar,fpr]
	lapply(rr, function(myfpr){
		mylist=getList(uscTbl, fprVsCutoffTbl,myvar, myfpr)
		if(myfpr==0.001){
			write.table(mylist,file=paste("interimData/usz_cohort/selectedsamples_",("blood")[USEBLOOD],("_uszBackg")[USEUSZBACKGROUND],"_",myvar," fpr",myfpr,".txt",sep=""), sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)
		}
		mylist
	})
})


ff5=lapply(myvars, function(myvar){
	mytbl=uscTbl[,c(c("Patient_or_Control_ID","target","mydate"),myvar),with=FALSE]
	setnames(mytbl,c("id","target","mydate","pheno"))
	#ff=mytbl[,list(mydate=mydate[1],target="averaged",pheno=(mean(pheno)*sqrt(3))),by=id]
	ff=mytbl[target!="NC_IgG",list(mydate=mydate[1],target="averagedSpikeRBD",pheno=(mean(pheno)*sqrt(2))),by=id]
	setkey(ff,mydate)
	ff[,mymonth:=floor(mydate/100)]
	SS=data.table(unique(ff[,mymonth]),c(1:length(unique(ff[,mymonth]))))
	ff[,monthRank:=match(ff[,mymonth],SS[,V1])]
	rr=fprVsCutoffTbl[method==myvar,fpr]
	ff3=lapply(rr, function(myfpr){
		mycutoff=fprVsCutoffTbl[fpr==myfpr & method==myvar,cutoff]
		ff2=ff[,list(fractionAboveCutoff=mean(pheno>=mycutoff),nrAboveCutoff=sum(pheno>=mycutoff),nrBelowCutoff=sum(pheno<mycutoff)),by=mymonth]
		ff2[,rankMonth:=c(1:length(mymonth))]
		ff2[,mymontFac:=factor(mymonth)]
		ff2[,fpr:=myfpr]
		return(ff2)
	})
	ff4=do.call("rbind",ff3)
	ff4[,myVar:=myvar]
	return(ff4)
})
ff6=do.call("rbind",ff5)
#            id mydate           target    pheno
# 1: J612290617 161229 averagedSpikeRBD 2.244706
ggplot(ff6[myVar=="neglog",], aes(y=fractionAboveCutoff,x=mymontFac, group=mymontFac)) + geom_bar(stat="identity") + facet_grid(fpr~.) #+ geom_hline(mycutoff)
ggsave(paste("Figs/usz_cohort/nCasesPerCutoff",("blood")[USEBLOOD],("_uszBackg")[USEUSZBACKGROUND],".pdf",sep=""),width=6,height=4)


getNormedTrueCasesEstimate=function(ff7){
	ff8=ff7[,list(above=sum(nrAboveCutoff),below=sum(nrBelowCutoff)),by=c("fpr","ispos")]
	ff8[,mytot:=above+below]
	ff8[,myfraction:=above/mytot]
	ff8pos=ff8[ispos==TRUE,]
	ff8neg=ff8[ispos==FALSE,]
	ff8pos[,negFraction:=ff8neg[,myfraction]]
	ff8pos[,expectedFalse:=negFraction*mytot]
	ff8pos[,esimated_cases_normed:=(above-expectedFalse)/mytot]
	ff8pos[,stratum:="overall"]
	return(ff8pos)
}

fromNonBloodFun=function(ff7){	
	ff8=ff7[,list(above=sum(nrAboveCutoff),below=sum(nrBelowCutoff)),by=c("fpr","ispos")]
	ff8[,mytot:=above+below]
	ff8[,myfraction:=above/mytot]
	ff8pos=ff8[ispos==TRUE,]
	ff8neg=ff8[ispos==FALSE,]
	ff8pos[,negFraction:=ff8neg[,myfraction]]
	ff8pos[,expectedFalse:=negFraction*mytot]
	ff8pos[,(above-expectedFalse)/mytot]
	ff9pos=getNormedTrueCasesEstimate(ff7)
}

if(USEBLOOD){
	ff6[,ispos:=mymonth>=2100]
	ff7=ff6[myVar=="neglog", ]
	ff9pos=getNormedTrueCasesEstimate(ff7)
	ff9pos[,category:="BSD"]
	return(ff9pos)
}else{
	ff6[,ispos:=mymonth>1900]
	ff7=ff6[myVar=="neglog",]
	ff9pos=getNormedTrueCasesEstimate(ff7)
	ff9pos=fromNonBloodFun(ff7)
	ff9pos[,category:="USZ"]
	ff9pos[,stratum:="overall"]
	return(ff9pos)
}
})

browser()
res=do.call("rbind",resl)
ggplot(res, aes(y=esimated_cases_normed,x=fpr,colour=category,group=category))+geom_point()+geom_line()+theme_bw() +ylab("prevalence estimate")
ggsave(paste("Figs/usz_cohort/estimated_true_positives_analysis",("_uszBackg")[USEUSZBACKGROUND],".pdf",sep=""),width=4,height=4)

})


