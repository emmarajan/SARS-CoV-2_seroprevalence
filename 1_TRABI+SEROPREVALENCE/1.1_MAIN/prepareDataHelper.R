
##############################################################
## parse data
##############################################################

require(data.table)
require(ggplot2)


getInHouseFilePaths=function(proteinName,addstr=""){
	allFileNames=list.files(paste("Data/SARS_CoV_2_share/",addstr,sep=""))
	allFileNames=allFileNames[!grepl("20200519",allFileNames)]
	allFileNames=allFileNames[grepl("IgG",allFileNames)]
	inds=grepl(proteinName,allFileNames)
	subset=allFileNames[inds]
	return(subset)
}

getReplicationFitPaths = function(nr1=TRUE){
	allFileNames=list.files(paste("Data/SARS_CoV_2_share/",sep=""))
	allFileNames=allFileNames[grepl("20200519",allFileNames)]
	if(nr1){
		allFileNames=allFileNames[grepl("9_1_",allFileNames)]
	}else{
		allFileNames=allFileNames[grepl("9_2_",allFileNames)]
	}
	out=paste("Data/SARS_CoV_2_share/",allFileNames,"/Analysis/fits_out.csv",sep="")
	return(out)
}

getInHouseFitPaths = function(proteinName,useImunology=FALSE,useKim=FALSE,useOxford=FALSE, getLayout=FALSE){
	if(useImunology + useKim + useOxford >1){
		stop("only one active allowed")
	}
	immustr=("20200503_Immunology_dataset/")[useImunology]
	kimstr=("20200506_KIM_reference_dataset/")[useKim]
	oxfordstr=("20200512_Oxford_cohort/")[useOxford]
	addstr=paste(immustr,kimstr,oxfordstr,sep="")
	paths = getInHouseFilePaths(proteinName,addstr)
	if(getLayout){
		mypath=paste("Data/SARS_CoV_2_share/",addstr,paths,"/Analysis/",sep="")
		ww=list.files(mypath)
		inds=grepl("hts.csv$",ww)
		combined=paste(mypath,ww,sep="")
		return(combined[inds])
		}else{
	return(paste("Data/SARS_CoV_2_share/",addstr,paths,"/Analysis/fits_out.csv",sep=""))
	}	
}



loadPath = function(){
	indsMat = kronecker(c(1:6),c(1:2))
	mats=lapply(c("Spike","NC"), function(x){
		paste(c(1:6),rep(x,6),sep="_")
	})
	res=do.call("c",mats)
	tbl=data.table(x=res)
	tbl[,antigen:=sub("\\d+_","",x)]
	tbl[,tablenr:=sub("_.+","",x)]
	tbl[,loadPath := paste("Data/Summary/",x,"-Table.csv",sep="")]
	out = lapply(c(1:nrow(tbl)), function(i){
		subtbl=tbl[i,]
		ww=fread(subtbl[,loadPath])
		antigen=subtbl[,antigen]
		ww[,antigen:=antigen]
		tablenr=subtbl[,tablenr]
		ww[,tablenr:=tablenr]
		return(ww)
	})
	outTbl = do.call("rbind",out)
}
##############################################################
## END: parse data
##############################################################
extractDate = function(dataTbl){
	dataTbl2=dataTbl[grepl("^J",Patient_or_Control_ID),]
	dataTbl2[,editedId:=sub("J0","20",Patient_or_Control_ID)]
	dataTbl2[,editedId:=sub("J","1",editedId)]
	dataTbl2[,samplenr:=as.numeric(sub("\\d{6}","",editedId))]
	dataTbl2[,mydate:=as.numeric(sub("(\\d{6}).+","\\1",editedId))]
	dataTbl2[,year:=floor(mydate/10000)]
	dataTbl2[,month:=floor((mydate - year*10000)/100)]
	dataTbl2[,day:= mydate %% 100]
	dataTbl2[mydate<190000,mydate:=180000]
	return(dataTbl2)
}

extractDateBDS = function(dataTbl){
	ids=dataTbl[,Patient_or_Control_ID]
	dateVec=rep(NA,length(ids))
	bdsnr=grepl("^BDS",ids)	
	dateVec[bdsnr]=180000
	dateVec[!bdsnr]=210000
	between_index=nchar(ids)==13
	tomap=ids[between_index]
	asnum=as.numeric(sub("BDS","",tomap))
	mydates=(floor(asnum/1E4) - 2E5)*100
	dateVec[between_index]=mydates
	newDate=dateVec
	dataTbl[,mydate:=newDate]
	return(dataTbl)
}


dataTbl=extractDate(loadPath())
##############################################################
## 
##############################################################
dataTbl[,inPandemic:=mydate>200200]
dataTbl=dataTbl[Error*1.5<corrlogEC50,]

require(reshape2)

dataTblNC=dataTbl[antigen=="NC",list(id=Patient_or_Control_ID, valNC=corrlogEC50,inPandemic)]
dataTblSpike=dataTbl[antigen=="Spike",list(id=Patient_or_Control_ID, valSpike=corrlogEC50)]
###
setkey(dataTblNC,id)
setkey(dataTblSpike,id)

getEnrichment = function(inPandemicSorted){
	cumsum(inPandemicSorted)/c(1:length(inPandemicSorted))
}


EE=merge(dataTblNC,dataTblSpike)
EE[,valNormed:=scale(valNC)+scale(valSpike)]
EE[,val:=valNC+valSpike]
toOrderBy = EE[,val]
enrichmentInPandemic=function(toOrderBy){
	matcher=order(toOrderBy)
	inPandemicSorted=rev(EE[matcher,inPandemic])
	getEnrichment(inPandemicSorted)
}

getTblEnrichmentInPandemic=function(toOrderBy){
	data.table(pos=c(1:length(toOrderBy)),enrVal=enrichmentInPandemic(toOrderBy))
}


# valEnr=getTblEnrichmentInPandemic(EE[,val])
# valNormedEnr=getTblEnrichmentInPandemic(EE[,valNormed])
# valNCEnr=getTblEnrichmentInPandemic(EE[,valNC])
# valSpikeEnr=getTblEnrichmentInPandemic(EE[,valSpike])

# valEnr[,type:="combined"]
# valNormedEnr[,type:="combinedNormed"]
# valNCEnr[,type:="NC"]
# valSpikeEnr[,type:="Spike"]

# XX=rbind(valEnr,valNormedEnr)
# XX=rbind(XX,valNCEnr)
# XX=rbind(XX,valSpikeEnr)

# ggplot(data=XX,aes(y=enrVal,x=sqrt(pos),colour=type,group=type))+geom_point() + theme_bw() + theme(legend.position=c(0.80,0.85)) + xlab("fraction of samples from Pandemic") + xlab("sample rank (sqrt)")
# ggsave("Figs/fractionInPandemic.pdf",width=4, height=4)

