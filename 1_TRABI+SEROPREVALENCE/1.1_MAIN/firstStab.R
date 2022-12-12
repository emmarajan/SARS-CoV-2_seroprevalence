

##############################################################
## parse data
##############################################################
require(data.table)
require(ggplot2)
source("Code/prepareDataHelper.R")
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


valEnr=getTblEnrichmentInPandemic(EE[,val])
valNormedEnr=getTblEnrichmentInPandemic(EE[,valNormed])
valNCEnr=getTblEnrichmentInPandemic(EE[,valNC])
valSpikeEnr=getTblEnrichmentInPandemic(EE[,valSpike])

valEnr[,type:="combined"]
valNormedEnr[,type:="combinedNormed"]
valNCEnr[,type:="NC"]
valSpikeEnr[,type:="Spike"]

XX=rbind(valEnr,valNormedEnr)
XX=rbind(XX,valNCEnr)
XX=rbind(XX,valSpikeEnr)

ggplot(data=XX,aes(y=enrVal,x=sqrt(pos),colour=type,group=type))+geom_point() + theme_bw() + theme(legend.position=c(0.80,0.85)) + xlab("fraction of samples from Pandemic") + xlab("sample rank (sqrt)")
ggsave("Figs/fractionInPandemic.pdf",width=4, height=4)


# TblNC=dataTbl[antigen=="NC",list(id=Patient_or_Control_ID, valNC=corrlogEC50,inPandemic)]
# dataTblSpike=dataTbl[antigen=="Spike",list(id=Patient_or_Control_ID, valSpike=corrlogEC50)]

# dataTblSpike=dataTbl[antigen=="Spike",]

# dataTblNC[,corrlogEC50]


# ggplot(dataTbl,aes(x=(corrlogEC50>2)*1,fill=inPandemic))+geom_histogram(alpha=0.6) +facet_grid(antigen~.)

# ggplot(dataTbl,aes(x=(corrlogEC50>2)*1,fill=inPandemic))+geom_density(alpha=0.6) +facet_grid(antigen~.)



# ggplot(dataTbl,aes(x=corrlogEC50,fill=paste(antigen,inPandemic)))


# dataTbl[antigen=="NC",sum(corrlogEC50>2),by=inPandemic]


# ggplot(dataTbl[corrlogEC50>2,],aes(x=corrlogEC50,fill=inPandemic))+geom_histogram(alpha=0.6) +facet_grid(antigen~.)

#  dataTbl[antigen=="NC" & corrlogEC50>2,by=inPandemic]

# #+geom_density(alpha=0.6) +facet_grid(antigen~.)

