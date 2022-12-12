require(data.table)
require(ggplot2)
source("Code/prepareDataHelper.R")

bdsClinical=list.files("Data/Clinical/")
bdsClinical=paste("Data/Clinical/",bdsClinical[grepl("csv$",bdsClinical)],sep="")

Luzern=bdsClinical[grepl("Luzern",bdsClinical)]
ZHBSD=bdsClinical[grepl("ZHBSD",bdsClinical)]

LuzernL=lapply(Luzern,function(x){
	ww=fread(x)
	ww2=ww[,list(fullnr=Entnahmenummer,nr=sub("([A-Z]+\\d{12}).+","\\1",Entnahmenummer),year=sub("\\t","",Geburtsjahr),sex=sub(".*(M|W).*","\\1",Geschlecht))]
	return(ww2)
})

luzernTbl=do.call("rbind",LuzernL)

zhL=lapply(ZHBSD,function(x){
	ww=fread(x)
	ww2=ww[,list(fullnr=Entnahmenummer,nr=sub("=?([A-Z]+\\d{12}).*","\\1",Entnahmenummer),year=sub(".*(\\d{4}).*","\\1",Geburtsdatum),sex=sub(".*(M|W).*","\\1",Geschlecht))]
	return(ww2)
})
zhTbl=do.call("rbind",zhL)

bdsTbl=rbind(luzernTbl,zhTbl)


bloodMeasured=fread("interimData/bds_cohort/cohort_combined_outTbl_Blood.txt")
bloodMeasured[,nr:=sub("([A-Z]+\\d{12}).+","\\1",Pseudo_ID)]

setkey(bloodMeasured,nr)
setkey(bdsTbl,nr)
mergedBDS=merge(bdsTbl,bloodMeasured)


uszMeasured=fread("interimData/bds_cohort/cohort_combined_outTbl_Above14.txt")

###############
descrTbl=fread("Data/20200524_descrp_stats.csv")
descrTbl[,Orgfa:=NULL]
descrTbl2=fread("interimData/withNegRemoved.txt")
descrTbl2[,PatientIDList:=Patient_or_Control_ID]
descrTbl2[,Orgfa:=bestfittingOrf]
setkey(descrTbl2,PatientIDList)
setkey(descrTbl,PatientIDList)
descrTbl=merge(descrTbl,descrTbl2)

print(dim(descrTbl))
setkey(descrTbl,Orgfa)


descrTbl[,Pseudo_ID:=PatientIDList]


setkey(descrTbl,Pseudo_ID)
setkey(uszMeasured,Pseudo_ID)

resTbl=merge(descrTbl,uszMeasured)
resTblForRepro=resTbl[,list(Pseudo_ID,posteriorProb)]
write.table(resTblForRepro, file="interimData/tableForClinicalRepro.txt", sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

resTbl[,age:=2020-YearBirth]


ggplot(resTbl[Gender=="M",], aes(x=age)) +geom_histogram(bins=30)
ggsave("Figs/uszAgeM.pdf")
ggplot(resTbl[Gender=="W",], aes(x=age)) +geom_histogram(bins=30)
ggsave("Figs/uszAgeF.pdf")


	
#resTbl[posteriorProb>0.5 & organization=="Notfall",]


orgTbl0=fread("Data/Organisationseinheit.csv")
orgTbl=fread("Data/OrganisationseinheitUpdate.csv")[,list(V1,V2,V3)]
setnames(orgTbl,c("org","organizationDeutsch","organization"))
orgTbl[,Orgfa:=as.numeric(sub("B","",org))]
setkey(resTbl,Orgfa)
setkey(orgTbl,Orgfa)
resTbl2=merge(resTbl,orgTbl)
resTbl3=resTbl2[,list(probM=sum(posteriorProb>0.5),nrOfSamples=length(age),organization=organization[1]),by=Orgfa]
setkey(resTbl3,nrOfSamples)

ggplot(resTbl2,aes(x=organization,y=posteriorProb,shape=Gender))+geom_point()


resTbl2[,myorg:=sub("/","",sub("\x8a","ae",organization))]
oo=names(rev(sort(table(resTbl2[,myorg]))))

resTbl2[,Organization:=factor(myorg,levels=oo)]

resTbl3=resTbl2[,list(Organization,posteriorProb)]
ggplot(resTbl3,aes(x=Organization)) + geom_bar( stat = "count") + scale_y_log10() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("Figs/clinicalTbl_cohort.pdf")

ggplot(resTbl3[posteriorProb>0.5,],aes(x=Organization)) + geom_bar( stat = "count") + scale_y_log10() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("Figs/clinicalTbl_cohort2.pdf")

resTbl3[,rr:=1]
resTbl4=resTbl3[,list(caseNr=length(rr)),by=list(Organization,posteriorProb>0.5)]


resTbl4=resTbl3[,list(`Nr of total Patients / 20`=(sum(posteriorProb>=0.0)-sum(posteriorProb>0.5)*20)/20,`Nr of corona Patients`=sum(posteriorProb>0.5)),by=Organization]
require(reshape2)
mm=as.data.table(melt(resTbl4,id.vars="Organization"))

myblue=rgb(59,109,156,maxColorValue= 255)
mylightblue=rgb(96,158,199,maxColorValue= 255)

ggplot(mm,aes(x=Organization,y=value,fill=variable)) + geom_bar( stat = "identity") +scale_fill_manual(values=c(myblue,mylightblue)) + theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position = c(0.8, 0.8))
ggsave("Figs/clinicalTbl_cohort3.pdf",width=4.5,height=4.5)


browser()


resTbl4_1=resTbl3[,list(caseNr=length(rr)/50),by=Organization]
resTbl4_2=resTbl3[posteriorProb>0.5,list(caseNr=length(rr)),by=Organization]


resTbl4_1[,dataType="nr of total Patients / 50"]
resTbl4_2[,dataType="nr of corona Patients"]

ggplot(resTbl4,aes(x=Organization,y=caseNr,fill=posteriorProb)) + geom_rect( stat = "identity") + scale_y_log10() +  theme(axis.text.x = element_text(angle = 90, hjust = 1))
# #ggsave("Figs/clinicalTbl_cohort2.pdf")
# myblue=rgb(59,109,156,maxColourValue= 255)
# mylightblue=rgb(96,158,199,maxColourValue= 255)

# R:96,G:158,B:199

# [16:23] Raphael Jacquat (Guest)
#     J005040742    0.998500
# J005061096    1.000000
# J005110942    0.999696
# J005130632    0.993047
# J005140441    0.726274





