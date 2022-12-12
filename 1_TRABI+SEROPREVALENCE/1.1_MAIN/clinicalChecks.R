require(data.table)

source("Code/prepareDataHelper.R")

tbl=fread("Data/20200530_descrp_stats_full_list_2.csv")
tbl[,Patient_or_Control_ID:=PatientIDList]

uszMeasured=fread("interimData/bds_cohort/cohort_combined_outTbl_Above14.txt")
uszMeasured[,Patient_or_Control_ID:=Pseudo_ID]
setkey(uszMeasured,Patient_or_Control_ID)
setkey(tbl,Patient_or_Control_ID)
tbl=merge(tbl,uszMeasured)




tbl=extractDate(tbl)
tbl[,yearUpd:=year+2000]
tbl[,totalDays:=yearUpd*365+ (month-1)*31 + day]

tbl[,dateFromEindat:=gsub("\\-","",sub("(.+) .+","\\1",Eindat))]
tbl[,dateFromEindatNum:=as.numeric(dateFromEindat)]

tbl[,yearEindat:=floor(dateFromEindatNum/10000)]
tbl[,monthEindat:=floor((dateFromEindatNum - yearEindat*10000)/100)]
tbl[,dayEindat:= dateFromEindatNum %% 100]


tbl[,totalDaysEindat:=yearEindat*365+ (monthEindat-1)*31 + dayEindat]

tbl[,diffDays:=totalDays-totalDaysEindat]
tbl2=tbl[diffDays>=0,]



myf=function(x,y){y[which.min(x)]}
#withNegRemoved=tbl2[,list(bestfittingOrf=myf(diffDays,Orgfa),bestfittingDep=myf(diffDays,Department)), by=Patient_or_Control_ID]
withNegRemoved=tbl2[,list(bestfittingOrf=myf(diffDays,Orgfa)), by=Patient_or_Control_ID]
myf=function(x,y){y[which.min(x)]}
#keepingNegRemoved=tbl[,list(bestfittingOrf=myf(abs(diffDays),Orgfa),bestfittingDep=myf(abs(diffDays),Department)), by=Patient_or_Control_ID]
keepingNegRemoved=tbl[,list(bestfittingOrf=myf(abs(diffDays),Orgfa)), by=Patient_or_Control_ID]


write.table(withNegRemoved,file="interimData/withNegRemoved.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)



aa=fread("Data/Organisationseinheit.csv")
setnames(aa,c("OrgfaStr","Department"))
aa[,Orgfa:=as.numeric(OrgfaStr)]
setkey(aa,Orgfa)
setkey(tbl,Orgfa)
tbl=merge(tbl,aa)





outTbl3=merge(withNegRemoved,uszMeasured)
outTbl4=merge(keepingNegRemoved,uszMeasured)

sort(table(outTbl4[,bestfittingDep]))

##outTbl[]
#outTbl2[bestfittingOrf!=(Orgfa),]





