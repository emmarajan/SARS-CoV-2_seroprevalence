require(data.table)
####################################################################
require(data.table)
require(ggplot2)


oxford=fread("interimData/oxford_cohort/oxford_cohort_inputData.txt")
oxford=oxford[category=="convalescent positive",]
kim=fread("interimData/kim_cohort/kim_cohort_RocData.R")
bds=fread("interimData/bds_cohort/bds_cohort_RocData.R")


setkey(kim, fpr)
setkey(oxford, fpr)
fprList=c(0, 0.001,0.005,0.01)
myfpr=fprList[1]

myfpr=0.001
wwl=lapply(fprList,function(myfpr){
	xl=lapply(unique(kim[,target]),function(mytarget){
			tail(kim[fpr<=myfpr & target==mytarget,],1)[,list(tpr,target,fpr=myfpr,cohort="kim")]
		})
	x=do.call("rbind",xl)
	yl=lapply(unique(oxford[,target]),function(mytarget){
			tail(oxford[fpr<=myfpr & target==mytarget,],1)[,list(tpr,target,fpr=myfpr,cohort="oxford")]
		})
	y=do.call("rbind",yl)
	zl=lapply(unique(bds[,target]),function(mytarget){
		tail(bds[fpr<=myfpr & target==mytarget,],1)[,list(tpr,target,fpr=myfpr,cohort="bds")]
	})
	z=do.call("rbind",zl)
	return(rbind(rbind(x,y),z))
})
ww=do.call("rbind",wwl)

ggplot(ww, aes(x=fpr,y=tpr, group=target,fill=target)) +geom_bar(colour="black",stat="identity",position="dodge")+ facet_grid(cohort~.) +theme_bw()
ggsave("Figs/compareRocs.pdf")

