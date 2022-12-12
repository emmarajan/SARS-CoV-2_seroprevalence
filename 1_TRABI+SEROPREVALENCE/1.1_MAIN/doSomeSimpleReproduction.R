
##############################################################
## parse data
##############################################################
require(data.table)
require(ggplot2)
source("Code/prepareDataHelper.R")

### loading summary data.
dataTbl=extractDate(loadPath())


###############################################################
## reproduce IC50 computation: (needs Kd and Ab concentration)
##
## note: thautan.py line 144: IC50 = (1 + 2 * self.kd_norm) / (2 * self.a0_norm)
## Kd => kd_norm
## Ab => a0_norm
###############################################################

myIC50 = dataTbl[,( 1 + 2 * KD) / (2 * Ab_conc)]

print("manually computed IC50. Difference should be very close to zero")
print(sum((myIC50 - dataTbl[,IC50])^2) / sum(myIC50^2))
print("looks ok")

###############################################################
## reproduce neglogEC50 computation: looks like IC50 and EC50 are used interchangeably
##
## neglogEC50 = -log10(IC50)
###############################################################

print("manually computed neglogEC50. Difference should be very close to zero")
 dataTbl[, diffneglogEC50:=(-log10(IC50) - neglogEC50)]
print("most close to zero")
print("those that aren't, are very negative i.e IC50 close to 0. Numeric instability? it doesn't really matter in practice.")
max(dataTbl[abs(diffneglogEC50)>1E-3,neglogEC50])


###############################################################
## i still don't quite understand the difference between neglogEC50 and corrlogEC50
##
## corrlogEC50 seems to provide a hard floor on low values (see second plot for example). 
## the floor is also given in the value neglogmaxfit. 
## however, in the data I have, changes only happen well below a reasonable positive signal. (see first plot) So it is probably not that 
## crucial at the moment.

###############################################################
dataTbl[,plot(corrlogEC50,neglogEC50)]
abline(0,1)
dataTbl[Plate_ID=="D200422d62254713031",plot(corrlogEC50,neglogEC50)]
abline(0,1)

shouldBeTrue = dataTbl[Plate_ID=="D200422d62254713031",min(corrlogEC50)]==dataTbl[Plate_ID=="D200422d62254713031",neglogmaxfit][1]
print(shouldBeTrue)
###############################################################
## Explore fit: show the distribution of signal and plot fit over actual data to see how the fitted curves
## look compared to normal data.
###############################################################

drawCurves = function(A0_norm, KD_norm, concentrations){
	sqrtTerm = sqrt((A0_norm * concentrations + KD_norm)^2 + 2 * (KD_norm - A0_norm * concentrations) + 1)                                                                            
    conc_bound = 0.5 * (A0_norm * concentrations + KD_norm + 1 - sqrtTerm)
    browser()
    return(conc_bound)
}

linearStretch = function(vals, baseline=0.2,plateau=2){
	(plateau - baseline)*vals + baseline
}


###############################################################
## explore fit.
## load some signal data for some plate.  

## base line and plateau seem to be hard coded.
## baseline=0.2, plateau=2
###############################################################
### loading example signal data.
signalTbl = fread("Data/SARS_CoV_2_share/20200422_1_Spike_IgG/Analysis/20200422_1_Spike_IgG_merged_data_from_hts.csv")

baseline=0.2
plateau=2

## merge signal data for some example plate with the fit data.
setkeyv(signalTbl,c("patient_or_control_id", "target"))
dataTbl[,patient_or_control_id:=Patient_or_Control_ID]
setkeyv(dataTbl,c("patient_or_control_id", "target"))
mergeTbl = merge(dataTbl,signalTbl)
mergeTblSpike = mergeTbl[antigen == "Spike",]

## pick sample with typical curve. (the one with  -log10(IC50) as close as possible to 2)
selected_id = mergeTblSpike[which.min(abs(neglogEC50-2)),patient_or_control_id]
examplePatientTbl = mergeTblSpike[ selected_id ==  patient_or_control_id,]
examplePatientTbl[,plot(log(concentration),signal)]

vals = examplePatientTbl[1,c(Ab_conc,KD)]
concs = examplePatientTbl[,concentration]
signal = examplePatientTbl[,signal]
predicted = drawCurves(vals[1], vals[2], concs)

plot(log10(concs),linearStretch(predicted))

newcons = c(concs,max(concs)*2,max(concs)*5,max(concs)*10)
predictedUpdated = drawCurves(vals[1], vals[2], newcons)

measuredTbl = data.table(logconc = log10(newcons), val = linearStretch(predictedUpdated),type="predicted")
predTbl = data.table(logconc = log10(concs), val = signal,type="measured")
plotTbl = rbind(measuredTbl, predTbl)
ggplot(plotTbl, aes(x = logconc, y = val, colour = type))+geom_point()
ggsave("Figs/")





x=2

sqrtTerm = sqrt((x + x)^2 + 2 * (x - x) + 1)                                                                            
conc_bound = 0.5 * (x + x + 1 - sqrtTerm)




