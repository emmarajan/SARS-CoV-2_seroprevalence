source("Code/evaluateHistoricalAssayCorrelation.R")


source("Code/bds_cohort_combined.R")


source("Code/pickSamplesForSecondaryScreen.R")



############################################################
############################################################

source("Code/evaluateHistoricalAssayCorrelation.R")


## out: Figs/usz_cohort/estimated_true_positives_analysis",("_uszBackg")[USEUSZBACKGROUND],".pdf
source("Code/define_cutoffs.R")


source("Code/pickSamplesForSecondaryScreen.R")

##output:
##
## Figs/Figs/kim_cohort/spike_vs_time.pdf
source("Code/prepare_kim_cohort_data.R")


source("Code/prepare_kim_cohort_plots.R")


source("Code/bds_cohort.R")



## out: Figs/compareRocs.pdf
source("Code/compareRocVals.R")


##ouput: 
##	
## Figs/KIMcommercialComparison.pdf:Figs/Preprint/Fig1C.pdf
## Figs/KIMinternalComparison.pdf:Figs/Preprint/Fig1B.pdf
source("Code/kimComparisonToComercial.R")


##ouput: 
##	
## Figs/cohort_corr_Prevalence_USZ.pdf:Figs/Preprint/Fig2A.pdf
## Figs/cohort_corr_Prevalence_Blood.pdf:Figs/Preprint/Fig2B.pdf
##	
## Figs/cohort_corr_Prevalence_USZ.pdf:Figs/Preprint/S1FigC.pdf
## Figs/cohort_corr_Prevalence_Blood.pdf:Figs/Preprint/S1FigD.pdf
##
## Figs/cohort_combined_RocData_Above14_edit.pdf:Figs/Preprint/S1FigE.pdf
## Figs/cohort_combined_RocData_Blood_edit.pdf:Figs/Preprint/S1FigF.pdf
##
## Figs/cohort_corr_Prevalence_Mix_USZ.pdf:Figs/Preprint/S1FigC
## Figs/cohort_corr_Prevalence_Mix_Blood.pdf:Figs/Preprint/S1FigD.pdf
source("Code/bds_cohort_combined.R")

##  Code/rescreenRes.R
##
##  output:
##  Figs/replication/replicationReplication.pdf:Figs/Preprint/S1FigG.pdf
##  Figs/SpikeReplication.pdf:Figs/Preprint/Fig2D.pdf
##	Figs/compareTopVals.pdf:Figs/Preprint/Fig4A.pdf
##  Figs/rbd_vs_sars.pdf:Figs/Preprint/Fig2E.pdf
source("Code/rescreenRes.R")

