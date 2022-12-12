########## EMMENEGGER ET AL. 2022, POST COVID-19 CONDITION ANALYSIS ##########

### PRESETS ---------------------------------------------------------------------- 

## Packages
library(openxlsx)
library(tidyverse)
library(lubridate)
library(ggpubr)
library(RColorBrewer)
library(table1)

## Presets
col <- brewer.pal(n=11, name="RdBu")
render.cont <- function(x) { with(stats.default(x), c("", "Mean (SD)" = sprintf("%0.1f (%0.1f)", MEAN, SD), "Median (IQR)" = sprintf("%0.1f (%0.1f to %0.1f)", MEDIAN, Q1, Q3), "Range" = sprintf("%0.0f to %0.0f", MIN, MAX))) }
render.cont.eq5d <- function(x) { with(stats.default(x), c("", "Mean (SD)" = sprintf("%0.2f (%0.2f)", MEAN, SD), "Median (IQR)" = sprintf("%0.2f (%0.2f to %0.2f)", MEDIAN, Q1, Q3), "Range" = sprintf("%0.2f to %0.2f", MIN, MAX))) }
render.cat <- function(x) { c("", sapply(stats.default(x), function(y) with(y, sprintf("%d (%0.1f%%)", FREQ, PCTnoNA)))) }


### IMPORT & PROCESS DATA ---------------------------------------------------------------------- 

## Read data
data <- read.csv("20220603_SARSCOV2_questionnaire_compact_reduced.csv", na.strings = c("NA", ""))

## Process data
data$cov2_survey_timestamp <- as.Date(data$cov2_survey_timestamp, format = "%d/%m/%Y %H:%M")

data$ever_infected <- recode_factor(ifelse(data$de_or_en == 0, data$q22, data$q1_engl_2), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$ever_infected_2l <- recode_factor(data$ever_infected, "Don't know" = "No")
data$infected_2020 <- recode_factor(ifelse(data$de_or_en == 0, data$q1, data$q1_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$first_date <- as.Date(ifelse(data$de_or_en == 0, data$q2, data$q2_engl), format = "%d/%m/%Y")
data$ever_reinfected <- recode_factor(ifelse(data$de_or_en == 0, data$q24, data$q1_engl_3), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$ever_reinfected_2l <- recode_factor(data$ever_reinfected, "Don't know" = "No")
data$last_infect_date <- as.Date(ifelse(data$de_or_en == 0, data$q25, data$q2_engl_2), format = "%d/%m/%Y")

data$hospital <- recode_factor(ifelse(data$de_or_en == 0, data$q3, data$q3_engl), "1" = "Yes", "0" = "No")
data$supp_oxygen <- recode_factor(ifelse(data$de_or_en == 0, data$q4, data$q4_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$icu <- recode_factor(ifelse(data$de_or_en == 0, data$q5, data$q5_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$mech_vent <- recode_factor(ifelse(data$de_or_en == 0, data$q6, data$q6_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")

data$symp_at_infection <- recode_factor(ifelse(data$de_or_en == 0, data$q7, data$q7_engl), "1" = "Yes", "0" = "No")
data$symp_fever <- recode_factor(ifelse(data$de_or_en == 0, data$q8___1, data$q8_engl___1), "1" = "Yes", "0" = "No")
data$symp_fatigue <- recode_factor(ifelse(data$de_or_en == 0, data$q8___2, data$q8_engl___2), "1" = "Yes", "0" = "No")
data$symp_headache <- recode_factor(ifelse(data$de_or_en == 0, data$q8___3, data$q8_engl___3), "1" = "Yes", "0" = "No")
data$symp_cough <- recode_factor(ifelse(data$de_or_en == 0, data$q8___4, data$q8_engl___4), "1" = "Yes", "0" = "No")
data$symp_dyspnea <- recode_factor(ifelse(data$de_or_en == 0, data$q8___5, data$q8_engl___5), "1" = "Yes", "0" = "No")
data$symp_chestpain <- recode_factor(ifelse(data$de_or_en == 0, data$q8___6, data$q8_engl___6), "1" = "Yes", "0" = "No")
data$symp_runnynose <- recode_factor(ifelse(data$de_or_en == 0, data$q8___7, data$q8_engl___7), "1" = "Yes", "0" = "No")
data$symp_irreyes <- recode_factor(ifelse(data$de_or_en == 0, data$q8___8, data$q8_engl___8), "1" = "Yes", "0" = "No")
data$symp_sorethroat <- recode_factor(ifelse(data$de_or_en == 0, data$q8___9, data$q8_engl___9), "1" = "Yes", "0" = "No")
data$symp_anosmia <- recode_factor(ifelse(data$de_or_en == 0, data$q8___10, data$q8_engl___10), "1" = "Yes", "0" = "No")
data$symp_diarrhea <- recode_factor(ifelse(data$de_or_en == 0, data$q8___11, data$q8_engl___11), "1" = "Yes", "0" = "No")
data$symp_nausea <- recode_factor(ifelse(data$de_or_en == 0, data$q8___12, data$q8_engl___12), "1" = "Yes", "0" = "No")
data$symp_abdominalpain <- recode_factor(ifelse(data$de_or_en == 0, data$q8___13, data$q8_engl___13), "1" = "Yes", "0" = "No")
data$symp_appetiteloss <- recode_factor(ifelse(data$de_or_en == 0, data$q8___14, data$q8_engl___14), "1" = "Yes", "0" = "No")
data$symp_myalgia <- recode_factor(ifelse(data$de_or_en == 0, data$q8___15, data$q8_engl___15), "1" = "Yes", "0" = "No")
data$symp_arthralgia <- recode_factor(ifelse(data$de_or_en == 0, data$q8___16, data$q8_engl___16), "1" = "Yes", "0" = "No")
data$symp_rash <- recode_factor(ifelse(data$de_or_en == 0, data$q8___17, data$q8_engl___17), "1" = "Yes", "0" = "No")
data$symp_tingling <- recode_factor(ifelse(data$de_or_en == 0, data$q8___18, data$q8_engl___18), "1" = "Yes", "0" = "No")
data$symp_other <- recode_factor(ifelse(data$de_or_en == 0, data$q8___19, data$q8_engl___19), "1" = "Yes", "0" = "No")
data$symp_other_txt <- ifelse(data$de_or_en == 0, data$q9, data$q9_engl)

# data$symp_start_date <- ifelse(data$de_or_en == 0, data$q10, data$q10_engl) ## question missing from german questionnaire
data$symp_duration <- recode_factor(ifelse(data$de_or_en == 0, data$q11, data$q11_engl), "1" = "<2 weeks", "2" = "2-4 weeks", "3" = "1-3 months", "4" = "3-6 months", "5" = ">6 months", "6" = "Ongoing")
data$diag_pneumonia <- recode_factor(ifelse(data$de_or_en == 0, data$q12, data$q12_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")

data$comorb_heartdisease <- recode_factor(ifelse(data$de_or_en == 0, data$q13___1, data$q13_engl___1), "1" = "Yes", "0" = "No")
data$comorb_stroke <- recode_factor(ifelse(data$de_or_en == 0, data$q13___2, data$q13_engl___2), "1" = "Yes", "0" = "No")
data$comorb_cancer <- recode_factor(ifelse(data$de_or_en == 0, data$q13___3, data$q13_engl___3), "1" = "Yes", "0" = "No")
data$comorb_depression <- recode_factor(ifelse(data$de_or_en == 0, data$q13___4, data$q13_engl___4), "1" = "Yes", "0" = "No")
data$comorb_diabetes <- recode_factor(ifelse(data$de_or_en == 0, data$q13___5, data$q13_engl___5), "1" = "Yes", "0" = "No")
data$comorb_hypertension <- recode_factor(ifelse(data$de_or_en == 0, data$q13___6, data$q13_engl___6), "1" = "Yes", "0" = "No")
data$comorb_arthritis <- recode_factor(ifelse(data$de_or_en == 0, data$q13___7, data$q13_engl___7), "1" = "Yes", "0" = "No")
data$comorb_osteoporosis <- recode_factor(ifelse(data$de_or_en == 0, data$q13___8, data$q13_engl___8), "1" = "Yes", "0" = "No")
data$comorb_asthma <- recode_factor(ifelse(data$de_or_en == 0, data$q13___9, data$q13_engl___9), "1" = "Yes", "0" = "No")
data$comorb_copd <- recode_factor(ifelse(data$de_or_en == 0, data$q13___10, data$q13_engl___10), "1" = "Yes", "0" = "No")
data$comorb_hepatitis <- recode_factor(ifelse(data$de_or_en == 0, data$q13___11, data$q13_engl___11), "1" = "Yes", "0" = "No")
data$comorb_ibd <- recode_factor(ifelse(data$de_or_en == 0, data$q13___12, data$q13_engl___12), "1" = "Yes", "0" = "No")
data$comorb_ckd <- recode_factor(ifelse(data$de_or_en == 0, data$q13___13, data$q13_engl___13), "1" = "Yes", "0" = "No")
data$comorb_pollenallergy <- recode_factor(ifelse(data$de_or_en == 0, data$q13___14, data$q13_engl___14), "1" = "Yes", "0" = "No")
data$comorb_foodallergy <- recode_factor(ifelse(data$de_or_en == 0, data$q13___15, data$q13_engl___15), "1" = "Yes", "0" = "No")
data$comorb_otherallergy <- recode_factor(ifelse(data$de_or_en == 0, data$q13___16, data$q13_engl___16), "1" = "Yes", "0" = "No")
data$comorb_autoimmune <- recode_factor(ifelse(data$de_or_en == 0, data$q13___17, data$q13_engl___17), "1" = "Yes", "0" = "No")
data$comorb_atopicderm <- recode_factor(ifelse(data$de_or_en == 0, data$q13___18, data$q13_engl___18), "1" = "Yes", "0" = "No")
data$comorb_migraine <- recode_factor(ifelse(data$de_or_en == 0, data$q13___19, data$q13_engl___19), "1" = "Yes", "0" = "No")
data$comorb_cysticfibrosis <- recode_factor(ifelse(data$de_or_en == 0, data$q13___20, data$q13_engl___20), "1" = "Yes", "0" = "No")
data$comorb_hiv <- recode_factor(ifelse(data$de_or_en == 0, data$q13___21, data$q13_engl___21), "1" = "Yes", "0" = "No")
data$comorb_immunsuppression <- recode_factor(ifelse(data$de_or_en == 0, data$q13___22, data$q13_engl___22), "1" = "Yes", "0" = "No")
data$comorb_other <- recode_factor(ifelse(data$de_or_en == 0, data$q13___23, data$q13_engl___23), "1" = "Yes", "0" = "No")

data$comorb_otherallergy_txt <- ifelse(data$de_or_en == 0, data$q14, data$q14_engl)
data$comorb_other_txt <- ifelse(data$de_or_en == 0, data$q15, data$q15_engl)
data$height <- ifelse(data$de_or_en == 0, data$q16, data$q16_engl)
data$weight <- ifelse(data$de_or_en == 0, data$q17, data$q17_engl)
data$smoking <- recode_factor(ifelse(data$de_or_en == 0, data$q18, data$q18_engl), "1" = "Yes, daily", "2" = "Yes, sometimes", "3" = "Ex-smoker", "4" = "Non-smoker")
data$alcohol <- recode_factor(ifelse(data$de_or_en == 0, data$q19, data$q19_engl), "1" = "Never", "2" = "Once per month or less", "3" = "2-3 times per month", "4" = "2-3 times per week", "5" = "4 times or more per week")
                              
data$recovery <- recode_factor(ifelse(data$de_or_en == 0, data$q20, data$q20_engl), "1" = "Recovered", "2" = "Better, but not fully recovered", "3" = "Neither better nor worse", "4" = "Worse")
data$lc_symp_fatigue <- recode_factor(ifelse(data$de_or_en == 0, data$q21___1, data$q21_engl___1), "1" = "Yes", "0" = "No")
data$lc_symp_sleeping <- recode_factor(ifelse(data$de_or_en == 0, data$q21___2, data$q21_engl___2), "1" = "Yes", "0" = "No")
data$lc_symp_headache <- recode_factor(ifelse(data$de_or_en == 0, data$q21___3, data$q21_engl___3), "1" = "Yes", "0" = "No")
data$lc_symp_appetiteloss <- recode_factor(ifelse(data$de_or_en == 0, data$q21___4, data$q21_engl___4), "1" = "Yes", "0" = "No")
data$lc_symp_concentration <- recode_factor(ifelse(data$de_or_en == 0, data$q21___5, data$q21_engl___5), "1" = "Yes", "0" = "No")
data$lc_symp_memory <- recode_factor(ifelse(data$de_or_en == 0, data$q21___6, data$q21_engl___6), "1" = "Yes", "0" = "No")
data$lc_symp_nervousness <- recode_factor(ifelse(data$de_or_en == 0, data$q21___7, data$q21_engl___7), "1" = "Yes", "0" = "No")
data$lc_symp_depression <- recode_factor(ifelse(data$de_or_en == 0, data$q21___8, data$q21_engl___8), "1" = "Yes", "0" = "No")
data$lc_symp_rash <- recode_factor(ifelse(data$de_or_en == 0, data$q21___9, data$q21_engl___9), "1" = "Yes", "0" = "No")
data$lc_symp_skinproblems <- recode_factor(ifelse(data$de_or_en == 0, data$q21___10, data$q21_engl___10), "1" = "Yes", "0" = "No")
data$lc_symp_gisymptoms <- recode_factor(ifelse(data$de_or_en == 0, data$q21___11, data$q21_engl___11), "1" = "Yes", "0" = "No")
data$lc_symp_hairloss <- recode_factor(ifelse(data$de_or_en == 0, data$q21___12, data$q21_engl___12), "1" = "Yes", "0" = "No")
data$lc_symp_reducedperformance <- recode_factor(ifelse(data$de_or_en == 0, data$q21___13, data$q21_engl___13), "1" = "Yes", "0" = "No")
data$lc_symp_palpitations <- recode_factor(ifelse(data$de_or_en == 0, data$q21___14, data$q21_engl___14), "1" = "Yes", "0" = "No")
data$lc_symp_dyspnea <- recode_factor(ifelse(data$de_or_en == 0, data$q21___15, data$q21_engl___15), "1" = "Yes", "0" = "No")
data$lc_symp_cough <- recode_factor(ifelse(data$de_or_en == 0, data$q21___16, data$q21_engl___16), "1" = "Yes", "0" = "No")
data$lc_symp_myalgia <- recode_factor(ifelse(data$de_or_en == 0, data$q21___17, data$q21_engl___17), "1" = "Yes", "0" = "No")
data$lc_symp_arthralgia <- recode_factor(ifelse(data$de_or_en == 0, data$q21___18, data$q21_engl___18), "1" = "Yes", "0" = "No")
data$lc_symp_fever <- recode_factor(ifelse(data$de_or_en == 0, data$q21___19, data$q21_engl___19), "1" = "Yes", "0" = "No")
data$lc_symp_tingling <- recode_factor(ifelse(data$de_or_en == 0, data$q21___20, data$q21_engl___20), "1" = "Yes", "0" = "No")
data$lc_symp_tremor <- recode_factor(ifelse(data$de_or_en == 0, data$q21___21, data$q21_engl___21), "1" = "Yes", "0" = "No")
data$lc_symp_hearing <- recode_factor(ifelse(data$de_or_en == 0, data$q21___22, data$q21_engl___22), "1" = "Yes", "0" = "No")
data$lc_symp_vertigo <- recode_factor(ifelse(data$de_or_en == 0, data$q21___23, data$q21_engl___23), "1" = "Yes", "0" = "No")
data$lc_symp_vision <- recode_factor(ifelse(data$de_or_en == 0, data$q21___24, data$q21_engl___24), "1" = "Yes", "0" = "No")

data$diag_lc <- recode_factor(ifelse(data$de_or_en == 0, data$t1, data$t1_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$diag_lung <- recode_factor(ifelse(data$de_or_en == 0, data$t2, data$t2_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$diag_heart <- recode_factor(ifelse(data$de_or_en == 0, data$t3, data$t3_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$diag_stroke <- recode_factor(ifelse(data$de_or_en == 0, data$t4, data$t4_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$diag_kidney <- recode_factor(ifelse(data$de_or_en == 0, data$t5, data$t5_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$diag_liver <- recode_factor(ifelse(data$de_or_en == 0, data$t6, data$t6_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$diag_gi <- recode_factor(ifelse(data$de_or_en == 0, data$t7, data$t7_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$diag_coagulation <- recode_factor(ifelse(data$de_or_en == 0, data$t8, data$t8_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$diag_immunesystem <- recode_factor(ifelse(data$de_or_en == 0, data$t9, data$t9_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$diag_bloodother <- recode_factor(ifelse(data$de_or_en == 0, data$t10, data$t10_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$diag_skin <- recode_factor(ifelse(data$de_or_en == 0, data$t11, data$t11_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$diag_thyroid <- recode_factor(ifelse(data$de_or_en == 0, data$t12, data$t12_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$diag_brain <- recode_factor(ifelse(data$de_or_en == 0, data$t13, data$t13_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$diag_nervefunction <- recode_factor(ifelse(data$de_or_en == 0, data$t14, data$t14_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$diag_other <- recode_factor(ifelse(data$de_or_en == 0, data$t15, data$t15_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")
data$diag_other_txt <- recode_factor(ifelse(data$de_or_en == 0, data$q23, data$q23_engl), "1" = "Yes", "0" = "No", "2" = "Don't know")

data$eq5d_mb <- ifelse(data$de_or_en == 0, data$eq5d_mb_5l_swi_ger, data$eq5d_mb_5l_uk_eng)
data$eq5d_sc <- ifelse(data$de_or_en == 0, data$eq5d_sc_5l_swi_ger, data$eq5d_sc_5l_uk_eng)
data$eq5d_ua <- ifelse(data$de_or_en == 0, data$eq5d_ua_5l_swi_ger, data$eq5d_ua_5l_uk_eng)
data$eq5d_pd <- ifelse(data$de_or_en == 0, data$eq5d_pd_5l_swi_ger, data$eq5d_pd_5l_uk_eng)
data$eq5d_ad <- ifelse(data$de_or_en == 0, data$eq5d_ad_5l_swi_ger, data$eq5d_ad_5l_uk_eng)
data$eq_vas <- ifelse(data$de_or_en == 0, data$eq5d5l_vas2_swi_ger, data$eq5d5l_vas2_uk_eng)

data$comments <- ifelse(data$de_or_en == 0, data$q26, data$q23_engl_2)

data$pcr_pos <- ifelse(data$pcr_pos == FALSE, "Yes", "No")

## Clean data
data$height[data$patient_id == "ZfJDiVyO3STP"] <- 173 # correct height data
data$first_date[data$patient_id == "cTnLDPiUSutO"] <- as_date("2020-10-01") # implausible date entered (date prior to start of pandemic, 2019-10-01), assumed one year off
data$first_date[data$patient_id == "v6smROrVNS27"] <- as_date("2022-03-16") # implausible date entered (date in future at completion of questionnaire, 2022-04-16), assumed one month off
data$first_date[data$patient_id == "zH9s8BYgOqcA"] <- as_date("2022-04-05") # implausible date entered (date in future at completion of questionnaire, 2022-04-16), assumed one month off
# data$serology_result_date[data$patient_id == "zGMniKeBmc0X"] <- as_date() # prepandemic sample? (was labelled as 18-00-00)

## Categorize seropositives & -negatives
data$serology_result <- factor(ifelse(data$posterior >= 0.5, "seropositive", "seronegative"))
data$serology_result_date <- as.Date(paste0("20", substr(data$month, 1, 2), "-", substr(data$month, 3, 4), "-", formatC(data$day, width = 2, format = "d", flag = "0")), format = "%Y-%m-%d")

## Define those reinfected prior to sampling
data$infect_prior_sampling <- factor(ifelse(!is.na(data$first_date) & data$first_date <= data$serology_result_date, "Yes", "No"), levels = c("Yes", "No"))

## Categorize known infection or seropositive
data$infect_or_seropos <- factor(ifelse(data$ever_infected_2l == "Yes" | data$serology_result == "seropositive", "Yes", "No"), levels = c("Yes", "No"))

## Categorize seropositive wihouth known infection as asymptomatic (alternative assumption)
data$symp_at_infection_alt <-  factor(case_when(data$ever_infected_2l == "Yes" | !is.na(data$symp_at_infection) ~ data$symp_at_infection, 
                                                data$ever_infected_2l == "No" & data$serology_result == "seropositive" & is.na(data$symp_at_infection) ~ factor(c("No"), levels = c("Yes", "No"))))

## Group by pandemic infection wave (or no infection)
cutoff_ss2020 <- as_date("2020-07-31")
cutoff_ws2021_22 <- as_date("2021-11-30")
cutoff_days <- 30

# Primary analysis: assumed any seropositive individual without reported first infection date was infected within prevalent wave, with those tested seropositive within the first 30 days assigned to the previous wave
data$infect_wave <- factor(case_when(data$ever_infected_2l == "No" & data$serology_result == "seronegative" ~ "Never infected",
                                     data$ever_infected_2l == "Yes" & data$first_date <= cutoff_ss2020 ~ "Spring/Summer 2020",
                                     data$ever_infected_2l == "Yes" & data$first_date > cutoff_ss2020 & data$first_date <= cutoff_ws2021_22 ~ "Fall/Winter 2020/2021",
                                     data$ever_infected_2l == "Yes" & data$first_date > cutoff_ws2021_22 ~ "Winter/Spring 2021/2022",
                                     data$ever_infected_2l == "No" & data$serology_result == "seropositive" & data$serology_result_date <= cutoff_ss2020 + cutoff_days ~ "Spring/Summer 2020",
                                     data$ever_infected_2l == "No" & data$serology_result == "seropositive" & data$serology_result_date > cutoff_ss2020 + cutoff_days & data$serology_result_date <= cutoff_ws2021_22 + cutoff_days ~ "Fall/Winter 2020/2021",
                                     data$ever_infected_2l == "No" & data$serology_result == "seropositive" & data$serology_result_date > cutoff_ws2021_22 + cutoff_days ~ "Winter/Spring 2021/2022",
                                     data$ever_infected_2l == "Yes" & data$serology_result == "seropositive" & is.na(data$first_date) & data$serology_result_date <= cutoff_ss2020 + cutoff_days ~ "Spring/Summer 2020",
                                     data$ever_infected_2l == "Yes" & data$serology_result == "seropositive" & is.na(data$first_date) & data$serology_result_date > cutoff_ss2020 + cutoff_days & data$serology_result_date <= cutoff_ws2021_22 + cutoff_days ~ "Fall/Winter 2020/2021",
                                     data$ever_infected_2l == "Yes" & data$serology_result == "seropositive" & is.na(data$first_date) & data$serology_result_date > cutoff_ws2021_22 + cutoff_days ~ "Winter/Spring 2021/2022",
                                     is.na(data$ever_infected_2l) & data$serology_result == "seropositive" & data$serology_result_date <= cutoff_ss2020 ~ "Spring/Summer 2020",
                                     is.na(data$ever_infected_2l) & data$serology_result == "seropositive" & data$serology_result_date > cutoff_ss2020 & data$serology_result_date <= cutoff_ws2021_22 ~ "Fall/Winter 2020/2021",
                                     is.na(data$ever_infected_2l) & data$serology_result == "seropositive" & data$serology_result_date > cutoff_ws2021_22 ~ "Winter/Spring 2021/2022",
                                     TRUE ~ NA_character_), levels = c("Never infected", "Spring/Summer 2020", "Fall/Winter 2020/2021", "Winter/Spring 2021/2022"))

# Sensitivity analysis: any seropositive individual without reported first infection date and non-unambiguous infection timing assigned NA
data$infect_wave_sens <- factor(case_when(data$ever_infected_2l == "No" & data$serology_result == "seronegative" ~ "Never infected", 
                                          data$ever_infected_2l == "Yes" & data$first_date <= cutoff_ss2020 ~ "Spring/Summer 2020", 
                                          data$ever_infected_2l == "Yes" & data$first_date > cutoff_ss2020 & data$first_date <= cutoff_ws2021_22 ~ "Fall/Winter 2020/2021", 
                                          data$ever_infected_2l == "Yes" & data$first_date > cutoff_ws2021_22 ~ "Winter/Spring 2021/2022", 
                                          data$ever_infected_2l == "No" & data$serology_result == "seropositive" & data$serology_result_date <= cutoff_ss2020 + cutoff_days ~ "Spring/Summer 2020", 
                                          data$ever_infected_2l == "No" & data$serology_result == "seropositive" & data$serology_result_date > cutoff_ss2020 + cutoff_days & data$serology_result_date <= cutoff_ws2021_22 + cutoff_days ~ NA_character_,
                                          data$ever_infected_2l == "No" & data$serology_result == "seropositive" & data$serology_result_date > cutoff_ws2021_22 + cutoff_days ~ NA_character_,
                                          data$ever_infected_2l == "Yes" & data$serology_result == "seropositive" & is.na(data$first_date) & data$serology_result_date <= cutoff_ss2020 + cutoff_days ~ "Spring/Summer 2020", 
                                          data$ever_infected_2l == "Yes" & data$serology_result == "seropositive" & is.na(data$first_date) & data$serology_result_date > cutoff_ss2020 + cutoff_days & data$serology_result_date <= cutoff_ws2021_22 + cutoff_days ~ NA_character_,
                                          data$ever_infected_2l == "Yes" & data$serology_result == "seropositive" & is.na(data$first_date) & data$serology_result_date > cutoff_ws2021_22 + cutoff_days ~ NA_character_,
                                          is.na(data$ever_infected_2l) & data$serology_result == "seropositive" & data$serology_result_date <= cutoff_ss2020 ~ "Spring/Summer 2020", 
                                          is.na(data$ever_infected_2l) & data$serology_result == "seropositive" & data$serology_result_date > cutoff_ss2020 & data$serology_result_date <= cutoff_ws2021_22 ~ NA_character_,
                                          is.na(data$ever_infected_2l) & data$serology_result == "seropositive" & data$serology_result_date > cutoff_ws2021_22 ~ NA_character_, 
                                          TRUE ~ NA_character_), levels = c("Never infected", "Spring/Summer 2020", "Fall/Winter 2020/2021", "Winter/Spring 2021/2022"))

# # manual check
# data %>% 
#   select(patient_id, ever_infected_2l, serology_result, infect_or_seropos, first_date, serology_result_date, infect_wave, infect_wave_sens, pcr_pos, days_since_PCR_pos) %>% 
#   # filter(infect_or_seropos == "Yes" & is.na(first_date)) %>% 
#   View()

# manually assign most likely infection timeframe (based on days_since_PCR_pos)
data$infect_wave_sens[data$patient_id == "NXItAP6l57Jb"] <- as.factor("Fall/Winter 2020/2021")
data$infect_wave_sens[data$patient_id == "YHDFNCPrxOnu"] <- as.factor("Fall/Winter 2020/2021")

## Finalize dataset
data <- data %>% 
  mutate(bmi = weight / (height/100)^2, 
         days_since_infection = cov2_survey_timestamp - first_date) %>% 
  select(patient_id, cov2_survey_timestamp, de_or_en, 
         sex, age, clinic, height, weight, bmi, 
         ever_infected, ever_infected_2l, infected_2020, first_date, infect_prior_sampling, days_since_infection, infect_wave, infect_wave_sens, infect_or_seropos, symp_at_infection_alt, ever_reinfected:comments, 
         pcr_pos, days_since_PCR_pos, 
         Spike_plogEC50, NC_plogEC50, RBD_plogEC50, posterior, serology_result, serology_result_date)


# ## Write dataset
# write.csv(data, "20221119_SARSCOV2_questionnaire_compact_processed.csv", row.names = FALSE)
# 
# ## Read dataset
# data <- read.csv("20221119_SARSCOV2_questionnaire_compact_processed.csv")


### POPULATION CHARACTERISTICS ---------------------------------------------------------------------- 

## Population characteristics
table1(~ age + sex + bmi + 
         comorb_heartdisease + comorb_stroke + comorb_cancer + comorb_diabetes + comorb_hypertension + 
         comorb_arthritis + comorb_osteoporosis + comorb_asthma + comorb_copd + comorb_hepatitis + comorb_ibd + comorb_ckd + 
         comorb_autoimmune + comorb_atopicderm + comorb_pollenallergy + comorb_foodallergy + comorb_otherallergy + 
         comorb_hiv + comorb_immunsuppression + 
         comorb_migraine + comorb_cysticfibrosis + 
         comorb_depression + comorb_other + 
         smoking + alcohol + 
         clinic,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = data)

## Average time since first infection
table1(~ days_since_infection
       | infect_prior_sampling,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = mutate(filter(data, ever_infected == "Yes"), days_since_infection = as.integer(days_since_infection)))

ggplot(data, aes(x = -days_since_infection)) + 
  geom_histogram(bins = 75) + 
  theme_bw()

ggplot(data, aes(x = first_date)) + 
  geom_histogram(bins = 75) + 
  theme_bw()

### SEROPOSITIVITY & REPORTED INFECTION ---------------------------------------------------------------------- 

## Seropositivity 
table(data$serology_result, useNA = "ifany")
prop.table(table(data$serology_result))

## Reported infection before serology testing
table(data$infect_prior_sampling, useNA = "ifany")
prop.table(table(data$infect_prior_sampling))

## Reported infection vs. seropositivity
table1(~ infect_prior_sampling + infected_2020 + ever_infected + ever_reinfected 
       | serology_result,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = data)

## Concordance reported infection before serology testing vs. seropositivity
# Cohen's Kappa
psych::cohen.kappa(table(recode_factor(data$infect_prior_sampling, "Yes" = "1", "No" = "0"), recode_factor(data$serology_result, "seropositive" = "1", "seronegative" = "0")))
# Percent concordance
sum(prop.table(table(data$infect_prior_sampling, data$serology_result))[1, 2], prop.table(table(data$infect_prior_sampling, data$serology_result))[2, 1])

### SEROPOSITIVITY & REPORTED SYMPTOMS ---------------------------------------------------------------------- 

## Reported symptoms at infection among those ever infected
table1(~ symp_at_infection,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = filter(data, ever_infected == "Yes"))

## Specific symptoms among all symptomatic infected
table1(~ symp_fever + symp_fatigue + symp_headache + symp_cough + symp_dyspnea + symp_chestpain + symp_runnynose + symp_irreyes + 
         symp_sorethroat + symp_anosmia + symp_diarrhea + symp_nausea + symp_abdominalpain + symp_appetiteloss + symp_myalgia + symp_arthralgia + 
         symp_rash + symp_tingling + symp_other, 
       render.continuous = render.cont, render.categorical = render.cat, 
       data = filter(data, ever_infected == "Yes" & symp_at_infection == "Yes"))

data %>% 
  filter(ever_infected == "Yes" & symp_at_infection == "Yes") %>% 
  pivot_longer(cols = symp_fever:symp_other, names_to = "symptom", values_to = "value") %>% 
  group_by(symptom, value) %>% 
  summarise(n = n()) %>% 
  mutate(ntot = sum(n), 
         freq = n / ntot) %>% 
  ungroup() %>% 
  filter(value == "Yes") %>% 
  mutate(symptom = fct_reorder(symptom, freq)) %>% 
  ggplot(aes(y = symptom, x = freq)) + 
  geom_bar(stat = "identity", width = 0.6, position = position_dodge(preserve = "single"), fill = col[c(11)]) + 
  theme_bw()

## Reported symptoms at infection among those infected before serology testing
table1(~ symp_at_infection
       | serology_result,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = filter(data, infect_prior_sampling == "Yes"))

## Specific symptoms among symptomatic infected before serology testing
table1(~ symp_fever + symp_fatigue + symp_headache + symp_cough + symp_dyspnea + symp_chestpain + symp_runnynose + symp_irreyes + 
         symp_sorethroat + symp_anosmia + symp_diarrhea + symp_nausea + symp_abdominalpain + symp_appetiteloss + symp_myalgia + symp_arthralgia + 
         symp_rash + symp_tingling + symp_other
       | serology_result,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = filter(data, infect_prior_sampling == "Yes" & symp_at_infection == "Yes"))

data %>% 
  filter(infect_prior_sampling == "Yes" & symp_at_infection == "Yes") %>% 
  pivot_longer(cols = symp_fever:symp_other, names_to = "symptom", values_to = "value") %>% 
  group_by(symptom, value) %>% 
  summarise(n = n()) %>% 
  mutate(ntot = sum(n), 
         freq = n / ntot) %>% 
  ungroup() %>% 
  filter(value == "Yes") %>% 
  mutate(symptom = fct_reorder(symptom, freq)) %>% 
  ggplot(aes(y = symptom, x = freq)) + 
  geom_bar(stat = "identity", width = 0.6, position = position_dodge(preserve = "single"), fill = col[c(11)]) + 
  theme_bw()

### SEVERITY OF INFECTION ---------------------------------------------------------------------- 

## Severity of infection among those ever infected, stratified by infection timing
table1(~ symp_at_infection + hospital + supp_oxygen + icu + mech_vent + diag_pneumonia
       | infect_wave,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = filter(data, ever_infected == "Yes"))

## Severity of infection among those ever infected, stratified by infection timing (sensitivity analysis)
table1(~ symp_at_infection + hospital + supp_oxygen + icu + mech_vent + diag_pneumonia
       | infect_wave_sens,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = filter(data, ever_infected == "Yes"))

### RECOVERY ---------------------------------------------------------------------- 

## Recovery by infection timing
table1(~ recovery
       | infect_wave,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = filter(data, ever_infected == "Yes"))

## Recovery by infection timing (sensitivity analysis)
table1(~ recovery
       | infect_wave_sens,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = filter(data, ever_infected == "Yes"))

## Symptom persistence among those with symptoms by infection timing
table1(~ symp_duration 
       | infect_wave,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = filter(data, ever_infected == "Yes" & symp_at_infection == "Yes"))

## Symptom persistence among those with symptoms by infection timing (sensitivity analysis)
table1(~ symp_duration 
       | infect_wave_sens,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = filter(data, ever_infected == "Yes" & symp_at_infection == "Yes"))

## Current symptoms by infection timing
table1(~ lc_symp_fatigue + lc_symp_sleeping + lc_symp_headache + lc_symp_appetiteloss + lc_symp_concentration + 
         lc_symp_memory + lc_symp_nervousness + lc_symp_depression + lc_symp_rash + lc_symp_skinproblems + 
         lc_symp_gisymptoms + lc_symp_hairloss + lc_symp_reducedperformance + lc_symp_palpitations + lc_symp_dyspnea + 
         lc_symp_cough + lc_symp_myalgia + lc_symp_arthralgia + lc_symp_fever + lc_symp_tingling + lc_symp_tremor + 
         lc_symp_hearing + lc_symp_vertigo + lc_symp_vision
       | infect_wave,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = filter(data, ever_infected == "Yes"))

data %>% 
  filter(ever_infected == "Yes") %>%
  pivot_longer(cols = lc_symp_fatigue:lc_symp_vision, names_to = "symptom", values_to = "value") %>% 
  group_by(infect_wave, symptom, value) %>% 
  summarise(n = n()) %>% 
  filter(!is.na(value)) %>% 
  mutate(ntot = sum(n), 
         freq = n / ntot) %>% 
  ungroup() %>% 
  filter(value == "Yes") %>% 
  mutate(symptom = fct_reorder(symptom, freq)) %>% 
  ggplot(aes(y = symptom, x = freq, fill = fct_rev(infect_wave))) + 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) + 
  scale_fill_manual(values = col[c(11,10,9)]) + 
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme_bw() + labs(fill = "Infection wave")

# sensitivity analysis
data %>% 
  filter(ever_infected == "Yes") %>%
  pivot_longer(cols = lc_symp_fatigue:lc_symp_vision, names_to = "symptom", values_to = "value") %>% 
  group_by(infect_wave_sens, symptom, value) %>% 
  summarise(n = n()) %>% 
  filter(!is.na(value)) %>% 
  mutate(ntot = sum(n), 
         freq = n / ntot) %>% 
  ungroup() %>% 
  filter(value == "Yes") %>% 
  mutate(symptom = fct_reorder(symptom, freq)) %>% 
  ggplot(aes(y = symptom, x = freq, fill = fct_rev(infect_wave_sens))) + 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) + 
  scale_fill_manual(values = col[c(11,10,9)]) + 
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme_bw() + labs(fill = "Infection wave")

## Diagnosis of long covid or associated complications by infection timing
table1(~ diag_lc + 
         diag_lung + diag_heart + diag_stroke + diag_kidney + diag_liver + diag_gi + 
         diag_coagulation + diag_immunesystem + diag_bloodother + diag_skin + diag_thyroid + 
         diag_brain + diag_nervefunction + diag_other
       | infect_wave,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = filter(data, ever_infected == "Yes"))

data %>% 
  filter(ever_infected == "Yes") %>%
  pivot_longer(cols = diag_lung:diag_other, names_to = "symptom", values_to = "value") %>% 
  group_by(infect_wave, symptom, value) %>% 
  summarise(n = n()) %>% 
  filter(!is.na(value)) %>%
  mutate(ntot = sum(n), 
         freq = n / ntot) %>% 
  ungroup() %>% 
  filter(value == "Yes") %>% 
  mutate(symptom = fct_reorder(symptom, freq)) %>% 
  ggplot(aes(y = symptom, x = freq, fill = fct_rev(infect_wave))) + 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) + 
  scale_fill_manual(values = col[c(11,10,9)]) + 
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme_bw() + labs(fill = "Infection wave")

# sensitivity analysis
data %>% 
  filter(ever_infected == "Yes") %>%
  pivot_longer(cols = diag_lung:diag_other, names_to = "symptom", values_to = "value") %>% 
  group_by(infect_wave_sens, symptom, value) %>% 
  summarise(n = n()) %>% 
  filter(!is.na(value)) %>%
  mutate(ntot = sum(n), 
         freq = n / ntot) %>% 
  ungroup() %>% 
  filter(value == "Yes") %>% 
  mutate(symptom = fct_reorder(symptom, freq)) %>% 
  ggplot(aes(y = symptom, x = freq, fill = fct_rev(infect_wave_sens))) + 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) + 
  scale_fill_manual(values = col[c(11,10,9)]) + 
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme_bw() + labs(fill = "Infection wave")

### SYMPTOMS & COMPLICATIONS - COMPARISONS WITH NON-INFECTED ---------------------------------------------------------------------- 

## Symptom prevalence compared to non-infected
table1(~ lc_symp_fatigue + lc_symp_sleeping + lc_symp_headache + lc_symp_appetiteloss + lc_symp_concentration + 
         lc_symp_memory + lc_symp_nervousness + lc_symp_depression + lc_symp_rash + lc_symp_skinproblems + 
         lc_symp_gisymptoms + lc_symp_hairloss + lc_symp_reducedperformance + lc_symp_palpitations + lc_symp_dyspnea + 
         lc_symp_cough + lc_symp_myalgia + lc_symp_arthralgia + lc_symp_fever + lc_symp_tingling + lc_symp_tremor + 
         lc_symp_hearing + lc_symp_vertigo + lc_symp_vision
       | ever_infected_2l,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = mutate(data, ever_infected_2l = ifelse(!is.na(ever_infected) & ever_infected == "Yes", "Yes", "No")))

data %>% 
  pivot_longer(cols = lc_symp_fatigue:lc_symp_vision, names_to = "symptom", values_to = "value") %>% 
  group_by(ever_infected_2l, symptom, value) %>% 
  summarise(n = n()) %>% 
  filter(!is.na(value)) %>% 
  mutate(ntot = sum(n), 
         freq = n / ntot) %>% 
  ungroup() %>% 
  filter(!is.na(ever_infected_2l), value == "Yes") %>% 
  mutate(symptom = fct_reorder(symptom, freq)) %>% 
  ggplot(aes(y = symptom, x = freq, fill = fct_rev(ever_infected_2l))) + 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) + 
  scale_fill_manual(values = col[c(9,11)]) + 
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme_bw() + labs(fill = "Ever infected")

broom::tidy(glm(fct_rev(lc_symp_fatigue) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_sleeping) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_reducedperformance) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_concentration) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_memory) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_arthralgia) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_dyspnea) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_skinproblems) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_myalgia) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_cough) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_rash) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_tingling) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_headache) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_vertigo) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_gisymptoms) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_palpitations) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_hearing) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_appetiteloss) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_tremor) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_nervousness) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_depression) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_hairloss) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_vision) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_fever) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)

## Symptom prevalence compared by infection status and infection wave
table1(~ lc_symp_fatigue + lc_symp_sleeping + lc_symp_headache + lc_symp_appetiteloss + lc_symp_concentration + 
         lc_symp_memory + lc_symp_nervousness + lc_symp_depression + lc_symp_rash + lc_symp_skinproblems + 
         lc_symp_gisymptoms + lc_symp_hairloss + lc_symp_reducedperformance + lc_symp_palpitations + lc_symp_dyspnea + 
         lc_symp_cough + lc_symp_myalgia + lc_symp_arthralgia + lc_symp_fever + lc_symp_tingling + lc_symp_tremor + 
         lc_symp_hearing + lc_symp_vertigo + lc_symp_vision
       | ever_infected_2l,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = mutate(data, ever_infected_2l = ifelse(!is.na(ever_infected) & ever_infected == "Yes", "Yes", "No")))

data %>% 
  pivot_longer(cols = lc_symp_fatigue:lc_symp_vision, names_to = "symptom", values_to = "value") %>% 
  group_by(infect_wave, symptom, value) %>% 
  summarise(n = n()) %>% 
  filter(!is.na(value)) %>% 
  mutate(ntot = sum(n), 
         freq = n / ntot) %>% 
  ungroup() %>% 
  filter(!is.na(infect_wave), value == "Yes") %>% 
  mutate(symptom = fct_reorder(symptom, freq)) %>% 
  ggplot(aes(y = symptom, x = freq, fill = fct_rev(infect_wave))) + 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) + 
  scale_fill_manual(values = col[c(11,10,9,12)]) + 
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme_bw() + labs(fill = "Infection wave")

broom::tidy(glm(fct_rev(lc_symp_fatigue) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_sleeping) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_reducedperformance) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_concentration) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_memory) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_arthralgia) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_dyspnea) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_skinproblems) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_myalgia) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_cough) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_rash) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_tingling) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_headache) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_vertigo) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_gisymptoms) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_palpitations) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_hearing) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_appetiteloss) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_tremor) ~ infect_wave + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_nervousness) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_depression) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_hairloss) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_vision) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_fever) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)

## Complications compared to non-infected
table1(~ diag_lc + 
         diag_lung + diag_heart + diag_stroke + diag_kidney + diag_liver + diag_gi + 
         diag_coagulation + diag_immunesystem + diag_bloodother + diag_skin + diag_thyroid + 
         diag_brain + diag_nervefunction + diag_other
       | ever_infected_2l,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = mutate(data, ever_infected_2l = ifelse(!is.na(ever_infected) & ever_infected == "Yes", "Yes", "No")))

data %>% 
  pivot_longer(cols = diag_lung:diag_other, names_to = "symptom", values_to = "value") %>% 
  group_by(ever_infected_2l, symptom, value) %>% 
  summarise(n = n()) %>% 
  filter(!is.na(value)) %>% 
  mutate(ntot = sum(n), 
         freq = n / ntot) %>% 
  ungroup() %>% 
  filter(!is.na(ever_infected_2l), value == "Yes") %>% 
  mutate(symptom = fct_reorder(symptom, freq)) %>% 
  ggplot(aes(y = symptom, x = freq, fill = fct_rev(ever_infected_2l))) + 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) + 
  scale_fill_manual(values = col[c(9,11)]) + 
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme_bw() + labs(fill = "Ever infected")

broom::tidy(glm(fct_rev(diag_lung) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_heart) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_stroke) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_kidney) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_liver) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_gi) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_coagulation) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_immunesystem) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_bloodother) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_skin) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_thyroid) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_brain) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_nervefunction) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_other) ~ ever_infected_2l + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)

## Complications compared to non-infected and by infection wave
table1(~ diag_lc + 
         diag_lung + diag_heart + diag_stroke + diag_kidney + diag_liver + diag_gi + 
         diag_coagulation + diag_immunesystem + diag_bloodother + diag_skin + diag_thyroid + 
         diag_brain + diag_nervefunction + diag_other
       | infect_wave,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = mutate(data, ever_infected_2l = ifelse(!is.na(ever_infected) & ever_infected == "Yes", "Yes", "No")))

data %>% 
  pivot_longer(cols = diag_lung:diag_other, names_to = "symptom", values_to = "value") %>% 
  group_by(infect_wave, symptom, value) %>% 
  summarise(n = n()) %>% 
  filter(!is.na(value)) %>% 
  mutate(ntot = sum(n), 
         freq = n / ntot) %>% 
  ungroup() %>% 
  filter(!is.na(infect_wave), value == "Yes") %>% 
  mutate(symptom = fct_reorder(symptom, freq)) %>% 
  ggplot(aes(y = symptom, x = freq, fill = fct_rev(infect_wave))) + 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) + 
  scale_fill_manual(values = col[c(11,10,9,12)]) + 
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme_bw() + labs(fill = "Infection wave")

broom::tidy(glm(fct_rev(diag_lung) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_heart) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_stroke) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_kidney) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_liver) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_gi) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_coagulation) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_immunesystem) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_bloodother) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_skin) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_thyroid) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_brain) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_nervefunction) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_other) ~ infect_wave + age + sex, data = data, family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)

## Calculate EQ-5D-5L score
library(eq5d)
data$eq5d_score <- eq5d(data.frame(MO = data$eq5d_mb, SC = data$eq5d_sc, UA = data$eq5d_ua, PD = data$eq5d_pd, AD = data$eq5d_ad), version = "5L", country = "Netherlands", type = "VT", ignore.invalid = TRUE)

## EQ-5D-5L scores compared to non-infected
table1(~ eq5d_score + eq_vas
       | ever_infected_2l,
       render.continuous = render.cont.eq5d, render.categorical = render.cat, 
       data = data)

broom::tidy(lm(eq5d_score ~ ever_infected_2l + age + sex, data = data), conf.int = TRUE, p.value = TRUE)
broom::tidy(lm(eq_vas ~ ever_infected_2l + age + sex, data = data), conf.int = TRUE, p.value = TRUE)

## EQ-5D-5L scores compared to non-infected and by infection wave
table1(~ eq5d_score + eq_vas
       | infect_wave,
       render.continuous = render.cont.eq5d, render.categorical = render.cat, 
       data = data)

data %>% 
  filter(!is.na(infect_wave)) %>% 
  ggplot(aes(y = eq5d_score, x = infect_wave)) + 
  geom_boxplot() + 
  theme_bw()

data %>% 
  filter(!is.na(infect_wave)) %>% 
  ggplot(aes(y = eq_vas, x = infect_wave)) + 
  geom_boxplot() + 
  theme_bw()

broom::tidy(lm(eq5d_score ~ infect_wave + age + sex, data = data), conf.int = TRUE, p.value = TRUE)
broom::tidy(lm(eq_vas ~ infect_wave + age + sex, data = data), conf.int = TRUE, p.value = TRUE)

### COMPARISON WITH KNOWN INFECTION OR SEROPOSITIVE ---------------------------------------------------------------------- 

## Symptom prevalence compared to seronegatives and never infected
table1(~ lc_symp_fatigue + lc_symp_sleeping + lc_symp_headache + lc_symp_appetiteloss + lc_symp_concentration + 
         lc_symp_memory + lc_symp_nervousness + lc_symp_depression + lc_symp_rash + lc_symp_skinproblems + 
         lc_symp_gisymptoms + lc_symp_hairloss + lc_symp_reducedperformance + lc_symp_palpitations + lc_symp_dyspnea + 
         lc_symp_cough + lc_symp_myalgia + lc_symp_arthralgia + lc_symp_fever + lc_symp_tingling + lc_symp_tremor + 
         lc_symp_hearing + lc_symp_vertigo + lc_symp_vision
       | infect_or_seropos,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = mutate(data, ever_infected_2l = ifelse(!is.na(ever_infected) & ever_infected == "Yes", "Yes", "No")))

data %>% 
  pivot_longer(cols = lc_symp_fatigue:lc_symp_vision, names_to = "symptom", values_to = "value") %>% 
  group_by(infect_or_seropos, symptom, value) %>% 
  summarise(n = n()) %>% 
  filter(!is.na(value)) %>% 
  mutate(ntot = sum(n), 
         freq = n / ntot) %>% 
  ungroup() %>% 
  filter(!is.na(infect_or_seropos), value == "Yes") %>% 
  mutate(symptom = fct_reorder(symptom, freq)) %>% 
  ggplot(aes(y = symptom, x = freq, fill = infect_or_seropos)) + 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) + 
  scale_fill_manual(values = col[c(9,11)]) + 
  theme_bw()

## Complications compared to seronegatives and never infected
table1(~ diag_lc + 
         diag_lung + diag_heart + diag_stroke + diag_kidney + diag_liver + diag_gi + 
         diag_coagulation + diag_immunesystem + diag_bloodother + diag_skin + diag_thyroid + 
         diag_brain + diag_nervefunction + diag_other
       | infect_or_seropos,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = mutate(data, ever_infected_2l = ifelse(!is.na(ever_infected) & ever_infected == "Yes", "Yes", "No")))

data %>% 
  pivot_longer(cols = diag_lung:diag_other, names_to = "symptom", values_to = "value") %>% 
  group_by(infect_or_seropos, symptom, value) %>% 
  summarise(n = n()) %>% 
  filter(!is.na(value)) %>% 
  mutate(ntot = sum(n), 
         freq = n / ntot) %>% 
  ungroup() %>% 
  filter(!is.na(infect_or_seropos), value == "Yes") %>% 
  mutate(symptom = fct_reorder(symptom, freq)) %>% 
  ggplot(aes(y = symptom, x = freq, fill = infect_or_seropos)) + 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) + 
  scale_fill_manual(values = col[c(9,11)]) + 
  theme_bw()

## EQ-5D-5L scores compared to seronegatives
table1(~ eq5d_score + eq_vas
       | infect_or_seropos,
       render.continuous = render.cont.eq5d, render.categorical = render.cat, 
       data = data)

### COMPARISON BETWEEN ASYMPTOMATIC INFECTIONS AND SYMPTOMATIC INFECTIONS ---------------------------------------------------------------------- 

## Overview symptom status among those ever infected or seropositive
table1(~ symp_at_infection + symp_at_infection_alt + ever_reinfected
       | infect_or_seropos, 
       render.continuous = render.cont, render.categorical = render.cat, 
       data = data)

## Symptom prevalence after symptomatic infection compared to asymptomatic infection
table1(~ lc_symp_fatigue + lc_symp_sleeping + lc_symp_headache + lc_symp_appetiteloss + lc_symp_concentration + 
         lc_symp_memory + lc_symp_nervousness + lc_symp_depression + lc_symp_rash + lc_symp_skinproblems + 
         lc_symp_gisymptoms + lc_symp_hairloss + lc_symp_reducedperformance + lc_symp_palpitations + lc_symp_dyspnea + 
         lc_symp_cough + lc_symp_myalgia + lc_symp_arthralgia + lc_symp_fever + lc_symp_tingling + lc_symp_tremor + 
         lc_symp_hearing + lc_symp_vertigo + lc_symp_vision
       | symp_at_infection_alt,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = filter(data, infect_or_seropos == "Yes"))

data %>% 
  filter(infect_or_seropos == "Yes") %>% 
  pivot_longer(cols = lc_symp_fatigue:lc_symp_vision, names_to = "symptom", values_to = "value") %>% 
  group_by(symp_at_infection_alt, symptom, value) %>%  # infect_wave
  summarise(n = n()) %>% 
  filter(!is.na(value)) %>% 
  mutate(ntot = sum(n), 
         freq = n / ntot) %>% 
  ungroup() %>% 
  filter(!is.na(symp_at_infection_alt), value == "Yes") %>% 
  mutate(symptom = fct_reorder(symptom, freq)) %>% 
  ggplot(aes(y = symptom, x = freq, fill = symp_at_infection_alt)) + 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) + 
  # facet_wrap(. ~ infect_wave) + 
  scale_fill_manual(values = col[c(9,11)]) + 
  theme_bw()

broom::tidy(glm(fct_rev(lc_symp_fatigue) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_sleeping) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_reducedperformance) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_concentration) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_memory) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_arthralgia) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_dyspnea) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_skinproblems) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_myalgia) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_cough) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_rash) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_tingling) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_headache) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_vertigo) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_gisymptoms) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_palpitations) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_hearing) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_appetiteloss) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_tremor) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_nervousness) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_depression) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_hairloss) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_vision) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(lc_symp_fever) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)

## Complications after symptomatic infection compared to asymptomatic infection
table1(~ diag_lc + 
         diag_lung + diag_heart + diag_stroke + diag_kidney + diag_liver + diag_gi + 
         diag_coagulation + diag_immunesystem + diag_bloodother + diag_skin + diag_thyroid + 
         diag_brain + diag_nervefunction + diag_other
       | symp_at_infection_alt,
       render.continuous = render.cont, render.categorical = render.cat, 
       data = filter(data, infect_or_seropos == "Yes"))

data %>% 
  filter(infect_or_seropos == "Yes") %>% 
  pivot_longer(cols = diag_lung:diag_other, names_to = "symptom", values_to = "value") %>% 
  group_by(symp_at_infection_alt, symptom, value) %>% 
  summarise(n = n()) %>% 
  filter(!is.na(value)) %>% 
  mutate(ntot = sum(n), 
         freq = n / ntot) %>% 
  ungroup() %>% 
  filter(!is.na(symp_at_infection_alt), value == "Yes") %>% 
  mutate(symptom = fct_reorder(symptom, freq)) %>% 
  ggplot(aes(y = symptom, x = freq, fill = symp_at_infection_alt)) + 
  geom_bar(stat = "identity", position = position_dodge(preserve = "single")) + 
  scale_fill_manual(values = col[c(9,11)]) + 
  theme_bw()

broom::tidy(glm(fct_rev(diag_lung) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_heart) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_stroke) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_kidney) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_liver) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_gi) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_coagulation) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_immunesystem) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_bloodother) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_skin) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_thyroid) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_brain) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_nervefunction) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)
broom::tidy(glm(fct_rev(diag_other) ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes"), family = "binomial"), conf.int = TRUE, p.value = TRUE, exp = TRUE)

## EQ-5D-5L scores after symptomatic infection compared to asymptomatic infection
table1(~ eq5d_score + eq_vas
       | symp_at_infection_alt,
       render.continuous = render.cont.eq5d, render.categorical = render.cat, 
       data = filter(data, infect_or_seropos == "Yes"))

broom::tidy(lm(eq5d_score ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes")), conf.int = TRUE, p.value = TRUE)
broom::tidy(lm(eq5d_score ~ symp_at_infection_alt + age + sex + ever_reinfected_2l, data = filter(data, infect_or_seropos == "Yes")), conf.int = TRUE, p.value = TRUE)
broom::tidy(lm(eq_vas ~ symp_at_infection_alt + age + sex, data = filter(data, infect_or_seropos == "Yes")), conf.int = TRUE, p.value = TRUE)
broom::tidy(lm(eq_vas ~ symp_at_infection_alt + age + sex + ever_reinfected_2l, data = filter(data, infect_or_seropos == "Yes")), conf.int = TRUE, p.value = TRUE)


