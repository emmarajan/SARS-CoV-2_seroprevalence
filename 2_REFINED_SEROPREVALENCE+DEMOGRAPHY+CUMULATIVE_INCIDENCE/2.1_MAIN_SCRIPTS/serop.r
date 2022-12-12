
# compute seroprevalence in different ways ----------------

## no adjustment
agg_serop = NULL
boot_tmp = boot_usz %>%
  group_by(month3,bootstrap_nr) %>%
  summarise(p=sum(posterior)/n()) %>%
  summarise(p_low=quantile(p,0.025),
            p_high=quantile(p,0.975))
agg_serop = usz %>%
  group_by(month3) %>%
  summarise(p=sum(posterior)/n(),
            sample_size=n()) %>%
  left_join(boot_tmp) %>%
  mutate(sample="USZ",
         type="raw") %>%
  bind_rows(agg_serop)

boot_tmp = boot_bds %>%
  group_by(month3,bootstrap_nr) %>%
  summarise(p=sum(posterior)/n()) %>%
  summarise(p_low=quantile(p,0.025),
            p_high=quantile(p,0.975))
agg_serop = bds %>%
  group_by(month3) %>%
  summarise(p=sum(posterior)/n(),
            sample_size=n()) %>%
  left_join(boot_tmp) %>%
  mutate(sample="BDS",
         type="raw",
         month3=factor(month3,levels=levels(agg_serop$month3))) %>%
  bind_rows(agg_serop)

## post-strat on age and sex
demog_tmp = demog_zh %>%
  group_by(age_group,sex) %>%
  summarise(n_pop=sum(n_pop)) %>%
  ungroup() %>%
  mutate(p_pop=n_pop/sum(n_pop))
boot_tmp = boot_usz %>%
  group_by(month3,age_group,sex,bootstrap_nr) %>%
  summarise(p=sum(posterior)/n()) %>%
  left_join(demog_tmp) %>%
  group_by(month3,bootstrap_nr) %>%
  summarise(p=sum(p*p_pop)) %>%
  summarise(p_low=quantile(p,0.025),
            p_high=quantile(p,0.975))
tmp = usz %>%
  filter(!is.na(age_group),!is.na(sex)) %>%
  group_by(month3) %>%
  summarise(sample_size=n())
agg_serop = usz %>%
  filter(!is.na(age_group),!is.na(sex)) %>%
  group_by(month3,age_group,sex) %>%
  summarise(p=sum(posterior)/n()) %>%
  left_join(demog_tmp) %>%
  group_by(month3) %>%
  summarise(p=sum(p*p_pop)) %>%
  left_join(boot_tmp) %>%
  mutate(sample="USZ",
         type="adj_age_sex") %>%
  left_join(tmp) %>%
  bind_rows(agg_serop) 


demog_tmp = demog_zh %>%
  filter(!(age_group %in% c("[0,18)","[80,110)"))) %>%
  mutate(age_group=factor(age_group)) %>%
  group_by(age_group,sex) %>%
  summarise(n_pop=sum(n_pop)) %>%
  ungroup() %>%
  mutate(p_pop=n_pop/sum(n_pop))
boot_tmp = boot_bds %>%
  filter(!is.na(age_group),!is.na(sex)) %>%
  filter(!(age_group %in% c("[0,18)","[80,110)"))) %>%
  group_by(month3,age_group,sex,bootstrap_nr) %>%
  summarise(p=sum(posterior)/n()) %>%
  left_join(demog_tmp) %>%
  group_by(month3,bootstrap_nr) %>%
  summarise(p=sum(p*p_pop)) %>%
  summarise(p_low=quantile(p,0.025),
            p_high=quantile(p,0.975))
tmp = bds %>%
  filter(!is.na(age_group),!is.na(sex)) %>%
  group_by(month3) %>%
  summarise(sample_size=n())
agg_serop = bds %>%
  group_by(month3,age_group,sex) %>%
  summarise(p=sum(posterior)/n()) %>%
  left_join(demog_tmp) %>%
  group_by(month3) %>%
  summarise(p=sum(p*p_pop)) %>%
  left_join(boot_tmp) %>%
  mutate(sample="BDS",
         type="adj_age_sex") %>%
  left_join(tmp) %>%
  bind_rows(agg_serop)


agg_serop = agg_serop %>%
  mutate(month3=factor(month3,levels=levels(usz$month3)),
         type=factor(type,levels=c("raw","adj_age_sex"),labels=c("No","On age and sex")))


agg_serop = arrange(agg_serop,sample,type,month3)



# icd-10 -----------------------

## import
icd10_raw = readxl::read_xlsx("shared_Julien/ICD_10_COVID_20201227.xlsx")

icd10_key = read_csv2("data/CIM10GM2021_CSV_S_FR_versionmtadonne_codes_20201222.csv",col_name=FALSE) %>%
  mutate(X11=if_else(is.na(X11),"",X11)) %>%
  mutate(icd10_text = paste(X10,X11)) %>%
  dplyr::select(icd10_code=X8,icd10_text) %>%
  mutate(icd10_code=if_else(icd10_code=="U690","U6900",icd10_code)) 

## wide format
icd10 = icd10_raw %>%
  pivot_longer(4:103,names_to="var",values_to="ICD10") %>%
  dplyr::select(patient_id=pat_research_id,ICD10) %>%
  filter(!(ICD10 %in% c("NULL","?",".","0","538","797","85"))) %>%
  mutate(present=TRUE) %>%
  distinct() %>%
  group_by(ICD10) %>%
  mutate(n=n()) %>%
  ungroup() %>%
  mutate(n=n/n()) %>%
  filter(n>0.001) %>%
  dplyr::select(-n) %>%
  arrange(ICD10) %>%
  pivot_wider(names_from=ICD10,values_from=present) %>%
  mutate_all(.funs = function(x) ifelse(is.na(x),FALSE,x))

## merge with posterior serology 
icd10 = usz %>%
  dplyr::select(patient_id,sex,age,posterior) %>%
  group_by(patient_id) %>%
  summarise(posterior=max(posterior),sex=max(sex),age=max(age)) %>%
  inner_join(icd10) 

# compute seroprevalence after removing likely hospitalisations caused by covid-19
icd_code_covid = icd10 %>%
  filter(J96.00 | U99.0) %>%
  pull(patient_id)
hospit_covid = usz %>%
  filter(clinic %in% c("Internal Medicine","Infectious Diseases and Hospital Hygiene")) %>%
  pull(patient_id)
id_covid = unique(c(icd_code_covid,hospit_covid))  


## no adjustment
agg_serop2 = NULL
boot_tmp = boot_usz %>%
  filter(!(patient_id %in% id_covid)) %>%
  group_by(month3,bootstrap_nr) %>%
  summarise(p=sum(posterior)/n()) %>%
  summarise(p_low=quantile(p,0.025),
            p_high=quantile(p,0.975))
agg_serop2 = usz %>%
  filter(!(patient_id %in% id_covid)) %>%
  group_by(month3) %>%
  summarise(p=sum(posterior)/n()) %>%
  left_join(boot_tmp) %>%
  mutate(sample="USZ",
         type="raw") %>%
  bind_rows(agg_serop2)

## post-strat on age and sex
demog_tmp = demog_zh %>%
  group_by(age_group,sex) %>%
  summarise(n_pop=sum(n_pop)) %>%
  ungroup() %>%
  mutate(p_pop=n_pop/sum(n_pop))
boot_tmp = boot_usz %>%
  filter(!(patient_id %in% id_covid)) %>%
  group_by(month3,age_group,sex,bootstrap_nr) %>%
  summarise(p=sum(posterior)/n()) %>%
  left_join(demog_tmp) %>%
  group_by(month3,bootstrap_nr) %>%
  summarise(p=sum(p*p_pop)) %>%
  summarise(p_low=quantile(p,0.025),
            p_high=quantile(p,0.975))
agg_serop2 = usz %>%
  filter(!(patient_id %in% id_covid)) %>%
  group_by(month3,age_group,sex) %>%
  summarise(p=sum(posterior)/n()) %>%
  left_join(demog_tmp) %>%
  group_by(month3) %>%
  summarise(p=sum(p*p_pop)) %>%
  left_join(boot_tmp) %>%
  mutate(sample="USZ",
         type="adj_age_sex") %>%
  bind_rows(agg_serop2)

agg_serop2 = agg_serop2 %>%
  mutate(month3=factor(month3,levels=levels(usz$month3)),
         type=factor(type,levels=c("raw","adj_age_sex"),labels=c("No","On age and sex")))



agg_serop2 = arrange(agg_serop2,sample,type,month3)
