library(tidyverse)
library(lubridate)
library(cowplot)
theme_set(theme_bw())
library(rstanarm)
options(mc.cores = parallel::detectCores())

# functions ---------------------
logit = function(x) log(x/(1-x))
inv.logit = function(x) exp(x)/(1+exp(x))
qsum = function(x) c(quantile(x,probs=c(0.5,0.025,0.975)))

# load --------------------------
usz = read_tsv("shared_Julien/usz_annot_qda.txt")
bds = read_tsv("shared_Julien/bds_annot_qda_OnlyCopandemic.txt") %>%
  select(-patient_id,-clinic,-known_pos,-known_neg,-days_since_PCR_pos,-pcr_pos)
bds_pre = read_tsv("shared_Julien/bds_annot_qda_untilMarch.txt") %>%
  filter(!is.na(Pseudo_ID)) %>%
  transmute(sample_id=Pseudo_ID,month=mymonth,Spike_plogEC50=as.numeric(Spike_IgG),NC_plogEC50=as.numeric(NC_IgG), RBD_plogEC50=as.numeric(RBD_IgG),posterior=posteriorProb) 
bds = bind_rows(bds_pre,bds)
boot_usz = read_tsv("shared_Julien/bootstrap_usz.tbl",n_max=nrow(usz)*100)
boot_bds = read_tsv("shared_Julien/bootstrap_bds.tbl",n_max=nrow(bds)*100)
demog_zh = read_csv2("shared_Julien/KANTON_ZUERICH_bevoelkerung_1jahresklassen.csv")


# demog data ------------------------
demog_zh = demog_zh %>%
  filter(JAHR==2019) %>%
  dplyr::select(sex=GESCHLECHT_CODE,age=ALTERSKLASSE_CODE,n_pop=ANZAHL_PERSONEN) %>%
  mutate(age_group=cut(age,c(0,18,seq(30,80,by=10),110),right=FALSE),
         sex=ifelse(sex==2,"female","male")) %>%
  group_by(sex,age_group) %>%
  summarise(n_pop=sum(n_pop)) %>%
  ungroup() %>%
  mutate(p_pop=n_pop/sum(n_pop))

# data management serology data ----------------------
## dates
usz_dup = usz %>%
  mutate(year=as.numeric(paste0("20",substr(month,1,2))),
         month2=substr(month,3,4),
         date=ymd(paste0(year,month2,day)),
         month3=factor(month,
                       levels=c("1800","1912","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012"),
                       labels=c("2018","2019","Jan. 2020","Feb. 2020","March 2020","April 2020","May 2020","June 2020","July 2020","Aug. 2020","Sept. 2020","Oct. 2020","Nov 2020","Dec. 2020"))) 
usz = usz_dup %>%
  distinct(patient_id,month3,.keep_all=TRUE)
bds = bds %>%
  mutate(year=as.numeric(paste0("20",substr(month,1,2))),
         month2=substr(month,3,4),
         date=ymd(paste0(year,month2,day)),
         month3=factor(month,
                       levels=c("1800","1912","2001","2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012"),
                       labels=c("2018","2019","Jan. 2020","Feb. 2020","March 2020","April 2020","May 2020","June 2020","July 2020","Aug. 2020","Sept. 2020","Oct. 2020","Nov 2020","Dec. 2020"))) %>%
  distinct(sample_id,.keep_all=TRUE)

## age groups
usz = usz %>%
  mutate(age_group=cut(age,c(0,18,seq(30,80,by=10),110),right=FALSE))
usz_dup = usz_dup %>%
  mutate(age_group=cut(age,c(0,18,seq(30,80,by=10),110),right=FALSE))
bds = bds %>%
  mutate(age_group=cut(age,c(0,18,seq(30,80,by=10)),right=FALSE))


## bootstrap
boot_usz = inner_join(boot_usz,usz,by=c("Pseudo_ID"="sample_id")) 
boot_bds = inner_join(boot_bds,bds,by=c("Pseudo_ID"="sample_id")) 

