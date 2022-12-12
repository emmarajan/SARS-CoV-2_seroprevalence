#' ---
#' title: "Seroprevalence Zürich - analysis report"
#' author: "Julien Riou"
#' date: "January 5, 2021"
#' output:
#'    html_document:
#'      code_folding : hide
#'      toc: true
#'      toc_float: true
#' ---

#+ results="hide", warnings="false", echo="false"
knitr::opts_chunk$set(warning = FALSE,message=FALSE)
suppressWarnings(suppressMessages(source("dm.r")))
suppressWarnings(suppressMessages(source("serop.r")))
# rmarkdown::render("report.r")

#+ fig.width=9, fig.height=3
g1 = usz %>%
  group_by(age_group) %>%
  summarise(n=n()) %>%
  mutate(p=n/sum(n)) %>%
  ggplot() +
  geom_col(aes(x=as.numeric(age_group),y=p,fill=age_group),colour="black") +
  scale_fill_brewer(type="seq",palette="Blues") +
  labs(x=NULL,y="Proportion",fill="Age group",title="Age distribution in USZ sample") +
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,.3)) +
  scale_x_continuous(expand=c(0,0),breaks=1:length(unique(usz$age_group)),labels=levels(usz$age_group),limits=c(0.5,8.5)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))
g2 = bds %>%
  group_by(age_group) %>%
  summarise(n=n()) %>%
  mutate(p=n/sum(n)) %>%
  ggplot() +
  geom_col(aes(x=as.numeric(age_group),y=p,fill=age_group),colour="black") +
  scale_fill_brewer(type="seq",palette="Reds") +
  labs(x=NULL,y="Proportion",fill="Age group",title="Age distribution in BDS sample") +
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,.3)) +
  scale_x_continuous(expand=c(0,0),breaks=1:length(levels(bds$age_group)),labels=levels(bds$age_group),limits=c(0.5,8.5)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))
plot_grid(g1,g2)

#+ fig.width=9, fig.height=3

demog_tmp = demog_zh %>%
  group_by(age_group) %>%
  summarise(n_pop=sum(n_pop)) %>%
  mutate(p_pop=n_pop/sum(n_pop))
g1 = usz %>%
  group_by(age_group) %>%
  summarise(n=n()) %>%
  mutate(p=n/sum(n)) %>%
  left_join(demog_tmp) %>%
  ggplot() +
  geom_col(aes(x=as.numeric(age_group),y=p,fill=age_group),colour="black") +
  geom_line(aes(x=as.numeric(age_group),y=p_pop),colour="black") +
  geom_point(aes(x=as.numeric(age_group),y=p_pop),colour="black",size=2) +
  scale_fill_brewer(type="seq",palette="Blues") +
  labs(x=NULL,y="Proportion",fill="Age group",title="Age distribution in USZ sample") +
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,.3)) +
  scale_x_continuous(expand=c(0,0),breaks=1:length(levels(usz$age_group)),labels=levels(usz$age_group)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))
demog_tmp = demog_zh %>%
  filter(age_group!="[80,110)") %>%
  mutate(age_group=factor(age_group,levels=levels(bds$age_group))) %>%
  group_by(age_group) %>%
  summarise(n_pop=sum(n_pop)) %>%
  mutate(p_pop=n_pop/sum(n_pop)) 
g2 = bds %>%
  group_by(age_group) %>%
  summarise(n=n()) %>%
  mutate(p=n/sum(n)) %>%
  left_join(demog_tmp) %>%
  ggplot() +
  geom_col(aes(x=as.numeric(age_group),y=p,fill=age_group),colour="black") +
  geom_line(aes(x=as.numeric(age_group),y=p_pop),colour="black") +
  geom_point(aes(x=as.numeric(age_group),y=p_pop),colour="black",size=2) +
  scale_fill_brewer(type="seq",palette="Reds") +
  labs(x=NULL,y="Proportion",fill="Age group",title="Age distribution in BDS sample") +
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,.3)) +
  scale_x_continuous(expand=c(0,0),breaks=1:length(levels(bds$age_group)),labels=levels(bds$age_group),limits=c(0.5,8.5)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))
plot_grid(g1,g2)


#' Evolution of age distribution with time in USZ and BDS samples.

#+ fig.width=9, fig.height=3
g1 = usz %>%
  group_by(month3,age_group) %>%
  summarise(n=n()) %>%
  mutate(p=n/sum(n)) %>%
  ggplot() +
  geom_area(aes(x=as.numeric(month3),y=p,color=age_group,fill=age_group),colour="black") +
  geom_vline(xintercept=2.5,linetype=2) +
  scale_fill_brewer(type="seq",palette="Blues") +
  labs(x=NULL,y="Proportion",fill="Age group",title="Age distribution in USZ sample") +
  scale_y_continuous(expand=c(0,0),labels=scales::percent) +
  scale_x_continuous(expand=c(0,0),breaks=1:length(unique(usz$month3)),labels=levels(usz$month3)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))
g2 = bds %>%
  group_by(month3,age_group) %>%
  summarise(n=n()) %>%
  mutate(p=n/sum(n)) %>%
  ggplot() +
  geom_area(aes(x=as.numeric(month3),y=p,color=age_group,fill=age_group),colour="black") +
  scale_fill_brewer(type="seq",palette="Reds") +
  labs(x=NULL,y="Proportion",fill="Age group",title="Age distribution in BDS sample") +
  scale_y_continuous(expand=c(0,0),labels=scales::percent) +
  scale_x_continuous(expand=c(0,0),breaks=as.numeric(unique(bds$month3)),labels=unique(bds$month3)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))
plot_grid(g1,g2)

#' ## Sex

#+ fig.width=9, fig.height=3
g1 = usz %>%
  group_by(month3,sex) %>%
  summarise(n=n()) %>%
  mutate(p=n/sum(n)) %>%
  ggplot() +
  geom_area(aes(x=as.numeric(month3),y=p,color=sex,fill=sex),colour="black") +
  geom_vline(xintercept=2.5,linetype=2) +
  geom_hline(yintercept=.5,linetype=3) + 
  scale_fill_brewer(type="seq",palette="Blues") +
  labs(x=NULL,y="Proportion",fill="Sex",title="Sex distribution in USZ sample") +
  scale_y_continuous(expand=c(0,0),labels=scales::percent) +
  scale_x_continuous(expand=c(0,0),breaks=1:length(unique(usz$month3)),labels=levels(usz$month3)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))

g2 = bds %>%
  group_by(month3,sex) %>%
  summarise(n=n()) %>%
  mutate(p=n/sum(n)) %>%
  ggplot() +
  geom_area(aes(x=as.numeric(month3),y=p,color=sex,fill=sex),colour="black") +
  scale_fill_brewer(type="seq",palette="Reds") +
  geom_hline(yintercept=.5,linetype=3) + 
  labs(x=NULL,y="Proportion",fill="Sex",title="Sex distribution in BDS sample") +
  scale_y_continuous(expand=c(0,0),labels=scales::percent) +
  scale_x_continuous(expand=c(0,0),breaks=as.numeric(unique(bds$month3)),labels=unique(bds$month3)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))
plot_grid(g1,g2)


#' ## Hospital ward

#+ fig.width=9, fig.height=4.5
tmp = usz %>%
  filter(clinic!="Oncology") %>%
  group_by(month3,clinic) %>%
  summarise(n=n()) %>%
  mutate(p=n/sum(n))  %>%
  group_by(clinic) %>%
  mutate(max_p=max(p)) %>%
  filter(max_p>0.05) %>%
  mutate(infectio=ifelse(clinic %in% c("Infectious Diseases and Hospital Hygiene","Cardiology"),TRUE,FALSE))
tmp_names = tmp %>%
  group_by(clinic) %>%
  summarise(n=sum(n)) %>%
  arrange(n)
tmp$clinic2 = factor(tmp$clinic,levels=tmp_names$clinic)

ggplot(tmp) +
  geom_area(aes(x=as.numeric(month3),y=p,fill=clinic2,colour=infectio)) +
  geom_vline(xintercept=2.5,linetype=2) +
  scale_fill_discrete() +
  scale_colour_manual(values=c("black","grey40"),guide=FALSE) +
  # scale_size_manual(values=c(0.5,1),guide=FALSE) +
  labs(x=NULL,y="Proportion",fill="Hospital ward",title="Hospital ward in USZ sample") +
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,1)) +
  scale_x_continuous(expand=c(0,0),breaks=1:length(unique(usz$month3)),labels=levels(usz$month3)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))

#+ fig.width=9, fig.height=3

g1 = 
  filter(agg_serop,sample=="USZ",type=="No") %>%
  ggplot() +
  geom_pointrange(aes(x=as.numeric(month3),y=p,ymin=p_low,ymax=p_high),color="#0000FF") +
  geom_vline(xintercept=2.5,linetype=2) +
  labs(x=NULL,y="Proportion",title="Seroprevalence in USZ sample") +
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,0.075)) +
  scale_x_continuous(expand=c(0,0),breaks=1:length(unique(usz$month3)),labels=levels(usz$month3),limits = c(0.5,14.5)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))
g2 = filter(agg_serop,sample=="BDS",type=="No") %>%
  ggplot() +
  geom_pointrange(aes(x=as.numeric(month3),y=p,ymin=p_low,ymax=p_high),color="firebrick3") +
  geom_vline(xintercept=2.5,linetype=2) +
  labs(x=NULL,y="Proportion",title="Seroprevalence in BDS sample") +
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,0.075)) +
  scale_x_continuous(expand=c(0,0),breaks=1:length(unique(usz$month3)),labels=levels(usz$month3),limits = c(0.5,14.5)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))
plot_grid(g1,g2)

#+ fig.width=9, fig.height=3
tmp1 = usz %>%
  group_by(age_group,sex) %>%
  summarise(p=sum(posterior)/n())
tmp1 = boot_usz %>%
  group_by(age_group,sex,bootstrap_nr) %>%
  summarise(serop1_usz=sum(posterior)/n()) %>%
  summarise(p_low=quantile(serop1_usz,0.025),
            p_high=quantile(serop1_usz,0.975)) %>%
  right_join(tmp1)
g1 = ggplot(tmp1) +
  geom_pointrange(aes(x=age_group,y=p,ymax=p_high,ymin=p_low,colour=sex),position=position_dodge(0.4)) +
  scale_colour_manual(values=c("#0000FF","#1E90FF")) +
  labs(x=NULL,y="Proportion",colour="Sex",title="Seroprevalence by age group and sex in USZ") +
  scale_y_continuous(expand=expand_scale(mult=c(0,0.05)),labels=scales::percent,limits=c(0,0.057)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))

tmp2 = bds %>%
  group_by(age_group,sex) %>%
  summarise(p=sum(posterior)/n())
tmp2 = boot_bds %>%
  group_by(age_group,sex,bootstrap_nr) %>%
  summarise(serop1_usz=sum(posterior)/n()) %>%
  summarise(p_low=quantile(serop1_usz,0.025),
            p_high=quantile(serop1_usz,0.975)) %>%
  right_join(tmp2)

g2 = ggplot(tmp2) +
  geom_pointrange(aes(x=age_group,y=p,ymax=p_high,ymin=p_low,colour=sex),position=position_dodge(0.4)) +
  scale_colour_manual(values=c("firebrick1","firebrick4")) +
  labs(x=NULL,y="Proportion",colour="Sex",title="Seroprevalence by age group and sex in BDS") +
  scale_y_continuous(expand=expand_scale(mult=c(0,0.05)),labels=scales::percent,limits=c(0,0.057)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))
plot_grid(g1,g2)

#+ fig.width=9, fig.height=3

g1 = 
  filter(agg_serop,sample=="USZ") %>%
  ggplot() +
  geom_pointrange(aes(x=as.numeric(month3),y=p,ymin=p_low,ymax=p_high,color=type),position=position_dodge(.4)) +
  geom_vline(xintercept=2.5,linetype=2) +
  labs(x=NULL,y="Proportion",title="Seroprevalence in USZ sample",colour="Post-stratification") +
  scale_colour_manual(values=c("#0000FF","#1E90FF")) +
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,0.075)) +
  scale_x_continuous(expand=c(0,0),breaks=1:length(unique(usz$month3)),labels=levels(usz$month3),limits = c(0.5,14.5)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))
g2 = filter(agg_serop,sample=="BDS") %>%
  ggplot() +
  geom_pointrange(aes(x=as.numeric(month3),y=p,ymin=p_low,ymax=p_high,color=type),position=position_dodge(.4)) +
  geom_vline(xintercept=2.5,linetype=2) +
  labs(x=NULL,y="Proportion",title="Seroprevalence in BDS sample",colour="Post-stratification") +
  scale_colour_manual(values=c("firebrick2","firebrick4")) +
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,0.075)) +
  scale_x_continuous(expand=c(0,0),breaks=1:length(unique(usz$month3)),labels=levels(usz$month3),limits = c(0.5,14.5)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))
plot_grid(g1,g2)

#+ fig.width=9, fig.height=3
tmp1 = usz %>%
  mutate(icd_code_covid=patient_id %in% icd_code_covid) %>%
  group_by(icd_code_covid) %>%
  summarise(p=sum(posterior)/n())
tmp1 = boot_usz %>%
  mutate(icd_code_covid=patient_id %in% icd_code_covid) %>%
  group_by(icd_code_covid,bootstrap_nr) %>%
  summarise(serop1_usz=sum(posterior)/n()) %>%
  summarise(p_low=quantile(serop1_usz,0.025),
            p_high=quantile(serop1_usz,0.975)) %>%
  right_join(tmp1)
tmp2 = usz %>%
  mutate(hospit_covid=patient_id %in% hospit_covid) %>%
  group_by(hospit_covid) %>%
  summarise(p=sum(posterior)/n())
tmp2 = boot_usz %>%
  mutate(hospit_covid=patient_id %in% hospit_covid) %>%
  group_by(hospit_covid,bootstrap_nr) %>%
  summarise(serop1_usz=sum(posterior)/n()) %>%
  summarise(p_low=quantile(serop1_usz,0.025),
            p_high=quantile(serop1_usz,0.975)) %>%
  right_join(tmp2)

g1 = ggplot(tmp1) +
  geom_pointrange(aes(x=icd_code_covid,y=p,ymax=p_high,ymin=p_low,colour=icd_code_covid),position=position_dodge(0.4)) +
  scale_colour_manual(values=c("#0000FF","#1E90FF")) +
  labs(x=NULL,y="Proportion",colour="ICD code\nJ96.00 or U99.0",title="Seroprevalence in USZ sample") +
  scale_y_continuous(expand=expand_scale(mult=c(0,0.05)),labels=scales::percent,limits=c(0,0.04)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))

g2 = ggplot(tmp2) +
  geom_pointrange(aes(x=hospit_covid,y=p,ymax=p_high,ymin=p_low,colour=hospit_covid),position=position_dodge(0.4)) +
  scale_colour_manual(values=c("#0000FF","#1E90FF")) +
  labs(x=NULL,y="Proportion",colour="Hospitalised in\nspecialist ward",title=" ") +
  scale_y_continuous(expand=expand_scale(mult=c(0,0.05)),labels=scales::percent,limits=c(0,0.04)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))


plot_grid(g1,g2,ncol=2)



#+ fig.width=9, fig.height=3
tmp = agg_serop2 %>%
  mutate(excluded="Yes") 
tmp = agg_serop %>%
  filter(sample=="USZ") %>%
  mutate(excluded="No (baseline)") %>%
  bind_rows(tmp)

ggplot(tmp) +
  geom_pointrange(aes(x=as.numeric(month3),y=p,ymin=p_low,ymax=p_high,color=type,shape=excluded,alpha=excluded),position=position_dodge(.4)) +
  labs(x=NULL,y="Proportion",title="Seroprevalence in USZ sample after removing patients hospitalized\nbecause of SARS-CoV-2 infection",colour="Post-stratification",shape="Exclusion of A and B",alpha="Exclusion of A and B") +
  scale_colour_manual(values=c("#0000FF","#1E90FF")) +
  scale_y_continuous(expand=c(0,0),labels=scales::percent,limits=c(0,0.075)) +
  scale_alpha_manual(values=c(.5,1)) +
  scale_x_continuous(expand=c(0,0),breaks=1:length(unique(usz$month3)),labels=levels(usz$month3),limits = c(0.5,14.5)) +
  theme(axis.text.x = element_text(angle=45,hjust=1))

#' ## Antibody waning

thres = 0.5
usz_repeat = usz_dup %>%
  arrange(patient_id,date) %>%
  group_by(patient_id) %>%
  mutate(dup=n(),
         max_post=max(posterior)) %>%
  filter(!is.na(date),dup>1,max_post>thres) %>%
  mutate(dup_rank=row_number(),
         rank_first_pos=if_else(posterior>thres,dup_rank,as.integer(0))) %>%
  filter(rank_first_pos>0) %>%
  group_by(patient_id) %>%
  mutate(dup=n(),
         dup_rank=row_number(),
         daydiff=as.numeric(date-first(date)),
         logit_post=logit(posterior),
         logit_Spike_plogEC50 = Spike_plogEC50 ,
         logit_NC_plogEC50 =NC_plogEC50,
         logit_RBD_plogEC50 = RBD_plogEC50,
         Spike_plogEC50 = inv.logit(Spike_plogEC50),
         NC_plogEC50 = inv.logit(NC_plogEC50),
         RBD_plogEC50 = inv.logit(RBD_plogEC50)) %>%
  filter(dup>1) 

g1 = stan_lmer(logit_Spike_plogEC50 ~ daydiff + (daydiff|patient_id),data=usz_repeat,iter=5000)
g2 = stan_lmer(logit_NC_plogEC50 ~ daydiff + (daydiff|patient_id),data=usz_repeat,iter=5000)
g3 = stan_lmer(logit_RBD_plogEC50 ~ daydiff + (daydiff|patient_id),data=usz_repeat,iter=5000)
g4 = stan_lmer(logit_post ~ daydiff + (daydiff|patient_id),data=usz_repeat,iter=5000)

daydiff_pred = tibble(daydiff=seq(0,max(usz_repeat$daydiff)))
p1 = tibble(Spike=t(apply(posterior_predict(g1,re.form = ~0,newdata=daydiff_pred),2,qsum)))
p2 = tibble(NC=t(apply(posterior_predict(g2,re.form = ~0,newdata=daydiff_pred),2,qsum)))
p3 = tibble(RBD=t(apply(posterior_predict(g3,re.form = ~0,newdata=daydiff_pred),2,qsum)))
p4 = tibble(Posterior=t(apply(posterior_predict(g4,re.form = ~0,newdata=daydiff_pred),2,qsum)))

ppred = bind_cols(daydiff_pred,p1,p2,p3,p4) %>%
  pivot_longer(2:5) %>% 
  mutate(name=factor(name,
                     levels=c("Spike","NC","RBD","Posterior")))

usz_repeat %>% 
  select(daydiff,starts_with("logit")) %>%
  pivot_longer(3:6) %>%
  mutate(name=factor(name,
                     levels=c("logit_Spike_plogEC50","logit_NC_plogEC50","logit_RBD_plogEC50","logit_post"),
                     labels=c("Spike","NC","RBD","Posterior"))) %>%
  ggplot() +
  geom_point(aes(x=daydiff,y=value,group=patient_id),alpha=.4) +
  geom_line(aes(x=daydiff,y=value,group=patient_id),alpha=.4) +
  geom_ribbon(data=ppred,aes(x=daydiff,ymin=value[,2],ymax=value[,3],fill=name),alpha=.3) +
  geom_line(data=ppred,aes(x=daydiff,y=value[,1],colour=name),size=1) +
  facet_wrap(~name,scales = "free",ncol=4) +
  scale_colour_discrete(guide=FALSE) +
  scale_fill_discrete(guide=FALSE)  +
  labs(x="Difference in days",y="Value")

summary(g1,pars="daydiff",digits=5,probs=c(0.025,0.975))
summary(g2,pars="daydiff",digits=5,probs=c(0.025,0.975))
summary(g3,pars="daydiff",digits=5,probs=c(0.025,0.975))
summary(g4,pars="daydiff",digits=5,probs=c(0.025,0.975))

#+ fig.width=9, fig.height=5

load("waning/post_samples6_2021_01_14_13_37_52.Rdata")

post_samples = samples_model6

day0 = ymd(20200216)
daycuts = ymd(c(20200315,20200419,20200517,20200614,20200719,20200816,20200913,20201018,20201115))
daymax = ymd(20201213)

month_data = agg_serop %>%
  filter(type=="No") %>%
  select(month3,p,sample,sample_size,p_low,p_high) %>%
  filter(!(month3 %in% c("2018","2019","Jan. 2020","Feb. 2020"))) %>%
  mutate(date=rep(c(daycuts,daymax),2)) %>%
  mutate(sample=factor(sample,levels=c("USZ","BDS"),labels=c("USZ sample","BDS sample")))

A2 = summary(samples_model5,pars="comp_A2")[[1]] %>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts2,
         type="Raw seroprevalence") 
A2 = summary(post_samples,pars="IAR")[[1]] %>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts2,
         type="Corrected for waning") %>%  
  filter(t<=289) %>% 
  bind_rows(A2) %>%
  mutate(type=factor(type,levels=c("Raw seroprevalence","Corrected for waning")),
         date=t+day0)


g1 = ggplot(A2) +
  geom_line(aes(x=date,y=`50%`,group=type,linetype=type)) +
  geom_ribbon(aes(x=date,ymin=`2.5%`,ymax=`97.5%`,fill=type),alpha=.3) +
  geom_pointrange(data=month_data,aes(x=date,y=p,ymin=p_low,ymax=p_high,colour=sample),position=position_dodge(6),size=.4) +
  geom_pointrange(data=filter(A2,type=="Infection attack rate",t%in%c(106,289)),aes(x=date,y=`50%`,ymin=`2.5%`,ymax=`97.5%`),colour="steelblue") +
  scale_y_continuous(labels=scales::percent) +
  scale_x_date(date_breaks="month",date_labels="%e %b") +
  labs(x="Date",y="Proportion",colour="Seroprevalence data",fill="Model prediction",linetype="Model prediction") +
  scale_fill_manual(values=c("grey40","steelblue")) +
  scale_colour_manual(values=c("#0000FF","firebrick2")) +
  theme(legend.position=c(.35,.65))


Re = summary(post_samples,pars="R_e")[[1]]%>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts2) %>%
  mutate(date=t+day0)
g2 = ggplot(Re) +
  geom_line(aes(x=date,y=`50%`)) +
  geom_hline(yintercept=1,linetype=2) +
  geom_ribbon(aes(x=date,ymin=`2.5%`,ymax=`97.5%`),fill="orange",alpha=.3) +
  labs(x="Date",y=expression(R[e])) +
  scale_x_date(date_breaks="2 month",date_labels="%e %b") 

post_w = tibble(value=extract(post_samples,pars=c("tau"))[[1]],
                type="Waning half-life (days)",
                prior=rexp(2000,50)) %>%
  mutate(value=log(2)/value,
         prior=log(2)/prior) %>%
  ggplot() +
  geom_density(aes(prior),fill="grey",n=2^16) +
  geom_density(aes(value),fill="pink",alpha=.8,n=2^19) +
  facet_wrap(~type) +
  scale_y_continuous(expand=expansion(c(0,.05))) +
  coord_cartesian(xlim=c(0,100),clip = TRUE) +
  labs(x="Value",y="Density")


post_ifr = tibble(value=extract(post_samples,pars=c("ifr"))[[1]],
                type="Infection-fatality ratio",
                prior=rbeta(2000,5,995)) %>%
  ggplot() +
  geom_density(aes(prior),fill="grey") +
  geom_density(aes(value),fill="pink",alpha=.8) +
  facet_wrap(~type) +
  scale_y_continuous(expand=expansion(c(0,0.05))) +
  scale_x_continuous(labels=scales::percent,breaks=c(0,.01,.02)) +
  coord_cartesian(xlim=c(0,.02)) +
  labs(x="Value",y="Density")

zh_cases = read_csv("shared_Julien/COVID19_Fallzahlen_Kanton_ZH_total.csv") %>%
  mutate(conf=ncumul_conf-lag(ncumul_conf,default = 0)) %>%
  select(date,conf)

ascertainment = summary(post_samples,pars="IAR")[[1]] %>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts2,
         type="Infection attack rate",
         date=day0+t,
         pop=list_usz_bds$popsize) %>%
  left_join(tibble(deaths=list_usz_bds$deaths,
                   t=list_usz_bds$ts))  %>%
  select(date,median=`50%`,lower_bound=`2.5%`,upper_bound=`97.5%`,pop,deaths) %>%
  left_join(zh_cases)  %>%
  mutate(conf=replace_na(conf,0),
         month=month(date)) %>%
  group_by(month) %>%
  summarise(median=max(median),
            lower_bound=max(lower_bound),
            upper_bound=max(upper_bound),
            death=sum(deaths,na.rm = TRUE),
            conf=sum(conf,na.rm = TRUE),
            date=max(date),
            pop=max(pop)) %>%
  mutate(median=((median-lag(median,default=0))*pop),
         lower_bound=((lower_bound-lag(lower_bound,default=0))*pop),
         upper_bound=((upper_bound-lag(upper_bound,default=0))*pop)) 


g3 = ggplot(ascertainment) +
  geom_line(aes(x=date,y=median/conf-1)) +
  geom_ribbon(aes(x=date,ymin=lower_bound/conf-1,ymax=upper_bound/conf-1),fill="seagreen",alpha=.3) +
  coord_cartesian(ylim=c(0,50)) +
  scale_y_continuous(expand=expansion(c(0,0.05))) +
  labs(x="Date",y="Hidden ratio") +
  scale_x_date(date_breaks="2 month",date_labels="%e %b") 


  

plot_grid(g1,
          plot_grid(g2,
                    g3,
                    plot_grid(post_w,post_ifr,labels=c("D","E")),
                    ncol=1,
                    labels=c("B","C","")),
          labels=c("A","")
          ) 

summary(post_samples,pars=c("IAR[106]","IAR[289]","tau","ifr"),probs=c(.5,.025,.975))[[1]] %>%
  as.data.frame() %>%
  rownames_to_column() %>%
  as_tibble() %>%
  mutate(variable=c("IAR (1 Jun)","IAR (1 Dec)","Waning rate","IFR"))%>%
  select(variable,median=`50%`,lower_boumd=`2.5%`,upper_bound=`97.5%`) 
  


#' # Appendix
#' 
#' All seroprevalence estimates:

as.data.frame(agg_serop)

#' All sereprevalence estimates after removing likely hospitalisations due to COVID-19
as.data.frame(agg_serop2)



save(list=ls(),file="all_objects.Rdata")