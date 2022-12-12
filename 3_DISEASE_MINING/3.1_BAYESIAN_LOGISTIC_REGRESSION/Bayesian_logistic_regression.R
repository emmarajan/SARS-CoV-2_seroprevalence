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

## merge with posterior serology 
icd10 = usz %>%
  dplyr::select(patient_id,sex,age,posterior) %>%
  group_by(patient_id) %>%
  summarise(posterior=max(posterior),sex=max(sex),age=max(age)) %>%
  inner_join(icd10)


## lasso
icd10_logit =  icd10 %>%
  dplyr::select(-patient_id) %>%
  mutate(posterior=ifelse(posterior==0,0+0.001,posterior),
         posterior=ifelse(posterior==1,1-0.001,posterior),
         posterior=logit(posterior),
         age=(age-mean(age))/sd(age)) 


# normal prior
m1_lin = stan_glm(posterior~.,
                  family=gaussian(),
                  prior=normal(0,10),
                  prior_intercept = normal(0, 10),
                  data=icd10_logit)
summary(m1_lin)
s1_lin = summary(m1_lin,probs = c(0.025,.5,0.975)) %>%
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>%
  filter(!(rowname %in% c("log-posterior","mean_PPD","(Intercept)","sigma"))) %>%
  mutate(rowname=gsub("TRUE","",rowname),
         OR=exp(`50%`),
         OR_low=exp(`2.5%`),
         OR_high=exp(`97.5%`),
         type="Uninformative prior") %>%
  arrange(-abs(`50%`))

# lasso prior with 1 df
m1_lasso = stan_glm(posterior~.,
                    family=gaussian(),
                    prior=lasso(df=1),
                    data=icd10_logit)
summary(m1_lasso)
s1_lasso = summary(m1_lasso,probs = c(0.025,.5,0.975)) %>%
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>%
  filter(!(rowname %in% c("log-posterior","mean_PPD","(Intercept)","sigma"))) %>%
  mutate(rowname=gsub("TRUE","",rowname),
         OR=exp(`50%`),
         OR_low=exp(`2.5%`),
         OR_high=exp(`97.5%`),
         type="LASSO") %>%
  arrange(-abs(`50%`))


# horseshoe prior (Piironen and Vehtari (2017))  

# which recommends setting the global_scale argument equal to the ratio of the expected
# number of non-zero coefficients to the expected number of zero coefficients, divided by
# the square root of the number of observations.
globalscale = (10/189)/sqrt(199)
m1_hs = stan_glm(posterior~.,
                 family=gaussian(),
                 prior=hs(global_scale=globalscale),
                 data=icd10_logit)
summary(m1_hs)
s1_hs = summary(m1_hs,probs = c(0.025,.5,0.975)) %>%
  as.data.frame() %>%
  rownames_to_column() %>% 
  as_tibble() %>%
  filter(!(rowname %in% c("log-posterior","mean_PPD","(Intercept)","sigma"))) %>%
  mutate(rowname=gsub("TRUE","",rowname),
         OR=exp(`50%`),
         OR_low=exp(`2.5%`),
         OR_high=exp(`97.5%`),
         type="Regularized horseshoe") %>%
  arrange(-abs(`50%`))


post_icd10 = bind_rows(s1_lasso,s1_hs,s1_lin) %>%
  mutate(rowname=factor(rowname,levels=s1_hs$rowname[order(s1_hs$`50%`)]),
         type=factor(type,levels=c("Uninformative prior","LASSO","Regularized horseshoe")))
saveRDS(post_icd10,file="post_icd10.rds")

#' 
post_icd10 = readRDS("post_icd10.rds")
ranks = filter(post_icd10,type=="Regularized horseshoe") %>%
  arrange(-`50%`) %>%
  mutate(rank=row_number()) %>%
  select(ICD10=rowname,rank)
post_icd10 %>%
  mutate(OR=paste0(sprintf("%.2f",OR)," (",sprintf("%.2f",OR_low),"-",sprintf("%.2f",OR_high),")")) %>%
  select(ICD10=rowname,OR,type) %>%
  pivot_wider(names_from=type,values_from=OR) %>%
  left_join(ranks) %>%
  arrange(rank) %>% 
  filter(rank %in% c(1:10,(max(rank)-9):max(rank))) %>%
  select(rank,ICD10,`Uninformative prior`,LASSO,`Regularized horseshoe`)

#+ fig.width=9, fig.height=18
post_icd10 %>%
  mutate(type=factor(type,labels=c("Unregularized","LASSO","Horseshoe"))) %>%
  ggplot() +
  geom_pointrange(aes(x=rowname,y=OR,ymin=OR_low,ymax=OR_high,colour=type)) +
  geom_hline(yintercept=1,linetype=2) +
  scale_colour_discrete(guide=FALSE) +
  scale_y_log10() +
  coord_flip() +
  facet_wrap(~type,ncol=3) +
  labs(y="Odds ratio (95% credible interval)",x="ICD 10 code",title="Association between ICD-10 and SARS-CoV-2 seropositivity in USZ sample.")



