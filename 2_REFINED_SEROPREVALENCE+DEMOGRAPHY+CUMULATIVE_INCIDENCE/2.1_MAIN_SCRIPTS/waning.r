library(tidyverse)
library(lubridate)
library(cowplot)
theme_set(theme_bw())
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
require(bayesplot)
require(posterior)

# setup ---------------------
source("dm.r")
source("serop.r")

# surveillance data
zh_deaths = read_csv("waning/COVID19_Fallzahlen_Kanton_ZH_altersklassen_geschlecht_dec.csv") %>%
  group_by(week_number=Week) %>%
  summarize(deaths=sum(NewDeaths)) %>%
  bind_rows(tibble(week_number=7:8,deaths=0)) %>%
  arrange(week_number)

# seroprevalence data
day0 = ymd(20200216)
daycuts = ymd(c(20200315,20200419,20200517,20200614,20200719,20200816,20200913,20201018,20201115))
daymax = ymd(20201213)
dayknots = ymd(c(as.character(day0),20200317,20200511,20200920,20201031,as.character(daymax)))
period = tibble(date=seq.Date(day0,daymax,1)) %>%
  mutate(day_number=row_number(),
         date_period=cut(date,breaks=c(day0-1,day0,daycuts,daymax),right=TRUE,labels=FALSE),
         knots=if_else(date %in% dayknots,1,0),
         weekday=weekdays(date),
         week_number=if_else(weekday=="Sunday",week(date),as.double(NA))) %>%
  left_join(zh_deaths) %>%
  ungroup()

period_ts = period  %>%
  filter(!is.na(deaths)) %>%
  mutate(rank_ts=row_number())

period_knots = period %>%
  filter(knots==1)

month_data = agg_serop %>%
  filter(type=="No") %>%
  select(month3,p,sample,sample_size) %>%
  pivot_wider(names_from=sample,values_from=c(p,sample_size)) %>%
  filter(!(month3 %in% c("2018","2019","Jan. 2020","Feb. 2020"))) %>%
  mutate(date=c(daycuts,daymax)) %>%
  left_join(period) %>%
  mutate(k_BDS=round(p_BDS*sample_size_BDS),
         k_USZ=round(p_USZ*sample_size_USZ)) %>%
  left_join(period_ts)

# Fom Linton et al
linton_mean = 20.2
linton_sd = 11.6
# Get lognormal parameters from mean and sd
get_par_lnorm = function(m,s) {
  mu = log(m) - 1/2*log((s/m)^2+1)
  sigma = sqrt(log((s/m)^2+1))
  return(list(mu=mu,sigma=sigma))
}
linton_pars = get_par_lnorm(linton_mean,linton_sd)
# Discretize
gamma = numeric(9)
ss = seq(1,64,by=7)
for(i in 1:9) {
  gamma[i] = plnorm(ss[i+1],linton_pars$mu,linton_pars$sigma)-plnorm(ss[i],linton_pars$mu,linton_pars$sigma)
}
gamma=gamma/sum(gamma)


# format data -------------------
list_usz_bds = list(
  N = nrow(period_ts),
  ts = period_ts$day_number, 
  popsize = sum(demog_zh$n_pop),
  num_usz = nrow(month_data),
  t_usz = month_data$rank_ts,
  k_usz = month_data$k_USZ,
  n_usz = month_data$sample_size_USZ,
  num_bds = nrow(month_data),
  t_bds = month_data$rank_ts,
  k_bds = month_data$k_BDS,
  n_bds = month_data$sample_size_BDS,
  num_knots=sum(period_knots$knots),
  knots=period_knots$day_number,
  D = length(period_ts$deaths),
  deaths = period_ts$deaths,
  
  N2 = nrow(period),
  ts2 = period$day_number,
  
  G=length(gamma),
  p_gamma=gamma,
  
  atol = 1e-4,
  rtol = 1e-4,
  max_num_steps = 5000,
  inference = 1
)
save(list_usz_bds,file="waning/list_usz_bds.Rdata")


# after adjusting on age sex

month_data2 = agg_serop %>%
  filter(type=="On age and sex") %>%
  select(month3,p,sample,sample_size) %>%
  pivot_wider(names_from=sample,values_from=c(p,sample_size)) %>%
  filter(!(month3 %in% c("2018","2019","Jan. 2020","Feb. 2020"))) %>%
  mutate(date=c(daycuts,daymax)) %>%
  left_join(period) %>%
  mutate(k_BDS=round(p_BDS*sample_size_BDS),
         k_USZ=round(p_USZ*sample_size_USZ)) %>%
  left_join(period_ts)

month_data_2_bds = filter(month_data2,!is.na(k_BDS))

list_usz_bds2 = list(
  N = nrow(period_ts),
  ts = period_ts$day_number, 
  popsize = sum(demog_zh$n_pop),
  num_usz = nrow(month_data2),
  t_usz = month_data2$rank_ts,
  k_usz = month_data2$k_USZ,
  n_usz = month_data2$sample_size_USZ,
  num_bds = nrow(month_data_2_bds),
  t_bds = month_data_2_bds$rank_ts,
  k_bds = month_data_2_bds$k_BDS,
  n_bds = month_data_2_bds$sample_size_BDS,
  num_knots=sum(period_knots$knots),
  knots=period_knots$day_number,
  D = length(period_ts$deaths),
  deaths = period_ts$deaths,
  
  N2 = nrow(period),
  ts2 = period$day_number,
  
  G=length(gamma),
  p_gamma=gamma,
  
  atol = 1e-4,
  rtol = 1e-4,
  max_num_steps = 5000,
  inference = 1
)
save(list_usz_bds2,file="waning/list_usz_bds2.Rdata")


# cluster -------------------------------
system("scp waning/waning_cluster.r UBELIX:/home/ubelix/ispm/jr18s506/projects/covid_serop_usz/waning/.")
system("scp waning/*.stan UBELIX:/home/ubelix/ispm/jr18s506/projects/covid_serop_usz/waning/.")
system("scp waning/list_usz_bds.Rdata UBELIX:/home/ubelix/ispm/jr18s506/projects/covid_serop_usz/waning/.")
system("scp waning/list_usz_bds2.Rdata UBELIX:/home/ubelix/ispm/jr18s506/projects/covid_serop_usz/waning/.")
system("scp waning/run_waning.sh UBELIX:/home/ubelix/ispm/jr18s506/projects/covid_serop_usz/waning/.")

system("scp UBELIX:/home/ubelix/ispm/jr18s506/projects/covid_serop_usz/waning/post_samples4_2021_01_08_17_30_46.Rdata waning/.")
system("scp UBELIX:/home/ubelix/ispm/jr18s506/projects/covid_serop_usz/waning/post_samples5_2021_01_08_18_37_37.Rdata waning/.")
system("scp UBELIX:/home/ubelix/ispm/jr18s506/projects/covid_serop_usz/waning/post_samples6_2021_01_10_10_15_22.Rdata waning/.")
system("scp UBELIX:/home/ubelix/ispm/jr18s506/projects/covid_serop_usz/waning/post_samples5bis_2021_01_11_00_49_54.Rdata waning/.")
system("scp UBELIX:/home/ubelix/ispm/jr18s506/projects/covid_serop_usz/waning/post_samples6_2021_01_14_13_37_52.Rdata waning/.")
system("scp UBELIX:/home/ubelix/ispm/jr18s506/projects/covid_serop_usz/waning/post_samples7_2021_01_12_07_47_51.Rdata waning/.")
system("scp UBELIX:/home/ubelix/ispm/jr18s506/projects/covid_serop_usz/waning/post_samples9_2021_01_13_18_15_00.Rdata waning/.")


# load results ---------------------------
# l = load("waning/post_samples4_2021_01_08_17_30_46.Rdata")
# l = load("waning/post_samples5_2021_01_08_18_37_37.Rdata")
# 
# load("waning/post_samples6_2021_01_10_10_15_22.Rdata")

load("waning/post_samples5bis_2021_01_11_00_49_54.Rdata")
load("waning/post_samples6_2021_01_14_13_37_52.Rdata")
load("waning/post_samples7_2021_01_12_07_47_51.Rdata")
load("waning/post_samples9_2021_01_13_18_15_00.Rdata")



# model 3 -------------------------------

# simulate from prior
# list_usz_bds$inference = 0
# samples_model3 = stan(file="waning/model3.stan",
#                       data=list_usz_bds,
#                       chains=1,iter=200,init=.5)

print(samples_model3,pars=c("beta_0","rho","I0_raw","I0","gamma","sigma","tau"))
stan_trace(samples_model3,pars=c("beta_0","rho","I0_raw","I0","gamma","sigma","tau"))

output_model3 = summary(samples_model3,pars="comp_R2")[[2]][,,4] %>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts2) %>%
  left_join(tibble(t=list_usz_bds$ts[list_usz_bds$t_usz],
                   p_usz=list_usz_bds$k_usz/list_usz_bds$n_usz,
                   p_bds=list_usz_bds$k_bds/list_usz_bds$n_bds))
ggplot(output_model3) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  geom_point(aes(x=t,y=p_usz),colour="seagreen") +
  geom_point(aes(x=t,y=p_bds),colour="firebrick") +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Week",y="Seroprevalence")

Re_model3 = summary(samples_model3,pars="R_e")[[1]]%>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts2)
ggplot(Re_model3) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  labs(x="Week",y="R_e")

inc_D_model3 = summary(samples_model3,pars="incidence_deaths")[[1]]%>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts,
         deaths=list_usz_bds$deaths)
ggplot(inc_D_model3) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  geom_point(aes(x=t,y=deaths)) +
  labs(x="Week",y="N") 



# model 4 -------------------------------

# simulate from prior
# list_usz_bds$inference = 0
# samples_model4 = stan(file="waning/model4.stan",
#                       data=list_usz_bds,
#                       chains=1,iter=200,init=.5)

check_hmc_diagnostics(samples_model4)
print(samples_model4,pars=c("beta_0","rho","I0_raw","I0","gamma","sigma","tau","waning_period"))
stan_trace(samples_model4,pars=c("beta_0","rho","I0_raw","I0","gamma","sigma","tau"))

output_model4 = summary(samples_model4,pars="comp_R2",probs=c(.5,.05,.95))[[1]]%>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts2) %>%
  left_join(tibble(t=list_usz_bds$ts[list_usz_bds$t_usz],
                   p_usz=list_usz_bds$k_usz/list_usz_bds$n_usz,
                   p_bds=list_usz_bds$k_bds/list_usz_bds$n_bds))
ggplot(output_model4) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`5%`,ymax=`95%`),alpha=.3) +
  geom_point(aes(x=t,y=p_usz),colour="seagreen") +
  geom_point(aes(x=t,y=p_bds),colour="firebrick") +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Week",y="Seroprevalence")


output_model4 = summary(samples_model4,pars="comp_R")[[1]] %>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts) %>%
  left_join(tibble(t=list_usz_bds$ts[list_usz_bds$t_usz],
                   p_usz=list_usz_bds$k_usz/list_usz_bds$n_usz,
                   p_bds=list_usz_bds$k_bds/list_usz_bds$n_bds))
ggplot(output_model4) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  geom_point(aes(x=t,y=p_usz),colour="seagreen") +
  geom_point(aes(x=t,y=p_bds),colour="firebrick") +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Week",y="Seroprevalence")

Re_model4 = summary(samples_model4,pars="R_e")[[1]]%>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts2)
ggplot(Re_model4) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`25%`,ymax=`75%`),alpha=.3) +
  labs(x="Week",y="R_e")

inc_D_model4 = summary(samples_model4,pars="inc_D")[[1]]%>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts,
         deaths=list_usz_bds$deaths)
ggplot(inc_D_model4) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  geom_point(aes(x=t,y=deaths)) +
  labs(x="Week",y="N") 





# model 5 -------------------------------

# simulate from prior
list_usz_bds$inference = 0
samples_model5 = stan(file="waning/model5.stan",
                      data=list_usz_bds,
                      chains=1,iter=200,init=.5)

print(samples_model5,pars=c("beta_0","rho","I0_raw","I0","tau","ifr","waning_period"))

output_model5 = summary(samples_model5,pars="comp_A2")[[1]] %>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts2) %>%
  left_join(tibble(t=list_usz_bds$ts[list_usz_bds$t_usz],
                   p_usz=list_usz_bds$k_usz/list_usz_bds$n_usz,
                   p_bds=list_usz_bds$k_bds/list_usz_bds$n_bds))
ggplot(output_model5) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  geom_point(aes(x=t,y=p_usz),colour="seagreen") +
  geom_point(aes(x=t,y=p_bds),colour="firebrick") +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Week",y="Seroprevalence")

Re_model5 = summary(samples_model5,pars="R_e")[[1]]%>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts2)
ggplot(Re_model5) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  labs(x="Week",y="R_e")

inc_D_model5 = summary(samples_model5,pars="inc_D")[[1]]%>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts) %>%
  left_join(tibble(deaths=list_usz_bds$deaths,
                   t=list_usz_bds$ts)) 
ggplot(inc_D_model5) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  geom_point(aes(x=t,y=deaths)) +
  labs(x="Week",y="N") 

sum(list_usz_bds$deaths)/.01
print(samples_model5,pars="incidence_infections[300]")


summary(samples_model5,pars="y2")[[1]] %>%
  as_tibble() %>%
  mutate(comp=rep(1:8,list_usz_bds$N2),
         t=rep(list_usz_bds$ts2,each=8)) %>%
  ggplot() +
  geom_line(aes(x=t,y=`50%`,colour=factor(comp))) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`,fill=factor(comp)),alpha=.3) +
  facet_grid(comp~.,scales="free")


# model 6 -------------------------------


# simulate from prior
list_usz_bds2$inference = 0
samples_model6 = stan(file="waning/model6.stan",
                      data=list_usz_bds2,
                      chains=1,iter=200,init=.5,control=list(max_treedepth=13))

print(samples_model6,pars=c("beta_0","rho","I0_raw","I0","tau","ifr","waning_period"))
stan_trace(samples_model6,pars=c("beta_0","rho","I0_raw","I0","tau","ifr","waning_period"))


output_model6 = summary(samples_model6,pars="comp_A")[[1]] %>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts) %>%
  left_join(tibble(t=list_usz_bds$ts[list_usz_bds$t_usz],
                   p_usz=list_usz_bds$k_usz/list_usz_bds$n_usz)) %>%
  left_join(tibble(t=list_usz_bds$ts[list_usz_bds$t_bds],
                   p_bds=list_usz_bds$k_bds/list_usz_bds$n_bds))
ggplot(output_model6) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  geom_point(aes(x=t,y=p_usz),colour="seagreen") +
  geom_point(aes(x=t,y=p_bds),colour="firebrick") +
  geom_vline(xintercept=list_usz_bds$knots,linetype=2) +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Week",y="Seroprevalence")

Re_model6 = summary(samples_model6,pars="R_e")[[1]]%>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts2)
ggplot(Re_model6) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  geom_vline(xintercept=list_usz_bds$knots,linetype=2) +
  labs(x="Week",y="R_e")

inc_D_model6 = summary(samples_model6,pars="inc_D")[[1]]%>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts) %>%
  left_join(tibble(deaths=list_usz_bds$deaths,
                   t=list_usz_bds$ts)) 
ggplot(inc_D_model6) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  geom_vline(xintercept=list_usz_bds$knots,linetype=2) +
  geom_point(aes(x=t,y=deaths)) +
  labs(x="Week",y="N") 

summary(samples_model6,pars="y2")[[1]] %>%
  as_tibble() %>%
  mutate(comp=rep(1:8,list_usz_bds$N2),
         t=rep(list_usz_bds$ts2,each=8)) %>%
  ggplot() +
  geom_line(aes(x=t,y=`50%`,colour=factor(comp))) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`,fill=factor(comp)),alpha=.3) +
  facet_grid(comp~.,scales="free")



# model 7 -------------------------------

# simulate from prior
list_usz_bds$inference = 0
samples_model7 = stan(file="waning/model7.stan",
                      data=list_usz_bds,
                      chains=1,iter=200,init=.5,control=list(max_treedepth=13))

print(samples_model7,pars=c("beta_0","rho","I0_raw","I0","tau","ifr","waning_period"))
stan_trace(samples_model7,pars=c("beta_0","rho","I0_raw","I0","tau","ifr","waning_period"))


output_model7 = summary(samples_model7,pars="comp_A")[[1]] %>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts) %>%
  left_join(tibble(t=list_usz_bds$ts[list_usz_bds$t_usz],
                   p_usz=list_usz_bds$k_usz/list_usz_bds$n_usz,
                   p_bds=list_usz_bds$k_bds/list_usz_bds$n_bds))
ggplot(output_model7) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  geom_point(aes(x=t,y=p_usz),colour="seagreen") +
  geom_point(aes(x=t,y=p_bds),colour="firebrick") +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Week",y="Seroprevalence")

Re_model7 = summary(samples_model7,pars="R_e")[[1]]%>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts2)
ggplot(Re_model7) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  labs(x="Week",y="R_e")

inc_D_model7 = summary(samples_model7,pars="comp_D")[[1]]%>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts) %>%
  left_join(tibble(deaths=list_usz_bds$deaths,
                   t=list_usz_bds$ts)) 
ggplot(inc_D_model7) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  # geom_point(aes(x=t,y=deaths)) +
  labs(x="Week",y="N") 


summary(samples_model6,pars="y2")[[1]] %>%
  as_tibble() %>%
  mutate(comp=rep(1:8,list_usz_bds$N2),
         t=rep(list_usz_bds$ts2,each=8)) %>%
  ggplot() +
  geom_line(aes(x=t,y=`50%`,colour=factor(comp))) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`,fill=factor(comp)),alpha=.3) +
  facet_grid(comp~.,scales="free")



# model 8 -------------------------------

# simulate from prior
list_usz_bds$inference = 0
samples_model8 = stan(file="waning/model8.stan",
                      data=list_usz_bds,
                      chains=1,iter=200,init=.5,control=list(max_treedepth=13))

print(samples_model8,pars=c("beta_0","rho","I0_raw","I0","tau","ifr","waning_period"))
stan_trace(samples_model8,pars=c("beta_0","rho","I0_raw","I0","tau","ifr","waning_period"))


output_model8 = summary(samples_model8,pars="comp_A")[[1]] %>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts) %>%
  left_join(tibble(t=list_usz_bds$ts[list_usz_bds$t_usz],
                   p_usz=list_usz_bds$k_usz/list_usz_bds$n_usz,
                   p_bds=list_usz_bds$k_bds/list_usz_bds$n_bds))
ggplot(output_model8) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  geom_point(aes(x=t,y=p_usz),colour="seagreen") +
  geom_point(aes(x=t,y=p_bds),colour="firebrick") +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Week",y="Seroprevalence")

Re_model8 = summary(samples_model8,pars="R_e")[[1]]%>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts2)
ggplot(Re_model8) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  labs(x="Week",y="R_e")

inc_D_model8 = summary(samples_model8,pars="incidence_deaths")[[1]]%>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts2) %>%
  left_join(tibble(deaths=list_usz_bds$deaths,
                   t=list_usz_bds$ts)) 
ggplot(inc_D_model8) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  geom_point(aes(x=t,y=deaths)) +
  labs(x="Week",y="N") 


# model 9 -------------------------------

# simulate from prior
list_usz_bds$inference = 0
samples_model9 = stan(file="waning/model9.stan",
                      data=list_usz_bds,
                      chains=1,iter=200,init=.5,control=list(max_treedepth=13))

print(samples_model9,pars=c("beta_0","rho","I0_raw","I0","tau","ifr","waning_period"))
stan_trace(samples_model9,pars=c("beta_0","rho","I0_raw","I0","tau","ifr","waning_period"))


output_model9 = summary(samples_model9,pars="comp_A")[[1]] %>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts) %>%
  left_join(tibble(t=list_usz_bds$ts[list_usz_bds$t_usz],
                   p_usz=list_usz_bds$k_usz/list_usz_bds$n_usz,
                   p_bds=list_usz_bds$k_bds/list_usz_bds$n_bds))
ggplot(output_model9) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  geom_point(aes(x=t,y=p_usz),colour="seagreen") +
  geom_point(aes(x=t,y=p_bds),colour="firebrick") +
  geom_vline(xintercept=list_usz_bds$knots,linetype=2) +
  scale_y_continuous(labels=scales::percent) +
  labs(x="Week",y="Seroprevalence")

Re_model9 = summary(samples_model9,pars="R_e")[[1]]%>%
  as_tibble() %>%
  mutate(t=list_usz_bds$ts2)
ggplot(Re_model9) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  geom_vline(xintercept=list_usz_bds$knots,linetype=2) +
  labs(x="Week",y="R_e")

inc_D_model9 = summary(samples_model9,pars="comp_diffM")[[1]]%>%
  as_tibble() %>%
  slice_head(n=length(list_usz_bds$ts)) %>%
  mutate(t=list_usz_bds$ts) %>%
  left_join(tibble(deaths=list_usz_bds$deaths,
                   t=list_usz_bds$ts)) 
ggplot(inc_D_model9) +
  geom_line(aes(x=t,y=`50%`)) +
  geom_ribbon(aes(x=t,ymin=`2.5%`,ymax=`97.5%`),alpha=.3) +
  geom_vline(xintercept=list_usz_bds$knots,linetype=2) +
  geom_point(aes(x=t,y=deaths)) +
  labs(x="Week",y="N") 
