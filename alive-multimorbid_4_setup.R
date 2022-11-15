
###########################################################################
#
# Project: Prevalence of multimorbidity
#
# Purpose: Prepare multimorbidity data for analysis
#
# Author: Jacqueline Rudolph
#
# Last Update: 14 Nov 2022
#
###########################################################################

packages <- c("tidyverse", "lubridate", "survival", "broom")
for (package in packages) {
  library(package, character.only=T)
}


# Read in data ------------------------------------------------------------

dat <- read_csv("./alive_multimorbid_imputed_visits.csv")

mort <- read_csv("./alive_mort.csv") 


# Manage data -------------------------------------------------------------

# Establish inclusion criteria:
  # Participant must have 2 visits, 2007-2019, that occur within 365 days of each other
  # Because of how outcomes are defined
  # Then, start clock at 2nd observed visit
dat2 <- dat %>% 
  group_by(id) %>% 
  mutate(visdate0 = lag(visdate),
         t_between = (visdate0 %--% visdate)/days(1),
         dead = 0) %>% 
  filter(!is.na(t_between))

ids <- dat2$id[!duplicated(dat2$id) & dat2$t_between>365]
dat3 <- dat2 %>% filter(!(id %in% ids))

# Censor someone when they had no visit for 365 days
dat4 <- dat3 %>% 
  mutate(censor = as.numeric(t_between>365),
         drop = as.numeric(cumsum(cumsum(censor))>=2),
         visdate = ifelse(censor==1, visdate0 + days(1), visdate),
         visdate = as.Date(visdate, origin="1970-01-01")) %>% 
  filter(drop==0) %>% 
  select(-c(t_between, n_vis, drop))

# Prepare mortality data
mort2 <- mort %>% 
  rename(visdate = dod) %>% 
  select(id, visdate) %>% 
  mutate(dead = 1) %>% 
  filter(id %in% unique(dat4$id),
         year(visdate)!=2020)

# Add mortality data
dat5 <- dat4 %>% 
  bind_rows(mort2) %>% 
  arrange(id, visdate) %>% 
  filter(lag(censor, default=0)==0) %>% 
  mutate(visdate0 = ifelse(dead==1, lag(visdate), visdate0),
         censor = ifelse(dead==1, 0, censor))

# Create extra record for individuals who dropped out and never returned
last <- filter(dat5, !duplicated(id, fromLast=T))
ids <- last$id[last$censor==0 & last$dead==0 & year(last$visdate)<2019]
drop <- tibble(id=ids, visdate=ymd("2019-12-31"))

dat6 <- dat5 %>% 
  bind_rows(drop) %>% 
  arrange(id, visdate) %>% 
  mutate(visdate = ifelse(is.na(censor), lag(visdate)+1, visdate),
         visdate = as.Date(visdate, origin="1970-01-01"),
         visdate0 = ifelse(is.na(censor), lag(visdate), visdate0),
         visdate0 = as.Date(visdate0, origin="1970-01-01"),
         dead = ifelse(is.na(censor), 0, dead),
         censor = ifelse(is.na(censor), 1, censor))

# Set up date ranges
  # To estimate state probabilities, clock will start on Jan 1, 2008
  # Exclude visits that occurred prior to 2008
  # If start of interval was before 2008, set to begin on Jan 1, 2008
dat7 <- dat6 %>% 
  mutate(days = (ymd("2008-01-01") %--% visdate)/days(1),
         days0 = (ymd("2008-01-01") %--% visdate0)/days(1),
         days0 = ifelse(days0<0, 0, days0),
         state = case_when(
           censor==1 ~ 0,
           dead==1 ~ 3,
           n_cond<2 ~ 1,
           T ~ 2)) %>% 
  filter(days>=0)


# Example survfit code ----------------------------------------------------

# Suppose we ignored the control piece for a second
# This is how I would specify the call to survfit

# If we don't specify the initial state, the model artificially creats a transition 
  # from (s0) to a participant's first state
dat8 <- dat7 %>% 
  mutate(state0 = lag(state),
         state0 = ifelse(is.na(state0), state, state0))

res <- tidy(survfit(Surv(days0, days, as.factor(state)) ~ 1, data=dat8, id=id, istate=state0)) %>% 
  mutate(time = as.Date(time, origin=ymd("2008-01-01"))) %>% 
  filter(state!="(s0)") 

res_wide <- res %>% 
  pivot_wider(id_cols="time", names_from="state", names_prefix="risk", values_from="estimate") %>% 
  mutate(no_multi=risk1,
         multi=risk1+risk2,
         death=risk1+risk2+risk3)

ggplot(data=res_wide, aes(x=time)) +
  geom_area(aes(y=death, fill="Death")) +
  geom_area(aes(y=multi, fill="Multimorbidity")) +
  geom_area(aes(y=no_multi, fill="No multimorbidity"))


