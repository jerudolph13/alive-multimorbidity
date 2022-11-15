
###########################################################################
#
# Project: Prevalence of multimorbidity
#
# Purpose: Combine imputed data sets
#
# Author: Jacqueline Rudolph
#
# Last Update: 14 Nov 2022
#
###########################################################################

library("tidyverse")
library("lubridate")


# Read in data ------------------------------------------------------------

base <- read_csv("../data/alive_use_impute.csv")

diab <- read_csv(file="../data/alive_diab_impute.csv")

hbp <- read_csv(file="../data/alive_hbp_impute.csv")

liver <- read_csv(file="../data/alive_liver_impute.csv")

lung <- read_csv(file="../data/alive_lung_impute.csv")

renal <- read_csv(file="../data/alive_renal_impute.csv")

dep <- read_csv(file="../data/alive_dep_impute.csv")

cancer <- read_csv("../data/alive_cancer.csv") %>% 
  select(id, dx_date) %>% 
  # Grab first cancers only
  filter(!duplicated(id))


# Merge data --------------------------------------------------------------

dat <- base %>% 
  left_join(hbp, by=c("id", "n_vis", "visdate")) %>% 
  left_join(diab, by=c("id", "n_vis", "visdate")) %>% 
  left_join(lung, by=c("id", "n_vis", "visdate")) %>% 
  left_join(renal, by=c("id", "n_vis", "visdate")) %>% 
  left_join(liver, by=c("id", "n_vis", "visdate")) %>% 
  left_join(dep, by=c("id", "n_vis", "visdate")) %>% 
  left_join(cancer, by="id")


# Create indicators of disease --------------------------------------------

dat2 <- dat %>% 
  # Diabetes: HbA1c > 6.5% or medication use
  mutate(hi_hba1c = as.numeric((hba1c>6.5)),

         # Obstructive lung disease: FEV/FVC <=0.70
         poor_lung_func = as.numeric((fev_fvc<=0.70)),
         
         # Kidney dysfunction: urine protein-creatinine ratio >200 or GFR <60
         renal_dysfunc = as.numeric((uprt_crt>200) | (gfr<60)),
         
         # Liver fibrosis: liver stiffness >=9.3
         fibrosis = as.numeric(fbscan>=9.3),
         
         # Liver cirrhosis: liver stiffness >=12.3
         cirrhosis = as.numeric(fbscan>=12.3),
         
         # Hypertension: Dias BP >=90, Sys BP >=140, or treatment in last 6 months
         hi_bp = as.numeric((bpdias>=90) | (bpsys>=140)),

         # Cancer: any diagnosis in yr or before
         cancer = ifelse(is.na(dx_date), 0, as.numeric(dx_date<=visdate)),
         
         # Depression: CESD>=23 or self-reported treatment in last 6 months
         cesd23 = as.numeric(cesdtot>=23))


# Calculate time between elevated measurements ----------------------------

hba1c <- dat2 %>% 
  filter(hi_hba1c==1) %>% 
  group_by(id) %>% 
  mutate(t_between = (lag(visdate) %--% visdate)/years(1),
         flag = t_between<1 & !is.na(t_between),
         diab_idx = cumsum(cumsum(flag))) %>% 
  filter(diab_idx==1) %>% 
  select(id, visdate, diab_idx) 

hbp <- dat2 %>% 
  filter(hi_bp==1) %>% 
  group_by(id) %>% 
  mutate(t_between = (lag(visdate) %--% visdate)/years(1),
         flag = t_between<1 & !is.na(t_between),
         hbp_idx = cumsum(cumsum(flag))) %>% 
  filter(hbp_idx==1) %>% 
  select(id, visdate, hbp_idx) 

lung_func <- dat2 %>% 
  filter(poor_lung_func==1) %>% 
  group_by(id) %>% 
  mutate(t_between = (lag(visdate) %--% visdate)/years(1),
         flag = t_between<1 & !is.na(t_between),
         lung_idx = cumsum(cumsum(flag))) %>% 
  filter(lung_idx==1) %>% 
  select(id, visdate, lung_idx) 

renal_func <- dat2 %>% 
  filter(renal_dysfunc==1) %>% 
  group_by(id) %>% 
  mutate(t_between = (lag(visdate) %--% visdate)/years(1),
         flag = t_between<1 & !is.na(t_between),
         renal_idx = cumsum(cumsum(flag))) %>% 
  filter(renal_idx==1) %>% 
  select(id, visdate, renal_idx) 

liver_fib <- dat2 %>% 
  filter(fibrosis==1) %>% 
  group_by(id) %>% 
  mutate(t_between = (lag(visdate) %--% visdate)/years(1),
         flag = t_between<1 & !is.na(t_between),
         liver_idx = cumsum(cumsum(flag))) %>% 
  filter(liver_idx==1) %>% 
  select(id, visdate, liver_idx) 

cesd <- dat2 %>% 
  filter(cesd23==1) %>% 
  group_by(id) %>% 
  mutate(t_between = (lag(visdate) %--% visdate)/years(1),
         flag = t_between<1 & !is.na(t_between),
         dep_idx = cumsum(cumsum(flag))) %>% 
  filter(dep_idx==1) %>% 
  select(id, visdate, dep_idx) 

dat3 <- dat2 %>% 
  left_join(hba1c, by=c("id", "visdate")) %>% 
  left_join(hbp, by=c("id", "visdate")) %>% 
  left_join(lung_func, by=c("id", "visdate")) %>% 
  left_join(renal_func, by=c("id", "visdate")) %>% 
  left_join(liver_fib, by=c("id", "visdate")) %>% 
  left_join(cesd, by=c("id", "visdate"))
  

# Define multimorbidity ---------------------------------------------------
# Someone flagged as having condition if 2 measurements within 1 year

dat4 <- dat3 %>%
  group_by(id) %>% 
  mutate(across(c(diab_idx, hbp_idx, lung_idx, renal_idx, liver_idx, dep_idx), ~ ifelse(is.na(.x), 0, .x)),
         diabetes = as.numeric((cumsum(diab_idx + diabtx6m))>=1),
         hypertension = as.numeric((cumsum(hbp_idx + hbptx6m))>=1),
         lung_disease = as.numeric(cumsum(lung_idx)>=1),
         renal_disease = as.numeric(cumsum(renal_idx)>=1),
         depression = as.numeric((cumsum(dep_idx + deptx6m))>=1),
         liver_disease = as.numeric(cumsum(liver_idx)>=1),
         
         # Multimorbidity
         n_cond = hypertension +
           diabetes +
           lung_disease + 
           renal_disease +
           liver_disease +
           depression +
           cancer +
           hiv) %>% 
  select(-c(yr, bpsys, bpdias, hba1c, fev_fvc, gfr, uprt_crt, fbscan, diab_idx, 
            hbp_idx, lung_idx, renal_idx, liver_idx, dep_idx, dx_date))


# Output data -------------------------------------------------------------

write_csv(dat4, file="../data/alive_multimorbid_imputed_visits.csv")



