
###################################################################################################################
#
# Project: Impute missing data in ALIVE
#
# Purpose: Combine imputed data sets
#
# Author: Jacqueline Rudolph
#
# Last Update: 12 Apr 2022
#
###################################################################################################################

library("tidyverse")
library("lubridate")

# Function to make data annual
make.annual <- function(dat, ...) {
  output <- data.frame()
  n <- 1
  for (i in list(...)) {
    dat2 <- select(dat, id, yr, !!i) %>% 
      group_by(id, yr) %>% 
      summarize(summ = ifelse(is.nan(mean(!!i, na.rm=T)), NA, mean(!!i, na.rm=T)))
    names(dat2) <- c("id", "yr", i)
    if (n==1) {
      output <- dat2
    } else {
      output <- merge(output, dat2)
    }
    n <- n+1
  }
  return(output)
}


# Read in data ------------------------------------------------------------

# Raw data
base <- read_csv("../data/alive_base.csv") %>% 
  select(-c(visdate, visit, bincome, bevcig, bevalc, bcurrpa, beduc))

socdem <- read_csv("../data/alive_socdem.csv") %>% 
  mutate(yr = year(visdate)) %>% 
  select(-c(visdate, enroll_period, visit, race)) 

hiv <- read_csv("../data/alive_hiv.csv") %>%
  mutate(yr = year(visdate)) %>% 
  select(-c(visdate, enroll_period, visit)) 

drug <- read_csv("../data/alive_drug.csv") %>% 
  mutate(yr = year(visdate)) %>% 
  select(-c(visdate, enroll_period, visit))

cancer <- read_csv("../data/alive_cancer.csv") %>% 
  select(id, dx_date) %>% 
  # Grab first cancers only
  mutate(first = as.numeric(!duplicated(id, fromLast=F))) %>% 
  filter(first==1) %>% 
  select(-first)

# Imputed data
bmi <- read_csv(file="../data/alive_bmi_impute.csv")

diab <- read_csv(file="../data/alive_diab_impute.csv")

hbp <- read_csv(file="../data/alive_hbp_impute.csv")

liver <- read_csv(file="../data/alive_liver_impute.csv")

lung <- read_csv(file="../data/alive_lung_impute.csv")

renal <- read_csv(file="../data/alive_renal_impute.csv")

dep <- read_csv(file="../data/alive_dep_impute.csv")


# Make data annual --------------------------------------------------------

# HIV data
hiv.annual <- make.annual(hiv, sym("hiv"), sym("hrt6m"), sym("hcvvis"), sym("hivvl")) %>% 
  mutate(across(c(hiv, hrt6m, hcvvis), ~ ceiling(.x)))

# Sociodemographic data
socdem.annual <- make.annual(socdem, sym("age"), sym("m0f1"), sym("black"), sym("anyhomels"), sym("inclt5k"),
                             sym("work")) %>% 
  mutate(across(c(black, m0f1, anyhomels, inclt5k, work), ~ ceiling(.x)))

# Drug use data
drug.annual <- make.annual(drug, sym("curuser"), sym("alcheavy"), sym("cigyn")) %>% 
  mutate(across(c(curuser, alcheavy, cigyn), ~ ceiling(.x)))


# Merge data --------------------------------------------------------------

# Merge raw data
annual <- left_join(hiv.annual, socdem.annual, by=c("id", "yr")) %>% 
  left_join(drug.annual, by=c("id", "yr")) %>% 
  left_join(base, by="id") %>% 
  left_join(cancer, by="id") %>% 
  filter(yr>2006 & yr<2020)

# Merge imputed data
impute <- bmi %>% 
  left_join(diab, by=c("id", "yr")) %>% 
  left_join(hbp, by=c("id", "yr")) %>% 
  left_join(liver, by=c("id", "yr")) %>% 
  left_join(lung, by=c("id", "yr")) %>% 
  left_join(renal, by=c("id", "yr")) %>% 
  left_join(dep, by=c("id", "yr"))

# Merge raw to imputed
alive.imputed <- annual %>% 
  left_join(impute, by=c("id", "yr"))


# Create indicators of disease --------------------------------------------

alive.imputed2 <- alive.imputed %>% 
  # Diabetes: HbA1c > 6.5% or medication use
  mutate(hi_hba1c = as.numeric((hba1c>6.5)),
         diab_idx = as.numeric((hba1c>6.5) | (diabtx6m==1)),
         
         # Obstructive lung disease: FEV/FVC <=0.70
         poor_lung_func = as.numeric((fev_fvc<=0.70)),
         
         # Kidney dysfunction: urine protein-creatinine ratio >200 or GFR <60
         renal_dysfunc = as.numeric((uprt_crt>200) | (gfr<60)),
         
         # Liver fibrosis: liver stiffness >=9.3
         liver_fibrosis = as.numeric(fbscan>=9.3),
         
         # Liver cirrhosis: liver stiffness >=12.3
         cirrhosis = as.numeric(fbscan>=12.3),
         
         # Obesity: BMI >=30
         obesity = as.numeric(bmimedh>=30),
         
         # Underweight: BMI <18.5
         underweight = as.numeric(bmimedh<18.5),
         
         # Hypertension: Dias BP >=90, Sys BP >=140, or treatment in last 6 months
         hi_bp = as.numeric((bpdias>=90) | (bpsys>=140)),
         hyper_idx = as.numeric((bpdias>=90) | (bpsys>=140) | (hbptx6m==1)),
         
         # Cancer: any diagnosis in yr or before
         cancer = ifelse(is.na(dx_date), 0, as.numeric(year(dx_date)<=yr)),
         
         # Depression: CESD>=23 or self-reported treatment in last 6 months
         depression = as.numeric(cesd23==1 | deptx6m==1))


# Define multimorbidity ---------------------------------------------------
# Permanent conditions: diabetes, hypertension, OLD with 2 consecutive measures of poor lung function, 
#                       renal disease if 2 consecutive measures of kidney dysfunction, cancer
# Otherwise, condition allowed to toggle on/off

multimorbid.imputed <- alive.imputed2 %>%
  group_by(id) %>% 
  mutate(
    # Define permanent conditions
    # Diabetes
    diabetes = as.numeric(cumsum(diab_idx)>=1),

    # Lung disease
    becomes_disease = ifelse(is.na(lag(poor_lung_func)), 0,
                             as.numeric(poor_lung_func==1 & lag(poor_lung_func)==1)),
    lung_disease = as.numeric(cumsum(becomes_disease)>=1),

    # Renal disease
    becomes_disease = ifelse(is.na(lag(renal_dysfunc)), 0,
                             as.numeric(renal_dysfunc==1 & lag(renal_dysfunc)==1)),
    renal_disease = as.numeric(cumsum(becomes_disease)>=1),

    # Liver disease
    liver_cirrhosis = as.numeric(cumsum(cirrhosis)>=1),

    # Hypertension
    hypertension = as.numeric(cumsum(hyper_idx)>=1),

    # Multimorbidity
    n_cond = hypertension +
             diabetes +
             as.numeric(poor_lung_func==1 | lung_disease==1) + 
             as.numeric(renal_dysfunc==1 | renal_disease==1) +
             as.numeric(liver_fibrosis==1) +
             as.numeric(depression==1) +
             as.numeric(cancer==1) +
             as.numeric(obesity==1 | underweight==1)) %>% 
  select(-becomes_disease)


# Output data -------------------------------------------------------------

write_csv(multimorbid.imputed, file="../data/alive_multimorbid_imputed.csv")



