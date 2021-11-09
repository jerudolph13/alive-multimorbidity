
###################################################################################################################
#
# Purpose: Manipulate data and deal with missingness
#
# Author: Jacqueline Rudolph
#
# Last Update: 21 Oct 2021
#
###################################################################################################################

packages <- c("tidyverse", "data.table", "zoo")
for (package in packages) {
  library(package, character.only=T)
}

# A function to get annual values
make_annual <- function(dat, ...) {
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

dat <- read_csv(file="../data/alive_comorbid.csv") %>% 
  mutate(yr = as.numeric(format(visdate, format="%Y")))

enroll <- dat %>% 
  select(id, enroll_period) %>% 
  filter(as.numeric(!duplicated(id, fromLast=F))==1)

hiv <- read_csv(file="../data/alive_hiv.csv") %>% 
  mutate(yr = as.numeric(format(visdate, format="%Y")))

socdem <- read_csv(file="../data/alive_socdem.csv") %>% 
  mutate(yr = as.numeric(format(visdate, format="%Y")))


# Summarize by year -------------------------------------------------------

# Comorbidity data
comorbid_annual <- make_annual(dat, sym("hba1c"), sym("fev_fvc"), sym("gfr"), sym("uprt_crt"), sym("fbscan"), 
                                    sym("bmimedh"), sym("bpdias"), sym("bpsys"), sym("hemo"), 
                                    sym("hbptx6m"), sym("diabtx6m")) %>% 
  mutate(hbptx6m = ceiling(hbptx6m),
         diabtx6m = ceiling(diabtx6m))

# HIV data
hiv_annual <- make_annual(hiv, sym("hiv")) %>% 
  mutate(hiv = ceiling(hiv))

# Sociodemographic data
socdem_annual <- make_annual(socdem, sym("age"), sym("m0f1"), sym("black")) %>% 
  mutate(age_cat = cut(age, breaks=c(0, 39, 49, 59, 69, 120), 
                       labels=c("<40", "40-49", "50-59", "60-69", "70+")),
         over60 = as.numeric(age>=60),
         black = ceiling(black),
         m0f1 = ceiling(m0f1))

# Merge data
comorbid_annual <- left_join(left_join(comorbid_annual, hiv_annual), socdem_annual, by=c("id", "yr"))
comorbid_annual <- left_join(comorbid_annual, enroll, by=c("id"))


# Create indicators of chronic disease ------------------------------------

comorbid_annual2 <- comorbid_annual %>% 
  # Diabetes: HbA1c > 6.5% or medication use
  mutate(diab_lab = as.numeric((hba1c>6.5) | (diabtx6m==1)),
         
  # Obstructive lung disease: FEV/FVC <=0.70
         lung_lab = as.numeric(fev_fvc<=0.70),
  
  # Kidney dysfunction: urine protein-creatinine ratio >200 or GFR <60
         renal_lab = as.numeric((uprt_crt>200) | (gfr<60)),
  
  # Liver fibrosis: liver stiffness >=9.3
         liver_lab = as.numeric(fbscan>=9.3),
  
  # Obesity: BMI >=30
         obese_lab = as.numeric(bmimedh>=30),
  # Hypertension: Dias BP >=90, Sys BP >=140, or treatment in last 6 onths
         hypertension_lab = as.numeric((bpdias>=90) | (bpsys>=140) | (hbptx6m==1)),
  
  # Anemia: hematocrit cut-points based on race, sex, and age
         anemia_lab = as.numeric((black==0 & m0f1==0 & over60==0 & hemo<13.7) | # non-Black men under 60
                                 (black==0 & m0f1==0 & over60==1 & hemo<13.2) | # non-Black men over 60
                                 (black==0 & m0f1==1 & hemo<12.2) |             # non-Black women
                                 (black==1 & m0f1==0 & over60==0 & hemo<12.9) | # Black men under 60
                                 (black==1 & m0f1==0 & over60==1 & hemo<12.7) | # Black men over 60 
                                 (black==1 & m0f1==1 & hemo<11.5))) %>%        # Black women
  select(-c(hba1c, diabtx6m, fev_fvc, uprt_crt, gfr, fbscan, bmimedh, bpdias, bpsys,
            hbptx6m, hemo))


# Examine data ------------------------------------------------------------

summary(select(comorbid_annual2, -c(id, yr, enroll_period)))

# Substantial missingness in diabetes, lung, renal, and liver data
# Persists once we remove 2006 and 2020 (some vars not collected in those years)
comorbid0719 <- filter(comorbid_annual2, yr %in% seq(2007, 2019, 1))
summary(select(comorbid0719, -c(id, yr, enroll_period)))

# Does amount of missingness vary over time and by HIV status?
summ_miss <- comorbid0719 %>% 
#  group_by(yr) %>% 
  summarize(miss_diab = mean(is.na(diab_lab)),
            miss_lung = mean(is.na(lung_lab)),
            miss_renal = mean(is.na(renal_lab)),
            miss_liver = mean(is.na(liver_lab)),
            miss_hyper = mean(is.na(hypertension_lab)),
            miss_obese = mean(is.na(obese_lab)),
            miss_hiv = mean(is.na(hiv)),
            n_tot = n())

# Visualize proportion missing
ggplot(data=summ_miss) +
  geom_bar(aes(x=yr, y=miss_obese), stat="identity") +
  facet_wrap(facets=vars(hiv))


# Deal with missingness ---------------------------------------------------

# Start with last object carried forward
comorbid0719_locf <- comorbid0719 %>% 
  group_by(id) %>% 
  mutate_all(~ na.locf(., na.rm=FALSE))

# Now do next observation carried backward
# Don't do this for anemia because it was structurally missing before 2014
comorbid0719_nocb <- comorbid0719_locf %>% 
  select(-anemia_lab) %>% 
  group_by(id) %>% 
  mutate_all(~ na.locf(., na.rm=FALSE, fromLast=TRUE))

# Now do next observation carried backward for anemia
anemia_nocb <- comorbid0719_locf %>% 
  select(id, yr, anemia_lab) %>%
  filter(yr>2013) %>% 
  group_by(id) %>% 
  mutate_all(~ na.locf(., na.rm=FALSE, fromLast=TRUE))

# Merge back together
comorbid0719_nocb <- left_join(comorbid0719_nocb, anemia_nocb, by=c("id", "yr"))

# How much missingness now?
summary(select(comorbid0719_nocb, -c(id, yr)))
  # Any missingness reflects a case of "always missing"
    # Most missingness in lung and liver disease data (~7%). 
    # Moderate missingness in diabetes and renal disease data (~2%).
    # Negligable missingess in hypertension and obesity data (~0.2%)

summ_miss <- comorbid0719_nocb %>% 
  ungroup() %>% 
  summarize(miss_diab = mean(is.na(diab_lab)),
            miss_lung = mean(is.na(lung_lab)),
            miss_renal = mean(is.na(renal_lab)),
            miss_liver = mean(is.na(liver_lab)),
            miss_hyper = mean(is.na(hypertension_lab)),
            miss_obese = mean(is.na(obese_lab)),
            miss_hiv = mean(is.na(hiv)),
            n_tot = n())
    

# Define multimorbidity ---------------------------------------------------
# Cumulative multimorbidity: once participant has disease, they keep it
# Year-specific multimorbidity: only consider disease value from that year

# 2007-2019: 6 disease definition (no anemia)
multimorbid6 <- comorbid0719_nocb %>%
  group_by(id) %>% 
  mutate(# Diabetes
    diab_lab = ceiling(diab_lab),
    cum_diab = as.numeric(cumsum(diab_lab)>=1),
    # Lung disease
    lung_lab = ceiling(lung_lab),
    cum_lung = as.numeric(cumsum(lung_lab)>=1),
    # Renal disease
    renal_lab = ceiling(renal_lab),
    cum_renal = as.numeric(cumsum(renal_lab)>=1),
    # Liver disease
    liver_lab = ceiling(liver_lab),
    cum_liver = as.numeric(cumsum(liver_lab)>=1),
    # Hypertension
    hypertension_lab = ceiling(hypertension_lab),
    cum_hyper = as.numeric(cumsum(hypertension_lab)>=1),
    # Obesity
    obese_lab = ceiling(obese_lab),
    cum_obese = as.numeric(cumsum(obese_lab)>=1),
    # Multimorbidity
    cum_multi = ifelse(is.na(cum_diab), 0, cum_diab) + ifelse(is.na(cum_lung), 0, cum_lung) +
                ifelse(is.na(cum_renal), 0, cum_renal) + ifelse(is.na(cum_liver), 0, cum_liver) +
                ifelse(is.na(cum_hyper), 0, cum_hyper) + ifelse(is.na(cum_obese), 0, cum_obese),
    yr_multi = ifelse(is.na(diab_lab), 0, diab_lab) + ifelse(is.na(lung_lab), 0, lung_lab) +
               ifelse(is.na(renal_lab), 0, renal_lab) + ifelse(is.na(liver_lab), 0, liver_lab) +
               ifelse(is.na(hypertension_lab), 0, hypertension_lab) + ifelse(is.na(obese_lab), 0, obese_lab))

# 2014-2019: 7 disease definition
multimorbid7 <- comorbid0719_nocb %>%
  filter(yr>2013) %>% 
  group_by(id) %>% 
  mutate(# Diabetes
    diab_lab = ceiling(diab_lab),
    cum_diab = as.numeric(cumsum(diab_lab)>=1),
    # Lung disease
    lung_lab = ceiling(lung_lab),
    cum_lung = as.numeric(cumsum(lung_lab)>=1),
    # Renal disease
    renal_lab = ceiling(renal_lab),
    cum_renal = as.numeric(cumsum(renal_lab)>=1),
    # Liver disease
    liver_lab = ceiling(liver_lab),
    cum_liver = as.numeric(cumsum(liver_lab)>=1),
    # Hypertension
    hypertension_lab = ceiling(hypertension_lab),
    cum_hyper = as.numeric(cumsum(hypertension_lab)>=1),
    # Obesity
    obese_lab = ceiling(obese_lab),
    cum_obese = as.numeric(cumsum(obese_lab)>=1),
    # Anemia
    anemia_lab = ceiling(anemia_lab),
    cum_anemia = as.numeric(cumsum(anemia_lab)>=1),
    # Multimorbidity
    cum_multi = ifelse(is.na(cum_diab), 0, cum_diab) + ifelse(is.na(cum_lung), 0, cum_lung) +
                ifelse(is.na(cum_renal), 0, cum_renal) + ifelse(is.na(cum_liver), 0, cum_liver) +
                ifelse(is.na(cum_hyper), 0, cum_hyper) + ifelse(is.na(cum_obese), 0, cum_obese) +
                ifelse(is.na(cum_anemia), 0, cum_anemia),
    yr_multi = ifelse(is.na(diab_lab), 0, diab_lab) + ifelse(is.na(lung_lab), 0, lung_lab) +
               ifelse(is.na(renal_lab), 0, renal_lab) + ifelse(is.na(liver_lab), 0, liver_lab) +
               ifelse(is.na(hypertension_lab), 0, hypertension_lab) + ifelse(is.na(obese_lab), 0, obese_lab) +
               ifelse(is.na(anemia_lab), 0, anemia_lab)) 


# Output data -------------------------------------------------------------

write_csv(multimorbid6, file="../data/alive_multimorbid6.csv")
write_csv(multimorbid7, file="../data/alive_multimorbid7.csv")


