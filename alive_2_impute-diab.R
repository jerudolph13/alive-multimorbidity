
###################################################################################################################
#
# Project: Impute missing data in ALIVE
#
# Purpose: Run random forest missing data model for diabetes
#
# Author: Jacqueline Rudolph
#
# Last Update: 12 Apr 2022
#
###################################################################################################################

# Goal: Use random forest to impute missing values of laboratory measurements
    # https://stat.ethz.ch/education/semesters/ss2012/ams/paper/missForest_1.2.pdf
# MI would be preference in associational/causal analyses but doesn't have theoretical support in 
# cluster analyses

packages <- c("tidyverse", "data.table", "lubridate", "tidyselect", "missForest", "doParallel")
for (package in packages) {
  library(package, character.only=T)
}

# Function to make data annual
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

base <- read_csv("../data/alive_base.csv") %>% 
  select(-c(visdate, visit, bincome, bevcig, bevalc, bcurrpa, beduc))

socdem <- read_csv("../data/alive_socdem.csv") %>% 
  mutate(yr = year(visdate)) %>% 
  select(-c(visdate, enroll_period, visit, race))

comorbid <- read_csv("../data/alive_comorbid.csv") %>% 
  mutate(yr = year(visdate)) %>% 
  select(-c(visdate, enroll_period, visit, hemo, ucreat, uprotn, fbscgrp, fbscge12_3))

drug <- read_csv("../data/alive_drug.csv") %>% 
  mutate(yr = year(visdate)) %>% 
  select(-c(visdate, enroll_period, visit))


# Make data annual --------------------------------------------------------

# Condition data
comorbid_annual <- make_annual(comorbid, sym("hba1c"), sym("diabtx6m"), sym("diabdx6m"), sym("evrdiab")) %>% 
  mutate(across(c(diabdx6m, diabtx6m, evrdiab), ~ ceiling(.x)))

# Sociodemographic data
socdem_annual <- make_annual(socdem, sym("age"), sym("m0f1"), sym("black"), sym("anyhomels"), sym("inclt5k"),
                             sym("work")) %>% 
  mutate(across(c(black, m0f1, anyhomels, inclt5k, work), ~ ceiling(.x)),
         age = round(age))

# Substance use data
drug_annual <- make_annual(drug, sym("curuser"), sym("alcheavy"), sym("cigyn")) %>% 
  mutate(across(c(curuser, alcheavy, cigyn), ~ ceiling(.x)))

# Merge data
annual <- left_join(left_join(comorbid_annual, drug_annual, by=c("id", "yr")), socdem_annual, by=c("id", "yr")) %>% 
  filter(yr>2006 & yr<2020) %>% 
  ungroup()


# Percent missing ---------------------------------------------------------

summary(filter(comorbid_annual, yr>2006 & yr<2020))


# Fill in missed visits ---------------------------------------------------

allvisits <- annual %>% expand(id, yr)
long <- left_join(allvisits, annual, by=c("id", "yr"))

# Fill in age for missed years
age_imp <- long %>% 
  filter(!is.na(age)) %>% 
  mutate(conv = yr-age,
         first = as.numeric(!duplicated(id))) %>% 
  filter(first==1) %>% 
  select(id, conv)

long2 <- left_join(long, age_imp, by="id") %>% 
  mutate(age=yr-conv) %>% 
  select(-conv)


# Make data sets wide -----------------------------------------------------

wide <- long2 %>% 
  pivot_wider(id_cols=id, names_from=yr, values_from=-c(yr, id))

wide2 <- right_join(base, wide, by="id")
  
wide3 <- wide2 %>% 
  mutate(enroll_period = factor(enroll_period, levels=c("1988-1992", 
                                                        "1994-1997", 
                                                        "1998-2001", 
                                                        "2005-2015", 
                                                        "2016-2018")), 
         across(c(beduchs, bempcurr, beverwk, bevinj, bhome1yr, bnevmar,
                  vars_select(names(wide), starts_with("diabdx6m")), 
                  vars_select(names(wide), starts_with("diabtx6m")),
                  vars_select(names(wide), starts_with("evrdiab")),
                  vars_select(names(wide), starts_with("m0f1")), 
                  vars_select(names(wide), starts_with("black")), 
                  vars_select(names(wide), starts_with("anyhomels")), 
                  vars_select(names(wide), starts_with("inclt5k")), 
                  vars_select(names(wide), starts_with("work")), 
                  vars_select(names(wide), starts_with("curuser")), 
                  vars_select(names(wide), starts_with("alcheavy")), 
                  vars_select(names(wide), starts_with("cigyn"))), ~ factor(.x, levels=c(0, 1))))

wide3 <- as.data.frame(wide3)
  

# Imputation --------------------------------------------------------------

set.seed(123)
registerDoParallel(cores=8)
impute <- missForest(select(wide3, -id), verbose=T, ntree=20, parallelize="forests")
wide.impute <- impute$ximp
error <- impute$OOBerror
#error <- data.frame(var = names(select(wide2, -id)), error = impute$OOBerror)

wide.impute2 <- bind_cols(id=wide3$id, wide.impute)


# Make data long ----------------------------------------------------------

hba1c <- wide.impute2 %>% 
  select(id, vars_select(names(wide), starts_with("hba1c"))) %>% 
  pivot_longer(cols=starts_with("hba1c"), names_to="yr", names_prefix="hba1c_", values_to="hba1c")
trt <- wide.impute2 %>% 
  select(id, vars_select(names(wide), starts_with("diabtx6m"))) %>% 
  pivot_longer(cols=starts_with("diabtx6m"), names_to="yr", names_prefix="diabtx6m_", values_to="diabtx6m")

diab <- left_join(hba1c, trt, by=c("id", "yr"))


# Export data -------------------------------------------------------------

write_csv(diab, file="../data/alive_diab_impute.csv")
