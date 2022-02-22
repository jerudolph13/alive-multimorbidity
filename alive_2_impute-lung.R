
###################################################################################################################
#
# Purpose: Impute missingness in the ALIVE data using random forest
#
# Author: Jacqueline Rudolph
#
# Last Update: 08 Dec 2021
#
###################################################################################################################

# Goal: Use random forest to impute missing values of laboratory measurements
    # https://stat.ethz.ch/education/semesters/ss2012/ams/paper/missForest_1.2.pdf
# MI would be preference in associational/causal analyses but doesn't have theoretical support in 
# cluster analyses
# Compare results using different imputation approaches?

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

hiv <- read_csv("../data/alive_hiv.csv") %>%
  mutate(yr = year(visdate)) %>% 
  select(-c(visdate, enroll_period, visit))


# Make data annual --------------------------------------------------------

comorbid_annual <- make_annual(comorbid, sym("fev_fvc"), sym("lungdx6m"), sym("lungtx6m"), sym("evrlungdx")) %>% 
  mutate(across(c(lungdx6m, lungtx6m, evrlungdx), ~ ceiling(.x)))

# HIV data
  # hivvl dropped from model because of poor performance
hiv_annual <- make_annual(hiv, sym("hiv"), sym("hrt6m"), sym("hcvvis")) %>% 
  mutate(across(c(hiv, hrt6m, hcvvis), ~ ceiling(.x)))

# Sociodemographic data
socdem_annual <- make_annual(socdem, sym("age"), sym("m0f1"), sym("black"), sym("anyhomels"), sym("inclt5k"),
                             sym("work")) %>% 
  mutate(across(c(black, m0f1, anyhomels, inclt5k, work), ~ ceiling(.x)),
         age = round(age))

# Merge data
annual <- left_join(left_join(comorbid_annual, hiv_annual, by=c("id", "yr")), socdem_annual, by=c("id", "yr")) %>% 
  filter(yr>2006 & yr<2020)


# Fill in missed visits ---------------------------------------------------

allvisits <- annual %>% expand(id, yr)
long <- left_join(allvisits, annual, by=c("id", "yr"))
  # rename(vl = hivvl)

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
  
wide2 <- wide2 %>% 
  mutate(enroll_period=factor(enroll_period), 
         across(c(beduchs, bempcurr, beverwk, bevinj, bhome1yr, bnevmar,
                  vars_select(names(wide), starts_with("m0f1")), vars_select(names(wide), starts_with("black")), 
                  vars_select(names(wide), starts_with("anyhomels")), vars_select(names(wide), starts_with("inclt5k")), 
                  vars_select(names(wide), starts_with("work")), vars_select(names(wide), starts_with("lungdx6m")), 
                  vars_select(names(wide), starts_with("lungtx6m")), vars_select(names(wide), starts_with("evrlungdx")),
                  vars_select(names(wide), starts_with("hiv")), vars_select(names(wide), starts_with("hrt6m")), 
                  vars_select(names(wide), starts_with("hcvvis"))), ~ factor(.x, levels=c(0, 1))))

wide2 <- as.data.frame(wide2)
  

# Imputation --------------------------------------------------------------

set.seed(123)
registerDoParallel(cores=3)
impute <- missForest(select(wide2, -id), verbose=T, ntree=20, parallelize="forests")
wide.impute <- impute$ximp
error <- impute$OOBerror
#error <- data.frame(var = names(select(wide2, -id)), error = impute$OOBerror)

wide.impute2 <- bind_cols(id=wide2$id, wide.impute)


# Make data long ----------------------------------------------------------

lung <- wide.impute2 %>% 
  select(id, vars_select(names(wide), starts_with("fev"))) %>% 
  pivot_longer(cols=starts_with("fev"), names_to="yr", names_prefix="fev_fvc_", values_to="fev_fvc")


# Export data -------------------------------------------------------------

write_csv(lung, file="../data/alive_lung_impute.csv")
