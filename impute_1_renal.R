
###################################################################################################################
#
# Purpose: Impute missingness in the ALIVE data using random forest
#
# Author: Jacqueline Rudolph
#
# Last Update: 19 Dec 2021
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

comorbid_annual <- make_annual(comorbid, sym("gfr"), sym("uprt_crt"), sym("renaldx6m"), sym("renaltx6m"), 
                               sym("evrrenal")) %>% 
  mutate(across(c(renaldx6m, renaltx6m, evrrenal), ~ ceiling(.x)))

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
  filter(yr>2006 & yr<2020) %>% 
  # See if log transforming uprt_crt leads to improved prediction accuracy
  mutate(uprt_crt = log(uprt_crt),
         gfr = log(gfr))


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

# Split long data set for GFR (2007-2010 & 2011-2019)
long0710 <- filter(long2, yr<2011)
  # Identify people with GFR<95: We will model missing records for these individuals
  gfrlt95 <- unique(long0710$id[long0710$gfr<95 & !is.na(long0710$gfr)]) 
  long0710_lt95 <- long0710[long0710$id %in% gfrlt95, ]
  # Identify people with GFR=95: We will set missing values to 95
  long0710_gt95 <- long0710[!(long0710$id %in% gfrlt95), ] %>% 
    mutate(gfr = 95)

long1119 <- filter(long2, yr>=2011)
  

# Make data sets wide -----------------------------------------------------

wide <- long2 %>% 
  pivot_wider(id_cols=id, names_from=yr, values_from=-c(yr, id))
wide0710_lt95 <- long0710_lt95 %>% 
  pivot_wider(id_cols=id, names_from=yr, values_from=-c(yr, id))
wide1119 <- long1119 %>% 
  pivot_wider(id_cols=id, names_from=yr, values_from=-c(yr, id))

wide2 <- right_join(base, wide, by="id")
wide0710_lt95 <- right_join(base, wide0710_lt95, by="id")
wide1119 <- right_join(base, wide1119, by="id")
  
wide2 <- wide2 %>% 
  mutate(enroll_period=factor(enroll_period), 
         across(c(beduchs, bempcurr, beverwk, bevinj, bhome1yr, bnevmar,
                  vars_select(names(wide), starts_with("m0f1")), vars_select(names(wide), starts_with("black")), 
                  vars_select(names(wide), starts_with("anyhomels")), vars_select(names(wide), starts_with("inclt5k")), 
                  vars_select(names(wide), starts_with("work")), vars_select(names(wide), starts_with("renaldx6m")), 
                  vars_select(names(wide), starts_with("renaltx6m")), vars_select(names(wide), starts_with("evrrenal")),
                  vars_select(names(wide), starts_with("hiv")), vars_select(names(wide), starts_with("hrt6m")), 
                  vars_select(names(wide), starts_with("hcvvis"))), ~ factor(.x, levels=c(0, 1))))
wide0710_lt95 <- wide0710_lt95 %>% 
  mutate(enroll_period=factor(enroll_period), 
         across(c(beduchs, bempcurr, beverwk, bevinj, bhome1yr, bnevmar,
                  vars_select(names(wide0710_lt95), starts_with("m0f1")), vars_select(names(wide0710_lt95), starts_with("black")), 
                  vars_select(names(wide0710_lt95), starts_with("anyhomels")), vars_select(names(wide0710_lt95), starts_with("inclt5k")), 
                  vars_select(names(wide0710_lt95), starts_with("work")), vars_select(names(wide0710_lt95), starts_with("renaldx6m")), 
                  vars_select(names(wide0710_lt95), starts_with("renaltx6m")), vars_select(names(wide0710_lt95), starts_with("evrrenal")),
                  vars_select(names(wide0710_lt95), starts_with("hiv")), vars_select(names(wide0710_lt95), starts_with("hrt6m")), 
                  vars_select(names(wide0710_lt95), starts_with("hcvvis"))), ~ factor(.x, levels=c(0, 1))))
wide1119 <- wide1119 %>% 
  mutate(enroll_period=factor(enroll_period), 
         across(c(beduchs, bempcurr, beverwk, bevinj, bhome1yr, bnevmar,
                  vars_select(names(wide1119), starts_with("m0f1")), vars_select(names(wide1119), starts_with("black")), 
                  vars_select(names(wide1119), starts_with("anyhomels")), vars_select(names(wide1119), starts_with("inclt5k")), 
                  vars_select(names(wide1119), starts_with("work")), vars_select(names(wide1119), starts_with("renaldx6m")), 
                  vars_select(names(wide1119), starts_with("renaltx6m")), vars_select(names(wide1119), starts_with("evrrenal")),
                  vars_select(names(wide1119), starts_with("hiv")), vars_select(names(wide1119), starts_with("hrt6m")), 
                  vars_select(names(wide1119), starts_with("hcvvis"))), ~ factor(.x, levels=c(0, 1))))

wide_uprt <- wide2 %>% 
  select(-c(vars_select(names(wide), starts_with("gfr"))))
wide_uprt <- as.data.frame(wide_uprt)

wide_gfr0710 <- wide0710_lt95 %>% 
  select(-c(vars_select(names(wide0710_lt95), starts_with("uprt"))))
wide_gfr0710 <- as.data.frame(wide_gfr0710)

wide_gfr1119 <- wide1119 %>% 
  select(-c(vars_select(names(wide1119), starts_with("uprt"))))
wide_gfr1119 <- as.data.frame(wide_gfr1119)

# Imputation --------------------------------------------------------------

# Impute uprt_crt
set.seed(123)
registerDoParallel(cores=3)
impute <- missForest(select(wide_uprt, -id), verbose=T, ntree=20, parallelize="forests")
wide.impute_uprt <- impute$ximp
error_uprt <- impute$OOBerror
#error <- data.frame(var = names(select(wide2, -id)), error = impute$OOBerror)

wide.impute_uprt2 <- bind_cols(id=wide2$id, wide.impute_uprt)

# Impute gfr 2007-2010
set.seed(123)
registerDoParallel(cores=3)
impute <- missForest(select(wide_gfr0710, -id), verbose=T, ntree=20, parallelize="forests")
wide.impute_gfr <- impute$ximp
error_gfr1 <- impute$OOBerror
#error_gfr <- data.frame(var = names(select(wide_gfr, -id)), error = impute$OOBerror)

wide.impute_gfr0710 <- bind_cols(id=wide_gfr0710$id, wide.impute_gfr)

# Impute gfr 2011-2019
set.seed(123)
registerDoParallel(cores=3)
impute <- missForest(select(wide_gfr1119, -id), verbose=T, ntree=20, parallelize="forests")
wide.impute_gfr <- impute$ximp
error_gfr2 <- impute$OOBerror
#error_gfr <- data.frame(var = names(select(wide_gfr, -id)), error = impute$OOBerror)

wide.impute_gfr1119 <- bind_cols(id=wide_gfr1119$id, wide.impute_gfr)


# Make data long ----------------------------------------------------------

uprt <- wide.impute_uprt2 %>% 
  select(id, vars_select(names(wide), starts_with("uprt"))) %>% 
  pivot_longer(cols=starts_with("uprt"), names_to="yr", names_prefix="uprt_crt_", values_to="uprt_crt") %>% 
  mutate(uprt_crt = exp(uprt_crt))

gfr1 <- wide.impute_gfr0710 %>% 
  select(id, vars_select(names(wide_gfr0710), starts_with("gfr"))) %>% 
  pivot_longer(cols=starts_with("gfr"), names_to="yr", names_prefix="gfr_", values_to="gfr")
gfr2 <- wide.impute_gfr1119 %>% 
  select(id, vars_select(names(wide_gfr1119), starts_with("gfr"))) %>% 
  pivot_longer(cols=starts_with("gfr"), names_to="yr", names_prefix="gfr_", values_to="gfr")
gfr3 <- long0710_gt95 %>% 
  select(id, yr, gfr) %>% 
  mutate(yr=as.character(yr))
gfr <- bind_rows(gfr1, gfr2, gfr3) %>% 
  arrange(id, yr) %>% 
  mutate(gfr = exp(gfr),
         gfr = ifelse((yr %in% c("2007", "2008", "2009", "2010")) & gfr>95, 95, gfr))

renal <- left_join(uprt, gfr, by=c("id", "yr"))


# Export data -------------------------------------------------------------

write_csv(renal, file="../data/alive_renal_impute.csv")
