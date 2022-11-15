
###################################################################################################################
#
# Project: Prevalence of multimorbidity
#
# Purpose: Run random forest missing data model for renal function
#
# Author: Jacqueline Rudolph
#
# Last Update: 29 Sep 2022
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


# Read in data ------------------------------------------------------------

base <- read_csv("../data/alive_base.csv") %>% 
  select(id, beduchs, enroll_period)

socdem <- read_csv("../data/alive_socdem.csv") %>% 
  filter(year(visdate) %in% 2007:2019) %>% 
  select(-c(enroll_period, visit, race))

comorbid <- read_csv("../data/alive_comorbid.csv") %>% 
  filter(year(visdate) %in% 2007:2019) %>% 
  select(id, visdate, gfr, uprt_crt, renaldx6m, renaltx6m, evrrenal)

drug <- read_csv("../data/alive_drug.csv") %>% 
  filter(year(visdate) %in% 2007:2019) %>% 
  select(id, visdate, alcheavy, alcmoderate, curuser, anntivwomj, cigyn)

hiv <- read_csv("../data/alive_hiv.csv") %>% 
  filter(year(visdate) %in% 2007:2019) %>% 
  select(id, visdate, hiv, hcvvis)


# Manage data -------------------------------------------------------------

# Create baseline variables
socdem.base <- socdem %>% 
  select(id, m0f1, black) %>% 
  filter(!duplicated(id))

base <- left_join(base, socdem.base, by="id")

# Merge longitudinal data
dat <- comorbid %>%
  left_join(select(socdem, -c(m0f1, black)), by=c("id", "visdate")) %>% 
  left_join(drug, by=c("id", "visdate")) %>% 
  left_join(hiv, by=c("id", "visdate")) %>% 
  group_by(id) %>% 
  mutate(n = 1,
         n_vis = as.character(cumsum(n)),
         yr = year(visdate)) %>% 
  ungroup()

# Split long data set for GFR (2007-2010 & 2011-2019)
long0710 <- filter(dat, yr<2011)
long1119 <- filter(dat, yr>=2011)
  

# Make data sets wide -----------------------------------------------------

wide <- dat %>% 
  pivot_wider(id_cols=id, names_from=n_vis, values_from=-c(id, n, n_vis, visdate))
wide0710 <- long0710 %>% 
  pivot_wider(id_cols=id, names_from=n_vis, values_from=-c(id, n, n_vis, visdate))
wide1119 <- long1119 %>% 
  pivot_wider(id_cols=id, names_from=n_vis, values_from=-c(id, n, n_vis, visdate))

wide2 <- right_join(base, wide, by="id")
wide0710 <- right_join(base, wide0710, by="id")
wide1119 <- right_join(base, wide1119, by="id")


# Urine protein/creatinine ------------------------------------------------

wide3 <- wide2 %>% 
  mutate(enroll_period = factor(enroll_period, levels=c("1988-1992", 
                                                        "1994-1997", 
                                                        "1998-2001", 
                                                        "2005-2015", 
                                                        "2016-2018")), 
         across(c(beduchs, m0f1, black,
                  vars_select(names(wide2), starts_with("renaldx6m")), 
                  vars_select(names(wide2), starts_with("renaltx6m")), 
                  vars_select(names(wide2), starts_with("evrrenal")),
                  vars_select(names(wide2), starts_with("anyhomels")), 
                  vars_select(names(wide2), starts_with("inclt5k")), 
                  vars_select(names(wide2), starts_with("work")), 
                  vars_select(names(wide2), starts_with("anntivwomj")), 
                  vars_select(names(wide2), starts_with("alcheavy")), 
                  vars_select(names(wide2), starts_with("alcmoderate")), 
                  vars_select(names(wide2), starts_with("cigyn")),
                  vars_select(names(wide2), starts_with("hiv")),
                  vars_select(names(wide2), starts_with("hcvvis"))), ~ factor(.x, levels=c(0, 1))),
         across(vars_select(names(wide2), starts_with("yr")), ~ factor(.x, levels=2007:2019)))

wide_uprt <- wide3 %>% 
  select(-c(vars_select(names(wide3), starts_with("gfr"))))
wide_uprt <- as.data.frame(wide_uprt)

# Impute uprt_crt
set.seed(123)
registerDoParallel(cores=8)
impute <- missForest(select(wide_uprt, -id), verbose=T, ntree=20, parallelize="forests")
wide.impute_uprt <- impute$ximp
error_uprt <- impute$OOBerror

wide.impute_uprt2 <- bind_cols(id=wide3$id, wide.impute_uprt)


# GFR, 2007-2010 ----------------------------------------------------------

wide0710 <- wide0710 %>% 
  mutate(enroll_period = factor(enroll_period, levels=c("1988-1992", 
                                                        "1994-1997", 
                                                        "1998-2001", 
                                                        "2005-2015", 
                                                        "2016-2018")), 
         across(c(beduchs, m0f1, black,
                  vars_select(names(wide0710), starts_with("renaldx6m")), 
                  vars_select(names(wide0710), starts_with("renaltx6m")), 
                  vars_select(names(wide0710), starts_with("evrrenal")),
                  vars_select(names(wide0710), starts_with("anyhomels")), 
                  vars_select(names(wide0710), starts_with("inclt5k")), 
                  vars_select(names(wide0710), starts_with("work")), 
                  vars_select(names(wide0710), starts_with("anntivwomj")), 
                  vars_select(names(wide0710), starts_with("alcheavy")), 
                  vars_select(names(wide0710), starts_with("alcmoderate")), 
                  vars_select(names(wide0710), starts_with("cigyn")),
                  vars_select(names(wide0710), starts_with("hiv")),
                  vars_select(names(wide0710), starts_with("hcvvis"))), ~ factor(.x, levels=c(0, 1))),
         across(vars_select(names(wide0710), starts_with("yr")), ~ factor(.x, levels=2007:2019)))
wide_gfr0710 <- wide0710 %>% 
  select(-c(vars_select(names(wide0710), starts_with("uprt"))))
wide_gfr0710 <- as.data.frame(wide_gfr0710)

# Impute gfr
set.seed(123)
registerDoParallel(cores=8)
impute <- missForest(select(wide_gfr0710, -id), verbose=T, ntree=20, parallelize="forests")
wide.impute_gfr <- impute$ximp
error_gfr1 <- impute$OOBerror

wide.impute_gfr0710 <- bind_cols(id=wide_gfr0710$id, wide.impute_gfr)


# GFR,  2011-2019 ---------------------------------------------------------

wide1119 <- wide1119 %>% 
  mutate(enroll_period = factor(enroll_period, levels=c("1988-1992", 
                                                        "1994-1997", 
                                                        "1998-2001", 
                                                        "2005-2015", 
                                                        "2016-2018")), 
         across(c(beduchs, m0f1, black,
                  vars_select(names(wide1119), starts_with("renaldx6m")), 
                  vars_select(names(wide1119), starts_with("renaltx6m")), 
                  vars_select(names(wide1119), starts_with("evrrenal")),
                  vars_select(names(wide1119), starts_with("anyhomels")), 
                  vars_select(names(wide1119), starts_with("inclt5k")), 
                  vars_select(names(wide1119), starts_with("work")), 
                  vars_select(names(wide1119), starts_with("anntivwomj")), 
                  vars_select(names(wide1119), starts_with("alcheavy")), 
                  vars_select(names(wide1119), starts_with("alcmoderate")), 
                  vars_select(names(wide1119), starts_with("cigyn")),
                  vars_select(names(wide1119), starts_with("hiv")),
                  vars_select(names(wide1119), starts_with("hcvvis"))), ~ factor(.x, levels=c(0, 1))),
         across(vars_select(names(wide1119), starts_with("yr")), ~ factor(.x, levels=2007:2019)))

wide_gfr1119 <- wide1119 %>% 
  select(-c(vars_select(names(wide1119), starts_with("uprt"))))
wide_gfr1119 <- as.data.frame(wide_gfr1119)

# Impute gfr 2011-2019
set.seed(123)
registerDoParallel(cores=8)
impute <- missForest(select(wide_gfr1119, -id), verbose=T, ntree=20, parallelize="forests")
wide.impute_gfr <- impute$ximp
error_gfr2 <- impute$OOBerror

wide.impute_gfr1119 <- bind_cols(id=wide_gfr1119$id, wide.impute_gfr)


# Make data long ----------------------------------------------------------

uprt <- wide.impute_uprt2 %>% 
  select(id, vars_select(names(wide3), starts_with("uprt"))) %>% 
  pivot_longer(cols=starts_with("uprt"), names_to="n_vis", names_prefix="uprt_crt_", values_to="uprt_crt") 
trt <- wide.impute_uprt2 %>% 
  select(id, vars_select(names(wide3), starts_with("renaltx6m"))) %>% 
  pivot_longer(cols=starts_with("renaltx6m"), names_to="n_vis", names_prefix="renaltx6m_", values_to="renaltx6m")

gfr1 <- wide.impute_gfr0710 %>% 
  select(id, vars_select(names(wide_gfr0710), starts_with("gfr"))) %>% 
  pivot_longer(cols=starts_with("gfr"), names_to="n_vis", names_prefix="gfr_", values_to="gfr")
gfr2 <- wide.impute_gfr1119 %>% 
  select(id, vars_select(names(wide_gfr1119), starts_with("gfr"))) %>% 
  pivot_longer(cols=starts_with("gfr"), names_to="n_vis", names_prefix="gfr_", values_to="gfr")

renal1 <- long0710 %>% 
  select(id, visdate, n_vis) %>% 
  left_join(gfr1, by=c("id", "n_vis")) %>% 
  # To match observed data
  mutate(gfr = ifelse(gfr>95, 95, gfr))
renal2 <- long1119 %>% 
  select(id, visdate, n_vis) %>% 
  left_join(gfr2, by=c("id", "n_vis"))
renal <- bind_rows(renal1, renal2) %>% 
  arrange(id, visdate) %>% 
  left_join(uprt, by=c("id", "n_vis")) %>% 
  left_join(trt, by=c("id", "n_vis"))


# Export data -------------------------------------------------------------

write_csv(renal, file="../data/alive_renal_impute.csv")
