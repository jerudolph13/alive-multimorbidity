
###################################################################################################################
#
# Project: Prevalence of multimorbidity
#
# Purpose: Run random forest missing data model for diabetes
#
# Author: Jacqueline Rudolph
#
# Last Update: 28 Sep 2022
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
  select(id, visdate, diabdx6m, diabtx6m, evrdiab, hba1c)

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

# Make data set wide
wide <- dat %>% 
  pivot_wider(id_cols=id, names_from=n_vis, values_from=-c(id, n, n_vis, visdate))

# Merge baseline data
wide2 <- right_join(base, wide, by="id")

# Prepare data for missForest  
wide3 <- wide2 %>% 
  mutate(enroll_period=factor(enroll_period), 
         across(c(beduchs, m0f1, black,
                  vars_select(names(wide), starts_with("diabdx6m")), 
                  vars_select(names(wide), starts_with("diabtx6m")),
                  vars_select(names(wide), starts_with("evrdiab")),
                  vars_select(names(wide), starts_with("anyhomels")), 
                  vars_select(names(wide), starts_with("inclt5k")), 
                  vars_select(names(wide), starts_with("work")), 
                  vars_select(names(wide), starts_with("curuser")), 
                  vars_select(names(wide), starts_with("anntivwomj")), 
                  vars_select(names(wide), starts_with("alcheavy")), 
                  vars_select(names(wide), starts_with("alcmoderate")), 
                  vars_select(names(wide), starts_with("cigyn")),
                  vars_select(names(wide), starts_with("hiv")),
                  vars_select(names(wide), starts_with("hcvvis"))), ~ factor(.x, levels=c(0, 1))),
         across(vars_select(names(wide), starts_with("yr")), ~ factor(.x, levels=2007:2019)))

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
  pivot_longer(cols=starts_with("hba1c"), names_to="n_vis", names_prefix="hba1c_", values_to="hba1c")
trt <- wide.impute2 %>% 
  select(id, vars_select(names(wide), starts_with("diabtx6m"))) %>% 
  pivot_longer(cols=starts_with("diabtx6m"), names_to="n_vis", names_prefix="diabtx6m_", values_to="diabtx6m")

# Keep imputed values from observed visits
diab <- dat %>% 
  select(id, visdate, n_vis) %>% 
  left_join(hba1c, by=c("id", "n_vis")) %>% 
  left_join(trt, by=c("id", "n_vis"))


# Export data -------------------------------------------------------------

write_csv(diab, file="../data/alive_diab_impute.csv")
