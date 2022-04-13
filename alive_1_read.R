
###################################################################################################################
#
# Purpose: Read in and organize ALIVE data
#
# Author: Jacqueline Rudolph
#
# Last Update: 04 Apr 2022
#
###################################################################################################################

library("tidyverse")
library("lubridate")

load(file="../../../raw_data/alive/jackiereq.rdata")
load(file="../../../raw_data/alive/jackiereq3.rdata")
load(file="../../../raw_data/alive/jackiereq4.rdata")
load(file="../../../raw_data/alive/jackiereq5.rdata")

dat <- tibble(jackiereq2) %>% 
  arrange(id, visdate)

enroll <- tibble(jackiereq3) %>% 
  arrange(id)

dat <- left_join(dat, enroll, by="id") %>% 
  mutate(days_study = (basdate %--% visdate)/ddays(1),
         yrs_study = (basdate %--% visdate)/dyears(1),
         yr_enroll = as.numeric(format(basdate, format="%Y")),
         enroll_period = cut(yr_enroll, breaks=c(1987, 1992, 1997, 2001, 2015, 2018),
                             labels=c("1988-1992","1994-1997","1998-2001","2005-2015","2016-2018")))


# Baseline data -----------------------------------------------------------

base <- dat %>% 
  mutate(firstvis = as.numeric(!duplicated(id, fromLast=F))) %>% 
  filter(firstvis==1) %>% 
  select(id, visit, visdate, enroll_period, bcurrpa, beduc, beduchs, bempcurr, bevalc, 
         bevcig, beverwk, bevinj, bhome1yr, bincome, bnevmar)
write_csv(base, file="../data/alive_base.csv")


# Sociodemographic data ---------------------------------------------------

socdem <- dat %>% 
  select(id, visit, visdate, enroll_period, age, m0f1, race, black, anyhomels, inclt5k, work)
write_csv(socdem, file="../data/alive_socdem.csv")


# Comorbidity data --------------------------------------------------------

comorbid <- dat %>% 
  mutate(gfr = ifelse(black==1, gfr_aa, gfr_naa)) %>% 
  select(id, visit, visdate, enroll_period,
         bmimedh,                                                       # Obesity
         lungdx6m, lungtx6m, copd, fev_fvc, maxfev, maxfvc, evrlungdx,  # Lung disease
         diabdx6m, diabtx6m, hba1c, evrdiab,                            # Diabetes
         hbpdx6m, hbptx6m, bpdias, bpsys, evrhbp,                       # High blood pressure
         hcholdx6m, hcholtx6m, evrhchol,                                # High cholesterol
         hrtpbdx6m, hrtpbtx6m, evrhrtpb,                                # Heart problems
         renaldx6m, renaltx6m, gfr, ucreat, uprotn, uprt_crt, evrrenal, # Renal disease
         strkedx6m, strketx6m, evrstrke,                                # Stroke
         seizrdx6m, seizrtx6m, evrseizr,                                # Seizure
         fbscan, fbscgrp, fbscge12_3,                                   # Liver fibrosis
         hemo)                                                          # Anemia
write_csv(comorbid, file="../data/alive_comorbid.csv")


# Mental health data ------------------------------------------------------

mental <- dat %>% 
  select(id, visit, visdate, enroll_period,  
         cesd23, cesdtot,                   # CESD score
         anxdx6m, anxtx6m, emedanx,         # Anxiety
         depdx6m, deptx6m, emeddep,         # Depression
         manddx6m, mandtx6m, emedmand,      # Manic depression
         schidx6m, schitx6m, emedschi,      # Schizophrenia
         ehspmntl)                          # Ever hospitalized for mental health
write_csv(mental, file="../data/alive_mental.csv")


# HIV/HCV data ------------------------------------------------------------

hiv <- dat %>% 
  select(id, visit, visdate, enroll_period, hiv, hivvl, hrt6m, hcvvis)
write_csv(hiv, file="../data/alive_hiv.csv")


# Health care data --------------------------------------------------------

care <- dat %>% 
  select(id, visit, visdate, enroll_period,
         primcare, samedoc, samedoc2, intrfsoc,                                       # General care info
         er6m, inpat6m, outpat6m, outer6m, ivamed6m,                                  # Health care visits
         insur6m, ihmo6m, imcaid6m, imcaidexpd, imcare6m, iother6m, ipriv6m, irwhite) # Insurance
write_csv(care, file="../data/alive_care.csv")


# Drug use data -----------------------------------------------------------

drug <- dat %>% 
  select(id, visit, visdate, enroll_period,
         acdrgtx, prscrbbup, prscrbmeth,                                           # Drug treatments
         alcheavy, alcmoderate, auditcat, audittot, auditgrp,                      # Alcohol use
         cigyn, cigpack,                                                           # Smoking
         curuser, ivstat, injcoc, injher, injothdrg, injpainklr, crystali, spdbal, # Injection drug use
         anntivwomj, crack, crystal, mj, smkher, snrtcoc, snrther, benzost,        # Non-injection drug use
         buprest, clonist, hallucgnntdoc, methst, oxycost, painkilrntdoc, percost, 
         sedativntdoc, stimultntdoc, tranqulntdoc) %>% 
  # Create indicators for any drug treatment and any non-injection drug use
  mutate(drgtrt = (acdrgtx==1) | (prscrbbup==1) | (prscrbmeth==1),
         drg_nonidu = rowSums(dat[, c("anntivwomj", "crack", "crystal", "mj", "smkher", "snrtcoc", 
                                       "snrther", "benzost", "buprest", "clonist", "hallucgnntdoc", 
                                       "methst", "oxycost", "painkilrntdoc", "percost", "sedativntdoc", 
                                       "stimultntdoc", "tranqulntdoc")] == 1, na.rm=T) > 0)
write_csv(drug, file="../data/alive_drug.csv")


# Quality of life data ----------------------------------------------------

qual <- dat %>% 
  select(id, visit, visdate, enroll_period, gnbdypn, gnhlthst, mosphs, mosmhs, jailge7d)
write_csv(qual, file="../data/alive_qual.csv")


# Mortality data ----------------------------------------------------------

deaths <- tibble(jackiereq5) %>% 
  filter(dead20f==1) %>% 
  mutate(dod = ymd(dthdate20f),
         yod = year(dod)) %>% 
  select(id, dod, yod)

write_csv(deaths, file="../data/alive_mort.csv")



