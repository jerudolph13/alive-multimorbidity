
###########################################################################
#
# Purpose: Read in and manipulate ALIVE cancer data
#
# Author: Jacqueline Rudolph
#
# Last Update: 26 Apr 2022
#
###########################################################################

library("tidyverse")
library("lubridate")


# Read in data ------------------------------------------------------------

load(file="../../../raw_data/alive/jackiereq6.rdata")
surv <- read_csv("../data/cancer_survival.csv")


# Diagnoses ---------------------------------------------------------------

cancer.dx <- tibble(jackiereq6) %>% 
  filter(MCR_DxDate!="NA") %>% 
  mutate(MCR_DxDate = ymd(MCR_DxDate),
         MCR_cancertype = ifelse(MCR_cancertype=="NA", NA, MCR_cancertype),
         MCR_Psite = ifelse(MCR_Psite=="NA", NA, MCR_Psite),
         
         # Combine stage data
         MCR_DerivedSS2000 = as.numeric(MCR_DerivedSS2000),
         
         MCR_SEERSumStg1977 = as.numeric(MCR_SEERSumStg1977),
         MCR_SEERSumStg1977 = ifelse(is.na(MCR_SEERSumStg1977), 99, MCR_SEERSumStg1977),
         
         MCR_SEERSumStg2000 = as.numeric(MCR_SEERSumStg2000),
         MCR_SEERSumStg2000 = ifelse(is.na(MCR_SEERSumStg2000) & !is.na(MCR_DerivedSS2000), 
                                     MCR_DerivedSS2000, MCR_SEERSumStg2000),
         MCR_SEERSumStg2000 = ifelse(is.na(MCR_SEERSumStg2000), 99, MCR_SEERSumStg2000),
         
         MCR_Summary_Stage_2018 = as.numeric(MCR_Summary_Stage_2018),
         MCR_Summary_Stage_2018 = ifelse(is.na(MCR_Summary_Stage_2018), 99, MCR_Summary_Stage_2018),
         
         MCR_DerivedSS2000 = as.numeric(MCR_DerivedSS2000),
         stage = ifelse(MCR_SEERSumStg2000==0 | MCR_SEERSumStg1977==0 | MCR_Summary_Stage_2018==0, "in situ",
                 ifelse(MCR_SEERSumStg2000==1 | MCR_SEERSumStg1977==1 | MCR_Summary_Stage_2018==1, "localized",
                 ifelse(MCR_SEERSumStg2000<=5 | MCR_SEERSumStg1977<=5 | MCR_Summary_Stage_2018<=5, "regional",
                 ifelse(MCR_SEERSumStg2000==7 | MCR_SEERSumStg1977==7 | MCR_Summary_Stage_2018==7, "distant",
                 ifelse(MCR_SEERSumStg2000==8 | MCR_SEERSumStg1977==8 | MCR_Summary_Stage_2018==8, "benign",
                 ifelse(MCR_SEERSumStg2000==9 | MCR_SEERSumStg1977==9 | MCR_Summary_Stage_2018==9, "unknown", NA))))))) %>% 
  rename(dx_date = MCR_DxDate, cancer_type = MCR_cancertype, psite = MCR_Psite) %>% 
  select(id, dx_date, cancer_type, psite, stage)


# Treatments --------------------------------------------------------------

cancer.tx <- tibble(jackiereq6) %>% 
  filter(MCR_DxDate!="NA") %>% 
  mutate(MCR_DxDate = ymd(MCR_DxDate),
         yr_BRM = as.numeric(substr(MCR_RxDateBRM, 1, 4)),
         yr_chemo = as.numeric(substr(MCR_RxDateChemo, 1, 4)),
         yr_horm = as.numeric(substr(MCR_RxDateHorm, 1, 4)),
         yr_mdsurg = as.numeric(substr(MCR_RxDateMostDefSurg, 1, 4)),
         yr_oth = as.numeric(substr(MCR_RxDateOth, 1, 4)),
         yr_rad = as.numeric(substr(MCR_RxDateRad, 1, 4)),
         yr_surg = as.numeric(substr(MCR_RxDateSurg, 1, 4)),
         yr_sys = as.numeric(substr(MCR_RxDateSystemic, 1, 4))) %>% 
  rowwise() %>% 
  mutate(yr_first_trt = min(c(yr_BRM, yr_chemo, yr_horm, yr_mdsurg, yr_oth, 
                              yr_rad, yr_surg, yr_sys), na.rm=T),
         yr_first_trt = ifelse(is.infinite(yr_first_trt), NA, yr_first_trt)) %>% 
  ungroup() %>% 
  rename(dx_date=MCR_DxDate) %>% 
  select(id, dx_date, yr_first_trt, yr_surg, yr_rad, yr_chemo)


# Combine data ------------------------------------------------------------

cancer <- bind_cols(cancer.dx, select(cancer.tx, -c("id", "dx_date"))) %>% 
  # Determine if standard treatment was received
  mutate(std_trt = case_when(
    # Regardless of stage, blood cancers that receive chemo/radiation are properly treated
    (cancer_type %in% c("chronlymphleuk", "chronmyelleuk", "acutemyelleuk")) & !is.na(yr_chemo) ~ 1,
    cancer_type=="hodgelymph" & !is.na(yr_chemo) ~ 1,
    cancer_type=="multmyeloma" & !is.na(yr_chemo) ~ 1,
    cancer_type=="nonhodgelymph" & !is.na(yr_chemo) ~ 1,
    
    # Now determine for all other types    
    cancer_type=="otherilldefined" | stage=="unknown" ~ 99, #case_when does not like NA
    stage=="distant" ~ 0,
    stage=="benign" ~ 1,
    cancer_type=="anus" & (!is.na(yr_chemo) & !is.na(yr_rad)) ~ 1,
    cancer_type=="AINgradeIII" ~ 1,
    cancer_type=="bladder" & stage=="local" & !is.na(yr_surg) ~ 1,
    cancer_type=="bladder" & stage=="regional" & (!is.na(yr_surg) & !is.na(yr_chemo)) ~ 1,
    cancer_type=="insitubadder" & !is.na(yr_surg) ~ 1,
    cancer_type=="brain" & (!is.na(yr_surg) | !is.na(yr_rad)) ~ 1,
    (cancer_type %in% c("fembreast", "malebreast")) & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="insitufembreast" & !is.na(yr_surg) ~ 1,
    cancer_type=="cervix" & (!is.na(yr_surg) | !is.na(yr_chemo)) ~ 1,
    cancer_type=="insitucervix" & !is.na(yr_surg) ~ 1,
    cancer_type=="CINgradeIII" & !is.na(yr_surg) ~ 1,
    cancer_type=="colonnotrect" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="esophagus" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="gallbladder" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="gumothermouth" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="kapsarc" ~ 1, # Check later if they got ART?
    cancer_type=="kidneyren" & !is.na(yr_surg) ~ 1,
    cancer_type=="larynx" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="liverandduct" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="lungbronch" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="mesothelio" & (!is.na(yr_surg) | !is.na(yr_chemo)) ~ 1,
    cancer_type=="oropharynx" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="othdig" & (!is.na(yr_surg) | !is.na(yr_chemo)) ~ 1,
    cancer_type=="otherbucphar" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="otherendoc" ~ 1,
    cancer_type=="othfemgen" & stage=="local" & !is.na(yr_surg) ~ 1,
    cancer_type=="othfemgen" & stage=="regional" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="ovary" & (!is.na(yr_surg) | !is.na(yr_chemo)) ~ 1,
    cancer_type=="pancreas" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="prostate" & stage=="local" & !is.na(yr_surg) ~ 1,
    cancer_type=="prostate" & stage=="regional" & (!is.na(yr_surg) | !is.na(yr_rad)) ~ 1,
    cancer_type=="rectandsig" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="skinmelan" & !is.na(yr_chemo) ~ 1,
    cancer_type=="smallintest" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="softtissue" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="stomach" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="testis" & !is.na(yr_surg) ~ 1,
    cancer_type=="thyroid" & !is.na(yr_surg) ~ 1,
    cancer_type=="tongue" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="pancreas" & (!is.na(yr_surg) | !is.na(yr_chemo) | !is.na(yr_rad)) ~ 1,
    cancer_type=="ureter" & !is.na(yr_surg) ~ 1,
    cancer_type=="VAINgradeIII" & !is.na(yr_surg) ~ 1,
    cancer_type=="VINgradeIII" & !is.na(yr_surg) ~ 1,
    TRUE ~ 0),
    std_trt = ifelse(std_trt==99, NA, std_trt))


# Survival ----------------------------------------------------------------

surv2 <- surv %>% 
  filter(race=="All races" & 
           ((sex=="Male and female" & !(site %in% c("Cervix Uteri", "Corpus and Uterus, NOS", "Ovary", "Prostate", "Testis"))) |
            (sex=="Female" & site %in% c("Cervix Uteri", "Corpus and Uterus, NOS", "Ovary")) |
            (sex=="Male" & site %in% c("Prostate", "Testis")))) %>%
  mutate(poor_prog = as.numeric(survival < 0.5)) %>% # Change threshold and see how results differ
  select(-c(sex, race))
surv2$stage <- recode(surv2$stage, "Unknown/unstaged" = "unknown",
                                   "Localized" = "localized",
                                   "Regional" = "regional",
                                   "Distant" = "distant")

cancer2 <- cancer %>% 
  mutate(site = case_when(
    cancer_type=="acutemyelleuk" ~ "Acute Myeloid Leukemia",
    cancer_type=="anus" ~ "Colon and Rectum",
    cancer_type=="bladder" ~ "Urinary Bladder",
    cancer_type=="brain" ~ "Brain and Other Nervous System",
    cancer_type=="cervix" ~ "Cervix Uteri",
    cancer_type=="chronlymphleuk" ~ "Chronic Lymphocytic Leukemia",
    cancer_type=="chronmyelleuk" ~ "Chronic Myeloid Leukemia",
    cancer_type %in% c("colonnotrect", "rectandsig") ~ "Colon and Rectum",
    cancer_type=="esophagus" ~ "Esophagus",
    cancer_type %in% c("fembreast", "malebreast") ~ "Breast",
    cancer_type=="gallbladder" ~ "Gallbladder",
    cancer_type %in% c("gumothermouth", "oropharynx", "otherbucphar", 
                       "tongue") ~ "Oral Cavity and Pharynx",
    cancer_type=="hodgelymph" ~ "Hodgkin Lymphoma",
    cancer_type=="kapsarc" ~ "Kaposi Sarcoma",
    cancer_type=="kidneyren" ~ "Kidney and Renal Pelvis",
    cancer_type=="larynx" ~ "Larynx",
    cancer_type=="liverandduct" ~ "Liver and Intrahepatic Bile Duct",
    cancer_type=="lungbronch" ~ "Lung and Bronchus",
    cancer_type=="mesothelio" ~ "Mesothelioma",
    cancer_type=="multmyeloma" ~ "Myeloma",
    cancer_type=="nonhodgelymph" ~ "Non-Hodgkin Lymphoma",
    cancer_type=="ovary" ~ "Ovary",
    cancer_type=="pancreas" ~ "Pancreas",
    cancer_type=="prostate" ~ "Prostate",
    cancer_type=="skinmelan" ~ "Melanoma of the Skin",
    cancer_type=="softtissue" ~ "Soft Tissue including Heart",
    cancer_type=="stomach" ~ "Stomach",
    cancer_type=="testis" ~ "Testis",
    cancer_type=="thyroid" ~ "Thyroid",
    cancer_type %in% c("AINgradeIII", "benigncnsendo", "benignothnerve",
                       "CINgradeIII", "insitubladder", "insitucervix",
                       "insitufembreast", "VAINgradeIII", "VINgradeIII",
                       "othdig", "otherendoc", "otherilldefined", 
                       "othfemgen", "smallintest", "ureter") ~ "NA")) %>% 
  left_join(surv2, by=c("site", "stage"))
  

# Output data -------------------------------------------------------------

write_csv(cancer2, file="../data/alive_cancer.csv")
