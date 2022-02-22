
###################################################################################################################
#
# Purpose: Explore ALIVE lab data (imputed data)
#
# Author: Jacqueline Rudolph
#
# Last Update: 19 Dec 2021
#
###################################################################################################################

packages <- c("tidyverse", "data.table", "RColorBrewer", "zoo", "ggpubr")
for (package in packages) {
  library(package, character.only=T)
}

source("./plot_thm.R")


# Read in data ------------------------------------------------------------

lab_annual <- read_csv("../data/alive_labs_imputed.csv")


# Diabetes (HbA1c) --------------------------------------------------------
# No data collected in 2006

# Distribution across all visits
summary(select(lab_annual, hba1c))

# Distribution by calendar year
summ_hba1c <- lab_annual %>% 
  group_by(yr) %>% 
  summarize(min = min(hba1c, na.rm=T),
            p5 = quantile(hba1c, probs=0.05, na.rm=T),
            p25 = quantile(hba1c, probs=0.25, na.rm=T),
            mean = mean(hba1c, na.rm=T),
            med = median(hba1c, na.rm=T),
            p75 = quantile(hba1c, probs=0.75, na.rm=T),
            p95 = quantile(hba1c, probs=0.95, na.rm=T),
            max = max(hba1c, na.rm=T))

jpeg(file="../figures/hba1c_imputed.jpg", height=5, width=5, units="in", res=300)  
ggplot(data=summ_hba1c, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="HbA1c (%)\n") +
  geom_ribbon(aes(ymin=p5, ymax=p95), alpha=0.3) +
  geom_ribbon(aes(ymin=p25, ymax=p75), alpha=0.5) +
  geom_line(aes(y=med), linetype="solid") +
  geom_hline(yintercept=6.5, linetype="dashed") +
  scale_x_continuous(expand=c(0,0), breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(limits=c(4, 9), breaks=seq(4, 9, 1))
dev.off()


# Obstructive lung disease (FEV/FVC) --------------------------------------
# fev_fvc, maxfev, maxfvc
# No data in 2006 or 2020

# Distribution across all visits
summary(select(lab_annual, fev_fvc))

# Distribution by calendar year
summ_fevfvc <- lab_annual %>% 
  group_by(yr) %>% 
  summarize(p5 = quantile(fev_fvc, probs=0.05, na.rm=T),
            p25 = quantile(fev_fvc, probs=0.25, na.rm=T),
            mean = mean(fev_fvc, na.rm=T),
            med = median(fev_fvc, na.rm=T),
            p75 = quantile(fev_fvc, probs=0.75, na.rm=T),
            p95 = quantile(fev_fvc, probs=0.95, na.rm=T))

jpeg(file="../figures/lung_imputed.jpg", height=5, width=5, units="in", res=300)
ggplot(data=summ_fevfvc, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="FEV/FVC\n") +
  geom_ribbon(aes(ymin=p5, ymax=p95), alpha=0.3) +
  geom_ribbon(aes(ymin=p25, ymax=p75), alpha=0.5) +
  geom_line(aes(y=med), linetype="solid") +
  geom_hline(yintercept=0.7, linetype="dashed") +
  scale_y_continuous(limits=c(0.55, 0.90), breaks=seq(0.6, 0.9, 0.1)) +
  scale_x_continuous(expand=c(0,0), breaks=seq(2007, 2019, 1)) 
dev.off()


# Renal disease -----------------------------------------------------------
# renaldx6m, renaltx6m, gfr, ucreat, uprotn, uprt_crt, evrrenal, renal_lab
# No data collected in 2006

# Distribution across all visits
summary(select(lab_annual, gfr, uprt_crt))

# Distribution by calendar year
summ_gfr <- lab_annual %>% 
  group_by(yr) %>% 
  summarize(p5 = quantile(gfr, probs=0.05, na.rm=T),
            p25 = quantile(gfr, probs=0.25, na.rm=T),
            mean = mean(gfr, na.rm=T),
            med = median(gfr, na.rm=T),
            p75 = quantile(gfr, probs=0.75, na.rm=T),
            p95 = quantile(gfr, probs=0.95, na.rm=T))

jpeg(file="../figures/gfr_imputed.jpg", height=5, width=5, units="in", res=300)
ggplot(data=summ_gfr, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y=bquote(atop("Glomerular Filtration Rate","(mL/min/1.73"~ m^3~")"))) +
  geom_ribbon(aes(ymin=p5, ymax=p95), alpha=0.3) +
  geom_ribbon(aes(ymin=p25, ymax=p75), alpha=0.5) +
  geom_line(aes(y=med), linetype="solid") +
  geom_hline(yintercept=60, linetype="dashed") +
  scale_y_continuous(expand=c(0,0), limits=c(35, 125), breaks=seq(40, 120, 20)) +
  scale_x_continuous(expand=c(0,0), breaks=seq(2007, 2019, 1))
dev.off()

summ_prot <- lab_annual %>% 
  group_by(yr) %>% 
  summarize(p5 = quantile(uprt_crt, probs=0.05, na.rm=T),
            p25 = quantile(uprt_crt, probs=0.25, na.rm=T),
            mean = mean(uprt_crt, na.rm=T),
            med = median(uprt_crt, na.rm=T),
            p75 = quantile(uprt_crt, probs=0.75, na.rm=T),
            max = max(uprt_crt, na.rm=T))

jpeg(file="../figures/uprot_imputed.jpg", height=5, width=5, units="in", res=300)
ggplot(data=summ_prot, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Urine Protein-Creatinine\n Concentration Ratio (mg/g)\n") +
  geom_ribbon(aes(ymin=p5, ymax=300), alpha=0.3) + # Max values are too large to show
  geom_ribbon(aes(ymin=p25, ymax=p75), alpha=0.5) +
  geom_line(aes(y=med), linetype="solid") +
  geom_hline(yintercept=200, linetype="dashed") +
  scale_x_continuous(expand=c(0,0), breaks=c(2007, 2009, 2011, 2013, 2015, 2017, 2019)) +
  scale_y_continuous(expand=c(0,0), limits=c(0, 300), breaks=seq(0, 300, 50))
dev.off()


# Liver fibrosis ----------------------------------------------------------
# fbscan, fbscgrp, fbscge12_3

# Distribution across all visits
summary(select(lab_annual, fbscan))

# Distribution by calendar year
summ_fbscan <- lab_annual %>% 
  group_by(yr) %>% 
  summarize(p5 = quantile(fbscan, probs=0.05, na.rm=T),
            p25 = quantile(fbscan, probs=0.25, na.rm=T),
            mean = mean(fbscan, na.rm=T),
            med = median(fbscan, na.rm=T),
            p75 = quantile(fbscan, probs=0.75, na.rm=T),
            p95 = quantile(fbscan, probs=0.95, na.rm=T))

jpeg(file="../figures/fbscan_imputed.jpg", height=5, width=5, units="in", res=300)
ggplot(data=summ_fbscan, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Liver Stiffness (kPA)\n") +
  geom_ribbon(aes(ymin=p5, ymax=p95), alpha=0.3) + 
  geom_ribbon(aes(ymin=p25, ymax=p75), alpha=0.5) +
  geom_line(aes(y=med), linetype="solid") +
  geom_hline(yintercept=9.3, linetype="dashed") +
  scale_x_continuous(expand=c(0,0), breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(expand=c(0,0), limits=c(3.5, 28), breaks=seq(4, 28, 2))
dev.off()


# Hypertension ------------------------------------------------------------
# bpdias, bpsys

# Distribution across all visits
summary(select(lab_annual, bpdias, bpsys))

# Distribution by calendar year
summ_bpdias <- lab_annual %>% 
  group_by(yr) %>% 
  summarize(p5 = quantile(bpdias, probs=0.05, na.rm=T),
            p25 = quantile(bpdias, probs=0.25, na.rm=T),
            mean = mean(bpdias, na.rm=T),
            med = median(bpdias, na.rm=T),
            p75 = quantile(bpdias, probs=0.75, na.rm=T),
            p95 = quantile(bpdias, probs=0.95, na.rm=T))

jpeg(file="../figures/bpdias_imputed.jpg", height=5, width=5, units="in", res=300)
ggplot(data=summ_bpdias, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Diastolic Blood Pressure (mm Hg)\n") +
  geom_ribbon(aes(ymin=p5, ymax=p95), alpha=0.3) +
  geom_ribbon(aes(ymin=p25, ymax=p75), alpha=0.5) +
  geom_line(aes(y=med), linetype="solid") +
  geom_hline(yintercept=90, linetype="dashed") +
  scale_x_continuous(expand=c(0,0), breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(expand=c(0,0), limits=c(60, 110), breaks=seq(65, 105, 5))
dev.off()

summ_bpsys <- lab_annual %>% 
  group_by(yr) %>% 
  summarize(p5 = quantile(bpsys, probs=0.05, na.rm=T),
            p25 = quantile(bpsys, probs=0.25, na.rm=T),
            mean = mean(bpsys, na.rm=T),
            med = median(bpsys, na.rm=T),
            p75 = quantile(bpsys, probs=0.75, na.rm=T),
            p95 = quantile(bpsys, probs=0.95, na.rm=T))

jpeg(file="../figures/bpsys_imputed.jpg", height=5, width=5, units="in", res=300)
ggplot(data=summ_bpsys, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Systolic Blood Pressure (mm Hg)\n") +
  geom_ribbon(aes(ymin=p5, ymax=p95), alpha=0.3) + 
  geom_ribbon(aes(ymin=p25, ymax=p75), alpha=0.5) +
  geom_line(aes(y=med), linetype="solid") +
  geom_hline(yintercept=140, linetype="dashed") +
  scale_x_continuous(expand=c(0,0), breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(limits=c(99, 170), breaks=seq(100, 170, 10))
dev.off()


# Obesity -----------------------------------------------------------------
# bmimedh

# Distribution across all visits
summary(select(lab_annual, bmimedh))

# Distribution by calendar year
summ_bmi <- lab_annual %>% 
  group_by(yr) %>% 
  summarize(p5 = quantile(bmimedh, probs=0.05, na.rm=T),
            p25 = quantile(bmimedh, probs=0.25, na.rm=T),
            mean = mean(bmimedh, na.rm=T),
            med = median(bmimedh, na.rm=T),
            p75 = quantile(bmimedh, probs=0.75, na.rm=T),
            p95 = quantile(bmimedh, probs=0.95, na.rm=T))

jpeg(file="../figures/bmi_imputed.jpg", height=5, width=5, units="in", res=300)
ggplot(data=summ_bmi, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y=bquote(atop("Body Mass Index (kg/"~m^2~")", ""))) +
  geom_ribbon(aes(ymin=p5, ymax=p95), alpha=0.3) +
  geom_ribbon(aes(ymin=p25, ymax=p75), alpha=0.5) +
  geom_line(aes(y=med), linetype="solid") +
  geom_hline(yintercept=30, linetype="dashed") +
  scale_x_continuous(expand=c(0,0), breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(limits=c(19, 41), breaks=seq(20, 40, 5))
dev.off()


