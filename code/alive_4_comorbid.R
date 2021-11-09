
###################################################################################################################
#
# Purpose: Explore ALIVE comorbidity data
#
# Author: Jacqueline Rudolph
#
# Last Update: 06 Oct 2021
#
###################################################################################################################

packages <- c("tidyverse", "data.table")
for (package in packages) {
  library(package, character.only=T)
}

source("./plot_thm.R")

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

dat0719 <- read_csv(file="../data/alive_multimorbid6.csv")
dat1419 <- read_csv(file="../data/alive_multimorbid7.csv")


# Diabetes ----------------------------------------------------------------

# Distribution across all visits
table(dat0719$diab_lab)
table(dat0719$cum_diab)

# Distribution by calendar year
summ_diab <- dat0719 %>%
  group_by(yr) %>% 
  summarize(n_cum = sum(cum_diab, na.rm=T),
            n_yr = sum(diab_lab, na.rm=T),
            n_tot = n(),
            prop_cum = n_cum/n_tot,
            prop_yr = n_yr/n_tot)

# Annual count figure
cols <- c("Total"="#a6bddb", "With diabetes"="#2b8cbe")
jpeg(file="../figures/single_diseases/diabetes_yr-n.jpg", height=5, width=8, units="in", res=300)  
ggplot(data=summ_diab, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Numer of participants") +
  geom_bar(aes(y=n_tot, fill="Total"), stat="identity") +
  geom_bar(aes(y=n_yr, fill="With diabetes"), stat="identity") +
  geom_text(aes(y=n_yr+30, label=ifelse(prop_yr==0, "", 
                                         paste0(round(prop_yr*100, digits=0),"%")))) +
  scale_fill_manual(name="", values=cols) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(expand=c(0,0))
dev.off()

# Cumulative count figure
cols <- c("Total"="#a6bddb", "With diabetes"="#2b8cbe")
jpeg(file="../figures/single_diseases/diabetes_cum-n.jpg", height=5, width=8, units="in", res=300)  
ggplot(data=summ_diab, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Numer of participants") +
  geom_bar(aes(y=n_tot, fill="Total"), stat="identity") +
  geom_bar(aes(y=n_cum, fill="With diabetes"), stat="identity") +
  geom_text(aes(y=n_cum+30, label=ifelse(prop_cum==0, "", 
                                                paste0(round(prop_cum*100, digits=0),"%")))) +
  scale_fill_manual(name="", values=cols) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(expand=c(0,0))
dev.off()


# Obstructive lung disease ------------------------------------------------

# Distribution across all visits
table(dat0719$lung_lab)
table(dat0719$cum_lung)

# Distribution by calendar year
summ_lung <- dat0719 %>%
  group_by(yr) %>% 
  summarize(n_cum = sum(cum_lung, na.rm=T),
            n_yr = sum(lung_lab, na.rm=T),
            n_tot = n(),
            prop_cum = n_cum/n_tot,
            prop_yr = n_yr/n_tot)

# Annual count figure
cols <- c("Total"="#a6bddb", "With lung disease"="#2b8cbe")
jpeg(file="../figures/single_diseases/lung_yr-n.jpg", height=5, width=8, units="in", res=300)  
ggplot(data=summ_lung, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Numer of participants") +
  geom_bar(aes(y=n_tot, fill="Total"), stat="identity") +
  geom_bar(aes(y=n_yr, fill="With lung disease"), stat="identity") +
  geom_text(aes(y=n_yr+30, label=ifelse(prop_yr==0, "", 
                                        paste0(round(prop_yr*100, digits=0),"%")))) +
  scale_fill_manual(name="", values=cols) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(expand=c(0,0))
dev.off()

# Cumulative count figure
cols <- c("Total"="#a6bddb", "With lung disease"="#2b8cbe")
jpeg(file="../figures/single_diseases/lung_cum-n.jpg", height=5, width=8, units="in", res=300)  
ggplot(data=summ_lung, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Numer of participants") +
  geom_bar(aes(y=n_tot, fill="Total"), stat="identity") +
  geom_bar(aes(y=n_cum, fill="With lung disease"), stat="identity") +
  geom_text(aes(y=n_cum+30, label=ifelse(prop_cum==0, "", 
                                          paste0(round(prop_cum*100, digits=0),"%")))) +
  scale_fill_manual(name="", values=cols) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(expand=c(0,0))
dev.off()


# Renal disease -----------------------------------------------------------

# Distribution across all visits
table(dat0719$renal_lab)
table(dat0719$cum_renal)

# Distribution by calendar year
summ_renal <- dat0719 %>%
  group_by(yr) %>% 
  summarize(n_cum = sum(cum_renal, na.rm=T),
            n_yr = sum(renal_lab, na.rm=T),
            n_tot = n(),
            prop_cum = n_cum/n_tot,
            prop_yr = n_yr/n_tot)

# Annual count figure
cols <- c("Total"="#a6bddb", "With renal disease"="#2b8cbe")
jpeg(file="../figures/single_diseases/renal_yr-n.jpg", height=5, width=8, units="in", res=300)  
ggplot(data=summ_renal, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Numer of participants") +
  geom_bar(aes(y=n_tot, fill="Total"), stat="identity") +
  geom_bar(aes(y=n_yr, fill="With renal disease"), stat="identity") +
  geom_text(aes(y=n_yr+30, label=ifelse(prop_yr==0, "", 
                                        paste0(round(prop_yr*100, digits=0),"%")))) +
  scale_fill_manual(name="", values=cols) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(expand=c(0,0))
dev.off()

# Cumulative count figure
cols <- c("Total"="#a6bddb", "With renal disease"="#2b8cbe")
jpeg(file="../figures/single_diseases/renal_cum-n.jpg", height=5, width=8, units="in", res=300)  
ggplot(data=summ_renal, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Numer of participants") +
  geom_bar(aes(y=n_tot, fill="Total"), stat="identity") +
  geom_bar(aes(y=n_cum, fill="With renal disease"), stat="identity") +
  geom_text(aes(y=n_cum+30, label=ifelse(prop_cum==0, "", 
                                           paste0(round(prop_cum*100, digits=0),"%")))) +
  scale_fill_manual(name="", values=cols) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(expand=c(0,0))
dev.off()


# Liver fibrosis ----------------------------------------------------------

# Distribution across all visits
table(dat0719$liver_lab)
table(dat0719$cum_liver)

# Distribution by calendar year
summ_liver <- dat0719 %>%
  group_by(yr) %>% 
  summarize(n_cum = sum(cum_liver, na.rm=T),
            n_yr = sum(liver_lab, na.rm=T),
            n_tot = n(),
            prop_cum = n_cum/n_tot,
            prop_yr = n_yr/n_tot)

# Annual count figure
cols <- c("Total"="#a6bddb", "With liver disease"="#2b8cbe")
jpeg(file="../figures/single_diseases/liver_yr-n.jpg", height=5, width=8, units="in", res=300)  
ggplot(data=summ_liver, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Numer of participants") +
  geom_bar(aes(y=n_tot, fill="Total"), stat="identity") +
  geom_bar(aes(y=n_yr, fill="With liver disease"), stat="identity") +
  geom_text(aes(y=n_yr+30, label=ifelse(prop_yr==0, "", 
                                        paste0(round(prop_yr*100, digits=0),"%")))) +
  scale_fill_manual(name="", values=cols) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(expand=c(0,0))
dev.off()

# Cumulative count figure
cols <- c("Total"="#a6bddb", "With liver disease"="#2b8cbe")
jpeg(file="../figures/single_diseases/liver_cum-n.jpg", height=5, width=8, units="in", res=300)  
ggplot(data=summ_liver, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Numer of participants") +
  geom_bar(aes(y=n_tot, fill="Total"), stat="identity") +
  geom_bar(aes(y=n_cum, fill="With liver disease"), stat="identity") +
  geom_text(aes(y=n_cum+30, label=ifelse(prop_cum==0, "", 
                                           paste0(round(prop_cum*100, digits=0),"%")))) +
  scale_fill_manual(name="", values=cols) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(expand=c(0,0))
dev.off()


# Hypertension ------------------------------------------------------------

# Distribution across all visits
table(dat0719$hypertension_lab)
table(dat0719$cum_hyper)

# Distribution by calendar year
summ_hyper <- dat0719 %>%
  group_by(yr) %>% 
  summarize(n_cum = sum(cum_hyper, na.rm=T),
            n_yr = sum(hypertension_lab, na.rm=T),
            n_tot = n(),
            prop_cum = n_cum/n_tot,
            prop_yr = n_yr/n_tot)

# Annual count figure
cols <- c("Total"="#a6bddb", "With hypertension"="#2b8cbe")
jpeg(file="../figures/single_diseases/hypertension_yr-n.jpg", height=5, width=8, units="in", res=300)  
ggplot(data=summ_hyper, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Numer of participants") +
  geom_bar(aes(y=n_tot, fill="Total"), stat="identity") +
  geom_bar(aes(y=n_yr, fill="With hypertension"), stat="identity") +
  geom_text(aes(y=n_yr+30, label=ifelse(prop_yr==0, "", 
                                        paste0(round(prop_yr*100, digits=0),"%")))) +
  scale_fill_manual(name="", values=cols) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(expand=c(0,0))
dev.off()

# Cumulative count figure
cols <- c("Total"="#a6bddb", "With hypertension"="#2b8cbe")
jpeg(file="../figures/single_diseases/hypertension_cum-n.jpg", height=5, width=8, units="in", res=300)  
ggplot(data=summ_hyper, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Numer of participants") +
  geom_bar(aes(y=n_tot, fill="Total"), stat="identity") +
  geom_bar(aes(y=n_cum, fill="With hypertension"), stat="identity") +
  geom_text(aes(y=n_cum+30, label=ifelse(prop_cum==0, "", 
                                           paste0(round(prop_cum*100, digits=0),"%")))) +
  scale_fill_manual(name="", values=cols) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(expand=c(0,0))
dev.off()


# Obesity -----------------------------------------------------------------

# Distribution across all visits
table(dat0719$obese_lab)
table(dat0719$cum_obese)

# Distribution by calendar year
summ_obese <- dat0719 %>%
  group_by(yr) %>% 
  summarize(n_cum = sum(cum_obese, na.rm=T),
            n_yr = sum(obese_lab, na.rm=T),
            n_tot = n(),
            prop_cum = n_cum/n_tot,
            prop_yr = n_yr/n_tot)

# Annual count figure
cols <- c("Total"="#a6bddb", "With obesity"="#2b8cbe")
jpeg(file="../figures/single_diseases/obesity_yr-n.jpg", height=5, width=8, units="in", res=300)  
ggplot(data=summ_obese, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Numer of participants") +
  geom_bar(aes(y=n_tot, fill="Total"), stat="identity") +
  geom_bar(aes(y=n_yr, fill="With obesity"), stat="identity") +
  geom_text(aes(y=n_yr+30, label=ifelse(prop_yr==0, "", 
                                        paste0(round(prop_yr*100, digits=0),"%")))) +
  scale_fill_manual(name="", values=cols) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(expand=c(0,0))
dev.off()

# Cumulative count figure
cols <- c("Total"="#a6bddb", "With obesity"="#2b8cbe")
jpeg(file="../figures/single_diseases/obesity_cum-n.jpg", height=5, width=8, units="in", res=300)  
ggplot(data=summ_obese, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Numer of participants") +
  geom_bar(aes(y=n_tot, fill="Total"), stat="identity") +
  geom_bar(aes(y=n_cum, fill="With obesity"), stat="identity") +
  geom_text(aes(y=n_cum+30, label=ifelse(prop_cum==0, "", 
                                           paste0(round(prop_cum*100, digits=0),"%")))) +
  scale_fill_manual(name="", values=cols) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  scale_y_continuous(expand=c(0,0))
dev.off()


# Anemia ------------------------------------------------------------------

# Distribution across all visits
table(dat1419$anemia_lab)
table(dat1419$cum_anemia)

# Distribution by calendar year
summ_anemia <- dat1419 %>%
  group_by(yr) %>% 
  summarize(n_cum = sum(cum_anemia, na.rm=T),
            n_yr = sum(anemia_lab, na.rm=T),
            n_tot = n(),
            prop_cum = n_cum/n_tot,
            prop_yr = n_yr/n_tot)

# Annual count figure
cols <- c("Total"="#a6bddb", "With anemia"="#2b8cbe")
jpeg(file="../figures/single_diseases/anemia_yr-n.jpg", height=5, width=5, units="in", res=300)  
ggplot(data=summ_anemia, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Numer of participants") +
  geom_bar(aes(y=n_tot, fill="Total"), stat="identity") +
  geom_bar(aes(y=n_yr, fill="With anemia"), stat="identity") +
  geom_text(aes(y=n_yr+30, label=ifelse(prop_yr==0, "", 
                                        paste0(round(prop_yr*100, digits=0),"%")))) +
  scale_fill_manual(name="", values=cols) +
  scale_x_continuous(breaks=seq(2014, 2019, 1)) +
  scale_y_continuous(expand=c(0,0))
dev.off()

# Cumulative count figure
cols <- c("Total"="#a6bddb", "With anemia"="#2b8cbe")
jpeg(file="../figures/single_diseases/anemia_cum-n.jpg", height=5, width=5, units="in", res=300)  
ggplot(data=summ_anemia, aes(x=yr)) + thm1 +
  labs(x="\nCalendar Year", y="Numer of participants") +
  geom_bar(aes(y=n_tot, fill="Total"), stat="identity") +
  geom_bar(aes(y=n_cum, fill="With anemia"), stat="identity") +
  geom_text(aes(y=n_cum+30, label=ifelse(prop_cum==0, "", 
                                            paste0(round(prop_cum*100, digits=0),"%")))) +
  scale_fill_manual(name="", values=cols) +
  scale_x_continuous(breaks=seq(2014, 2019, 1)) +
  scale_y_continuous(expand=c(0,0))
dev.off()
