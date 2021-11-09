
###################################################################################################################
#
# Purpose: Explore ALIVE multimorbidity data (6 disease definition)
#
# Author: Jacqueline Rudolph
#
# Last Update: 06 Oct 2021
#
###################################################################################################################

packages <- c("tidyverse", "data.table", "RColorBrewer", "zoo", "ggpubr")
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

multimorbid0719 <- read_csv(file="../data/alive_multimorbid6.csv")


# Distribution overall ----------------------------------------------------

# Distribution across all calendar years
ggplot(data=multimorbid0719, aes(x=cum_multi)) + thm2 +
  labs(y="Number of Participants\n", x="\nNumber of Comorbidities") +
  geom_bar() +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(0, 6, 1))

# Distribution by calendar year
jpeg(file="../figures/multimorbidity/multi6_n.jpg", height=5, width=8, units="in", res=300)  
ggplot(data=multimorbid0719, aes(fill=as.factor(cum_multi), x=yr)) + thm1 +
  labs(y="Number of Participants\n", x="\nYear", fill="Number of Comorbidities") +
  geom_bar(position="stack") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  guides(fill=guide_legend(ncol=7)) +
  scale_fill_brewer(palette="Blues")
dev.off()

jpeg(file="../figures/multimorbidity/multi6_prop.jpg", height=5, width=8, units="in", res=300)  
ggplot(data=multimorbid0719, aes(fill=as.factor(cum_multi), x=yr)) + thm1 +
  labs(y="Proportion of Participants\n", x="\nYear", fill="Number of Comorbidities") +
  geom_bar(position="fill") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  guides(fill=guide_legend(ncol=7)) +
  scale_fill_brewer(palette="Blues")
dev.off()


# Distribution by HIV status ----------------------------------------------

multimorbid0719$HIV <- recode(multimorbid0719$hiv, `1`="HIV Positive", `0`="HIV Negative")

jpeg(file="../figures/multimorbidity/multi6_n_hiv.jpg", height=5, width=12, units="in", res=300)  
ggplot(data=multimorbid0719, aes(fill=as.factor(cum_multi), x=yr)) + thm1 +
  labs(y="Number of Participants\n", x="\nYear", fill="Number of Comorbidities") +
  geom_bar(position="stack") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  guides(fill=guide_legend(ncol=7)) +
  scale_fill_brewer(palette="Blues") +
  facet_wrap(facets=vars(HIV), ncol=2)
dev.off()

jpeg(file="../figures/multimorbidity/multi6_prop_hiv.jpg", height=5, width=12, units="in", res=300)  
ggplot(data=multimorbid0719, aes(fill=as.factor(cum_multi), x=yr)) + thm1 +
  labs(y="Proportion of Participants\n", x="\nYear", fill="Number of Comorbidities") +
  geom_bar(position="fill") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  guides(fill=guide_legend(ncol=7)) +
  scale_fill_brewer(palette="Blues") +
  facet_wrap(facets=vars(HIV), ncol=2)
dev.off()

jpeg(file="../figures/multimorbidity/multi6_yr-prop_hiv.jpg", height=5, width=12, units="in", res=300)  
ggplot(data=multimorbid0719, aes(fill=as.factor(yr_multi), x=yr)) + thm1 +
  labs(y="Proportion of Participants\n", x="\nYear", fill="Number of Comorbidities") +
  geom_bar(position="fill") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  guides(fill=guide_legend(ncol=7)) +
  scale_fill_brewer(palette="Blues") +
  facet_wrap(facets=vars(HIV), ncol=2)
dev.off()


# Distribution by age categories ------------------------------------------

jpeg(file="../figures/multimorbidity/multi6_n_age.jpg", height=12, width=12, units="in", res=300)  
ggplot(data=multimorbid0719, aes(fill=as.factor(cum_multi), x=yr)) + thm1 +
  labs(y="Number of Participants\n", x="\nYear", fill="Number of Comorbidities") +
  geom_bar(position="stack") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  guides(fill=guide_legend(ncol=7)) +
  scale_fill_brewer(palette="Blues") +
  facet_wrap(facets=vars(age_cat), ncol=2)
dev.off()

jpeg(file="../figures/multimorbidity/multi6_prop_age.jpg", height=12, width=12, units="in", res=300)  
ggplot(data=multimorbid0719, aes(fill=as.factor(cum_multi), x=yr)) + thm1 +
  labs(y="Proportion of Participants\n", x="\nYear", fill="Number of Comorbidities") +
  geom_bar(position="fill") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  guides(fill=guide_legend(ncol=7)) +
  scale_fill_brewer(palette="Blues") +
  facet_wrap(facets=vars(age_cat), ncol=2)
dev.off()


# Distribution by HIV and age ---------------------------------------------

p1 <- ggplot(data=filter(multimorbid0719, hiv==1), aes(fill=as.factor(cum_multi), x=yr)) + thm1 +
        labs(y="Proportion of Participants\n", x="\nYear", fill="Number of Comorbidities", title="HIV Positive") +
        geom_bar(position="fill") +
        scale_y_continuous(expand=c(0,0)) +
        scale_x_continuous(breaks=seq(2007, 2019, 1)) +
        guides(fill=guide_legend(ncol=7)) +
        scale_fill_brewer(palette="Blues") +
        facet_wrap(facets=vars(age_cat), ncol=1)

p2 <- ggplot(data=filter(multimorbid0719, hiv==0), aes(fill=as.factor(cum_multi), x=yr)) + thm1 +
        labs(y="Proportion of Participants\n", x="\nYear", fill="Number of Comorbidities", title="HIV Negative") +
        geom_bar(position="fill") +
        scale_y_continuous(expand=c(0,0)) +
        scale_x_continuous(breaks=seq(2007, 2019, 1)) +
        guides(fill=guide_legend(ncol=7)) +
        scale_fill_brewer(palette="Blues") +
        facet_wrap(facets=vars(age_cat), ncol=1)

jpeg(file="../figures/multimorbidity/multi6_prop_hivage.jpg", height=25, width=16, units="in", res=300)  
ggarrange(p1, p2, ncol=2)
dev.off()

jpeg(file="../figures/multimorbidity/multi6_prop_hiv39.jpg", height=5, width=12, units="in", res=300)  
ggplot(data=filter(multimorbid0719, age_cat=="<40"), aes(fill=as.factor(cum_multi), x=yr)) + thm1 +
  labs(y="Proportion of Participants\n", x="\nYear", fill="Number of Comorbidities") +
  geom_bar(position="fill") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  guides(fill=guide_legend(ncol=7)) +
  scale_fill_brewer(palette="Blues") +
  facet_wrap(facets=vars(HIV), ncol=2)
dev.off()

jpeg(file="../figures/multimorbidity/multi6_prop_hiv49.jpg", height=5, width=12, units="in", res=300)  
ggplot(data=filter(multimorbid0719, age_cat=="40-49"), aes(fill=as.factor(cum_multi), x=yr)) + thm1 +
  labs(y="Proportion of Participants\n", x="\nYear", fill="Number of Comorbidities") +
  geom_bar(position="fill") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  guides(fill=guide_legend(ncol=7)) +
  scale_fill_brewer(palette="Blues") +
  facet_wrap(facets=vars(HIV), ncol=2)
dev.off()

jpeg(file="../figures/multimorbidity/multi6_prop_hiv59.jpg", height=5, width=12, units="in", res=300)  
ggplot(data=filter(multimorbid0719, age_cat=="50-59"), aes(fill=as.factor(cum_multi), x=yr)) + thm1 +
  labs(y="Proportion of Participants\n", x="\nYear", fill="Number of Comorbidities") +
  geom_bar(position="fill") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  guides(fill=guide_legend(ncol=7)) +
  scale_fill_brewer(palette="Blues") +
  facet_wrap(facets=vars(HIV), ncol=2)
dev.off()

jpeg(file="../figures/multimorbidity/multi6_prop_hiv69.jpg", height=5, width=12, units="in", res=300)  
ggplot(data=filter(multimorbid0719, age_cat=="60-69"), aes(fill=as.factor(cum_multi), x=yr)) + thm1 +
  labs(y="Proportion of Participants\n", x="\nYear", fill="Number of Comorbidities") +
  geom_bar(position="fill") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  guides(fill=guide_legend(ncol=7)) +
  scale_fill_brewer(palette="Blues") +
  facet_wrap(facets=vars(HIV), ncol=2)
dev.off()

jpeg(file="../figures/multimorbidity/multi6_prop_hiv70.jpg", height=5, width=12, units="in", res=300)  
ggplot(data=filter(multimorbid0719, age_cat=="70+"), aes(fill=as.factor(cum_multi), x=yr)) + thm1 +
  labs(y="Proportion of Participants\n", x="\nYear", fill="Number of Comorbidities") +
  geom_bar(position="fill") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  guides(fill=guide_legend(ncol=7)) +
  scale_fill_brewer(palette="Blues") +
  facet_wrap(facets=vars(HIV), ncol=2)
dev.off()


# Distribution by enrollment period ---------------------------------------

jpeg(file="../figures/multimorbidity/multi6_n_enroll.jpg", height=25, width=8, units="in", res=300)  
ggplot(data=multimorbid0719, aes(fill=as.factor(cum_multi), x=yr)) + thm1 +
  labs(y="Number of Participants\n", x="\nYear", fill="Number of Comorbidities") +
  geom_bar(position="stack") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  guides(fill=guide_legend(ncol=7)) +
  scale_fill_brewer(palette="Blues") +
  facet_wrap(facets=vars(enroll_period), ncol=1)
dev.off()

jpeg(file="../figures/multimorbidity/multi6_prop_enroll.jpg", height=25, width=8, units="in", res=300)  
ggplot(data=multimorbid0719, aes(fill=as.factor(cum_multi), x=yr)) + thm1 +
  labs(y="Proportion of Participants\n", x="\nYear", fill="Number of Comorbidities") +
  geom_bar(position="fill") +
  scale_y_continuous(expand=c(0,0)) +
  scale_x_continuous(breaks=seq(2007, 2019, 1)) +
  guides(fill=guide_legend(ncol=7)) +
  scale_fill_brewer(palette="Blues") +
  facet_wrap(facets=vars(enroll_period), ncol=1)
dev.off()

