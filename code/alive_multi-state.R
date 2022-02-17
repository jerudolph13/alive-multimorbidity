
###################################################################################################################
#
# Purpose: Trends in ALIVE multimorbidity over time
#
# Author: Jacqueline Rudolph
#
# Last Update: 17 Feb 2021
#
###################################################################################################################


packages <- c("tidyverse", "survival", "broom")
for (package in packages) {
  library(package, character.only=T)
}


# Read in data ------------------------------------------------------------

alive <- read_csv("../data/alive_comorbid_imputed.csv")


# Manipulate data ---------------------------------------------------------

alive2 <- alive %>% 
  group_by(id) %>% 
  mutate(
    # Indicator of multimorbidity
    multimorbid = as.numeric(n_cond>=2),
    
    # Indicator of drop out when someone misses a full year of visits
    drop = as.numeric(is.na(lag(yr)) | (yr - lag(yr)>=1)),
    
    # We censor them at the time they meet this definition (i.e., lag(yr)+1)
    yr = ifelse(drop==1, lag(yr) + 1, yr),
    
    # Build cumulative drop out to remove records after censoring
    cum_drop = cumsum(cumsum(cum_drop)),
    
    # Diseases under control
    n_cond_control = hypertension*(1 - hi_bp) + diabetes*(1 - hi_hba1c) + lung_disease*(1 - poor_lung_func) +
                      renal_disease*(1 - renal_dysfunc),
    prop_cond_control = n_cond_control/n_cond,
    all_control = as.numeric(prop_cond_control==1),
    most_control = as.numeric(prop_cond_control>0.5),
    
    # What state is a person in?
    #   0 = censored
    #   1 = no multimorbidity
    #   2 = controlled multimorbidity
    #   3 = uncontrolled multimorbidity
    state = ifelse(drop==1, 0,
            ifelse(multimorbid==0, 1,
            ifelse(multimorbid==1 & all_cond_control==1, 2 ,3)))) %>%
  # Remove any records that exist after a person missed >1 visit
  filter(cum_drop>1) %>% 
  select(-cum_drop) %>% 
  ungroup()


# Survival analysis -------------------------------------------------------

risk <- tidy(survfit(Surv(yr - 1, yr, as.factor(state)) ~ 1, data=alive2)) %>% 
  select(time, estimate, state)

# Transform data set:
  # One record for each year
  # Probability of each state (estimate) in that year


# Visualize ---------------------------------------------------------------

# Need to get cumulative probabilities to create stacked curves
# Order to accumulate:
  # 1. No multimorbidity
  # 2. No multimorbidity + Uncontrolled multimorbidity
  # 3. No multimorbidity + Uncontrolled multimorbidity + Controlled multimorbidity
  # 4. No multimorbidity + Uncontrolled multimorbidity + Controlled multimorbidity + Drop out

thm <- theme_classic() +
  theme(
    #Format text
    axis.title = element_text(family="Helvetica", size=16, color="black"),
    legend.text = element_text(family="Helvetica", size=16, color="black", margin=margin(t=0.25,b=0.25, unit="lines")),
    legend.title = element_text(family="Helvetica", size=16, color="black"),
    axis.text = element_text(family="Helvetica", size=16, color="black"),
    
    #Format axes
    axis.line = element_line(size=0.75),
    axis.ticks = element_line(size=0.75),
    
    #Format legend
    legend.title.align = 0.5,
    legend.position = c(0.85, 0.8),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA),
    legend.direction="vertical",
    legend.box.background = element_rect(colour = "black", size=0.75),
    
    #Add space around plot
    plot.margin = unit(c(2, 2, 2, 2), "lines")
  )

# We could alternatively use geom_area and plot in reverse order
# Create unique legend labels
cols<-c("No multimorbidity"="color1", 
        "Uncontrolled multimorbidity"="color2", 
        "Controlled multimorbidity"="color3",
        "Dropped out"="color4")

ggplot(data=risk2) + thm +
  labs(x="\n Calendar Year", y="Proportion of participants\n") +
  # Curve for 1
  geom_ribbon(aes(x=time, ymin=est_1_lo, ymax=est_1_hi, color="No multimorbidity")) +
  # Curve for 2
  geom_ribbon(aes(x=time, ymin=est_2_lo, ymax=est_2_hi, color="Uncontrolled multimorbidity")) +
  # Curve for 3
  geom_ribbon(aes(x=time, ymin=est_3_lo, ymax=est_3_hi, color="Controlled multimorbidity")) +
  # Curve for 4
  geom_ribbon(aes(x=time, ymin=est_4_lo, ymax=est_4_hi, color="Dropped out")) +
  scale_color_manual(name="State:", values=cols) +
  scale_x_continuous(expand=c(0, 0)) +
  scale_y_continious(expand=c(0, 0))

