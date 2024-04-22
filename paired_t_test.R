# Load packages 
library(stats)
library(rstatix)
library(TOSTER)
library(MOTE)
library(tidyverse)
library(janitor)

# Load data

data <- read_csv("replication_data.csv") %>%
  clean_names() %>%
  drop_na() 

paired_data <- data %>%
  select(id, aclr_single_leg_hop1_norm:con_single_leg_hop3_norm) %>%
  rowwise() %>%
  mutate(aclr_hop = mean(c(aclr_single_leg_hop1_norm, aclr_single_leg_hop2_norm,
                           aclr_single_leg_hop3_norm)),
         con_hop = mean(c(con_single_leg_hop1_norm, con_single_leg_hop2_norm,
                          con_single_leg_hop3_norm)))  %>% 
  mutate(difference = con_hop - aclr_hop) %>%
  as.data.frame()


# Prepare data ----------------

long_data <- paired_data %>% 
  select(id, con_hop,  aclr_hop) %>%
  pivot_longer(cols = c("con_hop", "aclr_hop"),
               names_to = "group",
               values_to = "normalised_hop_distance") 

## Descriptives -------------------------------------

summary_rep <- long_data %>%
  group_by(group) %>%
  summarise(count = n(),
            mean = mean(normalised_hop_distance),
            sd = sd(normalised_hop_distance)) %>%
  mutate(mean_diff = mean(paired_data$difference), 
         sd_diff = sd(paired_data$difference)
  )
summary_rep

## Resolving assumptions  ------------------------------------
### Distribution ---------------------------------------

ggplot(long_data, aes(normalised_hop_distance)) +
  geom_histogram(color="black", fill="white", 
                 bins = 10) +
  facet_wrap(~ group,
          labeller = label_both)

ggplot(long_data, aes(group, normalised_hop_distance, color = group)) +
  geom_boxplot(show.legend = FALSE) +
  theme_minimal()

  ggplot(long_data, aes(group, normalised_hop_distance, color = group)) +  
  geom_violin(fill = "light gray") +
  geom_boxplot(width = .07,
               fill = "white") +
  geom_jitter(position = position_jitter(0.21)) +
  stat_summary(fun = mean,
               geom = "point",
               shape = 12,
               color = "black",
               size = 5) +
  theme_bw()

### Outliers on difference score -----------------------------------

paired_data %>%
  identify_outliers(difference)

### Normality ----------------------------------------------------------
  
paired_data %>% 
    shapiro_test(difference) 
  
# Paired t-test ---------------------------------------------------

long_data$group <- as.factor(long_data$group)
  
# R compares conditions alphabetically, I am reordering here to match the original study
  
  long_data$group <- forcats::fct_relevel(long_data$group, "con_hop", "aclr_hop")
  
  
results <- t.test(normalised_hop_distance ~ group, long_data, 
                  alternative = "two.sided", paired = TRUE, conf.level = 0.95) %>%
  tidy()
results

# Analyse the replication ------

## Calculate replication ES ------

rep_dz <- d.dep.t.diff(mdiff = summary_rep$mean_diff[1], sddiff = summary_rep$sd_diff[1], 
                       n = summary_rep$count[1], a = 0.05)
rep_dz

## Original values ------

orig_values <- data.frame(
  ori_pval = 0.0001,
  N = 15,
  m1 = 197.7,
  sd1 = 26.1,
  m2 = 183.9,
  sd2 = 26.1,
  d = 0.53
)

## Calculate Original ES  ------

quantile = 1 - (orig_values$ori_pval/2)# for two-tailed

orig_tval <- qt(quantile, df = 14) 

orig_dz <- d.dep.t.diff.t(t = orig_tval, n = 15, a = 0.05)
orig_dz

# the original study reported d = 0.53 which seems to be Cohen's dav

orig_dav <- d.dep.t.avg(m1=orig_values$m1, m2=orig_values$m2, 
                        sd1=orig_values$sd1, sd2=orig_values$sd2, 
                        n = 15, a = 0.05)
orig_dav

## Z-test  --------

rep_test <- compare_smd(
  smd1 = orig_dz$d,
  n1 = orig_values$N,
  smd2 = rep_dz$d,
  n2 = summary_rep$count[1],
  paired = TRUE,
  alternative = "greater")
rep_test

