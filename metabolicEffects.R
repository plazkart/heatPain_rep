#get data

library(dplyr)
library(tidyverse)
library(tidyr)

# 1. Search differences between different time points in dynamic
#get data
dats <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\spectraDynamic_sm.csv')
#bold-corrected dAta
dats <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\spectra-TP6-SM-BC.csv')

subjectNames <- paste0("sub_", 1:32)
dats <- dats %>% mutate(subNames = subjectNames) 
#exclude subjects QA
excludedSubjects <- c(1:6, 8, 9, 10, 11, 13,17, 20, 22, 23,26, 28, 32 )

dats <- dats[-excludedSubjects,]
#pivot column
dats <- dats %>% pivot_longer(cols = "act_NAA_tp_01_sm":"sham_Glx_tp_06_sm",
                              names_to = "time",
                              values_to = "MetValue") %>% select(subNames, time, MetValue)
dats <- dats %>% separate_wider_delim(time, delim = "_tp_", names = c("condition", "time"))
dats <- dats %>% separate_wider_delim(condition, delim = "_", names = c("condition", "met"))
dats$time = gsub('_sm','',dats$time)
dats$time <- as.numeric(dats$time)

#Get normilised to Cr data
dats_Cr <- dats %>% pivot_wider(names_from = met, values_from = MetValue) %>% 
  mutate(NAA_Cr = NAA/Cr,Glx_Cr = Glx/Cr ) %>% 
  select(subNames, condition, time, NAA_Cr, Glx_Cr) %>% 
  pivot_longer(cols = c("NAA_Cr", "Glx_Cr"), names_to = "met", values_to = "MetValue")

################################################################################
## Statistical analysis
#check for normality
normalityStatistics <- dats %>% filter(met == "Glx", condition == "act") %>%
  select(time, MetValue) %>%
  group_by(time) %>%
  summarise_all(.funs = funs(statistic = t.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
normalityStatistics
normalityStatistics <- dats %>% filter(met == "Glx", condition == "sham") %>%
  select(time, MetValue) %>%
  group_by(time) %>%
  summarise_all(.funs = funs(statistic = t.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
normalityStatistics
# For Glx both conditons are normally distributed for each time point
# For Cr both statistics are normally distributed for each time point
# For NAA statistics in act condition are not normally distributed for a few points

normalityStatistics <- dats_Cr %>% filter(met == "Glx_Cr", condition == "act") %>%
  select(time, MetValue) %>%
  group_by(time) %>%
  summarise_all(.funs = funs(statistic = t.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
normalityStatistics
normalityStatistics <- dats_Cr %>% filter(met == "Glx_Cr", condition == "sham") %>%
  select(time, MetValue) %>%
  group_by(time) %>%
  summarise_all(.funs = funs(statistic = t.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
normalityStatistics

# For Glx/Cr only Act group is normally distributed for each time point
# For NAA/Cr only Sham group is normally distributed for each time point

### T-test
datsWider <- dats %>% filter(met == "Glx", condition == "act") %>% select(subNames, time, MetValue)
getTest <- function(df, x){
  df <- df %>% filter(time == 1 | time == x)
  t.test(MetValue ~ time, data = df, paired = TRUE)
}
tests <- lapply(2:6, function(x) getTest(datsWider, x))
results <- data.frame(
  lag = 1:5, 
  xsquared = sapply(tests, "[[", "statistic"), 
  pvalue = sapply(tests, "[[", "p.value"),
  meanDiference = sapply(tests, "[[", "estimate"),
  stdDiff = sapply(tests, "[[", "stderr")
)
view(results)

#For Cr normilised
datsWider <- dats_Cr %>% filter(met == "Glx_Cr", condition == "act") %>% select(subNames, time, MetValue)
getTest <- function(df, x){
  df <- df %>% filter(time == 1 | time == x)
  t.test(MetValue ~ time, data = df, paired = TRUE)
}
tests <- lapply(2:6, function(x) getTest(datsWider, x))
results <- data.frame(
  lag = 1:5, 
  xsquared = sapply(tests, "[[", "statistic"), 
  pvalue = sapply(tests, "[[", "p.value")
)
view(results)

# Second point for Glx/Cr has trend to increase p = 0.06
# Third point for NAA/Cr increasing p = 0.045

## Glx differs from 1st point in act: 2nd point 0.02
## Glx differs from 1st point in sham: no differences

################################################################################
#Plot the dynamics
####
dats %>% filter(met == 'Glx') %>% group_by(subNames, condition) %>% 
  mutate(MetValue_1 = MetValue[time == 1]) %>%
  mutate(Z = (MetValue-MetValue_1)) %>% ungroup() %>%
  group_by(condition, time) %>% summarise(meanZ = mean(Z), stdZ = sd(Z)) %>%
  ggplot(aes(x = time, y = meanZ, color = condition)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymax = meanZ+stdZ, ymin = meanZ-stdZ)) +
  facet_wrap(~ condition)

#get normilised to Cr data
dats_Cr <- dats %>% pivot_wider(names_from = met, values_from = LWH) %>% 
  mutate(NAA_Cr = NAA/Cr,Glx_Cr = Glx/Cr ) %>% 
  select(subNames, condition, time, NAA_Cr, Glx_Cr)
dats_Cr %>% select(-NAA_Cr) %>% group_by(subNames, condition) %>%
  mutate(Z = (Glx_Cr-mean(Glx_Cr))/sqrt(var(Glx_Cr))) %>% ungroup() %>%
  group_by(condition, time) %>% summarise(meanZ = mean(Z)) %>%
  ggplot(aes(x = time, y = meanZ, color = condition)) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~ condition)