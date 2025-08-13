# 1. Analysis of metabolic concentrations 

library(dplyr)
library(tidyverse)
library(tidyr)
library(readxl)

datsMets <- read_excel(r"(C:\Users\Science\YandexDisk\Work\data\fMRS-hp\results\Отдельные анализы\ALL_analisys.xlsx)")

excludedSubjects <- c(1:3, 6, 8,  10, 23, 25, 14, 19, 5, 17, 29 ) 
dats <- datsMets[-excludedSubjects,]

dats <- dats %>% pivot_longer(cols = "act_NAA":"sham_Glx",
                              names_to = "time",
                              values_to = "MetValue") %>% select(subj, time, MetValue)
dats <- dats %>% separate_wider_delim(time, delim = "_", names = c("condition", "met"))


#Get normilised to Cr data (UNEDITED)
dats_Cr <- dats %>% pivot_wider(names_from = met, values_from = MetValue) %>% 
  mutate(NAA_Cr = NAA/Cr,Glx_Cr = Glx/Cr ) %>% 
  select(subNames, condition, time, NAA_Cr, Glx_Cr) %>% 
  pivot_longer(cols = c("NAA_Cr", "Glx_Cr"), names_to = "met", values_to = "MetValue")

################################################################################
## Statistical analysis
#check for normality
normalityStatistics <- dats %>% filter(met == "Glx") %>%
  select(subj, condition, MetValue) %>%
  group_by(condition) %>%
  summarise_all(.funs = funs(statistic = t.test(.)$p.value, 
                             p.value = shapiro.test(.)$p.value))

normalityStatistics <- dats  %>%
  group_by(met, condition) %>% select(MetValue) %>%
  summarise_all(.funs = funs(statistic = t.test(.)$statistics, 
                             p.value = shapiro.test(.)$p.value))
# Glx is non-normal, Cr and NAA - OK

### Descriptive statistics

df_stat <- dats %>% group_by(condition, met) %>% summarize(
  count = n(),
  mean = mean(MetValue, na.rm = TRUE), 
  sd = sd(MetValue, na.rm = TRUE)) 
head(df_stat)



# compute the difference
d <- with(datsWider, 
          datsWider %>% filter(condition == "act") %>% select(MetValue) -
            datsWider %>% filter(condition == "sham") %>% select(MetValue))
# Shapiro-Wilk normality test for the differences
shapiro.test(d$MetValue) # => p-value = 0.6141

datsWider <- dats %>% filter(met == "Cr") %>% select(subj, condition, MetValue)
res <- t.test(MetValue ~ condition, data = datsWider, paired = TRUE)
res

datsWider <- dats %>% filter(met == "NAA") %>% select(subj, condition, MetValue)
res <- t.test(MetValue ~ condition, data = datsWider, paired = TRUE)
res


# S3 method for default
datsWider <- dats %>% filter(met == "Glx") %>% select(subj, condition, MetValue)
res <- wilcox.test(MetValue ~ condition,  data = datsWider)