# This scripts are dedicated for statistic quantification of 
# Glx dynamic during short heat pain stimulation.

# 1. Search differences between different time points in dynamic
#get data
dats <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\GLX_normilised.csv')

sham_df <- dats[1:13, ]
bc_df <- dats[-26:-1, ]

#check for normality
apply(sham_df,2,shapiro.test)
apply(bc_df,2,shapiro.test)

#make statistic miracle

# 2. Open SHAM condition results
dats <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\all_all_sham.csv')

# Check normality of concentratons and see historam
apply(dats$Glx_AC,2,shapiro.test)

#3. Process BOLD-MRS data
library(dplyr)
library(tidyverse)
library(tidyr)

dats <- read.csv('C:\\Users\\user\\YandexDisk\\Work\\data\\fMRS-hp\\results\\BOLD_MRS.csv')
columnsActSham <- cbind(paste0("act", "_tp_0", 1:6, "_LWCr"), paste0("sham", "_tp_0", 1:6, "_LWCr"))
dats <- dats[, columnsActSham]
subjectNames <- paste0("sub_", 1:32)
dats <- dats %>% mutate(subNames = subjectNames) 
row.names(dats) <- subjectNames

#3.1 QA check for outliers
# library
library(ggplot2)

# basic scatterplot
dats %>% filter(condition == 'sham') %>% group_by(subNames) %>% mutate(LW_diff=LW/mean(LW)) %>%
ungroup() %>% group_by(time) %>%
ggplot(aes(x=timeD, y=LW_diff, colour=subNames)) + 
  geom_point()


excludedSubjects <- c(1:3, 11)
dats <- dats[-excludedSubjects,]

dats <- dats %>% pivot_longer(cols = "act_tp_01_LWCr":"sham_tp_06_LWCr",
                              names_to = "time",
                              values_to = "LW")

dats <- dats %>% 
  mutate(condition = case_when(
    str_detect(time, "act") ~ "act",
    str_detect(time, "sham") ~ "sham",
    TRUE ~ time))

dats <- dats %>% 
  mutate(timeD = case_when(
    str_detect(time, "01") ~ '1',
    str_detect(time, "02") ~ '2',
    str_detect(time, "03") ~ '3',
    str_detect(time, "04") ~ '4',
    str_detect(time, "05") ~ '5',
    str_detect(time, "06") ~ '6',
    TRUE ~ time))

dats <- dats %>% select(-time)

df <- data.frame(dats)
library(afex)
# Преобразуем переменные в факторы
df$subName <- as.factor(df$subName)
df$condition <- as.factor(df$condition)
df$time <- as.factor(df$timeD)

# Выполняем Repeated Measures ANOVA с ковариатой
res <- aov_ez(
  id = "subName",
  dv = "LW",
  within = c("time", "condition"),
  data = df,
  type = 3
)

# Выводим результаты
print(res)


mod <- lm(LW ~ subNames+timeD+condition, data=dats)
library(lme4)
lmer(LW ~ subNames + condition + timeD, data = dats)