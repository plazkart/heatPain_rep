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

dats <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\BOLD_MRS.csv')
columnsActSham <- cbind(paste0("act", "_tp_0", 1:6, "_LWCr"), paste0("sham", "_tp_0", 1:6, "_LWCr"))
#dats <- dats[, columnsActSham]
subjectNames <- paste0("sub_", 1:32)
dats <- dats %>% mutate(subNames = subjectNames) 
row.names(dats) <- subjectNames

#exclusion and pivoting
excludedSubjects <- c(1:3, 11)
dats <- dats[-excludedSubjects,]

dats <- dats %>% pivot_longer(cols = "act_tp_01_LWCr":"sham_tp_06_HNAA",
                              names_to = "time",
                              values_to = "LWH")

dats <- dats %>% separate_wider_delim(time, delim = "_tp_", names = c("condition", "time"))
dats <- dats %>% separate_wider_delim(time, delim = "_", names = c("time", "Value"))


dats <- dats %>% 
  mutate(ValueLWH = case_when(
    str_detect(Value, "LW") ~ "LW",
    str_detect(Value, "H") ~ "H"))

dats <- dats %>% 
  mutate(Metabolite = case_when(
    str_detect(Value, "Cr") ~ "Cr",
    str_detect(Value, "NAA") ~ "NAA"))

dats <- dats %>% select(-'Value')

#3.1 QA check for outliers
# library
library(ggplot2)


# basic scatterplot
dats %>% filter(condition == 'act', ValueLWH == 'LW') %>% group_by(subNames, Metabolite) %>% mutate(LW_diff=LWH/mean(LWH)) %>%
ungroup() %>% group_by(time, Metabolite) %>% summarise(average_LW = mean(LW_diff, na.rm = TRUE)) %>%
ggplot(aes(x=time, y=average_LW, colour = Metabolite)) + 
  geom_point() +
  geom_line()


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