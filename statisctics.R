# This scripts are dedicated for statistic quantification of 
# Glx dynamic during short heat pain stimulation.
library(dplyr)
library(tidyverse)
library(tidyr)

# 1. Search differences between different time points in dynamic
#get data
dats <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\spectraDynamic_sm.csv')
#bold-corrected dAta
dats <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\spectra-TP6-SM-BC.csv')
# dynamics dAta in ACT state
dats <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\spectra-TP6-SM-BC.csv')

subjectNames <- paste0("sub_", 1:32)
dats <- dats %>% mutate(subNames = subjectNames) 
#exzclude subjects
excludedSubjects <- c(1:3, 11, 4, 5, 6, 8, 19, 20, 32, 22 )
#excludeaccording to QA
excludedSubjects <- c(1:3, 9, 11, 13, 26, 28 )

dats <- dats[-excludedSubjects,]
#pivot column
dats <- dats %>% pivot_longer(cols = "act_NAA_tp_01_sm_bc":"act_Glx_tp_06_sm_bc",
                              names_to = "time",
                              values_to = "MetValue") %>% select(subNames, time, MetValue)
dats <- dats %>% separate_wider_delim(time, delim = "_tp_", names = c("condition", "time"))
dats <- dats %>% separate_wider_delim(condition, delim = "_", names = c("condition", "met"))
dats$time = gsub('_sm_bc','',dats$time)
dats$time <- as.numeric(dats$time)
#Plot the dynamics
dats %>% filter(met == 'Glx') %>% group_by(subNames, condition) %>%
   mutate(Z = (MetValue-mean(MetValue))/sqrt(var(MetValue))) %>% ungroup() %>%
  group_by(condition, time) %>% summarise(meanZ = mean(Z)) %>%
  ggplot(aes(x = time, y = meanZ, color = condition)) +
  geom_point(size = 3) +
  geom_line() +
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

###########################
#Dynamics un ACT state
###########################

# dynamics dAta in ACT state
dats <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\spectra-TP6.csv')

subjectNames <- paste0("sub_", 1:32)
dats <- dats %>% mutate(subNames = subjectNames) 
#excludeaccording to QA
excludedSubjects <- c(1:3, 9, 11, 13, 26, 28 )
dats <- dats[-excludedSubjects,]
# excluson according to low fMRI activation
excludedSubjects <- c(4, 5, 6, 8, 19, 20, 32, 17, 22)
dats <- dats[-excludedSubjects,]
#pivot column
dats <- dats %>% pivot_longer(cols = "act_NAA_tp_01_sm_bc":"act_Glx_tp_06_sm",
                              names_to = "time",
                              values_to = "MetValue") %>% select(subNames, time, MetValue)
dats <- dats %>% separate_wider_delim(time, delim = "_tp_", names = c("condition", "time"))
dats <- dats %>% separate_wider_delim(condition, delim = "_", names = c("condition", "met"))
dats <- dats %>% separate_wider_delim(time, delim = "_s", names = c("time", "boldCorrected"))
dats$time <- as.numeric(dats$time)

dats <- dats %>% 
  mutate(bc = case_when(
    str_detect(boldCorrected, "m_bc") ~ "1",
    str_detect(boldCorrected, "m") ~ "0"))
dats <- dats %>% select(-c(condition, boldCorrected))
dats$bc <- as.numeric(dats$bc)

#Plot the mean or median dynamics
dats %>% filter(met == 'Glx') %>% group_by(subNames, bc) %>%
  mutate(Z = (MetValue-mean(MetValue))/sqrt(var(MetValue))) %>% ungroup() %>%
  group_by(bc, time) %>% summarise(meanZ = median(Z)) %>%
  ggplot(aes(x = time, y = meanZ, color = bc)) +
  geom_point(size = 3) +
  geom_line() +
  facet_wrap(~ bc)

#Plot the every point in dynamics for each subject
dats %>% filter(met == 'Glx') %>% group_by(subNames, bc) %>%
  mutate(Z = (MetValue-mean(MetValue))/sqrt(var(MetValue))) %>% ungroup() %>%
  ggplot(aes(x = time, y = Z, color = subNames)) +
  geom_point(size = 3) +
  facet_wrap(~ bc)

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

## Statisical analisys
#check for normality
normalityStatistics <- dats %>% select(met, time, bc, MetValue) %>%
  filter(met == "Glx") %>%
  group_by(met, time, bc) %>%
  summarise_all(.funs = funs(statistic = t.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))

datsWider <- dats %>% filter(met == "Glx", bc == 1) %>% select(subNames, time, MetValue) %>%
  filter(time == 1 | time== 2)
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
#find the difference
###################
# END of the part
###################

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
#for temporally smoothed data
columnsActSham <- cbind(paste0("act", "_tp_0", 1:6, "_sm_LWCr"), paste0("sham", "_tp_0", 1:6, "_sm_HNAA"))
#dats <- dats[, columnsActSham]
subjectNames <- paste0("sub_", 1:32)
dats <- dats %>% mutate(subNames = subjectNames) 
row.names(dats) <- subjectNames

#exclusion and pivoting
#excludedSubjects <- c(1:3, 11)
excludedSubjects <- c(1:3, 11, 4, 5, 6, 8, 19, 20, 32, 22 )
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

dats <- dats[-(which(dats$time %in% "06")),]

#3.1 QA check for outliers
# library
library(ggplot2)


# basic scatterplot
plotDats <- dats %>% filter(ValueLWH == 'LW') %>% 
  group_by(subNames, condition, Metabolite) %>% 
  mutate(LW_diff=LWH/mean(LWH)) %>%
  ungroup() %>% group_by(time, Metabolite, condition) %>% 
  summarise(average_LW = mean(LW_diff, na.rm = TRUE))
plotDats %>% pivot_wider(names_from = Metabolite, values_from = average_LW) %>%
  ggplot(aes(x=time)) + 
  geom_point(aes(y=Cr, colour = 'coral'), size = 3) +
  geom_line(aes(y=Cr, colour = 'coral', group = 1)) +
  geom_point(aes(y=NAA), colour = 'darkcyan', size = 3) +
  geom_line(aes(y=NAA,colour = 'darkcyan', group = 1)) +
  scale_fill_discrete(name = "Value", labels = c("LW", "H")) +
  facet_wrap(~ condition)



plotDats <- dats %>% filter(Metabolite == 'Cr') %>% 
  group_by(subNames, condition, ValueLWH) %>% 
  mutate(LW_diff=LWH/mean(LWH)) %>%
ungroup() %>% group_by(time, ValueLWH, condition) %>% 
  summarise(average_LW = mean(LW_diff, na.rm = TRUE)) 

plotDats %>% pivot_wider(names_from = ValueLWH, values_from = average_LW) %>%
  ggplot(aes(x=time)) + 
  geom_point(aes(y=LW, colour = 'coral'), size = 3) +
  geom_line(aes(y=LW, colour = 'coral', group = 1)) +
  geom_point(aes(y=H), colour = 'darkcyan', size = 3) +
  geom_line(aes(y=H,colour = 'darkcyan', group = 1)) +
  scale_fill_discrete(name = "Value", labels = c("LW", "H")) +
  facet_wrap(~ condition)

plotDats %>% pivot_wider(names_from = Metabolite, values_from = average_LW) %>%
  ggplot(aes(x=time)) + 
  geom_point(aes(y=Cr, colour = 'coral'), size = 3) +
  geom_line(aes(y=Cr, colour = 'coral', group = 1)) +
  geom_point(aes(y=NAA), colour = 'darkcyan', size = 3) +
  geom_line(aes(y=NAA,colour = 'darkcyan', group = 1)) +
  facet_wrap(~ condition)

plotDats <- dats %>% filter(ValueLWH == 'H') %>% group_by(subNames, Metabolite, condition) %>% mutate(LW_diff=LWH/mean(LWH)) %>%
  ungroup() %>% group_by(time, Metabolite, condition) %>% summarise(average_LW = mean(LW_diff, na.rm = TRUE)) 

plotDats %>% pivot_wider(names_from = Metabolite, values_from = average_LW) %>%
  ggplot(aes(x=time)) + 
  geom_point(aes(y=Cr, colour = 'coral'), size = 3) +
  geom_line(aes(y=Cr, colour = 'coral', group = 1)) +
  geom_point(aes(y=NAA), colour = 'darkcyan', size = 3) +
  geom_line(aes(y=NAA,colour = 'darkcyan', group = 1)) +
  facet_wrap(~ condition)  

plotDats <- dats %>% filter(ValueLWH == 'LW', Metabolite == 'NAA') %>%
  select(subNames, time, LWH, condition)
plotDats$time <- as.numeric(factor(plotDats$time))
plotDats <- plotDats[-(which(plotDats$subNames %in% "sub_28")),]
plotDats <- plotDats[-(which(plotDats$subNames %in% "sub_9")),]

plotDats %>% group_by(condition, subNames) %>% mutate(LWH = LWH-mean(LWH)) %>%
  ggplot(aes(x = time, y = LWH, color = subNames)) +
  geom_line() +
  facet_wrap(~condition)

#3.2 Comparison of real MRS linewidthes and theoretically dervied data
#get LW and H
dats <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\BOLD_MRS.csv')
#get theoretical data
dats2 <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\BOLD_MRS_hrf.csv')
subjectNames <- paste0("sub_", 1:32)
dats2 <- dats2 %>% mutate(subNames = subjectNames) 

excludedSubjects <- c(1:3, 11, 4, 5, 6, 8, 19, 20, 32, 22 )

dats2 <- dats2[-excludedSubjects,]
dats2 <- dats2 %>% select(-"X06")

dats2 <- dats2 %>% pivot_longer(cols = "X01":"X05",
                                names_to = "time",
                                values_to = "LWH")

dats2$time <-gsub("X","",as.character(dats2$time))

dats <- dats %>% mutate(subNames = subjectNames) 
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

dats <- dats[-(which(dats$time %in% "06")),]



#Quantification and figure analysis
dats <- dats %>% filter(condition=="sham") %>% 
  group_by(subNames, ValueLWH, Metabolite) %>% mutate(normilisedLWH = LWH/mean(LWH)-1) %>%
  ungroup() 

#dats$normilisedLWH[dats$ValueLWH=='H'] <-  dats$normilisedLWH[dats$ValueLWH=='H']*-1

mergedDats <- merge(dats, dats2, by = c("time", "subNames"))
# load ggplot2
library(ggplot2)

# A basic scatterplot with color depending on Species
mergedDats %>% 
  ggplot(aes(x=normilisedLWH, y = LWH.y)) + 
  geom_point() 

#Compare three groups
mergedDats %>% 
  mutate(actLevel = case_when(
    LWH.y>0.4 ~ 3,
    (LWH.y>0.1 && LWH.y<0.4) ~ 2,
    LWH.y<0.1 ~ 1,
    TRUE ~ 2)) %>%
  group_by(actLevel) %>% summarise(mean(normilisedLWH))


# Statistics

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