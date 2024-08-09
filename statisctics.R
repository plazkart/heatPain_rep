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