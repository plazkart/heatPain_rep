#get data

library(dplyr)
library(tidyverse)
library(tidyr)

pathMain <- 'C:\\Users\\Admin\\YandexDisk\\Work\\data\\fMRS-hp\\results\\'
# 1. Search differences between different time points in dynamic
#get data
datsInit <- read.csv(paste(pathMain, 'spectraDynamic_sm.csv', sep = ''))
#bold-corrected dAta
datsInit <- read.csv(paste(pathMain, 'spectra-TP6-SM-BC.csv', sep = ''))

subjectNames <- paste0("sub_", 1:32)
dats <- datsInit %>% mutate(subNames = subjectNames) 
#exclude subjects QA
excludedSubjects <- c(1:6, 8, 9, 10, 11, 13, 17, 20, 22, 23,26, 28, 32 )
#excluded subjects for sham condition
excludedSubjects <- c(1:3, 24, 26, 28 )

dats <- dats[-excludedSubjects,]
#pivot column
dats <- dats %>% pivot_longer(cols = "act_NAA_tp_01_sm":"sham_Glx_tp_06_sm",
                              names_to = "time",
                              values_to = "MetValue") %>% select(subNames, time, MetValue)
dats <- dats %>% separate_wider_delim(time, delim = "_tp_", names = c("condition", "time"))
dats <- dats %>% separate_wider_delim(condition, delim = "_", names = c("condition", "met"))
dats$time = gsub('_sm','',dats$time)
dats$time <- as.numeric(dats$time)

#for BOLD-corrected DATA
dats <- dats %>% pivot_longer(cols = "act_NAA_tp_01_sm_bc":"act_Glx_tp_06_sm_bc",
                              names_to = "time",
                              values_to = "MetValue") %>% select(subNames, time, MetValue)
dats <- dats %>% separate_wider_delim(time, delim = "_tp_", names = c("condition", "time"))
dats <- dats %>% separate_wider_delim(condition, delim = "_", names = c("condition", "met"))
dats$time = gsub('_sm_bc','',dats$time)
dats$time <- as.numeric(dats$time)


#exclude 6th data point
dats <- dats[-(which(dats$time %in% 6)),]

#Get normilised to Cr data
dats_Cr <- dats %>% pivot_wider(names_from = met, values_from = MetValue) %>% 
  mutate(NAA_Cr = NAA/Cr,Glx_Cr = Glx/Cr ) %>% 
  select(subNames, condition, time, NAA_Cr, Glx_Cr) %>% 
  pivot_longer(cols = c("NAA_Cr", "Glx_Cr"), names_to = "met", values_to = "MetValue")

################################################################################
## Statistical analysis
#check for normality
normalityStatistics <- dats %>% filter(met == "Cr", condition == "act") %>%
  select(time, MetValue) %>%
  group_by(time) %>%
  summarise_all(.funs = funs(statistic = t.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
normalityStatistics
normalityStatistics <- dats %>% filter(met == "NAA", condition == "sham") %>%
  select(time, MetValue) %>%
  group_by(time) %>%
  summarise_all(.funs = funs(statistic = t.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
normalityStatistics
# For Glx both conditons are normally distributed for each time point
# For Cr both statistics are normally distributed for each time point
# For NAA statistics in act condition are not normally distributed for a few points

normalityStatistics <- dats_Cr %>% filter(met == "NAA_Cr", condition == "act") %>%
  select(time, MetValue) %>%
  group_by(time) %>%
  summarise_all(.funs = funs(statistic = t.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
normalityStatistics$p.value
normalityStatistics <- dats_Cr %>% filter(met == "NAA_Cr", condition == "sham") %>%
  select(time, MetValue) %>%
  group_by(time) %>%
  summarise_all(.funs = funs(statistic = t.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
normalityStatistics$p.value

# For Glx/Cr only Act group is normally distributed for each time point
# For NAA/Cr only Sham group is normally distributed for each time point

### T-test
datsWider <- dats %>% filter(met == "NAA", condition == "act") %>% select(subNames, time, MetValue)
results <- pairwise.t.test(datsWider$MetValue, datsWider$time,
                           p.adjust.method = 'none', paired = TRUE, FUN = mean)
p.adjust(results$p.value[,1], method = 'BH')

getTest <- function(df, x){
  df <- df %>% filter(time == 1, time == x)
  pairwise.t.test(MetValue, time, data = df, p. adjust.methods = 'none',  paired = TRUE)
}
#set 2:5 if there is 5 data points and lag = 1:4
tests <- lapply(2:5, function(x) getTest(datsWider, x))
results <- data.frame(
  lag = 1:4, 
  xsquared = sapply(tests, "[[", "statistic"), 
  pvalue = sapply(tests, "[[", "p.value"),
  meanDiference = sapply(tests, "[[", "estimate"),
  stdDiff = sapply(tests, "[[", "stderr")
)
results <- results %>% select(pvalue, meanDiference, stdDiff)
results$meanDiference <- results$meanDiference*-1

p.adjust(results$pvalue, method = 'BH')

#For Cr normilised
datsWider <- dats_Cr %>% filter(met == "NAA_Cr", condition == "act") %>% select(subNames, time, MetValue)
getTest <- function(df, x){
  df <- df %>% filter(time == 1 | time == x)
  t.test(MetValue ~ time, data = df, paired = TRUE)
}
#set 2:5 if there is 5 data points and lag = 1:4
tests <- lapply(2:5, function(x) getTest(datsWider, x))
results <- data.frame(
  lag = 1:4, 
  xsquared = sapply(tests, "[[", "statistic"), 
  pvalue = sapply(tests, "[[", "p.value"),
  meanDiference = sapply(tests, "[[", "estimate"),
  stdDiff = sapply(tests, "[[", "stderr")
)
results <- results %>% select(pvalue, meanDiference, stdDiff)
results$meanDiference <- results$meanDiference*-1
view(results)
p.adjust(results$pvalue, method = 'BY')

# Second point for Glx/Cr has trend to increase p = 0.06
# Third point for NAA/Cr increasing p = 0.045

## Glx differs from 1st point in act: 2nd point 0.02
## Glx differs from 1st point in sham: no differences

#################### repeated measures two-way ANOVA 
# Statistics

df <- data.frame(dats %>% filter(met == 'Glx') %>% select(subNames, condition, time, MetValue))
library(afex)
# Преобразуем переменные в факторы
df$subNames <- as.factor(df$subNames)
df$condition <- as.factor(df$condition)
df$time <- as.factor(df$time)

# Выполняем Repeated Measures ANOVA с ковариатой
res <- aov_ez(
  id = "subNames",
  dv = "MetValue",
  within = c("time", "condition"),
  data = df,
  type = 3
)

# Выводим результаты
print(res)

##################################
#Get values for a TABLE data
outTableData <- dats %>% filter(met == 'Glx') %>% group_by(condition, time) %>% 
  summarise(meanValue = mean(ValueLWH), stdZ = sd(ValueLWH))
outTableData <- dats_Cr %>% filter(met == 'Glx_Cr') %>% group_by(condition, time) %>% 
  summarise(meanValue = mean(MetValue), stdZ = sd(MetValue))

################################################################################
#Plot the dynamics
###Define colors 
rhg_cols <- c("#F28E2B", "#76B7B2", "#EDC948", "#F27314", "#F8A31B", 
              "#E2C59F", "#B6C5CC", "#8E9CA3", "#556670", "#000000")
####
metName <- 'NAA'
dats %>% filter(met == metName) %>% group_by(subNames, condition) %>%
  mutate(MetValue_1 = MetValue[time == 1]) %>%
  mutate(Z = (MetValue-MetValue_1)) %>% ungroup() %>%
  ggplot(aes(x = time*2-2, y = Z)) +
  geom_jitter(size = 2,  color = "#F28E2B", alpha  = 0.7, width = 0.1) +
  stat_summary(fun.y =mean, geom="line", size= 1, color = "#8E9CA3") +
  stat_summary(fun.y =mean, geom="point", shape=20, size=5, color = "#8E9CA3") +
  stat_summary(fun.data =mean_se, geom="errorbar", width = 0.15, size = 1.0, alpha = 0.9, color = "#8E9CA3",
               linetype = "solid",position=position_dodge(width=0.3)) +
  scale_colour_manual(values =  rhg_cols)+
  facet_grid(vars(condition))+
  xlab('time, s') + # for the x axis label
  ylab(paste('change', metName,', mM')) + # for the y axis label
  theme_minimal()


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


####experimental
dats %>% filter(met == 'Glx', condition == 'act') %>% group_by(time) %>% 
  mutate(Z = (MetValue-mean(MetValue))/sd(MetValue)) %>% ungroup() %>%
  ggplot(aes(x = time, y = Z)) +
  geom_jitter(size = 2,  color = "#F28E2B", alpha  = 0.7, width = 0.1) +
  stat_summary(fun.y =mean, geom="line", size= 1, color = "#8E9CA3") +
  stat_summary(fun.y =mean, geom="point", shape=20, size=5, color = "#8E9CA3") +
  stat_summary(fun.data =mean_se, geom="errorbar", width = 0.3, size = 1.5, alpha = 0.9, color = "#8E9CA3",
               linetype = "solid",position=position_dodge(width=0.3)) +
  scale_colour_manual(values =  rhg_cols)+
  theme_minimal()