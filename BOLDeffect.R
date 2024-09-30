
library(dplyr)
library(tidyverse)
library(tidyr)
#get LW and H
dats <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\BOLD_MRS.csv')

subjectNames <- paste0("sub_", 1:32)
dats <- dats %>% mutate(subNames = subjectNames) 

excludedSubjects <- c(1:6, 8, 9, 10, 11, 13,17, 20, 22, 23,26, 28, 32 )
#excluded subjects for sham condition
excludedSubjects <- c(1:3, 24, 26, 28 )

dats <- dats[-excludedSubjects,]
dats <- dats %>% pivot_longer(cols = "act_tp_01_sm_LWCr":"sham_tp_06_sm_HNAA",
                              names_to = "time",
                              values_to = "LWH")

dats <- dats %>% select(subNames, time, LWH )
dats <- dats %>% separate_wider_delim(time, delim = "_tp_", names = c("condition", "time"))
dats <- dats %>% separate_wider_delim(time, delim = "_sm_", names = c("time", "Value"))

dats$time <- as.numeric(dats$time)

dats <- dats %>% 
  mutate(ValueLWH = case_when(
    str_detect(Value, "LW") ~ "LW",
    str_detect(Value, "H") ~ "H"))

dats <- dats %>% 
  mutate(met = case_when(
    str_detect(Value, "Cr") ~ "Cr",
    str_detect(Value, "NAA") ~ "NAA"))

dats <- dats %>% select(-'Value')
dats <- dats[-(which(dats$time %in% 6)),]


################################################################################
### Statistical tests
normalityStatistics <- dats %>% filter(met == "NAA", ValueLWH == "LW", condition == "sham") %>%
  select(time, LWH) %>%
  group_by(time) %>%
  summarise_all(.funs = funs(statistic = t.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
normalityStatistics$p.value
normalityStatistics <- dats %>% filter(met == "NAA", ValueLWH == "H") %>%
  select(time, LWH, condition) %>%
  group_by(time, condition) %>%
  summarise_all(.funs = funs(statistic = t.test(.)$statistic, 
                             p.value = shapiro.test(.)$p.value))
normalityStatistics$p.value

#### T-test
datsWider <- dats %>% filter(met == "Cr", condition == "sham",  ValueLWH == "LW") %>% 
  select(subNames, time, LWH)
getTest <- function(df, x){
  df <- df %>% filter(time == 1 | time == x)
  t.test(LWH ~ time, data = df, paired = TRUE)
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

p.adjust(results$pvalue)

##################################
#Get values for a TABLE data
outTableData <- dats %>% filter(met == 'Cr', ValueLWH == "LW") %>% group_by(condition, time) %>% 
  summarise(meanValue = mean(LWH), stdZ = sd(LWH))
outTableData <- dats_Cr %>% filter(met == 'Glx_Cr') %>% group_by(condition, time) %>% 
  summarise(meanValue = mean(MetValue), stdZ = sd(MetValue))



######################################################
#Plot the data


dats %>% filter(condition=="sham", ValueLWH == "LW", met == "Cr") %>% 
  ggplot(aes(LWH)) +
  geom_histogram()


###############################
# Analysis of dependence of concentration and Linewidth changes

excludedSubjects <- c(1:6, 8, 9, 10, 11, 13,17, 20, 22, 23,26, 28, 32 )
subjectNames <- paste0("sub_", 1:32)
#get LW and H
datsLWH <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\BOLD_MRS_sm.csv')
datsLWH <- datsLWH %>% mutate(subNames = subjectNames)
datsLWH <- datsLWH[-excludedSubjects,]
#get data
datsMet1 <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\spectraDynamic_sm.csv')
datsMet1 <- datsMet1 %>% mutate(subNames = subjectNames)
datsMet1 <- datsMet1[-excludedSubjects,]
#bold-corrected dAta
datsMet2 <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\spectra-TP6-SM-BC.csv')
datsMet2 <- datsMet2 %>% mutate(subNames = subjectNames)
datsMet2 <- datsMet2[-excludedSubjects,]


datsMet <- merge(datsMet1, datsMet2, by = "subNames")
dats <- merge(datsLWH, datsMet, by = "subNames")

dats <- dats %>% pivot_longer(cols = "act_tp_01_sm_LWCr":"act_Glx_tp_06_sm_bc",
                              names_to = "time",
                              values_to = "value")
dats <-  dats %>% select(subNames, time, value)

dats <- dats %>% 
  mutate(valueType = case_when(
    str_detect(time, "LW") ~ "LW",
    str_detect(time, "H") ~ "H",
    str_detect(time, "NAA_") ~ "conc",
    str_detect(time, "Glx_") ~ "conc", 
    str_detect(time, "Cr_") ~ "conc"))

dats <- dats[complete.cases(dats), ]

dats <- dats %>% 
  mutate(met = case_when(
    str_detect(time, "Cr") ~ "Cr",
    str_detect(time, "NAA") ~ "NAA",
    str_detect(time, "Glx") ~ "Glx"))

dats <- dats %>% 
  mutate(t = case_when(
    str_detect(time, "01") ~ "1",
    str_detect(time, "02") ~ "2",
    str_detect(time, "03") ~ "3",
    str_detect(time, "04") ~ "4",
    str_detect(time, "05") ~ "5",
    str_detect(time, "06") ~ "6"))

dats <- dats %>% 
  mutate(bc = case_when(
    str_detect(time, "bc") ~ 1,
    .default = 0
    ))

dats <- dats %>% 
  mutate(cond = case_when(
    str_detect(time, "act") ~ "act",
    str_detect(time, "sham") ~ "sham"
  ))
dats <- dats %>%  select(-time)

### get LW change
tempDats0 <- dats %>% filter(valueType == "LW", cond == "act", met == "NAA") %>% group_by(subNames) %>%
  mutate(LWchange = value - mean(value)) %>% select(subNames, LWchange)

### get conc change
tempDats <- dats %>% filter(valueType == "conc", cond == "act", met == "Glx") %>%
  pivot_wider(names_from = bc, values_from = value) %>% mutate(GlxDiff = )
colnames(tempDats) <- c("subNames",  "valueType", "met", "t",    "cond", "bc0", "bc1")
tempDats <- tempDats %>% mutate(metChange = bc1 - bc0) %>% select(subNames, metChange)

tempDats$LWChange <- tempDats0$LWchange

### Plot the dependence
tempDats %>% ggplot(aes(x = metChange, y = LWChange)) +
  geom_point()
