#get data

library(dplyr)
library(tidyverse)
library(tidyr)

batch_preprocessing <- function(bold, Cr, group_number, metName, conditonName) {
  if (bold==0) {
  #get data
  datsInit <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\spectraDynamic_sm.csv')
  } else {
  #bold-corrected dAta
  datsInit <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\spectra-TP6-SM-BC.csv')
  }
  
  subjectNames <- paste0("sub_", 1:32)
  dats <- datsInit %>% mutate(subNames = subjectNames) 
  #exclude subjects QA
  #excludedSubjects <- c(1:6, 8, 9, 10, 11, 13, 17, 20, 22, 23,26, 28, 32 )
  excludedSubjects <- c(1:6, 9, 11, 13, 26, 28, 22 ) # Due to the bad quality MRS
  dats <- dats[-excludedSubjects,]
  
  if (bold==0) {
  #pivot column
  dats <- dats %>% pivot_longer(cols = "act_NAA_tp_01_sm":"sham_Glx_tp_06_sm",
                                names_to = "time",
                                values_to = "MetValue") %>% select(subNames, time, MetValue)
  dats <- dats %>% separate_wider_delim(time, delim = "_tp_", names = c("condition", "time"))
  dats <- dats %>% separate_wider_delim(condition, delim = "_", names = c("condition", "met"))
  dats$time = gsub('_sm','',dats$time)
  dats$time <- as.numeric(dats$time) }
  else {
    dats <- dats %>% pivot_longer(cols = "act_NAA_tp_01_sm_bc":"act_Glx_tp_06_sm_bc",
                                  names_to = "time",
                                  values_to = "MetValue") %>% select(subNames, time, MetValue)
    dats <- dats %>% separate_wider_delim(time, delim = "_tp_", names = c("condition", "time"))
    dats <- dats %>% separate_wider_delim(condition, delim = "_", names = c("condition", "met"))
    dats$time = gsub('_sm_bc','',dats$time)
    dats$time <- as.numeric(dats$time)
  }
  
  # groups_1 <- c(17, 22, 20, 8, 32, 19, 23, 29, 31, 15)
  # groups_2 <- c(14, 7, 24, 30, 10, 21, 25, 12, 27, 16, 18)
  
  groups_1 <- c(17, 20, 8, 32, 19, 23, 29, 31, 15, 14)
  groups_2 <- c(7, 24, 30, 10, 21, 25, 12, 27, 16, 18)

  #exclude 6th data point
  dats <- dats[-(which(dats$time %in% 6)),]
  
  #add column for grouping
  dats$group = 1
  sub_groups_2 <- paste0("sub_", groups_2)  
  if (group_number==1) {
    sub_groups <- paste0("sub_", groups_1)
  } else {
    sub_groups <- paste0("sub_", groups_2)
  }
  
  if (Cr==1) {
    #Get normilised to Cr data
    dats_Cr <- dats %>% pivot_wider(names_from = met, values_from = MetValue) %>% 
      mutate(NAA_Cr = NAA/Cr,Glx_Cr = Glx/Cr ) %>% 
      select(subNames, condition, time, NAA_Cr, Glx_Cr) %>% 
      pivot_longer(cols = c("NAA_Cr", "Glx_Cr"), names_to = "met", values_to = "MetValue")
    # Creatine relative case
    dats_groups <- dats_Cr[(which(dats_Cr$subNames %in% sub_groups)),]
    metName <- paste0(metName, '_Cr')
    
  } else {
    # get needed group
    dats_groups <- dats[(which(dats$subNames %in% sub_groups)),]
  }
  dats[(which(dats$subNames %in% sub_groups_2)), 'group'] = 2
  return(dats)
  
}  

dats_groups <- batch_preprocessing(0, 0, 2, 'Glx', 'act')
dats_group <- dats_groups %>% filter(group == 2)
res <- batch_stattest(dats_group, 0, 0, 2, 'Glx', 'act')

batch_stattest <- function(dats_groups, bold, Cr, group_number, metName, conditonName) {
  datsWider <- dats_groups %>% filter(met == metName, condition == conditonName) %>% select(subNames, time, MetValue)
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
  mean_values <- dats_groups %>% filter(met == metName, condition == conditonName) %>% select(subNames, time, MetValue)  %>% 
    group_by(time) %>% summarise(mu = mean(MetValue))
  results$mean_values <- mean_values$mu[1]
  results$pvalue <- p.adjust(results$pvalue, method = 'BH')
  view(results)
  
  return(results)
}
################################################################################
# Test most activated vs weak activated persons


### T-test
# datsWider <- dats_groups %>% filter(met == "Glx", condition == "act") %>% select(subNames, time, MetValue)
# getTest <- function(df, x){
#   df <- df %>% filter(time == 1 | time == x)
#   t.test(MetValue ~ time, data = df, paired = TRUE)
# }
# #set 2:5 if there is 5 data points and lag = 1:4
# tests <- lapply(2:5, function(x) getTest(datsWider, x))
# results <- data.frame(
#   lag = 1:4, 
#   xsquared = sapply(tests, "[[", "statistic"), 
#   pvalue = sapply(tests, "[[", "p.value"),
#   meanDiference = sapply(tests, "[[", "estimate"),
#   stdDiff = sapply(tests, "[[", "stderr")
# )
# results <- results %>% select(pvalue, meanDiference, stdDiff)
# results$meanDiference <- results$meanDiference*-1
# 
# p.adjust(results$pvalue, method = 'BH')

## 0.195317291 0.006202118 0.060027264 0.149780352 In act GROUP 2!!!





# 
# dats_groups %>% filter(met == 'Glx', condition == 'act') %>% group_by(subNames) %>% 
#   +     mutate(MetValue_1 = MetValue[time == 1]) %>%
#   +     mutate(Z = (MetValue-MetValue_1)) %>% ungroup() %>% select(subNames, time, Z) %>% 
#   +     group_by(time) %>% summarise(mean = mean(Z), sd =std(Z))