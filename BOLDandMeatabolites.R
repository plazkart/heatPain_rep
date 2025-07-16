library(dplyr)
library(tidyverse)
library(tidyr)

subjectNames <- paste0("sub_", 1:32)
#excludedSubjects <- c(1:6, 8, 9, 11, 13,17, 20, 22,26, 28, 32 )
excludedSubjects <- c(1:6, 9, 11, 13, 26, 28 ) # Only due to the bad quality MRS
#get metabolic data
datsMet <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\spectraAll.csv')
datsMet <- datsMet %>% mutate(subNames = subjectNames)
datsMet <- datsMet[-excludedSubjects,]
# get metabolic dynamic data
#get data
datsInit <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\spectraDynamic_sm.csv')
#bold-corrected dAta
datsInit <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\spectra-TP6-SM-BC.csv')
datsMetDyn <- datsInit %>% mutate(subNames = subjectNames)
#get BOLD-fMRI data
datsBOLD <- read.csv('C:\\Users\\Science\\YandexDisk\\Work\\data\\fMRS-hp\\results\\fMRI_metrics.csv')
datsBOLD <- datsBOLD %>% mutate(subNames = subjectNames)
datsBOLD <- datsBOLD[-excludedSubjects,]

dats <- merge(datsMetDyn, datsBOLD, by = 'subNames')

## Correlation BOLD-metrics with Metabolites concentrations 
cor(dats[,2:7], dats[,24:27])
## Results
#no correlation detected

## Correlation BOLD-metrics with Metabolites concentrations change
# Define the components to vary
metabolites <- c("Glx", "NAA", "Cr")
numbers <- sprintf("%02d", 1:5)  # Formats numbers as two-digit strings

# Generate all combinations
particular_columns <- unlist(lapply(metabolites, function(m) {
  paste0("act_", m, "_tp_", numbers, "_sm")
}))


diff_dats <- dats %>% select(which(colnames(dats) %in% particular_columns))
diff_dats <-  diff_dats %>% mutate(diff_dats[, 1:5] - diff_dats[,1],
                                   diff_dats[, 6:10] - diff_dats[,6],
                                   diff_dats[, 11:15] - diff_dats[,11])
diff_dats <- diff_dats[, -c(1, 6, 11)]

cor( dats[,38:41], diff_dats)

### Get values for statistics
dats %>% select(c(2:7)) %>% apply(2, mean)
dats %>% select(c(2:7)) %>% apply(2, sd)


