library(dplyr)
library(tidyverse)
library(tidyr)

subjectNames <- paste0("sub_", 1:32)
#excludedSubjects <- c(1:6, 8, 9, 11, 13,17, 20, 22,26, 28, 32 )
#excludedSubjects <- c(1:6, 9, 11, 13, 26, 28 ) # Only due to the bad quality MRS
excludedSubjects <- c(1:6, 9, 11, 26, 28, 22, 15 , 17 ) # Due to the bad quality MRS (aug 2025)

groups_1 <- c( 20, 8, 32, 19, 23, 29, 31, 14, 13)
groups_2 <- c(7, 24, 30, 10, 21, 25, 12, 27, 16, 18)
#get metabolic data
main_path = 'B:\\YandexDisk'
datsMet <- read.csv(paste(main_path,'\\Work\\data\\fMRS-hp\\results\\spectraAll.csv', sep = ""))
datsMet <- datsMet %>% mutate(subNames = subjectNames)
datsMet <- datsMet[-excludedSubjects,]
# get metabolic dynamic data
bold <- 0
 if (bold==0) {
  #get data
  datsInit <- read.csv(paste(main_path, '\\Work\\data\\fMRS-hp\\results\\spectraDynamic_sm.csv', sep = ""))
  } else {
  #bold-corrected dAta
  datsInit <- read.csv(paste(main_path, '\\Work\\data\\fMRS-hp\\results\\spectra-TP6-SM-BC.csv', sep = ""))
  }
datsMetDyn <- datsInit %>% mutate(subNames = subjectNames)
datsMetDyn <- datsMetDyn[-excludedSubjects,]

#add column for grouping
datsMetDyn$group = 1
sub_groups_2 <- paste0("sub_", groups_2)
datsMetDyn[(which(datsMetDyn$subNames %in% sub_groups_2)), 'group'] = 2

#get BOLD-fMRI data
datsBOLD <- read.csv(paste(main_path,'\\Work\\data\\fMRS-hp\\results\\fMRI_metrics.csv', sep = ""))
datsBOLD <- datsBOLD %>% mutate(subNames = subjectNames)
datsBOLD <- datsBOLD[-excludedSubjects,]


dats <- merge(datsMetDyn, datsBOLD, by = 'subNames')

## Correlation BOLD-metrics with Metabolites concentrations
dats_grouped <- dats %>% filter(group == 1)
cor(dats[,2:7], dats[,39:42])
## Results
#no correlation detected

dats_grouped <- dats %>% filter(group == 2)
## Correlation BOLD-metrics with Metabolites concentrations change
# Define the components to vary
metabolites <- c("Glx", "NAA", "Cr")
numbers <- sprintf("%02d", 1:5)  # Formats numbers as two-digit strings

# Generate all combinations
particular_columns <- unlist(lapply(metabolites, function(m) {
  paste0("sham_", m, "_tp_", numbers, "_sm")
}))


diff_dats <- dats_grouped %>% select(which(colnames(dats) %in% particular_columns))
diff_dats <-  diff_dats %>% mutate(diff_dats[, 1:5] - diff_dats[,1],
                                   diff_dats[, 6:10] - diff_dats[,6],
                                   diff_dats[, 11:15] - diff_dats[,11])
diff_dats <- diff_dats[, -c(1, 6, 11)]

correlation_results <- cor( dats_grouped[,39:42], diff_dats)

### Get values for statistics
dats %>% select(c(2:7)) %>% apply(2, mean)
dats %>% select(c(2:7)) %>% apply(2, sd)

