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