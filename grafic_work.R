################################################################################
#Plot the dynamics
###Define colors 
rhg_cols <- c("#F28E2B", "#76B7B2", "#EDC948", "#F27314", "#F8A31B", 
              "#E2C59F", "#B6C5CC", "#8E9CA3", "#556670", "#000000")
####

dats_groups %>% filter(met == 'Glx', condition == 'act') %>% group_by(subNames) %>% 
  mutate(MetValue_1 = MetValue[time == 1]) %>%
  mutate(Z = (MetValue-MetValue_1)/MetValue_1) %>% ungroup() %>% 
  ggplot(aes(x = time*2-2, y = Z, color = group)) +
  geom_jitter(size = 2, alpha  = 0.7, width = 0.1) +
  stat_summary(fun.y =mean, geom="line", size= 1) +
  stat_summary(fun.y =mean, geom="point", shape=20, size=5, color = "#8E9CA3") +
  stat_summary(fun.data =mean_se, geom="errorbar", width = 0.15, size = 1, alpha = 0.9, color = "#8E9CA3",
               linetype = "solid",position=position_dodge(width=0.3)) +
  scale_colour_manual(values =  rhg_cols)+
  xlab("Time, s") + ylab("change of [Glx], mM") +
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

##### FIND expected MRS BOLD CURVE
dats_hrf_mrs <- read.csv('C:/Users/Science/YandexDisk/Work/data/fMRS-hp/results/BOLD_MRS_hrf.csv')
subjectNames <- paste0("sub_", 1:32)
dats <- dats_hrf_mrs %>% mutate(subNames = subjectNames) 
#exclude subjects QA
#excludedSubjects <- c(1:6, 8, 9, 10, 11, 13, 17, 20, 22, 23,26, 28, 32 )
excludedSubjects <- c(1:6, 9, 11, 13, 26, 28, 22 ) # Due to the bad quality MRS
dats <- dats[-excludedSubjects,]

dats <- dats %>% pivot_longer(cols = "X01":"X06",
                              names_to = "time",
                              values_to = "BOLD")
dats$time = gsub('X','',dats$time)
dats$time <- as.numeric(dats$time)
dats <- dats[-(which(dats$time %in% 6)),]

#add column for grouping
groups_2 <- c(7, 24, 30, 10, 21, 25, 12, 27, 16, 18)
dats$group = 1
sub_groups_2 <- paste0("sub_", groups_2)  

additional_curve <- dats %>% 
  group_by(subNames) %>% mutate(BOLD_1 = BOLD[time == 1]) %>% ungroup %>%
  group_by(time) %>%
  mutate(Z = (BOLD-BOLD_1)/BOLD_1) %>%
  summarise(mean_Z = mean(Z, na.rm = TRUE))

additional_curve <- dats %>%
  group_by(subNames) %>% mutate(BOLD_1 = BOLD[time == 1]) %>% ungroup %>%
  group_by(time) %>%
  mutate(Z = (BOLD-BOLD_1)) %>%
  group_by(time) %>%
  summarise(mean_Z = mean(Z, na.rm = TRUE)/2.5)

# 3. Scale the additional data to match the primary y-axis range
scaled_data <- dats_groups %>% 
  filter(met == 'Glx', condition == 'act') %>% 
  group_by(subNames) %>% 
  mutate(MetValue_1 = MetValue[time == 1]) %>%
  mutate(Z = (MetValue-MetValue_1)/MetValue_1)

y1_range <- range(scaled_data$Z, na.rm = TRUE)
y2_range <- range(additional_curve$mean_Z, na.rm = TRUE)
scale_factor <- diff(y1_range) / diff(y2_range)

# Use that for relative change
# mutate(Z = (MetValue-MetValue_1)/MetValue_1) %>% 
  
dats_groups %>% 
  filter(met == 'Glx', condition == 'act') %>% 
  group_by(subNames) %>% 
  mutate(MetValue_1 = MetValue[time == 1]) %>%
  mutate(Z = (MetValue-MetValue_1)/MetValue_1) %>% 
  ungroup() %>% mutate(group = as.factor(group))  %>%
  ggplot(aes(x = time*2-2, y = Z)) +  # Added group aesthetic
  geom_jitter(aes(color = group), size = 2, alpha = 0.7, width = 0.1) +
  scale_colour_manual(values = c("#D55E00", "#0072B2")) +
  stat_summary(aes(color = group, group = group), fun = mean, geom = "line", size = 1) +
  stat_summary(aes(color = group, group = group), fun = mean, geom = "point", shape = 20, size = 5) +
  stat_summary(aes(color = group, group = group), fun.data = mean_se, geom = "errorbar", width = 0.15, size = 1, 
               alpha = 0.9, linetype = "solid") +
  xlab("Time, s") + ylab("change of [Glx], mM") +
  theme_minimal() + geom_line(data = additional_curve, 
            aes(x = time*2-2, y = mean_Z), 
            color = "#EDC948", size = 2, linetype = "dashed") +
  scale_y_continuous(
    sec.axis = sec_axis(
      ~ (. ) *2.5,  # Inverse scaling
      name = "Additional Metric (units)"  # Adjust label
    )
  )

