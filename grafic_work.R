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
dats_hrf_mrs <- read.csv("C:\Users\Science\YandexDisk\Work\data\fMRS-hp\results\BOLD_MRS_hrf.csv")


dats_groups %>% 
  filter(met == 'Glx', condition == 'act') %>% 
  group_by(subNames) %>% 
  mutate(MetValue_1 = MetValue[time == 1]) %>%
  mutate(Z = (MetValue-MetValue_1)/MetValue_1) %>% 
  ungroup() %>% 
  ggplot(aes(x = time*2-2, y = Z, color = group, group = group)) +  # Added group aesthetic
  geom_jitter(size = 2, alpha = 0.7, width = 0.1) +
  stat_summary(fun = mean, geom = "line", size = 1) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 5) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.15, size = 1, 
               alpha = 0.9, linetype = "solid") +
  xlab("Time, s") + ylab("change of [Glx], mM") +
  theme_minimal()
