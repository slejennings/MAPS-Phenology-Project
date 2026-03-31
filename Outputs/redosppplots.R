# set up divergent palette
colorspace::diverging_hcl(n = 13, h = c(245, 125), c = c(30, 55), l = c(15, 95), power = c(0.75, 0.85), 
                          register = "custom_bluegreen" )
swatchplot(custom = diverging_hcl(n = 13, "custom_bluegreen" )) # examine palette



tstats_plot <- left_join(tstats, codetospec)

for(i in unique(tstats_plot$SPEC)){
  
  tstats_sub <- tstats_plot %>%
    filter(SPEC == i)
  
  name <- unique(tstats_sub$COMMONNAME)


##### Change in Breeding Phenology
FY <- tstats_sub %>%
  ggplot(aes(x = FY_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-5,5), name="t-statistic", rev=F, breaks=c(-4,-2,0,2,4)) + 
  geom_vline(aes(xintercept = 0), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in breeding phenology", 
       y = "Number of models",
       title = "Bird Breeding Phenology") +
  coord_cartesian(xlim=c(-6,6)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  theme(legend.key.size = unit(.4, 'cm'))

##### Change in Temp Anomalies - DW

TA_dw <- tstats_sub %>%
  ggplot(aes(x = tempanom_DW_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-6,6), name="t-statistic", rev=F, breaks=c(-6,-3,0,3,6)) + 
  geom_vline(aes(xintercept = 0), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in temperature anomalies", 
       y = "Number of models",
       title = "Temperature Anomaly: Decision Window") +
  coord_cartesian(xlim=c(-6,6)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  theme(legend.key.size = unit(.4, 'cm'))

TA_dw

##### Change in Temp Anomalies - LW

TA_lw <- tstats_sub %>%
  ggplot(aes(x = tempanom_LW_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-8,8), name="t-statistic", rev=F, breaks=c(-8,-4,0,4,8)) + 
  geom_vline(aes(xintercept = 0), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in temperature anomalies", 
       y = "Number of models",
       title = "Temperature Anomaly: Long-term Window") +
  coord_cartesian(xlim=c(-6,8)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  theme(legend.key.size = unit(.4, 'cm'))

##### Change in Cumulative (total) Precipitation - DW

TP_dw <- tstats_sub %>%
  ggplot(aes(x = prcp_DW_total_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-6,6), name="t-statistic", rev=F, breaks=c(-6,-3,0,3,6)) + 
  geom_vline(aes(xintercept = 0), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in cumulative precipitation", 
       y = "Number of models",
       title = "Cumulative Precipitation: Decision Window") +
  coord_cartesian(xlim=c(-6,6)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  theme(legend.key.size = unit(.4, 'cm'))

##### Change in Cumulative (total) Precipitation - LW

TP_lw <- tstats_sub %>%
  ggplot(aes(x = prcp_LW_total_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-5,5), name="t-statistic", rev=F, breaks=c(-4,-2,0,2,4)) + 
  geom_vline(aes(xintercept = 0), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in cumulative precipitation", 
       y = "Number of models",
       title = "Cumulative Precipitation: Long-term Window") +
  coord_cartesian(xlim=c(-6,6)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  theme(legend.key.size = unit(.4, 'cm'))

##### Change in Precipitation Variability - DW
PV_dw <- tstats_sub %>%
  ggplot(aes(x = prcp_DW_cov_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-6,6), name="t-statistic", rev=F, breaks=c(-6,-3,0,3,6)) + 
  geom_vline(aes(xintercept = 0), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in coefficient of variation \n for precipitation", 
       y = "Number of models",
       title = "Precipitation Variability: Decision Window") +
  coord_cartesian(xlim=c(-6,6)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  theme(legend.key.size = unit(.4, 'cm'))

##### Change in Precipitation Variability - LW

PV_lw <- tstats_sub %>%
  ggplot(aes(x = prcp_LW_cov_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-8,8), name="t-statistic", rev=F, breaks=c(-8,-4,0,4,8)) + 
  geom_vline(aes(xintercept = 0), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in coefficient of variation \n for precipitation", 
       y = "Number of models",
       title = "Precipitation Variability: Long-term window") +
  coord_cartesian(xlim=c(-8,8)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  theme(legend.key.size = unit(.4, 'cm'))

##### Change in Light

NL <- tstats_sub %>%
  ggplot(aes(x = light_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-18,18), name="t-statistic", rev=F, breaks=c(-14,-7,0,7,14)) + 
  geom_vline(aes(xintercept = 0), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in nighttime light", 
       y = "Number of models",
       title = "Light Pollution") +
  coord_cartesian(xlim=c(-18,18), ylim=c(0,8)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9))
  theme(legend.key.size = unit(.4, 'cm'))

# combine plots
# 1st row
R1 <-wrap_elements((FY + plot_spacer() + NL) + plot_layout(widths = c(0.48, 0.04, 0.48)) + 
                       plot_annotation(tag_levels = list(c('a', 'b'))) & theme(plot.tag = element_text(size = 14, face ="bold")))
# 2nd row
R2 <-wrap_elements((TA_dw + plot_spacer() + TA_lw) + plot_layout(widths = c(0.48, 0.04, 0.48)) + 
                       plot_annotation(tag_levels = list(c('c', 'd'))) & theme(plot.tag = element_text(size = 14, face ="bold")))
# third row
R3 <-wrap_elements((TP_dw + plot_spacer() + TP_lw) + plot_layout(widths = c(0.48, 0.04, 0.48)) + 
                      plot_annotation(tag_levels = list(c('e', 'f'))) & theme(plot.tag = element_text(size = 14, face ="bold")))
# fourth row
R4<-wrap_elements((PV_dw + plot_spacer() + PV_lw) + plot_layout(widths = c(0.48, 0.04, 0.48)) + 
                       plot_annotation(tag_levels = list(c('g', 'h'))) & theme(plot.tag = element_text(size = 14, face ="bold")))

spp_tstat_hist <- R1 / R2 / R3 / R4 +
    plot_annotation(title = name,
                    theme = theme(plot.title = element_text(size = 14, face="bold", hjust=0.5))) 

ggsave(spp_tstat_hist, filename = paste(i, "_changestat.png", sep = "_"), path = here("Figures/Spp change stats"), width=26, height=32, units = "cm")

}
