
# Make histograms of t-stats for breeding phenology, light and climate variables
# All species and stations pooled

# load packages
library(tidyverse)
library(here)
library(patchwork)
library(colorspace)

# import files 

# t-statistics for change in bird breeding phenology and environmental variables
tstats <- readRDS(here("Outputs", "combined_t_stats.rds"))

##### Change in Breeding Phenology

min(tstats$FY_tstat)
max(tstats$FY_tstat)

FYchangestat <- tstats %>%
  ggplot(aes(x = FY_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-6,6), name="t-statistic", rev=F) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
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
  annotate("segment", x=3, y=50, xend=5, yend=50, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-3, y=50, xend=-5, yend= 50,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 4, y=60, label = "Delaying", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -4, y=60, label = "Advancing", color = "#00294D", size=4)


FYchangestat


##### Change in Temp Anomalies - DW

min(tstats$tempanom_DW_tstat)
max(tstats$tempanom_DW_tstat)

TAchangestat_dw <- tstats %>%
  ggplot(aes(x = tempanom_DW_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-6,6), name="t-statistic", rev=F) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in temperature anomalies - decision", 
       y = "Number of models",
       title = "Temperature Anomaly - Decision") +
  coord_cartesian(xlim=c(-6,6)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=3, y=50, xend=5, yend=50, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-3, y=50, xend=-5, yend= 50,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 4, y=60, label = "Warming", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -4, y=60, label = "Cooling", color = "#00294D", size=4)

TAchangestat_dw

##### Change in Temp Anomalies - LW

min(tstats$tempanom_LW_tstat)
max(tstats$tempanom_LW_tstat)

TAchangestat_lw <- tstats %>%
  ggplot(aes(x = tempanom_LW_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-6,8), name="t-statistic", rev=F) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in temperature anomalies - long-term", 
       y = "Number of models",
       title = "Temperature Anomaly - Long-term") +
  coord_cartesian(xlim=c(-6,8)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=3, y=50, xend=5, yend=50, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-3, y=50, xend=-5, yend= 50,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 4, y=60, label = "Warming", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -4, y=60, label = "Cooling", color = "#00294D", size=4)

TAchangestat_lw

##### Change in Total Precipitation - DW

min(tstats$prcp_DW_total_tstat)
max(tstats$prcp_DW_total_tstat)

TPchangestat_dw <- tstats %>%
  ggplot(aes(x = prcp_DW_total_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-6,6), name="t-statistic", rev=F) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in total precipitation - decision", 
       y = "Number of models",
       title = "Total Precipitation - Decision") +
  coord_cartesian(xlim=c(-6,6)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=3, y=50, xend=5, yend=50, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-3, y=50, xend=-5, yend= 50,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 4, y=60, label = "Wetter", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -4, y=60, label = "Drier", color = "#00294D", size=4)

TPchangestat_dw

##### Change in Total Precipitation - LW

min(tstats$prcp_LW_total_tstat)
max(tstats$prcp_LW_total_tstat)

TPchangestat_lw <- tstats %>%
  ggplot(aes(x = prcp_LW_total_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-6,6), name="t-statistic", rev=F) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in total precipitation - long-term", 
       y = "Number of models",
       title = "Total Precipitation - Long-term") +
  coord_cartesian(xlim=c(-6,6)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=3, y=50, xend=5, yend=50, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-3, y=50, xend=-5, yend= 50,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 4, y=60, label = "Wetter", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -4, y=60, label = "Drier", color = "#00294D", size=4)

TPchangestat_lw


##### Change in Precipitation Variability - DW

min(tstats$prcp_DW_cov_tstat)
max(tstats$prcp_DW_cov_tstat)


PVchangestat_dw <- tstats %>%
  ggplot(aes(x = prcp_DW_cov_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-6,6), name="t-statistic", rev=F) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in precipitation variability - decision", 
       y = "Number of models",
       title = "CV of Precipitation - Decision") +
  coord_cartesian(xlim=c(-6,6)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=3, y=50, xend=5, yend=50, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-3, y=50, xend=-5, yend= 50,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 4, y=60, label = "More Variable", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -4, y=60, label = "Less Variable", color = "#00294D", size=4)

PVchangestat_dw

##### Change in Precipitation Variability - LW

min(tstats$prcp_LW_cov_tstat)
max(tstats$prcp_LW_cov_tstat)

PVchangestat_lw <- tstats %>%
  ggplot(aes(x = prcp_LW_cov_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-8,8), name="t-statistic", rev=F) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in precipitation variability - long-term", 
       y = "Number of models",
       title = "CV of Precipitation - Long-term") +
  coord_cartesian(xlim=c(-8,8)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=4, y=50, xend=6, yend=50, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-4, y=50, xend=-6, yend= 50,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 5, y=60, label = "More Variable", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -5, y=60, label = "Less Variable", color = "#00294D", size=4)

PVchangestat_lw


##### Change in Light

min(tstats$light_tstat)
max(tstats$light_tstat)

NLchangestat <- tstats %>%
  ggplot(aes(x = light_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-18,18), name="t-statistic", rev=F) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in nighttime light", 
       y = "Number of models",
       title = "Light Pollution") +
  coord_cartesian(xlim=c(-18,18)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=5, y=150, xend=10, yend=150, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-5, y=150, xend=-10, yend= 150,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 7, y=175, label = "Brighter", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -7, y=175, label = "Darker", color = "#00294D", size=4)

NLchangestat

tstat_hist <- (FYchangestat + NLchangestat) / (TAchangestat_dw + TAchangestat_lw) /
  (TPchangestat_dw + TPchangestat_lw) / (PVchangestat_dw + PVchangestat_lw)

ggsave(tstat_hist, filename = "TStatHistograms.pdf", path = here("Figures"), width=25, height=35, units = "cm")
