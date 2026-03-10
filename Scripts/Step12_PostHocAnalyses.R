###### MAPS Phenology #######
### Script name: Step12_PosthocAnalyses.R
## Author(s): CDF, SLJ

########## Objective/Description of Script #####################
# run posthoc analysis to compare mean change in temperature experienced by each species with their STI
# we examine this relationship using change in temp for both the decision and long-term window
#################################################################

#### Setup ####

# load packages
library(tidyverse)
library(here)
library(ggeffects)
library(nlme)
library(performance)
library(broom.mixed)

# import files

# t-statistics for change in bird breeding phenology and environmental variables
tstats <- readRDS(here("Outputs", "combined_t_stats.rds"))

# species temperature indices
STI <- read.csv(here("Data", "species_STI.csv"), header=T)

################################################################################
### Organize data and run post-hoc models
################################################################################

# calculate the mean t-statistic (and standard error) for temperature anomaly change across all stations in which each species occurs
# do this for both the decision and the long-term window
combT_summary <- tstats %>%
  group_by(SPEC) %>%
  summarise(
    n_decision = sum(!is.na(tempanom_DW_tstat)),
    mean_decision = mean(tempanom_DW_tstat, na.rm = TRUE),
    se_decision = sd(tempanom_DW_tstat, na.rm = TRUE) / sqrt(n_decision),
    n_long = sum(!is.na(tempanom_LW_tstat)),
    mean_long = mean(tempanom_LW_tstat, na.rm = TRUE),
    se_long = sd(tempanom_LW_tstat, na.rm = TRUE) / sqrt(n_long),
    .groups = "drop"
  )

# combine these data with STI for species
combDAT <- left_join(combT_summary, STI, by ="SPEC")

# examine the relationship between STI and the mean change in temperature anomalies for each species 
# first run the decision model with the SE used as a weighting function
m1 <- gls(mean_decision ~ STI, data = combDAT, weights = varFixed(~1/sqrt(se_decision)), method = "ML")

summary(m1)
confint(m1)
check_model(m1)

# next run the long-term model with SE as a weighting function
m2 <- gls(mean_long ~ STI, data = combDAT, weights = varFixed(~1/sqrt(se_long)), method = "ML")

summary(m2)
confint(m2)
check_model(m2)

################################################################################
### Make plots
################################################################################

# get the predicted effect from the decision window model
eff_STI_dw <- plot(ggeffects::predict_response(m1, terms =c("STI")), colors = "darkseagreen")

# add data, labels, nice formatting to plot
STI_dw_plot <- eff_STI_dw +
  geom_point(data = combDAT, aes(x = STI, y = mean_decision), color = "darkseagreen", size = 3.1, pch = 19) +
  geom_errorbar(aes(ymin = (mean_decision-se_decision), ymax = (mean_decision+se_decision), x=STI),
                color="darkseagreen", data = combDAT,inherit.aes = FALSE, lwd=.9) +
  labs(title= "", x = "Species temperature index (STI)", y = "Mean temperature change statistic")+
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), 
        axis.title.x = element_text(color = "black", size = 12), axis.title.y = element_text(color = "black", size = 12)) 
  

STI_dw_plot


# get the predicted effect from the decision window model
eff_STI_lw <- plot(ggeffects::predict_response(m2, terms =c("STI")), colors = "#1768B3")

# add data, labels, nice formatting to plot
STI_lw_plot <- eff_STI_lw +
  geom_point(data = combDAT, aes(x = STI, y = mean_long), color = "#1768B3", size = 3.1, pch = 19) +
  geom_errorbar(aes(ymin = (mean_long-se_long), ymax = (mean_long+se_long), x=STI),
                color="#1768B3", data = combDAT,inherit.aes = FALSE, lwd=.9) +
  labs(title= "", x = "Species temperature index (STI)", y = "Mean temperature change statistic")+
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), 
        axis.title.x = element_text(color = "black", size = 12), axis.title.y = element_text(color = "black", size = 12)) 


STI_lw_plot

posthoc_plot <- (STI_dw_plot + plot_spacer() + STI_lw_plot) + plot_layout(widths = c(0.49, 0.01, 0.49)) + # add a small space between the two plots
  plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 14, face ="bold"))

ggsave(posthoc_plot, filename = "Fig6_PosthocPlot.pdf", path = here("Figures"), width=20, height=10, units = "cm")