###### MAPS Phenology #######
### Script name: Step14_PosthocAnalyses.R
## Author(s): XXX (removed for peer review)

########## Objective/Description of Script #####################
# run posthoc analysis to compare mean change in temperature experienced by each species with their STI
# we examine this relationship using change in temp for both the decision and long-term window
# create Figure 5
#################################################################

#### Setup ####

# load packages
library(tidyverse)
library(here)
library(ggeffects)
library(nlme)
library(performance)
library(broom.mixed)
library(patchwork)

# import files

# t-statistics for change in bird breeding phenology and environmental variables
tstats <- readRDS(here("Outputs", "combined_t_stats.rds"))

# species temperature indices
STI <- read.csv(here("Data", "species_STI.csv"), header=T)

### list of MAPS stations
stations <- readRDS(here("Outputs", "STA_finallist.rds"))

################################################################################
### Organize data and run post-hoc models
################################################################################

# calculate the mean t-statistic (and standard error) for temperature anomaly change across all stations in which each species occurs
# do this for both the decision and the long-term window
combT_summary <- tstats %>%
  group_by(SPEC) %>%
  summarise(
    n_decision = sum(!is.na(tempanom_DW_tstat)),
    mean_decision = mean(tempanom_DW_tstat),
    se_decision = sd(tempanom_DW_tstat) / sqrt(n_decision),
    n_long = sum(!is.na(tempanom_LW_tstat)),
    mean_long = mean(tempanom_LW_tstat),
    se_long = sd(tempanom_LW_tstat) / sqrt(n_long),
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
eff_STI_dw <- plot(ggeffects::predict_response(m1, terms =c("STI")), colors = "#6ab43e")

# add data, labels, nice formatting to plot
STI_dw_plot <- eff_STI_dw +
  geom_point(data = combDAT, aes(x = STI, y = mean_decision), color = "#6ab43e", size = 2.5, pch = 19) +
  geom_errorbar(aes(ymin = (mean_decision-se_decision), ymax = (mean_decision+se_decision), x=STI),
                color="#6ab43e", data = combDAT,inherit.aes = FALSE, lwd=.9) +
  labs(title= "", x = "Species temperature index (STI)", y = "Mean temperature change statistic")+
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), 
        axis.title.x = element_text(color = "black", size = 12, margin = margin(t=0.3, unit="cm")), # add space between axis title and axis labels
        axis.title.y = element_text(color = "black", size = 12, margin = margin(r=0.3, unit="cm"))) + # add space between axis title and axis labels
  geom_hline(yintercept=0, linetype="dashed", color = "gray", linewidth=.6)

STI_dw_plot


# get the predicted effect from the decision window model 
eff_STI_lw <- plot(ggeffects::predict_response(m2, terms =c("STI")), colors = "#1768B3")

# add data, labels, nice formatting to plot
STI_lw_plot <- eff_STI_lw +
  geom_point(data = combDAT, aes(x = STI, y = mean_long), color = "#1768B3", size = 2.5, pch = 19) +
  geom_errorbar(aes(ymin = (mean_long-se_long), ymax = (mean_long+se_long), x=STI),
                color="#1768B3", data = combDAT,inherit.aes = FALSE, lwd=.9) +
  labs(title= "", x = "Species temperature index (STI)", y = "Mean temperature change statistic")+
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), 
        axis.title.x = element_text(color = "black", size = 12, margin = margin(t=0.3, unit="cm")), # add space between axis title and axis labels
        axis.title.y = element_text(color = "black", size = 12, margin = margin(r=0.3, unit="cm"))) + # add space between axis title and axis labels
  geom_hline(yintercept=0, linetype="dashed", color = "gray", linewidth=.6)


STI_lw_plot

####################################################
# Examine how the STI of species relates to the latitude of the stations where they breed

head(stations)
summary(stations)
stations$STA <- as.factor(stations$STA) # convert STA to factor

# combine t-stats with station coordinates
tstat_sts <- left_join(tstats, stations)

# find the mean latitude of the MAPS stations where each species breeds
spp_lat_STI <- tstat_sts %>% group_by(SPEC) %>%
  summarize(mean_lat = mean(DECLAT),
            n_lat = n(),
            se_lat = sd(DECLAT)/sqrt(n_lat)) %>%
  left_join(., STI)

# model relationship between species' STI and median latitude of the MAPS stations where they breed
STI_lat_m <- gls(mean_lat ~ STI, data = spp_lat_STI, weights = varFixed(~1/sqrt(se_lat)), method = "ML")
summary(STI_lat_m)
confint(STI_lat_m)

# get the predicted effect from the model
eff_STI_lat <- plot(ggeffects::predict_response(STI_lat_m, terms =c("STI")), colors = "#4F2691")

# add data, labels, nice formatting to plot
(STI_lat_plot <- eff_STI_lat +
    geom_point(data = spp_lat_STI, aes(x = STI, y = mean_lat), color = "#4F2691", size = 2.5, pch = 19) +
    geom_errorbar(aes(ymin = (mean_lat-se_lat), ymax = (mean_lat+se_lat), x=STI),
                  color="#4F2691", data = spp_lat_STI, inherit.aes = FALSE, lwd=.9) +
    labs(title= "", x= "Species temperature index (STI)", y="Mean latitude of MAPS stations \n where each species breeds")+
    theme_classic() +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), 
          axis.title.x = element_text(color = "black", size = 12, margin = margin(t=0.3, unit="cm")), # add space between axis title and axis labels
          axis.title.y = element_text(color = "black", size = 12, margin = margin(r=0.3, unit="cm")))) # add space between axis title and axis labels

####################################################
# find the mean change in temperature anomaly for the decision and long-term window at each MAPS stations
# this averages all the species at a station
# reminder: each species has it's own 30-day and 60-day decision window 
# this means that the range of dates encompassed by the DW and LW differ between species
# because of these differences, averaging is probably more appropriate for the long-term window where there is more likely to be overlap in the dates covered

STA_tempanom <- tstats %>% 
  group_by(STA) %>%
  summarize(mean_lw_tempanom = mean(tempanom_LW_tstat), # mean change in temp anomaly at each station for LW
            n_lw_tempanom = n(),
            se_lw_tempanom = sd(tempanom_LW_tstat)/n_lw_tempanom) %>% 
  mutate(se_lw_tempanom = ifelse(is.na(se_lw_tempanom), 0, se_lw_tempanom)) %>% # some stations only have one species, so we must replace NAs for SE with 0
  left_join(., stations) # add the coordinates of the stations

# model relationship between species' STI and median latitude of the MAPS stations where they breed
# for the long-term window
STA_tempanom_lw_m <- gls(mean_lw_tempanom ~ DECLAT, data = STA_tempanom, weights = varFixed(~1/sqrt(se_lw_tempanom)), method = "ML")
summary(STA_tempanom_lw_m) 
confint(STA_tempanom_lw_m)
# positive values for ave change in temp anomalies are warming whereas negative values are cooling
# stations at lower latitudes (more southern) are experiencing greater warming than stations at higher latitudes

# make a plot for the long-term window relationship
# get the predicted effect from the model
eff_STA_tempanom_lw <- plot(ggeffects::predict_response(STA_tempanom_lw_m, terms =c("DECLAT")), colors = "#077B82")

# add data, labels, nice formatting to plot

# for stations with only one species, we need an NA for the errorbars to keep ggplot happy
# create a modified data frame for plotting
STA_tempanom_plot <- STA_tempanom %>%
  rowwise() %>%
  mutate(
    ymin_bar = case_when(
    n_lw_tempanom == 1 ~ NA,
    n_lw_tempanom >1 ~ mean_lw_tempanom-se_lw_tempanom),
    ymax_bar = case_when(
      n_lw_tempanom == 1 ~ NA,
      n_lw_tempanom >1 ~ mean_lw_tempanom+se_lw_tempanom),
    )

(TA_tempanom_lw_plot <- eff_STA_tempanom_lw +
    geom_point(data = STA_tempanom_plot, aes(x = DECLAT, y = mean_lw_tempanom), color = "#077B82", size = 2, pch = 19) +
    geom_errorbar(aes(ymin = ymin_bar, ymax = ymax_bar, x= DECLAT),
                  color="#077B82", data = STA_tempanom_plot, inherit.aes = FALSE, lwd=.8, na.rm=T) +
    labs(title= "", y = "Mean temperature change statistic \n for the long-term window", x ="MAPS station latitude")+
    theme_classic() +
    geom_hline(yintercept=0, linetype="dashed", color = "gray", linewidth=.6) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), 
          axis.title.x = element_text(color = "black", size = 12, margin = margin(t=0.3, unit="cm")), # add space between axis title and axis labels
          axis.title.y = element_text(color = "black", size = 12, margin = margin(r=0.3, unit="cm")))) # add space between axis title and axis labels

#### combine plots

posthoc_R1 <- wrap_elements((STI_dw_plot + plot_spacer() + STI_lw_plot) + plot_layout(widths = c(0.49, 0.01, 0.49)) + # add a small space between the two plots
  plot_annotation(tag_levels = list(c('(a)', '(b)'))) & theme(plot.tag = element_text(size = 14, face ="bold")))

posthoc_R2 <- wrap_elements((STI_lat_plot + plot_spacer() + TA_tempanom_lw_plot) + plot_layout(widths = c(0.49, 0.01, 0.49)) + # add a small space between the two plots
  plot_annotation(tag_levels = list(c('(c)', '(d)'))) & theme(plot.tag = element_text(size = 14, face ="bold")))

(posthoc_plot <- posthoc_R1/posthoc_R2)

ggsave(posthoc_plot, filename = "Fig5_PosthocPlot.pdf", path = here("Figures"), width=22, height=20, units = "cm", device=cairo_pdf)
#ggsave(posthoc_plot, filename = "Fig5_PosthocPlot.png", path = here("Figures"), width=22, height=20, units = "cm")

