###### MAPS Phenology #######
### Script name: Step8_SpatialSpeciesModels.R
## Author(s): XXX (removed for peer review)

########## Objective/Description of Script #####################
# the goal of this script is to run models for each species examining the relationship between change in breeding phenology and change in environmental variables (light, climate)
# we want to obtain one model for each species for the decision window and another for the long-term window
# because MAPS stations are distributed across North America, we need to account for the spatial relationships in the data (i.e., the distance between MAPS stations) 
# to do this, we the spmodel package to fit models with a Matérn spatial covariance matrix
# in the subsequent script, we will compare the models to find the model specifications that work best for species across the long-term and decision window
#################################################################

#################################################################
#### Setup ####
#################################################################

# load packages
library(tidyverse)
library(here)
library(gridExtra)
library(patchwork)
library(colorspace)
library(spmodel)
library(sf)

# import data 

# t-statistics for bird phenology and environmental change
tstats <- readRDS(here("Outputs", "combined_t_stats.rds"))


# MAPS stations and their latitude/longitude
stations <- readRDS(here("Outputs", "STA_finallist.rds")) %>%
  mutate(LatR = round(DECLAT, 3), # round coordinates to 3 digits
         LongR = round(DECLNG, 3),
         STA = as.factor(STA)) %>%
  select(STA, LatR, LongR)

# join t-stats and station info
# set STA and SPEC as factors
sppmodels_dat <- left_join(tstats, stations) %>%
  mutate(STA = as.factor(STA),
         SPEC = as.factor(SPEC))

summary(sppmodels_dat)

# examine samples sizes for each species
# every species should have data from at least 30 stations
sample_check <- sppmodels_dat %>%
  group_by(SPEC) %>%
  count()


####### Decision Window Models ################

# approach: fit model for each species accounting for spatial autocorrelation
# convert from NAD83 (geographic coordinates) to Albers equal area (projected coordinates)
sppmodels_dat_sf <- st_as_sf(sppmodels_dat, coords = c(12:11), crs="EPSG:4269") %>% # expects longitude followed by latitude
  st_transform(., crs = "EPSG:9822") 


# group data by species, apply a spatial linear model with matern covariance to each species
# tidy the model using broom.mixed package
dwsppmods_matern <- sppmodels_dat_sf %>%
  group_by(SPEC) %>%
  group_map(
    ~broom.mixed::tidy(splm(data=.x, 
                            FY_tstat ~ scale(tempanom_DW_tstat) + scale(prcp_DW_total_tstat) + 
                              scale(prcp_DW_cov_tstat) + scale(light_tstat),
                            spcov_type = "matern"), conf.int=T, conf.level=0.95)) %>%
  setNames(unique(sort(sppmodels_dat_sf$SPEC))) %>%
  bind_rows(., .id="Species")


# Generate model diagnostic plots for each species and save those to a pdf
# we will make a df containing the species-specific models
# this is similar to what we did above but we are not tidying the output
dw.mods.df <- sppmodels_dat_sf %>%
  group_by(SPEC) %>%
  nest() %>%
  mutate(model = map(data,
                     ~splm(data=., FY_tstat ~ scale(tempanom_DW_tstat) + scale(prcp_DW_total_tstat) + 
                             scale(prcp_DW_cov_tstat) + scale(light_tstat),
                           spcov_type = "matern")))

# find variance inflation factor for models/variables and identify any with VIF >=3 indicating multicollinearity
dw.mods.VIF <- sppmodels_dat_sf %>%
  group_by(SPEC) %>%
  nest() %>%
  mutate(model = purrr::map(data,
                     ~splm(data=., FY_tstat ~ scale(tempanom_DW_tstat) + scale(prcp_DW_total_tstat) + 
                             scale(prcp_DW_cov_tstat) + scale(light_tstat),
                           spcov_type = "matern")),
         VIF = map(model, ~car::vif(.))) %>%
  select(-data, -model) %>%
  unnest(VIF) %>%
  mutate(variable = paste(rep(c("tempanom_DW_tstat", "prcp_DW_total_tstat", "prcp_DW_cov_tstat", "light_tstat"
                                ), length.out=n()))) %>%
  filter(VIF >= 3)


# get residual plots for all models
save_residual_plots_dw <- function(SPEC, model) {
  pdf(here("Figures/DWModelPlots", paste0(SPEC, " DW Species Diagnostic Plots", ".pdf")))
  par(mfrow=c(2,2), ask = F)
  plot(model, which = c(1, 2, 3, 6))
  dev.off()
}

map2(dw.mods.df$SPEC, dw.mods.df$model, save_residual_plots_dw)


################# Long-term window models ############################

# group data by species, apply a spatial linear model with matern covariance to each species
# tidy the model using broom.mixed package
lwsppmods_matern <- sppmodels_dat_sf %>%
  group_by(SPEC) %>%
  group_map(
    ~broom.mixed::tidy(splm(data=.x, 
                            FY_tstat ~ scale(tempanom_LW_tstat) + scale(prcp_LW_total_tstat) + 
                              scale(prcp_LW_cov_tstat) + scale(light_tstat),
                            spcov_type = "matern"), conf.int=T, conf.level=0.95)) %>%
  setNames(unique(sort(sppmodels_dat_sf$SPEC))) %>%
  bind_rows(., .id="Species")


# Generate model diagnostic plots for each species and save those to a pdf
# we will make a df containing the species-specific models
# this is similar to what we did above but we are not tidying the output
lw.mods.df <- sppmodels_dat_sf %>%
  group_by(SPEC) %>%
  nest() %>%
  mutate(model = map(data,
                     ~splm(data=., FY_tstat ~ scale(tempanom_LW_tstat) + scale(prcp_LW_total_tstat) + 
                             scale(prcp_LW_cov_tstat) + scale(light_tstat),
                           spcov_type = "matern")))

# find variance inflation factor for models/variables and identify any with VIF >=3 indicating multicollinearity
lw.mods.VIF <- sppmodels_dat_sf %>%
  group_by(SPEC) %>%
  nest() %>%
  mutate(model = map(data,
                     ~splm(data=., FY_tstat ~ scale(tempanom_LW_tstat) + scale(prcp_LW_total_tstat) + 
                             scale(prcp_LW_cov_tstat) + scale(light_tstat),
                           spcov_type = "matern")),
         VIF = map(model, ~car::vif(.))) %>%
  select(-data, -model) %>%
  unnest(VIF) %>%
  mutate(variable = paste(rep(c("tempanom_lw_tstat", "prcp_lw_total_tstat", "prcp_lw_cov_tstat", "light_tstat"
  ), length.out=n()))) %>%
  filter(VIF >= 3)


# get residual plots for all models
save_residual_plots_lw <- function(SPEC, model) {
  pdf(here("Figures/LWModelPlots", paste0(SPEC, " LW Species Diagnostic Plots", ".pdf")))
  par(mfrow=c(2,2), ask = F)
  plot(model, which = c(1, 2, 3, 6))
  dev.off()
}

map2(lw.mods.df$SPEC, lw.mods.df$model, save_residual_plots_lw)

# save and export models and model summaries
saveRDS(lwsppmods_matern, here("Models", "lw_model_summaries.rds"))
saveRDS(dwsppmods_matern, here("Models", "dw_model_summaries.rds"))
write.csv(lwsppmods_matern, here("Models", "lw_model_summaries.csv"))
write.csv(dwsppmods_matern, here("Models", "dw_model_summaries.csv"))
saveRDS(lw.mods.df, here("Models", "lw_models.rds"))
saveRDS(dw.mods.df, here("Models", "dw_models.rds"))

