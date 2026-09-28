###### MAPS Phenology #######
### Script name: Step9_CompareDW&LWModels.R
## Author(s): XXX (removed for peer review)

########## Objective/Description of Script #####################
# this script compares the long-term and decision window model for each species using AICc
#################################################################

#### Setup ####
# load packages
library(tidyverse)
library(here)
library(spmodel)

# import decision window models 

dw_sppmodels <- readRDS(here("Models/Model Outputs", "dw_models.rds")) %>%
  select(SPEC, model) %>% 
  rename(dw_model = model)

# import long-term window models

lw_sppmodels <- readRDS(here("Models/Model Outputs", "lw_models.rds")) %>%
  select(SPEC, model) %>%
  rename(lw_model = model)

# put the two models together into a single df
models_df <- left_join(dw_sppmodels, lw_sppmodels)

# get AICc for each model
models_AICc <- models_df %>%
  mutate(AICc_dw = purrr::map(dw_model, ~spmodel::AICc(.)),
         AICc_lw = purrr::map(lw_model, ~spmodel::AICc(.))) %>%
  select(-dw_model, -lw_model)


# compare models
models_AICc_diff <- models_AICc %>%
  rowwise() %>%
  mutate(AICc_diff = round(abs(AICc_dw-AICc_lw),1)) %>%
  filter(AICc_diff >=4) %>%
  mutate(preferredmod = ifelse(AICc_dw < AICc_lw, "decision", "long-term"))

# conclusions of model selection:
# requiring difference of at least 4 to identify a difference in model performance
# ORJU, WIWA, HOWR - long-term
# SOSP, AMRO - decision

