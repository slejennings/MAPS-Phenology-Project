###### MAPS Phenology #######
### Author(s): SLJ, LP
### Script name: Step7_CalculatingTStatistics.R

########## Objective/Description of Script #####################
# the goal of this script is to measure the change over time in breeding phenology, climate, and light
# we get values for each species/station combination
# for the climate variables, we measure change during the 30-day decision and the long-term windows
# we use the t-statistic from a linear model as our measure of change
# this variables takes into account the strength of the estimate (the slope) and the standard error associated with it
#################################################################

#### Setup #####
# load packages
library(tidyverse)
library(here)
library(broom)

# import data
env_FY_data <- readRDS(here("Outputs", "combined_climate_light_FY.rds"))

# set STA and SPEC as factors
env_FY_data <- env_FY_data %>%
  mutate(STA = as.factor(STA),
         SPEC = as.factor(SPEC))

colnames(env_FY_data)

##########################
### Get t-statistics for change in bird phenology
# using date of first hatch year capture for each species-station-year
# use column "doy_1stcapture" as the response variable, which is the ordinal date

FY_tstat <- env_FY_data %>%
  group_by(STA, SPEC) %>%
  nest() %>% # make a nested dataframe with one df for each station-species combination
  mutate(model= purrr::map(data, ~glm(doy_1stcapture ~ year, family = gaussian, data =.)), # run a model on each nested df
         tidy_m = purrr::map(model, broom::tidy)) %>% # tidy the model
  unnest(tidy_m) %>% # unnest the tidied models
  filter(term == "year") %>% # keep just the model coefficients for year (drop the intercept)
  select(STA, SPEC, statistic) %>% # identify columns to keep
  rename(FY_tstat = statistic)

#################################
### Get t-statistics for change in temperature anomalies 
# do this for the decision window and long-term window

tempanom_tstat <- env_FY_data %>%
  group_by(STA, SPEC) %>%
  nest() %>%
  mutate(DW_model= purrr::map(data, ~glm(tempanom_DW ~ year, family = gaussian, data =.)), # model for the decision window
         LW_model= purrr::map(data, ~glm(tempanom_LW ~ year, family = gaussian, data =.)), # model for the long-term window
         DW_tidy = purrr::map(DW_model, broom::tidy),
         LW_tidy = purrr::map(LW_model, broom::tidy)) %>%
  unnest(DW_tidy) %>% # unnest each set of tidied models one at a time. Starting with the decision window models
  filter(term == "year") %>%
  select(STA, SPEC, statistic, LW_tidy) %>%
  rename(tempanom_DW_tstat = statistic) %>% 
  unnest(LW_tidy) %>% # next, unnest the tidied long-term window models and add those to the df
  filter(term == "year") %>%
  select(STA, SPEC, tempanom_DW_tstat, statistic) %>%
  rename(tempanom_LW_tstat = statistic)


################################
### Get t-statistics for change in total precipitation 
# do this for the decision and long-term windows
# using a similar approach to the previous step

prcp_total_tstat <- env_FY_data %>%
  group_by(STA, SPEC) %>%
  nest() %>%
  mutate(DW_model= purrr::map(data, ~glm(prcp_DW_total ~ year, family = gaussian, data =.)),
         LW_model= purrr::map(data, ~glm(prcp_LW_total ~ year, family = gaussian, data =.)),
         DW_tidy = purrr::map(DW_model, broom::tidy),
         LW_tidy = purrr::map(LW_model, broom::tidy)) %>%
  unnest(DW_tidy) %>%
  filter(term == "year") %>%
  select(STA, SPEC, statistic, LW_tidy) %>%
  rename(prcp_DW_total_tstat = statistic) %>%
  unnest(LW_tidy) %>%
  filter(term == "year") %>%
  select(STA, SPEC, prcp_DW_total_tstat, statistic) %>%
  rename(prcp_LW_total_tstat = statistic)

################################
#### Get t-statistics for change in coefficient of variation for precipitation
# do this for the decision and long-term windows
# using a similar approach to the two previous steps

prcp_cov_tstat <- env_FY_data %>%
  group_by(STA, SPEC) %>%
  nest() %>%
  mutate(DW_model= purrr::map(data, ~glm(prcp_DW_cov ~ year, family = gaussian, data =.)),
         LW_model= purrr::map(data, ~glm(prcp_LW_cov ~ year, family = gaussian, data =.)),
         DW_tidy = purrr::map(DW_model, broom::tidy),
         LW_tidy = purrr::map(LW_model, broom::tidy)) %>%
  unnest(DW_tidy) %>%
  filter(term == "year") %>%
  select(STA, SPEC, statistic, LW_tidy) %>%
  rename(prcp_DW_cov_tstat = statistic) %>%
  unnest(LW_tidy) %>%
  filter(term == "year") %>%
  select(STA, SPEC, prcp_DW_cov_tstat, statistic) %>%
  rename(prcp_LW_cov_tstat = statistic)


#########################################
#### Get t-statistics for change in nighttime light 
# nighttime light data was only available as annual values, so there is no decision window or long term window
# group by STA to get one t-stat for each station

light_tstat <- env_FY_data %>%
  group_by(STA) %>%
  nest() %>%
  mutate(model= purrr::map(data, ~glm(light ~ year, family = gaussian, data =.)),
         m_tidy = purrr::map(model, broom::tidy)) %>%
  unnest(m_tidy) %>%
  filter(term == "year") %>%
  select(STA, statistic) %>%
  mutate(statistic = coalesce(statistic, 0)) %>% # change NAs to 0. These are stations with no change in light over time
  rename(light_tstat = statistic) 
  
################################################
#### Combine data frames containing t-statistics
combined_t_stats <- left_join(FY_tstat, tempanom_tstat, by=c("STA", "SPEC")) %>%
  left_join(., prcp_total_tstat, by=c("STA", "SPEC")) %>% 
  left_join(., prcp_cov_tstat, by=c("STA", "SPEC")) %>%
  left_join(., light_tstat, by="STA")

#View(combined_t_stats)
  
# saving data frame
saveRDS(combined_t_stats, here("Outputs", "combined_t_stats.rds"))
