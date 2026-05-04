###### MAPS Phenology #######
### Author(s): XXX (removed for peer review)
### Script name: Step6_CalculatingClimateVariables.R

########## Objective/Description of Script #####################
# we have two types of climate data: monthly temperature anomalies from Berkeley Earth, daily precipitation from Daymet 
# this script uses these data to obtain 3 climate measure for both the long-term and decision windows for each species/station/year
# the measures we calculate are: average temperature anomaly, total precipitation, and the coefficient of variation for precipitation
# at the end of the script, we combined the climate measures for the two windows with annual light values and the breeding phenology data
#################################################################

#### Setup ####
# load packages
library(tidyverse)
library(here)
library(terra)
library(sf)
library(ncdf4)

### Load data

### decision and long-term window dates
# this also contains the earliest capture of a first year bird for each species/station/year
HY_timeline <- readRDS(here("Outputs", "HY_timeline.rds"))

### precipitation
# there are 27 files that need to be imported and combined
prcp_csv <- list.files(path= here("Data/Daymet Prcp Files"), pattern = "*.csv", full.names = TRUE)

# read all files into a list of data frames and combine them by rows
prcp <- prcp_csv %>%
  purrr::map_dfr(read_csv) %>%
  select(ID, Date, DAYMET_004_prcp) %>% # simplify columns
  rename(STA = ID, # rename some columns to facilitate joining with other data frames
         prcp = DAYMET_004_prcp) %>%
  mutate(year = year(Date)) # add a column with year

### temperature anomalies
ras.global <- rast(here("Data", "North_America_TAVG_Gridded_0p25deg.nc"))

### list of MAPS stations
stations <- readRDS(here("Outputs", "STA_finallist.rds"))

### annual light values for stations
light <- readRDS(here("Outputs", "light.rds"))

#######################################################################
###################### PRECIPITATION DATA ############################# 

# daily precipitation data for North America
# data were obtained from Daymet and accessed using the NASA AppEEARS website by extracting daily values to the MAPS station coordinates
# we use these data to obtain total precipitation and coefficient of variation (COV) of precipitation during the two windows
# the goal is to obtain total and COV values for the decision window (DW) and long-term window (LW) for each species-station-year combination

head(prcp)

prcp_count <- prcp %>% group_by(STA, year) %>% count()
# there should be 365 per STA and year for all the days of the year
# for leaps years, there should also be 365 days - Daymet leaves out Dec 31 during leap years
# for 1991, we only pulled Oct 1 to Dec 31 (92 days)
# for 2018, we only pulled Jan 1 to Oct 1 (274 days)
# for more info see: https://daymet.ornl.gov/single-pixel-tool-guide

# simplify daymet data
# group data by station and make into a nested data frame
prcp_nest <- prcp %>%
  select(STA, Date, prcp) %>% 
  rename(date = Date) %>%
  group_by(STA) %>%
  nest()

# join the species-station-year window dates to the nested prcp data
prcp_dates <- left_join(HY_timeline, prcp_nest) %>%
  rename(prcp_data = data)

# calculate total prcp and coefficient of variation for prcp for the decision and long-term window at each species-station-year
prcp_DW_LW <- prcp_dates %>%
  mutate(prcp_DW = purrr::pmap(list(prcp_data, decision_start, decision_end), # use pmap because we need to use 3 columns of the nested data frame. Specify them here and put them into a list
                         \(prcp, start_dw, end_dw) # specify how to refer to the columns/variables we want to use
                         filter(prcp, date >= start_dw & date <= end_dw)), # we want to create a data frame for each species-station-year that contains the daily prcp for the decision window
         prcp_DW_total = purrr::map(prcp_DW, ~sum(.x$prcp)), # using the data frame made in the previous step, find the total/cumulative prcp across the decision window
         prcp_DW_cov = purrr::map(prcp_DW, ~sd(.x$prcp)/mean(.x$prcp)), # also find the coefficient of variation for decision window prcp
         prcp_LW = purrr::pmap(list(prcp_data, RY_start, decision_start), # repeat a similar process for the long-term window. The RY start is the start date and the beginning of the DW is the end date for this window
                        \(prcp, start_lw, end_lw)
                        filter(prcp, date >= start_lw & date <= end_lw)),
         prcp_LW_total = purrr::map(prcp_LW, ~sum(.x$prcp)), # using the data frame made in the previous step, find the total/cumulative prcp across the decision window
         prcp_LW_cov = purrr::map(prcp_LW, ~sd(.x$prcp)/mean(.x$prcp))) %>% # also find the coefficient of variation for decision window prcp
  select(STA, SPEC, year, prcp_DW_total, prcp_DW_cov, prcp_LW_total, prcp_LW_cov) %>% # select columns to keep
  mutate(prcp_DW_cov = ifelse(prcp_DW_cov=="NaN", 0, prcp_DW_cov), # where total prcp for the DW is zero produces NA for cov. Replace these NAs with zeros
         prcp_LW_cov = ifelse(prcp_LW_cov=="NaN", 0, prcp_LW_cov)) %>%  # where total prcp for the LW is zero produces NA for cov. Replace these NAs with zeros 
  mutate(across(c(prcp_DW_total, prcp_DW_cov, prcp_LW_total, prcp_LW_cov), as.numeric)) # set columns as numeric

#View(prcp_DW_LW)

# examine the correlation between decision window and long-term precipitation
cor(prcp_DW_LW$prcp_DW_total, prcp_DW_LW$prcp_LW_total)

(prcp_corplot <- prcp_DW_LW %>%
    ggplot(., aes(x=prcp_DW_total, y=prcp_LW_total)) +
    geom_point() +
    geom_smooth(method="lm") +
    labs(x="total precip - decision window", y= "total precip - long-term window") +
    theme_classic())

#######################################################################
################ TEMPERATURE ANOMALY DATA ############################# 

##### Extract anomaly values #####
# first, we need to  raster and extract values for the MAPS stations from the temperature anomaly raster
# we will be using ras.global 

# make a vector of MAPS stations
sta <- stations$STA

# store projection details for temperature anomalies as an object called crs_TA
crs_TA <- crs(ras.global) # EPSG:4326 is the projection (+proj=longlat +datum=WGS84 +no_defs)

# make the station locations into spatial points and set CRS as NAD83
sta_coord <- stations  %>% dplyr::select(DECLNG, DECLAT) %>% # keep columns DECLNG and DECLAT
  rename(lon=DECLNG, lat = DECLAT) # make a df of coordinates with columns lon and lat

sta_points <- vect(sta_coord, crs ="EPSG:4269") # make the df into a vector and specify CRS as NAD83

TA_points_prj <- sta_points %>% project(crs_TA) # project points to CRS of the temp anomaly layer

# use the extract function from the terra package with the projected points to get temperature anomaly values for each station
TA_pts_extract <- terra::extract(ras.global,TA_points_prj)

# pivot data and add month and year information
# data includes monthly values that span from 1850 through Feb 2025 (at the time of acquisition)
# each column called temperature_x refers to a specific month and year. temperature_1 is Jan of 1850 and temperature_2102 is Feb of 2025
tempanom_pivot <- TA_pts_extract %>%
  select(starts_with("temperature")) %>%
  mutate(STA = sta) %>% # add station info 
  pivot_longer(., cols = starts_with("temperature"),
               names_to = "timepoint",
               values_to = "TAVGanom") %>% 
  group_by(STA) %>%
  mutate(month = c(rep(1:12, times = 175), 1, 2), # add column with month 
         year = c(rep(1850:2024, each = 12), rep(2025, each =2))) # add column with year
  
# filter the data to keep only relevant years for MAPS
TAVGanom <- tempanom_pivot %>%
  select(STA, TAVGanom, month, year) %>%
  filter(year >= 1991 & 
           year <= 2018)

# confirm there are no years/months with missing temperature anomaly values
sum(is.na(TAVGanom$TAVGanom))

saveRDS(TAVGanom, here("Outputs", "TAVGanomaly.rds"))

###############################################
##### Calculate decision and long-term window #####
# get temperature anomalies for the 30 day decision window and 60 day long term windows
# calculating average anomaly across the windows for each station/species/year combo
# temp anomaly data is monthly
# we need to identify the months and number of days for each month within each window
# then we can calculate the weighted average temp anomaly for each window

# define a function that takes two inputs
date_df <- function(x, y) {
  tibble( # create a tibble
    dates = seq(x, y, by="days"), # make a column "dates" that is a sequence of dates from x to y
    month = month(dates) # add a column that contains the month associated with each date
  ) %>%
    group_by(month) %>% # group the tibble by month
    summarize(
      num_days = n() # count the number of days within each month
    )
}

# apply the function created above to the two columns of dates that reflect the start and end of each window
# decision window: spans decision_start to decision_end
# long-term window: spans long_start_temp to decision_start
HY_timeline_windows <- HY_timeline %>%
  mutate(DW_dates = purrr::map2(decision_start, decision_end, date_df),
         LW_dates = purrr::map2(long_start_temp, decision_start, date_df))

# find the average temperature anomaly for the decision window within each station-species-year         
tempanom_DW <- HY_timeline_windows %>% 
  select(STA, SPEC, year, DW_dates) %>% # keep only relevant columns
  unnest(DW_dates) %>% # unnest the DW_dates data frame created in the previous step
  left_join(., TAVGanom) %>% # join the temp anomaly data
  rowwise() %>% 
  mutate(TAVGanom_prop = num_days*TAVGanom) %>% # multiply the monthly temp anomaly value by the number of days for that month within the window
  group_by(STA, SPEC, year) %>% # group the data by STA, SPEC and year
  summarize(tempanom_DW = sum(TAVGanom_prop)/30) # sum the values in the TAVGanom_prop column and divide by 30 days for the DW

tempanom_LW <- HY_timeline_windows %>% 
  select(STA, SPEC, year, LW_dates) %>% # keep only relevant columns
  unnest(LW_dates) %>% # unnest the LW_dates data frame created in the previous step
  left_join(., TAVGanom) %>% # join the temp anomaly data
  rowwise() %>% 
  mutate(TAVGanom_prop = num_days*TAVGanom) %>% # multiply the monthly temp anomaly value by the number of days for that month within the window
  group_by(STA, SPEC, year) %>% # group the data by STA, SPEC and year
  summarize(tempanom_LW = sum(TAVGanom_prop)/60) # sum the values in the TAVGanom_prop column and divide by 60 days for the LW
  
# examine the correlation between decision window and long-term window temp anomalies
cor(tempanom_DW$tempanom_DW, tempanom_LW$tempanom_LW)

(temp_corplot <- left_join(tempanom_DW,tempanom_LW) %>%
  ggplot(., aes(x=tempanom_DW, y=tempanom_LW)) +
  geom_point() +
  geom_smooth(method="lm") +
  labs(x="average temp anomaly - decision window", y= "average temp anomaly - long-term window") +
  theme_classic())

##################################################################
##### Combine Environmental Data with Bird Phenology Data #####

# combine the earliest hatch year capture date with light values, temp anomalies for DW & LW and precip for DW and LW (total and COV)

combined_climate_light_FY <- light %>%
  mutate(STA = as.integer(STA)) %>% # STA is stored as an integer in all the other dfs, so convert it here to enable join
  left_join(HY_timeline, .) %>% # join hatch year capture data to light
  select(STA, SPEC, year, date_1stcapture, doy_1stcapture, light) %>% # simplify to drop columns that are no longer needed
  left_join(., prcp_DW_LW) %>% # join prcp DW and LW values
  left_join(., tempanom_DW) %>% # join temp anomaly values for DW and LW
  left_join(., tempanom_LW)

# export
saveRDS(combined_climate_light_FY, here("Outputs", "combined_climate_light_FY.rds"))
