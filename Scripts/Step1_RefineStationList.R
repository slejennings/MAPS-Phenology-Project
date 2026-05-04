###### MAPS Phenology #######
### Author(s): XXX (removed for peer review)
### Script name: Step1_RefineStationList.R

########## Objective/Description of Script #####################
# the goal of this script is to remove banding stations that do not meet basic requirements for our analysis:
# a long time series
# fairly consistent monitoring effort
# geographic coordinates available
#################################################################

#### Setup ####
# load packages
library(tidyverse)
library(here)

# Import list of MAPS stations
stationdat <- read.csv(here("Data", "MAPS_STATION_location_and_operations.csv"), header=T)
head(stationdat)
str(stationdat)

# Import bird capture data 
capturedat <- read.csv(here("Data", "MAPS_BANDING_capture_data.csv"), header=T)
head(capturedat)
str(capturedat)

#############################################################################################
#### Reduce station list to only those with more than 10 years of data ####

# use capture records to find the number of years of data at each Station (STA)
# filter to keep only stations with 10 or more years
NumYrsOp <- capturedat %>%
  group_by(STA) %>% # group by STA which is a unique identifier for station
  summarize(NumYrsOp = n_distinct(year)) %>% # find the number of years for each station
  filter(NumYrsOp >=10) # reduce the list to stations with 10 or more years of data

nrow(NumYrsOp) # 394 stations

#############################################################################################
#### Identify stations with more than 10 years of data that have lags in their operation ####

NonConsecutive <- capturedat %>%
  group_by(STA) %>% # group by STA which is a unique identifier for station
  arrange(year) %>%
  distinct(year) %>%
  summarize(NumYrsOp = n_distinct(year), # find the number of years for each station
            MinYr = min(year), # the earliest year
            MaxYr = max(year), # the most recent year
            YrList = str_c(year, collapse=" ")) %>%  # print list of the years with data
  filter(NumYrsOp >= 10) %>% # restrict to stations with 10 or more years
  mutate(TimeRange = ((MaxYr-MinYr)+1), # find number of years between earliest (min) and most recent (max)
         Consecutive = ifelse(TimeRange==NumYrsOp, "yes", "no")) %>% # identify stations with non consecutive data 
  filter(Consecutive=="no") %>% # filter to keep only stations with non consecutive data
  mutate(YrsDiff = TimeRange - NumYrsOp) %>% arrange(desc(YrsDiff)) # find the difference between the range of years and the number of years of operation 

print(NonConsecutive, n=20)
nrow(NonConsecutive) # 82 stations out of 394 with broken time series

# How much of a lag (in years) exists for each station?
Lag <- capturedat %>%
  group_by(STA) %>% # group by STA which is a unique identifier for station
  arrange(year) %>%
  distinct(year) %>%
  mutate(yrlag = year - lag(year)) %>% # calculates the lag between the consecutive rows of year
  arrange(STA, year) %>%
  filter(yrlag > 5) # look at everywhere that the lag is more than 5 years
# Note: this is all the stations, not restricted to 10 or more years of data

# combine lag and non consecutive lists to get a list of which stations with 10 or more years of data may be problematic
STA_10yr_Lag <- inner_join(NonConsecutive, Lag) %>%
  select(STA, NumYrsOp, year, yrlag, YrList) %>%
  arrange(yrlag)

print(STA_10yr_Lag, n=Inf)
# 12 rows but 11 unique stations
# SLJ and CDF decided what to do with each of these station by manually examining its time series
# The decisions were as follows:
# 15591 - Drop. Has two 6 year gaps
# 15650 - Drop. Doesn't have 10 yrs after removing the years before gap
# 16658 and 16659 - Keep but remove 2015, 2016
# 16696 - Keep 2007 - 2016 but drop 2001
# 13319 - Drop. Seven year gap in middle of time series
# 11109 and 11110 - Keep. Remove 1992 and 1993
# 15595 - Drop. Gap in middle of time series
# 16609 - Drop. Gap in middle of time series
# 15510 - Keep. Remove 2015-2018 and keep 1992-2001

# Make a data frame of the stations that are being kept but have certain years that need to be removed
# these will all still have 10 or more years of data after dropping these years
remove <- data.frame(STA = as.integer(c(rep("16658", 2), rep("16659", 2),"16696", rep("11109",2), rep("11110",2), rep("15510",4))),
                     year = as.integer(c(rep(c("2015", "2016"), 2), "2001", rep(c("1992","1993"),2), "2015", "2016", "2017", "2018")))
head(remove)

# make a list of stations and the number of years they operated after making the changes for lags
STAYrsOp <- capturedat %>%
  filter(!STA %in% c("15591", "15650", "13319", "15595", "16609")) %>% # remove stations to be dropped
  anti_join(., remove) %>% # drop certain years for stations with large time lags
  group_by(STA) %>% # group by STA which is a unique identifier for station
  summarize(NumYrsOp = n_distinct(year)) %>% # find the number of years for each station
  filter(NumYrsOp >=10) # reduce the list to stations with 10 or more years of data

nrow(STAYrsOp) # 389 stations  
print(STAYrsOp, n=20)  

###################################################################################
####  Remove Stations that are missing or have erroneous coordinates #####

# using stationdat data frame to accomplish this step
colnames(stationdat)
sort(unique(stationdat$STATE)) # identify all the US states and Canadian provinces represented in the data

STA_finallist <- left_join(STAYrsOp, stationdat) %>%
  filter(DECLAT<=90) %>% # keep all records where latitude is less than or equal to 90
  filter(DECLNG<=180) %>% # keep all records where longitude is less than or equal to 180
  select(STA, NumYrsOp, STATE, DECLAT, DECLNG)
nrow(STA_finallist) # 388 stations remain

# export the list for use in subsequent scripts
saveRDS(STA_finallist, here("Outputs", "STA_finallist.rds"))

