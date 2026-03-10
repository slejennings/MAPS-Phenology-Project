###### MAPS Phenology #######
### Script name: Step4_IdentifySpeciesStationsforModels.R
## Author(s): SLJ

########## Objective/Description of Script #####################
# The goal of this script is to examine remaining capture records
# Decide on the criteria to keep a species/station combination for subsequent models
# Determine the long-term and decision windows for each species/station combination
#################################################################

#### Setup ####
# load packages
library(tidyverse)
library(here)
library(ggforce)

# import filtered capture data
capturedat_final <- readRDS(here("Outputs", "capturedat_final.rds"))

# import list of retain stations, years and IPs
STA_Yr_IP_final <- readRDS(here("Outputs", "STA_Yr_IP_final.rds"))

# make a new df that contains the number of years of data at each station
STAYrsOp <- STA_Yr_IP_final %>% 
  group_by(STA) %>%
  count() %>%
  rename(NumYrsOp = n)

range(STAYrsOp$NumYrsOp) # perform a quick check to make sure all stations have between 10 and 27 yrs
mean(STAYrsOp$NumYrsOp) # average station timeseries is 14.6 yrs

####################################################################

# Find counts of adults and first years for each station, species, and year
allsppcounts <- capturedat_final %>% 
  group_nest(STA, SPEC, AgeCat) %>% # put each station, species and age category into its own tibble
  mutate(
  yrcount = map(data,. %>% group_by(year) %>% summarize(NumObs=n())) # within each tibble, group by the year and find the number observed
  ) %>%
  unnest(yrcount) %>% # unnest to see results in a flat tibble
  arrange(STA, SPEC)

# Reorganize the data frame slightly 
allsppcounts_age <- allsppcounts %>% 
  dplyr::select(-data) %>% # remove the data
  pivot_wider(., names_from=AgeCat, values_from = NumObs) %>% # make wider to put Adult counts and First year counts into their own columns
  mutate(Adult=replace_na(Adult,0), # replace all NAs in column Adult with 0
         FirstYear = replace_na(FirstYear, 0)) # replace all NAs in column FirstYear with 0

# check there are no NAs left in these two columns (Adult and FirstYear):         
sum(is.na(allsppcounts_age$Adult)) # should equal 0
sum(is.na(allsppcounts_age$FirstYear)) # should equal 0

# Expand the data to get every combination of years/species at each station
allsppcounts_expand <- allsppcounts_age %>% 
  group_by(STA) %>% # group by station ID
  tidyr::expand(SPEC, year) # expand the data to make sure all combinations of species and years are present

# did this work as expected?
allsppcounts_expand %>% group_by(STA) %>% summarize(count=n_distinct(SPEC)) # do we have different number of species at each station? yes
# next, confirm that SppYrs==NumYrsOp in allspp_expand
equal.test <- allsppcounts_expand %>%
  group_by(STA, SPEC) %>% # group by the station and species
  summarize(SppYrs = n()) %>% # count the number of years associated with each species
  inner_join(., STAYrsOp, by="STA") %>% # add the number of years of operation for each station
  mutate(check = if_else(SppYrs==NumYrsOp, "ok", # compare number of years for a species with number of years at that station
                         if_else(SppYrs < NumYrsOp, "less", "more"))) # identify rows where they are not equal. Mark if SppYrs is either Less than or More than NumOpYrs

unique(equal.test$check) # All should be "ok" if this worked

# Join adult and first year counts with new expanded data and make all rows with missing values have a count of zero
allsppcounts_combine <- allsppcounts_expand %>%
  left_join(., allsppcounts_age) %>% # left_join will keep all rows in allspp_expand and add rows that match in allspp_age
  mutate(Adult=replace_na(Adult,0), # replace all NAs in column Adult with 0
         FirstYear = replace_na(FirstYear, 0)) # replace all NAs in column FirstYear with 0
print(allsppcounts_combine, n=50)

# summarize data for each station and species combination
# make columns that contain:
# NumYrs_ZeroFY = number of years where no first years were captured at the station
# NumYrs_AnyFY = number of years where at least 1 first year bird was captured at the station -> needed to get earliest capture date
# TotalFY = the number of first year birds caught across all the station's years of operation
# NumYrsOp = total number of years of data for the station

FY_summary <- allsppcounts_combine %>% 
  group_by(STA, SPEC) %>% # 
  summarize(NumYrs_ZeroFY = sum(FirstYear==0), # number of years where no FY birds captured
            NumYrs_AnyFY = sum(FirstYear >0), # number of years where one or more FY bird captured
            TotalFY = sum(FirstYear)) %>% # total number of FY birds of each species at each station across all years the station operated
            inner_join(.,STAYrsOp, by="STA") # add number years operation for each station
#View(FY_summary)  

# we don't want to measure change in phenology at stations where many of the years had no captures of first year birds
# this step applies a filter to keep only species/station combinations where there were 10 or more operating years where first year birds were caught
FY_filter <- FY_summary %>%
  filter(NumYrs_AnyFY >= 10)

# this prints the number of stations for each species based on the thresholds set in the previous two steps
spp_summary <- FY_filter %>%
  group_by(SPEC) %>%
  count() %>%
  arrange(desc(n))
#View(spp_summary)

# keep all species with 30 or more models (aka 30 or more stations)
spp_summary_30 <- spp_summary %>%
  filter(n >=30)

spp30 <- spp_summary_30 %>% pull(SPEC) # get vector of species with 30 or more models

# use the vector of species to get a vector of species/station combinations to examine
SPECSTA_spp30 <- FY_filter %>%
  filter(SPEC %in% spp30) %>%
  mutate(SPEC_STA = paste(SPEC, STA, sep ="_")) %>%
  pull(SPEC_STA)

# examine each species/station time series and calculate time lags where no first year birds were caught
lags <- allsppcounts_combine %>% 
  mutate(SPEC_STA = paste(SPEC, STA, sep ="_")) %>%
  filter(SPEC_STA %in% SPECSTA_spp30) %>%
  group_by(SPEC, STA) %>% # group by SPEC and STA 
  filter(FirstYear >= 1) %>%
  arrange(SPEC, STA, year) %>%
  mutate(yrlag = year - lag(year)) # calculates the lag between the consecutive rows of year

# restrict to species/stations combinations where the lags are 5 or more years
# we will examine these visually using plots
lags5ormore <- lags %>%
  filter(yrlag >=5) %>%
  select(SPEC, STA, SPEC_STA, yrlag)
nrow(lags5ormore) # there are 106 to examine

# make vector of species/stations that we want to plot
STASPEClags <- lags5ormore %>% pull(SPEC_STA)

# make plots
# creating a bar chart for each species/station where time (year) is on x-axis and number of FY caught is on y-axis
# we will examine these plots to locate large time stretches where there are no bars (aka no FY birds)
plots <- allsppcounts_combine %>% 
  mutate(SPEC_STA = paste(SPEC, STA, sep ="_")) %>%
  filter(SPEC_STA %in% STASPEClags) %>%
  ggplot(aes(x = year, y=FirstYear)) +
  geom_col() +
  labs(x = "years", y = "# of first year captures") +
  theme_classic()

# print on multiple pages with one plot per species/station
plots + facet_wrap_paginate(~SPEC_STA, scales="free", nrow=3, ncol=3, page = 1) # view page 1 to check this is working as expected

# export all the plots in a single pdf to examine
pdf(here("Figures", "SPEC_STA_lag_plots.pdf"), onefile = TRUE)
for(i in 1:11){ 
  print(plots + facet_wrap_paginate(~SPEC_STA, scales="free", nrow =3 , ncol = 3, page = i))
}
dev.off()
# after examining these plots, we decided to keep all stations and species


# combine identified species/stations/years with the actual capture data
# keep only earliest FY birds for each species/station/year:
FirstFY_spp30 <- capturedat_final %>% 
  mutate(SPEC_STA = paste(SPEC, STA, sep ="_"),
         date_capture = as.Date(DATE, format = "%Y-%m-%d")) %>% # convert capture date to date format %>%
  filter(SPEC_STA %in% SPECSTA_spp30) %>%
  filter(AgeCat == "FirstYear") %>%
  group_by(STA, SPEC, year) %>%
  arrange(IP, date_capture) %>%
  slice_head() %>% # keep the first row in each group, which will be the earliest date
  rename(date_1stcapture = date_capture) %>% # rename the column as it now reflects the date of the earliest capture
  mutate(doy_1stcapture = yday(date_1stcapture)) # add column that converts date to day of year (doy or julian date)


# get mean earliest capture date for each species and station
mean1stcapture_spp30 <- FirstFY_spp30 %>%
  group_by(STA, SPEC) %>% # group data by station and species
  summarize(mean_doy_1stcapture = mean(doy_1stcapture)) %>% # find the mean julian of capture date for each species/station combination
  mutate(across(starts_with("mean"), ~round(.x, 0))) # round all columns that start with "mean" to nearest integer

# join the date ranges identified in the previous step back to the larger data frame
FirstFY_final <- left_join(FirstFY_spp30, mean1stcapture_spp30) %>%
  select(LOC, STATION, STA, SPEC, SPEC_STA, year, AgeCat:mean_doy_1stcapture)
# there are a lot of unnecessary columns here
# I have simplified this but any of the dropped columns could be added back in during this step if needed

# import breeding life history values from "The Birder's Handbook" by Ehrlich et al. 1988
BreedingTimeline <- read_csv(here("Data", "BreedingTimeline.csv"))

# keep only the average delay column, which will be used to estimate nest initiation date
BreedingTimeline <- BreedingTimeline %>%
  select(SPEC, DelayAvg) %>%
  mutate(DelayAvg = round(DelayAvg))

# calculate 30-day decision window and 60-day long-term window for temperature
# using day of year
window_doy <- FirstFY_final %>%
  left_join(BreedingTimeline, by = "SPEC") %>%
  mutate(SPEC_STA = paste(SPEC, STA, sep = "_")) %>%
  rowwise() %>%
  mutate(
    decision_end_doy = mean_doy_1stcapture - DelayAvg, # end of the decision window
    decision_start_doy = decision_end_doy - 29, # beginning of the decision window
    long_start_temp_doy = decision_start_doy - 59 # beginning of the long-term window for temperature
  )

# get window start and end points as dates rather than day of year
# as a final step, assign Oct 1 of the previous year as the beginning of the rain year
HY_timeline <-
  window_doy %>%
  mutate(
    origindate = as.Date(paste0(year, "-01-01")), # set the first day of the year for each row
    decision_end = as.Date(decision_end_doy-1, origin = origindate), # we need to subtract 1 because Jan 1 (origin date) is day 1
    decision_start = as.Date(decision_start_doy-1, origin = origindate),
    long_start_temp = as.Date(long_start_temp_doy-1, origin = origindate)
  ) %>%
  ungroup() %>%
  mutate(RY_start = as.Date(paste0(year-1, "-10-01")), # rain year start in date format
         RY_start_doy = yday(RY_start)) # rain year start as ordinal date
  
# the range of dates for each window are as follows:
# decision window: decision_start to decision_end
# long-term window for temp: long_start_temp to decision_start
# long-term window for rain: RY_start to decision_start

# perform a final check to confirm the windows are the correct length
# add 1 to the interval between dates to include both the start and end date in the window
HY_timeline_check <- 
  HY_timeline %>%
  rowwise() %>%
  mutate(intervalDW = (decision_end - decision_start) + 1,
         intervalLW = (decision_start - long_start_temp) + 1)
  

# export for use in subsequent scripts
saveRDS(HY_timeline, here("Outputs", "HY_timeline.rds"))

