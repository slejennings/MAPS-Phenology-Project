###### MAPS Phenology #######
### Script name: Step2_StandardizeStationEffort.R
## Author(s): SLJ

########## Objective/Description of Script #####################
# a goal of the MAPS program is that stations perform consistent/standardized monitoring effort within and across years
# however, this goal is not always met
# the objective of this script is to find the most consistent time series at each station 
# this includes finding the combination of years and IPs within each year at each station to retain that maximize the time series while meeting the requirement for standardized/consistent effort
# we drop stations and/or years within stations and/or IPs (within years) that have highly inconsistent effort
#################################################################

#### Setup ####
# load packages
library(tidyverse)
library(here)

# import MAPS effort data
# 2 data frames
effort <- read.csv(here("Data", "MAPS_EFFORT_net_open_and_close_times.csv"), header=T)
head(effort)
str(effort)

# some of the more recent years of data (2016-2018) were missing from MAPS_EFFORT_net_open_and_close_times even though they had capture records
# we reached out to MAPS staff to obtain the missing information for these stations and years
# import them here
missingeffort <-read.csv(here("Data", "MAPS_missing_effort.csv"), header=T)
str(missingeffort)

# import station list from previous script
STA_finallist <- readRDS(here("Outputs", "STA_finallist.rds"))
nrow(STA_finallist) # 388 stations

# get vector of STA
STAlistfinal <- STA_finallist %>% pull(STA) 

# format DATE and DELAY in both the effort data frames to ensure they match and can be successfully joined
missingeffort <- missingeffort %>% mutate(DATE=mdy(DATE), DELAY=as.character(DELAY))
effort <- effort %>% mutate(DATE=as_date(DATE))

# combine effort and reduce to only stations of interest
effort_all <- bind_rows(effort, missingeffort) %>%
  filter(STA %in% STAlistfinal)

#############################################

# station effort is distributed across IP (intended periods)
# examine how consistent effort is at each STA across all the years where it has data
# we will do this by looking at the IPs that were monitored

# identify the IPs at each station in each year
STA_Yr_IP_list <- effort_all %>% 
  select(STA, YEAR, IP) %>%
  arrange(STA, YEAR, IP) %>%
  distinct() %>%
  group_by(STA, YEAR) %>%
  summarize(IPlist_length = n(),
            IPlist = str_c(IP, collapse = " "),
            min_IP = min(IP),
            max_IP = max(IP))


#View(STA_Yr_IP_list)
# a quick scan of this shows that some STA (stations) are pretty consistent year to year, but others are not

# find the number of years at each station for each combination of IPs
IPlist_count <- STA_Yr_IP_list %>%
  group_by(STA, IPlist) %>%
  count() %>%
  rename(IPlist_count = n)

# get the number of years of data at each station
STA_Yr_count <- STA_Yr_IP_list %>%
  group_by(STA) %>%
  count() %>%
  rename(STAYr_count = n)

# combine the data frames from the previous two steps
IP_summary <- left_join(IPlist_count, STA_Yr_count) 
View(IP_summary) # scroll through this to see the variation across years

# do any stations have the same IP list across all their years of operation?
IP_consistent <- IP_summary %>%
  rowwise() %>%
  mutate(proportion = (IPlist_count/STAYr_count)) %>%
  filter(proportion == 1)

nrow(IP_consistent) # 80 stations

#############################################

# all the remaining stations (n = 308) have some degree of inconsistencies in the IPs when they monitored across years

# identify these more complex stations and put them into a data frame
STAfilter <- IP_consistent %>% pull(STA) # get vector of all the "consistent" stations

IP_complex <- IP_summary %>%
  filter(! STA %in% STAfilter) %>% # remove all the stations previously identified as having consistent monitoring effort
  rowwise() %>%
  mutate(proportion = (IPlist_count/STAYr_count))

length(unique(IP_complex$STA)) # this is number of stations with more "complicated" monitoring effort

#View(IP_complex) 


# add the length of each IP list back into the data frame
# e.g., the IP list of 2 3 4 5 6 7 8 9 contains 8 unique IP periods, so the IPlist_length = 8

IP_complex_length <- STA_Yr_IP_list %>% # the info on IPlist_length is contained in STA_Yr_list
  select(STA, IPlist, IPlist_length) %>%
  distinct() %>%
  left_join(IP_complex, .)

# for each station, identify the list of IPs that is the most common (has largest value for IPlist_count)
max_IPlistcount <- IP_complex_length %>%
  group_by(STA) %>%
  arrange(desc(IPlist_count)) %>% # arrange in descending order based on IPlist_count
  slice(1) %>% # use slice(1) to keep the row with largest IPlist_count for each station
  select(STA, IPlist, IPlist_count, IPlist_length) %>%
  rename(IPlist_1st = IPlist, IPlist_count_1st = IPlist_count, IPlist_length_1st = IPlist_length) # rename some columns


#  identify the second most common list of IPs at each station
second_IPlistcount <- IP_complex_length %>%
  group_by(STA) %>%
  arrange(desc(IPlist_count)) %>% # arrange in descending order based on IPlist_count
  slice(2) %>% # keep only the row that is second
  select(STA, IPlist, IPlist_count, IPlist_length) %>%
  rename(IPlist_2nd = IPlist, IPlist_count_2nd = IPlist_count, IPlist_length_2nd = IPlist_length) # rename some columns

# compare the most common and 2nd most common IP list at each STA
# calculate how much their IP lists differ in length
# we are interested in stations where the 1st and 2nd most common IP list differ by only one IP
compare_IPlistcount <- left_join(max_IPlistcount, second_IPlistcount) %>%
  rowwise() %>%
  mutate(length_diff = IPlist_length_1st - IPlist_length_2nd) %>%
  select(-IPlist_length_1st, - IPlist_length_2nd)

range(compare_IPlistcount$length_diff)
# we can see that some of the second most common IP list are longer and some are shorter than the most common combination of IPs
# negative values indicate the second most common IP list contains more IPs and is longer in length
hist(compare_IPlistcount$length_diff)
compare_IPlistcount %>% group_by(length_diff) %>% count()

#############################################
# Identify how the IP lists differ
# aka: which IP(s) are in one list but not the other?
# we only want to modify stations where the 1st and 2nd most common IP list differ by one IP

# First, identify how the IP lists differ for all situations where the most common IP list at a station is longer than the 2nd most common
# To harmonize the 1st and 2nd most common IP list we would need to delete data from the most common list to make it the same length/combination of IPs as the 2nd most common

# Note: we are filtering for stations where the length of the two lists only differs by one
# However, the unique IPs in each list can be different by more than one
# This is a rare situation, but occurs a handful of times
# We add an extra column called IPlist_IPdiff to find and remove these situations

IPlist_diff_one <- compare_IPlistcount %>%
  filter(length_diff == 1) %>%
  rowwise() %>%
  mutate(strdiff=list(setdiff(unlist(strsplit(as.character(IPlist_1st), " ")), 
                              unlist(strsplit(as.character(IPlist_2nd), " "))))) %>%
  unnest(c(STA, strdiff)) %>%
  group_by(STA) %>%
  summarize(IPlistdiff1 =  str_c(strdiff, collapse = " "), # identify IP(s) to delete from the most common list
            IPlist_IPdiff  = n()) %>%
  filter(IPlist_IPdiff == 1) %>% # remove any rare situations where unique IPs that differ between the two lists is actually more than 1
  mutate(IPlistdiff2 = "NA") # add column of NAs to facilitate binding rows in subsequent step


# next, identify how the IP lists differ for all for all situations where the 2nd most common IP list at a station is longer than the most common
# to harmonize the 1st and 2nd most common IP list we would need to delete data from the second most common list to make it the same same length/combination of IPs as the most common
# we will again create and use the IPlist_IPdiff column as we did above 
IPlist_diff_negativeone <- compare_IPlistcount %>%
  filter(length_diff == -1) %>%
  rowwise() %>%
  mutate(strdiff=list(setdiff(unlist(strsplit(as.character(IPlist_2nd), " ")), 
                              unlist(strsplit(as.character(IPlist_1st), " "))))) %>%
  unnest(c(STA, strdiff)) %>%
  group_by(STA) %>%
  summarize(IPlistdiff2 =  str_c(strdiff, collapse = " "), # identify IP(s) to delete from second most common list
            IPlist_IPdiff  = n()) %>%
  filter(IPlist_IPdiff == 1) %>% # remove any rare situations where unique IPs that differ between the two lists is actually more than 1
  mutate(IPlistdiff1 = "NA") # add column of NAs to facilitate binding rows in subsequent step

# combine the lists
IPlist_diff_complete <- bind_rows(IPlist_diff_one, IPlist_diff_negativeone) %>%
  select(-IPlist_IPdiff) # remove this column as it is no longer needed

# join to other relevant info from compare_IPlistcount
compare_IPlist_complete <- left_join(IPlist_diff_complete, compare_IPlistcount) %>%
  rowwise() %>%
  mutate(combined_IPcount = IPlist_count_1st + IPlist_count_2nd) %>% # calculate number of years of data if able to use the two IP lists
  left_join(., STA_Yr_count) %>% # add total years of data at each station
  rowwise() %>%
  mutate(years_deleted = STAYr_count - combined_IPcount, # add column with number of years deleted
         prop_retained = combined_IPcount/STAYr_count) %>% # add column with proportion of yrs of data retained using these two IP lists
  select(STA, IPlist_1st, IPlist_count_1st, IPlist_2nd, IPlist_count_2nd, IPlistdiff1, IPlistdiff2, length_diff,
         combined_IPcount, STAYr_count, years_deleted, prop_retained) %>% # reorder columns 
  filter(combined_IPcount >= 10) # keep only stations that would have 10 or more years after harmonizing the 1st and 2nd most common IP lists
  
#View(compare_IPlist_complete)
# the IPlistdiff1 and IPlistdiff2 columns identify which IPs would need to be removed
# IPlistdiff1 are IPs to be dropped from the most common list
# IPlistdiff2 are IPs to be dropped from the 2nd most common list

#############################################
# For STA where length_diff = -1
# these are locations where the 2nd most common list of IPs is longer (contains more IPs) than the most common
# for these years, we will delete the extra IPs from the years that have IPlist_2nd
# and then combine all years with IPlist_1st and IPlist_2nd

IPlist_lengthdiffneg1 <- compare_IPlist_complete %>% 
  filter(length_diff == -1) %>%
  select(STA, IPlist_1st, IPlist_count_1st, IPlist_2nd, IPlist_count_2nd, IPlistdiff1, IPlistdiff2)

nrow(IPlist_lengthdiffneg1) # 55 stations

# For STA where length_diff = 1
# these are locations where the most common list of IPs is longer (contains more IPs) than the 2nd most common
# for these, we need to find a balance between deleting IPs and/or deleting years

# if the number of years of data with the 1st (most) common IP list is not equal to 10 or more, then we need to delete IPs and combine years
# identify these stations
IPlist_lengthdiff1_combine <- compare_IPlist_complete %>% 
  filter(length_diff == 1) %>% # keep stations with length_diff == 1
  filter(IPlist_count_1st < 10) %>% # keep stations where IPlist_1st has fewer than 10 years of data
  select(STA, IPlist_1st, IPlist_count_1st, IPlist_2nd, IPlist_count_2nd, IPlistdiff1, IPlistdiff2)

nrow(IPlist_lengthdiff1_combine) # 34 stations

# For stations where the station's 1st (most) common IP list has 10 or more years of data, then consider two options:
# 1) keeping only years with the most common IP list
# 2) combining years with 1st and 2nd most common IP list after deleting IPs
# we will examine the total number of IPs retained under each scenario:

IPlist_lengthdiff1_compare <- left_join(max_IPlistcount, second_IPlistcount) %>% # we need to add the length of each IP list back into data frame
  left_join(compare_IPlist_complete, .) %>% 
  filter(length_diff == 1) %>% # filter to only on stations with length_diff == 1
  filter(IPlist_count_1st >= 10) %>% # filter to stations where IPlist_1st has 10 or more years of data
  rowwise() %>% # 
  mutate(total_IPs_only1st = IPlist_count_1st*IPlist_length_1st, # find the total # of IPs across all years with list_1st
          total_IPs_combine =  # this will find the total number of IPs if we delete IPs to combine years with 1st and 2nd most common list
           (((IPlist_count_1st*IPlist_length_1st) - IPlist_count_1st)  # find the total IPs across all years with list_1st, then subtract one IP from each year to mimic the IP deletion that would need to occur
          + (IPlist_count_2nd*IPlist_length_2nd)), # add the total # of IPs across all years with list_2nd
        which_to_keep = if_else(total_IPs_only1st >= total_IPs_combine, "only_1st", "combine_1st_2nd"))

IPlist_lengthdiff1_compare %>% group_by(which_to_keep) %>% count()


# identify STA where we will keep only years with 1st (most) common IP list
IPlist_lengthdiff1_1stonly <- IPlist_lengthdiff1_compare %>% 
  filter(which_to_keep == "only_1st") %>%
  select(STA, IPlist_1st, IPlist_count_1st)

nrow(IPlist_lengthdiff1_1stonly) # 81 stations

# identify STA where we will drop one IP from years with 1st (most) common IP list and combine with years with 2nd most common IP list
IPlist_lengthdiff1_combine1st2nd <- IPlist_lengthdiff1_compare %>% 
  filter(which_to_keep == "combine_1st_2nd") %>%
  select(STA, IPlist_1st, IPlist_count_1st, IPlist_2nd, IPlist_count_2nd, IPlistdiff1, IPlistdiff2)

nrow(IPlist_lengthdiff1_combine1st2nd) # 17 stations

#############################################
# Next, return to stations where length_diff is not equal to -1 or 1
# we do not want to delete IPs from these stations, but we can delete years to standardize effort
# find the stations that have more than 10 years for IPlist_1st

STAlengthdiff1 <- compare_IPlist_complete %>% pull(STA)

IPlist_lengthdiff_not1 <- compare_IPlistcount %>%
  filter(!STA %in% STAlengthdiff1) %>%
  filter(IPlist_count_1st >= 10) %>%
  select(-IPlist_2nd, -IPlist_count_2nd, -length_diff)

nrow(IPlist_lengthdiff_not1) # 44 stations

#############################################
# Put together all the stations and IPlists that remain
# This requires a few steps

# First, get all STAs that need no further processing
IPlist_noprocessing <- IP_consistent %>%
  select(STA, IPlist, IPlist_count)

nrow(IPlist_noprocessing) # 80 stations

# Second, get STAs where we need to delete years
IPlist_lengthdiff_not1 # 44 stations
IPlist_lengthdiff1_1stonly # 81 stations

IPlist_deleteYrs <- bind_rows(IPlist_lengthdiff_not1, IPlist_lengthdiff1_1stonly) %>%
  rename(IPlist = IPlist_1st, IPlist_count = IPlist_count_1st)

nrow(IPlist_deleteYrs) # 125 stations

# Finally, get STAs where we need to delete IPs and combine years
IPlist_lengthdiffneg1 # 55 stations
IPlist_lengthdiff1_combine # 34 stations
IPlist_lengthdiff1_combine1st2nd # 17 stations

IPlist_deleteIPs <- bind_rows(IPlist_lengthdiffneg1, IPlist_lengthdiff1_combine, IPlist_lengthdiff1_combine1st2nd) %>% # combine 3 data frames
  rowwise() %>% 
  mutate(IPlist_keep = if_else(IPlistdiff2>0, IPlist_1st, IPlist_2nd), # identify which IPlist to keep for each station
         IPlist_totalyrs = IPlist_count_1st + IPlist_count_2nd ) %>% # find total number of years at each station after combining IPlists
  pivot_longer(cols = c(IPlist_1st, IPlist_2nd), names_to = "IPlist_type", values_to = "IPlist") %>% # make longer
  mutate(IPtype = if_else(IPlist_type == "IPlist_1st", "1st", "2nd")) %>% # add better labels to each row
  rowwise() %>% 
  mutate(IPlist_count = if_else(IPtype == "1st", IPlist_count_1st, IPlist_count_2nd)) %>% # add number of years for each IPlist
  select(-IPlist_count_1st, -IPlist_count_2nd, -IPlist_type, -IPlistdiff1, -IPlistdiff2) %>% # remove 
  left_join(., STA_Yr_IP_list)  %>% # join with information about years
  select(STA, YEAR, IPlist, IPtype, IPlist_count, IPlist_length, IPlist_keep, IPlist_totalyrs)
         
length(unique(IPlist_deleteIPs$STA)) # 106 stations

# Put all STAs and IPlists together in one data frame
# add the YEAR associated with each IPlist
STA_Yr_IP_step1 <- bind_rows(IPlist_noprocessing, IPlist_deleteYrs) %>% 
  left_join(., STA_Yr_IP_list) %>% # combine with STA_Yr_IP_list to associate years with the retain IPlists at each station
  select(-min_IP, -max_IP, -IPlist_length)

# for stations that need IPs deleting, we are keeping all the Years with 1st and 2nd most common IP list 
# but only the IPlist that needs to be used on all years
STA_Yr_IP_step2 <- IPlist_deleteIPs %>%
  select(STA, YEAR, IPlist_keep, IPlist_totalyrs) %>%
  rename(IPlist = IPlist_keep, 
         IPlist_count = IPlist_totalyrs)

STA_Yr_IP_keep <- bind_rows(STA_Yr_IP_step1, STA_Yr_IP_step2)
STA_Yr_IP_keep
length(unique(STA_Yr_IP_keep$STA)) # 311 stations

# As a final filtering step, we need to check again for lags in the timeseries at each station
# our requirement is no lags of 5 or more years

# find all stations with either no lags or lags less than or equal to 5 yrs
lags_lessthan5 <- STA_Yr_IP_keep %>%
  arrange(STA, YEAR) %>%
  group_by(STA) %>% # group by STA which is a unique identifier for station
  mutate(yrlag = YEAR - lag(YEAR)) %>% # calculates the lag between the consecutive rows of year
  filter(yrlag <= 5) # look at everywhere that the lag is less than or equal to 5 years
# all of these stations are good to go
# they have 10 or more years and no lags longer than 5 years

# identify stations with lags that are longer than 5 yrs
lags_morethan5 <- STA_Yr_IP_keep %>%
  arrange(STA, YEAR) %>%
  group_by(STA) %>% # group by STA which is a unique identifier for station
  mutate(yrlag = YEAR - lag(YEAR)) %>% # calculates the lag between the consecutive rows of year
  filter(yrlag > 5) # look at everywhere that the lag is more than 5 years
# a few stations have lags longer than 5 years

STA_lags <- lags_morethan5 %>% pull(STA) # put those stations into a vector

# use the vector to better examine what is going on at these stations
STA_lags_examine <- STA_Yr_IP_keep %>%
  filter(STA %in% STA_lags) %>%
  arrange(STA, YEAR) %>%
  group_by(STA) %>% # group by STA which is a unique identifier for station
  mutate(yrlag = YEAR - lag(YEAR)) %>%
  arrange(STA, YEAR)
length(unique(STA_lags_examine$STA)) # 9 stations

# decide what to do for each station:
# STA 11109 - delete 1992, 1993 and retain 2001 to 2010 (10 yrs)
# STA 11110 - delete 1992, 1993 and retain 2001 to 2010 (10 yrs)
# STA 11935 - delete 1992, retain 1998 to 2018 (with some small lags, 14 yrs)
# STA 11936 - delete, can't make a 10 yr time series
# STA 11948 - delete 2012, 2013, retain 1996 to 2006 (11 yrs)
# STA 15510 - delete 2015, 2016, 2017, retain 1992 to 2001 (10 yrs)
# STA 16658 - delete 2015, 2016, retain 1996 to 2009 (12 yrs)
# STA 16659 - delete 2015, 2016, retain 1996 to 2009 (12 yrs)
# STA 16696 - delete, can't make a 10 yr time series (9 is max possible)

# put these decisions into a data frame
STA_remove <- data.frame(STA = as.integer(c(rep("11109", 2), rep("11110", 2), "11935", rep("11948", 2), rep("15510", 3), rep("16658",2), rep("16659",2))),
                         YEAR = as.integer(c(rep(c("1992", "1993"), 2), "1992", "2012", "2013", "2015", "2016", "2017", rep(c("2015", "2016"), 2))))
print(STA_remove)

# make a final list of stations, with retained IP lists and associated years
STA_Yr_IP_final <- STA_Yr_IP_keep %>%
  filter(!STA %in% c("11935", "16696")) %>% # remove STA 11935 and 16696 which need to be dropped due to not enough years without long lags
  anti_join(., STA_remove)  # remove the STAs and years identified above that need to be dropped to deal with lags


nrow(STA_Yr_IP_final)
length(unique(STA_Yr_IP_final$STA)) # 309 stations

# export this list
saveRDS(STA_Yr_IP_final, here("Outputs", "STA_Yr_IP_final.rds"))

#############################################
# complete a final check of any remaining stations to ensure we have kept as many as possible

STA_clean <- STA_Yr_IP_final %>% pull(STA) # get vector of stations that are being kept

IP_complex_finalcheck <- IP_complex %>%
  filter(!STA %in% STA_clean) %>% # remove stations that have already been selected to keep
  filter(!STA %in% c("11935", "16696")) # remove stations that have long lags identified in previous step

#View(IP_complex_finalcheck)
length(unique(IP_complex_finalcheck$STA)) # 77 stations
# none of these are able to be kept 
