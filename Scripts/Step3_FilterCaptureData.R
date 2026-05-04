###### MAPS Phenology #######
### Script name: Step3_FilterCaptureData.R
## Author(s): XXX (removed for peer review)

########## Objective/Description of Script #####################
# This script filters the bird capture data to retain only records from the stations, years and IPs identified in Step3
# It also takes multiple steps to filter the capture data to:
# remove recaptured individual
# classify age (identify hatch year birds)
# remove species and records identified by MAPs staff as potentially unreliable
# retain only locally breeding species at each station
#################################################################

#### Setup ####
# load packages
library(tidyverse)
library(here)

# import data
# list of retained stations with their years and IPs
STA_Yr_IP_final <- readRDS(here("Outputs", "STA_Yr_IP_final.rds"))

# bird capture data 
capturedat <- read.csv(here("Data", "MAPS_BANDING_capture_data.csv"), header=T)
head(capturedat)
str(capturedat)
nrow(capturedat)

###############################################################
# Retain only capture records from STAs/YEARs/IPs identified in the previous script/step

# we need to split the column IPlist and make a longer df before using it to filter capture data
# we want one row for each STA/YEAR/IP combination
head(STA_Yr_IP_final)
str(STA_Yr_IP_final)

STA_Yr_IP_long <- STA_Yr_IP_final %>%
  select(STA, IPlist, YEAR) %>%
  separate_longer_delim(IPlist, " ") %>% # make a new row for each character in IPlist
  mutate(IP = as.integer(IPlist)) %>% # format at integer to faciliate join in next step
  rename(year=YEAR)
head(STA_Yr_IP_long)

# keep only capture records that are associated with the stations, years and IPs we want to keep
capture_filter <- left_join(STA_Yr_IP_long, capturedat)

nrow(capturedat) - nrow(capture_filter) # number of rows we have removed
nrow(capture_filter) # this is the number of capture records that remain at this stage

###############################################################
# Resolve some species names

# Examine species that are in the data
unique(sort(capture_filter$SPEC))

# look at PSFL, COFL, and WEFL, which should all be classified under WEFL
flycatchers <-capture_filter %>% filter(SPEC %in% c("WEFL", "PSFL", "COFL"))
nrow(flycatchers) # there are 7897 rows
flycatchers %>% count(SPEC) # most were identified as PSFL

# make the change to label all as WEFL
capturespp <- capture_filter %>% mutate(SPEC = case_match(SPEC, c("PSFL", "COFL") ~ "WEFL", .default = SPEC))

# do we still see PSFL and COFL in the list of species names?
sort(unique(capturespp$SPEC)) # no

# confirm that all flycatchers are now listed as WEFL
capturespp %>% filter(SPEC=="WEFL") %>% nrow() == nrow(flycatchers) # should equal TRUE

###############################################################
# Remove records from individuals that were recaptured at a station within a year
# Remove records not labeled as suitable for productivity and survivorship analyses

# examine entries using column "C" which denotes recaptured individuals
recap <- capturespp %>% filter(C=="R") %>% arrange(BAND) # R is the code for recaptured birds

# look at values for BAND, which denotes the band number
unique(recap$BAND) # some recaptured birds have no band number (coded as "")
recap_noband <- recap %>% filter(BAND=="") # keep only entries that are missing band numbers
# which recaptured species are missing band numbers?
unique(recap_noband$SPEC)
# multiple species of hummingbirds (CAHU, ALHU, COHU, ANHU, BBIH, BTAH) and Northern bobwhite (NOBO)

# Next, examine data using BAND (band number) as an alternative to using column C
band <- capturespp %>% 
  filter(BAND!="") %>% # remove entries that are missing band numbers
  group_by(BAND) %>% # group the data by band number
  filter(n()>1) %>% # retain all band numbers that occur more than once (birds that are recaptured)
  arrange(BAND, year) # sort data frame by band number and year
#View(band) # look at results

# By examining "band", we can see that recaptured birds can occur in several ways
# First way: an individual is recaptured within a breeding season.
# We need to remove instances where individual birds were recaptured within a particular year/breeding season
# Second way: recaptured across multiple breeding seasons
# this is less of an issue for our analysis
# so we do not want to remove all instances where column C == R

# Make 3 separate data frames to keep track of different types of captures
# data frame 1: entries with no band number
no_band <- capturespp %>%
  filter(BAND=="")

# data frame 2: banded individuals that were not recaptured (band number only occurs once)
one_band <- capturespp %>%
  filter(BAND!="") %>% # remove entries that are missing band numbers
  group_by(BAND) %>% # group the data by band number
  filter(n()==1) # keep all band numbers that occur only once
# confirm there are no band numbers with duplicate entries:
sum(duplicated(one_band[,8])) # BAND is the 8th column. This should equal 0 if there are no duplicates

# data frame 3: banded individuals that were recaptured, but modified to only include one entry for each year an individual was sighted/caught
band_distinct <- capturespp %>%
  filter(BAND!="") %>% # remove entries that are missing band numbers
  group_by(BAND) %>%
  filter(n()>1) %>% # retain all band numbers that occur more than once (birds that are recaptured)
  arrange(BAND, year, DATE) %>% # sort data
  group_by(BAND, year) %>% # group data
  slice(1) # keep only the earliest date/capture for each band within each year

# check this worked as expected
# if it did, we should only have 1 record for each BAND in each year
band_check <- band_distinct %>%
  group_by(BAND, year) %>%
  count()
range(band_check$n) # all values of n = 1. This is what we we want!

# combine 3 types back into a single data frame
combine <- bind_rows(no_band, one_band, band_distinct)
nrow(capturespp) - nrow(combine) # the above steps removed 157809 rows

# Finally, use column "N" for as another cleaning step
# this column has the value of "-" for records that can be used in productivity and abundance analyses
unique(combine$N) # "G" for gallinaceous bird, "H" for hummingbird, "U" for unbanded should be removed
dat <- combine %>% filter(N=="-") # keep only records coded as "-"
head(dat)
nrow(combine)-nrow(dat) # removed 59149 rows


###############################################################
# Retain only locally breeding birds
# Classify all remaining captures as either first year or adults

# Keep only capture records from locally breeding species
# use column BRSTAT
unique(dat$BRSTAT)
# we want to keep records with B: Breeder, U: Usual Breeder, and O: Occasional Breeder
dat_local <- dat %>%
  filter(BRSTAT %in% c("B", "U", "O"))

# classify birds caught based on their age
capturedat_final <- dat_local %>%
  filter(AGE != "0") %>% # remove birds with unknown age
  mutate(AgeCat = if_else(AGE== 4 | AGE == 2, "FirstYear", "Adult")) # make all birds "FirstYear" or "Adult"

nrow(capturedat_final) # number of capture records remaining
unique(capturedat_final$SPEC) # list of species codes that remain

# export as .rds for use in subsequent scripts
saveRDS(capturedat_final, here("Outputs", "capturedat_final.rds"))
