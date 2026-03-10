###### MAPS Phenology #######
### Script name: Step5_LightCalculation.R
## Author(s): SLJ

########## Objective/Description of Script #####################
# extract annual nighttime light values for each MAPS station for 1992 through 2018
# we are using data files from Li et al. 2020 (full reference in readme file)
# we downloaded version 7 of the files, which covered 1992-2021
# link to folder where we downloaded the files: https://doi.org/10.6084/m9.figshare.9828827.v7
# these data layers are being regularly updated to add additional years, so more recent versions may be available
#################################################################

#### Setup ####
# load packages
library(sf)
library(terra)
library(tidyverse)
library(here)
library(broom)

### load MAPS station data 
stations <- readRDS(here("Outputs", "STA_finallist.rds"))


#############################################

### Import Nighttime Light Rasters

# import all .tif for nighttime light. Annual values
nightlight_rasters <- list.files(path= here("Data/Nighttime Light Files/Other years"), pattern='tif$', full.names = T)
# there should be 24 files/years

# combine rasters into a stack
NTL_1992to18 <- terra::rast(nightlight_rasters)

# 2009 to 2011 have a slightly different extent
NTL_diffextent <- list.files(path = here("Data/Nighttime Light Files/2009 to 2011"), pattern='tif$', full.names = T)
# there should be 3 files
NTL_2009to11 <-terra::rast(NTL_diffextent)

ext(NTL_1992to18)
ext(NTL_2009to11)
# different extents

# change extent of 2009 to 2011 to match other years
ext(NTL_2009to11) <- ext(NTL_1992to18)

# put all years together into a single raster
stackNL <- c(NTL_1992to18, NTL_2009to11)

# store projection details for Nighttime Lights rasters as an object called crs_NL
crs_NL <- crs(stackNL) # EPSG:4326 is the projection (+proj=longlat +datum=WGS84 +no_defs)

# make the station locations into spatial points and set CRS as NAD83
st_coord <- stations  %>% dplyr::select(DECLNG, DECLAT) %>% # keep columns DECLNG and DECLAT
  rename(lon=DECLNG, lat = DECLAT) # make a df of coordinates with columns lon and lat
points <- vect(st_coord, crs ="EPSG:4269") # make the df into a vector and specify CRS as NAD83
points_prj <- points %>% project(crs_NL) # project points to WGS84 to match the light layers

### Get nighttime light values for stations

# use the extract function with the stack of light rasters and the points to get nighttime light values for each station
pts_extract <- terra::extract(stackNL, points_prj) %>%
  mutate(ID = row_number()) # add a unique identifier

### Combine light values with station information and pivot
stlist <- stations %>% 
  dplyr::select(STA, DECLNG, DECLAT) %>% # keeping some identifier columns for the stations
  rename (Long=DECLNG, Lat = DECLAT) %>% # rename columns to Long and Lat
  mutate(ID = row_number())
  
sts_extract <- left_join(stlist, pts_extract, by = "ID") 

sts_extract_pivot <- sts_extract %>% 
  pivot_longer(
    cols = `Harmonized_DN_NTL_1992_calDMSP`:`Harmonized_DN_NTL_2011_calDMSP`, # includes 1992 - 2018, but the columns are not in order
    names_to = "year",
    values_to = "light"
  ) %>%
  mutate(year = as.numeric(parse_number(year)), 
         STA = as.character(STA)) %>%
  select(-ID)


# format data for export
light <- sts_extract_pivot %>%
  select(STA, year, light)

# confirm each station has 27 years of light values
light_check <- light %>% group_by(STA) %>% count()

# export data
saveRDS(light, here("Outputs/light.rds"))
