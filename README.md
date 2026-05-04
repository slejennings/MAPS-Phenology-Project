# MAPS Phenology README

############################################################################
### OVERVIEW
############################################################################

This repository contains data and code for the analyses in the manuscript "It's about time: thermal niche and visual light sensitivity explain breeding phenological change due to climate and light pollution". The analyses are organized into an R Project, which will reproduce the results, tables, and figures presented in the manuscript. This Read Me file describes the required software and the organization of the R Project and the associated files.

Important: certain files are too large to sync with this GitHub repository and need to be downloaded from the accompanying Dryad folder. This applies to any file marked with an asterisk (*) at the end of the filename. Link to Dryad folder: XXXXXX (will be added upon acceptance).


############################################################################
### CORRESPONDENCE
############################################################################

Please direct questions about the data, analysis, and results to:

Removed for peer review

############################################################################
### SOFTWARE
############################################################################

Rstudio 2025.05.0+496
R version 4.4.3 (2025-02-28)

R Packages and versions:
ape (5.8-1), bayesplot (1.11.1), brms (2.22.0), broom (1.0.7), broom.mixed (0.2.9.6), colorspace (2.1-1), geiger (2.0.11), ggeffects (2.2.0), ggforce (0.4.2), ggtree (3.15.0), gridExtra (2.3), here (1.0.1), loo (2.8.0), ncdf4 (1.23), nlme (3.1-167), patchwork (1.3.0), performance (0.15.1), phylolm (2.6.5), posterior (1.6.1), rstan (2.32.6), sf (1.0-19), terra (1.8-29), tidyverse (2.0.0)


############################################################################
### R PROJECT FOLDER STRUCTURE
############################################################################

The R Project contains folders that are organized using the following organizational structure:

Data
- Nighttime Light Files
    - 2009 to 2011
    - Other years
- Daymet Prcp Files

Scripts

Outputs

Models
- Model Plots
- Model Outputs

Figures
- Figure S2


############################################################################
### DESCRIPTIONS OF R PROJECT FOLDERS AND THEIR CONTENTS
############################################################################

### DATA

**Structure of Data Folder (contains multiple subfolders):**

Data
- Nighttime Light
    - 2009 to 2011
    - Other years
- Daymet Prcp Files

-------------------------------------------------------------------------------------------------------------------
**Folder**: Data 

**Folder Contents:** 12 files, file format(s): .csv, .tre, .nc, .pdf

**Folder Description:** contains the datafiles that were used in the analysis, mostly stored as .csv. There are two subfolders: Nighttime Light and Daymet Prcp Files.

**Note:** files with * at end of name are too large to sync with the GitHub repository and need to be downloaded from the accompanying Dryad folder

**List of Files:**

***File name:*** BirdCodetoSpecies.csv

***File Description:*** lists the 4-letter species codes, common names, and scientific names for species in this analysis

***Columns:***
  - SPEC: four-letter species code
  - COMMONNAME: common name of species
  - SCINAME: scientific name of species

-------------------------------------------------------------------------------------------------------------------
***File name:*** BreedingTimeline.csv

***File Description:*** contains life history values for 20 species in this analysis taken from
Ehrlich et al. 1988.

***Columns:***
  - SPEC: four-letter species code
  - NestBuilding_days: number of days to build nest, may be a range
  - Laying_eggs: number of days to lay eggs, may be a range
  - Incubations_days: number of days spent incubating eggs, may be a range
  - Fledging_days: number days raising chicks to fledging, may be a range
  - DelayLow: minimum number of days from nest initiation to fledging for each species
  - DelayHigh:  maximum number of days from nest initiation to fledging for each species
  - DelayAvg: average number of days from nest initiation to fledging for each species

-------------------------------------------------------------------------------------------------------------------
***File name:*** Jetz_ConsensusPhy.tre

***File Description:*** consensus phylogeny based on 10,000 for 9,993 bird species based on Jetz et al. 2012

-------------------------------------------------------------------------------------------------------------------
***File name:*** MAPS-database_codes-and-structures.pdf

***File Description:*** Explains the MAPS program and the contents and structure of all of the subsequent files in this folder that have MAPS in their filename

-------------------------------------------------------------------------------------------------------------------
***File name:*** MAPS_BANDING_capture_data.csv*

***File Description:*** Contains all the banding records for captured birds

***Columns:*** refer to “MAPS-database_codes-and-structures.pdf” for a complete description of the columns and data within this file

-------------------------------------------------------------------------------------------------------------------
***File name:*** MAPS_EFFORT_net_open_and_close_times.csv*

***File Description:*** Contains information about each station's operation, including the time each mist net was opened and closed

***Columns:*** refer to “MAPS-database_codes-and-structures.pdf” for a complete description	 
of the columns and data within this file

-------------------------------------------------------------------------------------------------------------------
***File name:*** MAPS_missing_effort.csv

***File Description:*** Certain stations/years were missing effort information in MAPS_EFFORT_net_open_and_close_times.csv. This file was provided to us by MAPS staff to fill in the missing values. Its structure is similar to MAPS_EFFORT_net_open_and_close_times.csv

***Columns:*** refer to “MAPS-database_codes-and-structures.pdf” for a complete description
of the columns and data within this file

-------------------------------------------------------------------------------------------------------------------
***File name:*** MAPS_STATION_location_and_operations.csv

***File Description:*** Contains information about each MAPS station, including geographic coordinates and years of operation

***Columns:*** refer to “MAPS-database_codes-and-structures.pdf” for a complete description of the columns and data within this file

-------------------------------------------------------------------------------------------------------------------
***File name:*** North_America_TAVG_Gridded_Op25deg.nc*

***File Description:*** high resolution monthly gridded temperature anomaly data from 1850 to present for North America with 0.25 X 0.25 resolution. Obtained from Berkeley Earth (Rohde & Hausfather 2020).

-------------------------------------------------------------------------------------------------------------------
***File name:*** species_eye.csv

***File Description:*** eye dimensions for the species in this study

***Columns:***
  - Species: scientific name
  - SPEC: four-letter species code
  - Tree_name: scientific name as listed in phylogenetic tree
  - C.T: ratio to the corneal diameter to transverse diameter of the eye. This is a proxy for dim light vision
  - Source: source for each species’ C.T values
  - Museum Collection IDs: for species added in this study, this lists the museum and collection IDs for the specimens used

-------------------------------------------------------------------------------------------------------------------
***File name:*** species_STI.csv

***File Description:*** Species temperature indices. These values reflect the mean temperature within each species' summer (breeding) range in degrees celsius. Obtained from Lehikoinon et al. 2021

***Columns:***
  - SPEC: four-letter species code
  - Species: scientific name
  - STI: species temperature index

-------------------------------------------------------------------------------------------------------------------
***File name:*** species_traits.csv

***File Description:*** traits for 20 species in this study from Sheard et al. 2020 including hand wing index (HWI), body mass and the annual precipitation within a species’ range

***Columns:***
  - SPEC: four-letter species code
  - Species_name: scientific name
  - HWI: hand wing index 
  - Body_mass_log: natural log of body mass in grams
  - AnnualPrecip: annual precipitation in mm within each species’ geographic range
-------------------------------------------------------------------------------------------------------------------

**Folder: Daymet Prcp Files**

**Folder Contents:** subfolder within Data, 27 files, file format(s): .csv

**Folder Description:** daily precipitation values for 270 MAPS banding stations for 1991 through 2018. Values were obtained from Daymet and requested via the Application for Extracting and Exploring Analysis Ready Samples (AppEEARS). This service provides the daily data from the nearest 1 km x 1km Daymet grid cell for each of the geographic points provided, in this case the coordinates for the MAPS stations. There is one file for each year that spans from October of the earlier year to October of the later year. The years included in each file are listed in the filename.

**List of Files:**

MAPS-Phenology-1991-to-1992-prcp-DAYMET-004-results.csv

MAPS-Phenology-1992-to-1993-prcp-DAYMET-004-results.csv

MAPS-Phenology-1993-to-1994-prcp-DAYMET-004-results.csv

MAPS-Phenology-1994-to-1995-prcp-DAYMET-004-results.csv

MAPS-Phenology-1995-to-1996-prcp-DAYMET-004-results.csv

MAPS-Phenology-1996-to-1997-prcp-DAYMET-004-results.csv

MAPS-Phenology-1997-to-1998-prcp-DAYMET-004-results.csv

MAPS-Phenology-1998-to-1999-prcp-DAYMET-004-results.csv

MAPS-Phenology-1999-to-2000-prcp-DAYMET-004-results.csv

MAPS-Phenology-2001-to-2002-prcp-DAYMET-004-results.csv

MAPS-Phenology-2002-to-2003-prcp-DAYMET-004-results.csv

MAPS-Phenology-2003-to-2004-prcp-DAYMET-004-results.csv

MAPS-Phenology-2004-to-2005-prcp-DAYMET-004-results.csv

MAPS-Phenology-2005-to-2006-prcp-DAYMET-004-results.csv

MAPS-Phenology-2006-to-2007-prcp-DAYMET-004-results.csv

MAPS-Phenology-2007-to-2008-prcp-DAYMET-004-results.csv

MAPS-Phenology-2008-to-2009-prcp-DAYMET-004-results.csv

MAPS-Phenology-2009-to-2010-prcp-DAYMET-004-results.csv

MAPS-Phenology-2010-to-2011-prcp-DAYMET-004-results.csv

MAPS-Phenology-2011-to-2012-prcp-DAYMET-004-results.csv

MAPS-Phenology-2012-to-2013-prcp-DAYMET-004-results.csv

MAPS-Phenology-2013-to-2014-prcp-DAYMET-004-results.csv

MAPS-Phenology-2014-to-2015-prcp-DAYMET-004-results.csv

MAPS-Phenology-2015-to-2016-prcp-DAYMET-004-results.csv

MAPS-Phenology-2016-to-2017-prcp-DAYMET-004-results.csv

MAPS-Phenology-2017-to-2018-prcp-DAYMET-004-results.csv

-------------------------------------------------------------------------------------------------------------------

**Folder: Nighttime Light Files**

**Folder Contents:** subfolder within Data, contains 2 folders, all files stored within subfolders

**Folder Description:** contains geospatial layers with information on annual average nighttime light for North America. There is one file per year. Files are organized into two folders as the files from certain years (2009, 2010, and 2011) differ slightly in their spatial extent from the other years. Files were obtained from Li et al. 2020. We used version 7 of the files (they are updated regularly with additional years) and downloaded 1992 through 2018.

**Note:** all files within these folders are too large to sync with the GitHub repository and need to be downloaded from the accompanying Dryad folder

-------------------------------------------------------------------------------------------------------------------

**Folder: 2009 to 2011**

**Folder Contents:** subfolder within Nighttime Light Files, 3 files, file format(s): .tif

**Folder Description:** contains geospatial data files with average annual light level values for North America for 2009, 2010 and 2011

**Note:** all files within this folders are too large to sync with the GitHub repository and need to be downloaded from the accompanying Dryad folder

**List of Files:**

Harmonized_DN_NTL_2009_calDMSP.tif*

Harmonized_DN_NTL_2010_calDMSP.tif*

Harmonized_DN_NTL_2011_calDMSP.tif*

-------------------------------------------------------------------------------------------------------------------

**Folder: Other years**

**Folder Contents:** subfolder within Nighttime Light Files, 24 files, file format(s): .tif

**Folder Description:** contains geospatial data files with average annual light level values for North America for 1992-2008, 2012-2018.

**Note:** all files within this folders are too large to sync with the GitHub repository and need to be downloaded from the accompanying Dryad folder

**List of Files:**

Harmonized_DN_NTL_1992_calDMSP.tif*

Harmonized_DN_NTL_1993_calDMSP.tif*

Harmonized_DN_NTL_1994_calDMSP.tif* 

Harmonized_DN_NTL_1995_calDMSP.tif* 

Harmonized_DN_NTL_1996_calDMSP.tif*

Harmonized_DN_NTL_1997_calDMSP.tif*

Harmonized_DN_NTL_1998_calDMSP.tif*

Harmonized_DN_NTL_1999_calDMSP.tif*

Harmonized_DN_NTL_2000_calDMSP.tif*

Harmonized_DN_NTL_2001_calDMSP.tif*

Harmonized_DN_NTL_2002_calDMSP.tif*

Harmonized_DN_NTL_2003_calDMSP.tif*

Harmonized_DN_NTL_2004_calDMSP.tif*

Harmonized_DN_NTL_2005_calDMSP.tif*

Harmonized_DN_NTL_2006_calDMSP.tif*

Harmonized_DN_NTL_2007_calDMSP.tif*

Harmonized_DN_NTL_2008_calDMSP.tif*

Harmonized_DN_NTL_2012_calDMSP.tif*

Harmonized_DN_NTL_2013_calDMSP.tif*

Harmonized_DN_NTL_2014_simVIIRS.tif*

Harmonized_DN_NTL_2015_simVIIRS.tif*

Harmonized_DN_NTL_2016_simVIIRS.tif*

Harmonized_DN_NTL_2017_simVIIRS.tif*

Harmonized_DN_NTL_2018_simVIIRS.tif*

############################################################################

### SCRIPTS

**Folder Contents:** 14 files, file format(s): .R

**Folder Description:** contains R scripts to analyze the data, to generate the findings in the manuscript, and to create the tables and results figures. The scripts are sequential, should be run in numerical order, and are labeled accordingly. At the beginning of each R script is a section titled “Objectives/Description of the Script” which details the main tasks accomplished by each file. We refer you to this section in lieu of providing an additional description of each file here.

**List of Files:**

Step1_RefineStationList.R

Step2_StandardizeStationEffort.R

Step3_FilterCaptureData.R

Step4_IdentifySpeciesStationsforModels

Step5_LightCalculation.R

Step6_CalculatingClimateVariables.R

Step7_CalculatingTStatistics.R

Step8_BayesianSpeciesModels.R

Step9_CompareModels.R

Step10_ModelValidation&Estimates.R

Step11_DWLWModelComparison.R

Step12_SpeciesFigures.R 

Step13_TraitAnalysis.R

Step14_PostHocAnalysis.R


############################################################################

### OUTPUTS 

**Folder Contents:** 8 files, file format(s): all .rds

**Folder Description:** Contains various files stored as .rds that contain organized data that is being moved between R scripts and used in subsequent analysis steps. All of these objects are generated using the scripts, so we have not described them individually here.

**List of Files:**

STA_finallist.rds

STA_Yr_IP_final.rds

capturedat_final.rds

HY_timeline.rds

combined_t_stats.rds

TAVGanomoly.rds

elpd_SE.rds

waic_SE.rds

############################################################################

### MODELS 

**Models Folder Structure (contains multiple subfolders):**

Models
  - Model Outputs
  - Model Plots
    
-------------------------------------------------------------------------------------------------------------------
**Folder:** Models

**Folder Contents:** 2 folders, all files stored within subfolders

**Folder Description:** contains models, their corresponding outputs, and diagnostic plots associated with the models organized into two subfolders. All files are generated by the accompanying R scripts. 

-------------------------------------------------------------------------------------------------------------------
**Folder: Model Outputs**

**Folder Contents:** subfolder within Models, 19 files, files format(s): .rds, .csv, .xlsx

**Folder Description:** contains the models and their corresponding outputs that are produced by the R scripts. As all files are generated by the accompanying code, we have not provided individual descriptions. Files are exported and saved as .rds, .csv or .xlsx files. All .xlsx files contain model summaries that have been formatted for inclusion in the accompanying publication. The .csv files contain the same information as the .xlsx files but have not been formatted for publication (they are the direct outputs from R). All .rds files contain saved versions of the models or their summary tables. 

**Note:** files with * at end of their name are too large to sync with the GitHub repository and need to be downloaded from the accompanying Dryad folder

**List of Files:**

dw_sppmodels.rds*

dw_sppmodels_k20.rds*

dw_sppmodels_noGP.rds*

lw_sppmodels.rds*

lw_sppmodels_k20.rds*

lw_sppmodels_noGP.rds*

dw_model_summaries.csv

dw_model_summaries.rds

lw_model_summaries.csv

lw_model_summaries.rds

dw_elpd_SE.rds

lw_elpd_SE.rds

traitmodels.csv

TraitModels_Formatted.xlsx

dwmodels_loo.rds*

lwmodels_loo.rds*

DWLW_elpd_waic.csv

DWLW_elpc_waic_Formatted.xlsx

DW&LWModelOutputs_Formatted.xlsx

-------------------------------------------------------------------------------------------------------------------

**Folder: Model Plots**

**Folder Contents:** subfolder within Models, 4 files, file format(s): .pdf

**Folder Description:** contains model diagnostic plots that are generated by the R scripts and exported as .pdf files.

**List of Files:**

dw_tracedensityplots.pdf

lw_tracedensityplots.pdf

dw_ppc_densityplots.pdf

lw_pcc_densityplots.pdf

############################################################################

### FIGURES

**Figures Folder Structure (contains 1 subfolder):**

Figures
  - Figure S2

-------------------------------------------------------------------------------------------------------------------
**Folder: Figures**

**Folder Contents:** 6 files, file format(s): .pdf. There is one subfolder called Figure S2

**Folder Description:** contains copies of the figures produced and included in the accompanying publication. Files are numbered to match the order of figures in the manuscript.

**List of Files:**

Fig2_ChangeStat.pdf

Fig3_TreePlusHeatMap.pdf

Fig4_TraitPanel.pdf

Fig5_PosthocPlot.pdf

FigS1_TStatHistograms.pdf

SPEC_STA_lag_plots.pdf

-------------------------------------------------------------------------------------------------------------------
**Folder: Figure S2**

**Folder Contents:** 20 files, file format(s): .pdf

**Folder Description:** Figure S2 includes one panel of plots for each of the 20 species. Each .pdf file in this folder is for one species. The four letter code in the filename reflects the species and the common name for each is printed at the top of each panel (visible upon opening the files).

**List of Files:**

AMRO_FigS2.pdf

AUWA_FigS2.pdf

BCCH_FigS2.pdf

BEWR_FigS2.pdf

CARW_FigS2.pdf

COYE_FigS2.pdf

DOWO_FigS2.pdf

GRCA_FigS2.pdf

HOWR_FigS2.pdf

LISP_FigS2.pdf

MGWA_FigS2.pdf

NOCA_FigS2.pdf

ORJU_FigS2.pdf

OVEN_FigS2.pdf

SOSP_FigS2.pdf

SWTH_FigS2.pdf

TUTI_FigS2.pdf

WIWA_FigS2.pdf

WOTH_FigS2.pdf

YEWA_FigS2.pdf

############################################################################
### REFERENCES
############################################################################

AppEEARS Team. 2026. Application for Extracting and Exploring Analysis Ready Samples (AppEEARS). Ver. 3.111. NASA EOSDIS Land Processes Distributed Active Archive Center (LP DAAC), USGS/Earth Resources Observation and Science (EROS) Center, Sioux Falls, South Dakota, USA. https://appeears.earthdatacloud.nasa.gov

Ehrlich, PR, Dobkin, DS, & Wheye, D. 1988. The birder’s handbook: A field guide to the natural history of North American birds. Simon & Schuster

Jennings, SL, EM Moylan, CD Francis. In revision. Avian traits link divergent population responses to a changing world.

Jetz, W, GH Thomas, JB Joy, K Hartmann, and AO Mooers. 2012. The Global Diversity of Birds in Space and Time. Nature 491 (7424): 444–48. https://doi.org/10.1038/nature11631

Lehikoinen, A, Lindström, A Santangeli, PM Sirkiä, L Brotons, V Devictor, J Elts, RPB Foppen, H Heldbjerg, S Herrange, M Herremans, MAR Hudson, F Jiguet, A Johnston, R Lorrilliere, EL Marjakangas, NL Michel, CM Moshøj, R Nellis, JV Paquet, AC Smith, T Szép, C van Turnhout. 2021. Wintering bird communities are tracking climate change faster than breeding communities. Journal of Animal Ecology 90:1085-1095. https://doi.org/10.1111/1365-2656.13433

Li, X, Zhou, Y, Zhao, M. 2020. A harmonized global nighttime light dataset 1992-2018. Scientific data 7: 168. https://doi.org/10.1038/s41597-020-0510-y
  - Note: We used version 7 of the associated data files: https://doi.org/10.6084/m9.figshare.9828827.v7

Rohde, RA and Z Hausfather. 2020. The Berkeley Earth Land/Ocean Temperature Record, Earth Syst. Sci. Data, 12, 3469–3479. https://doi.org/10.5194/essd-12-3469-2020

Ritland, SM. 1982. The Allometry of the Vertebrate Eye. The University of Chicago. https://www.proquest.com/docview/303088065?pq-origsite=gscholar&fromopenview=true&sourcetype=Dissertations%20&%20Theses.

Senzaki, M, JR Barber, JN Phillips, NH Carter, CB Cooper, MA Ditmer, KM Fristrup, CJW McClure, DJ Mennitt, LP Tyrell, J Vukomanovic, AA Wilson, CD Francis. 2020. Sensory pollutants alter bird phenology and fitness across a continent. Nature 587(7835), 605–610. https://doi.org/10.1038/s41586-020-2903-7 

Sheard, C, MHC Neate-Clegg, N Alioravainen, SEI Jones, C Vincent, HEA MacGregor, TP Bregman, S Claramunt, JA Tobias. 2020. Ecological Drivers of Global Gradients in Avian Dispersal Inferred from Wing Morphology. Nature Communications 11:2463. https://doi.org/10.1038/s41467-020-16313-6

Thornton, MM, R Shrestha, Y Wei, PE Thornton, and S-C Kao. “Daymet: Daily Surface Weather Data on a 1-Km Grid for North America, Version 4 R1.” ORNL Distributed Active Archive Center, January 1, 2022. doi:10.3334/ORNLDAAC/2129

Wilson, AA, MA Ditmer, JR Barber, NH Carter, ET Miller, LP Tyrrell, CD Francis. 2021. Artificial Night Light and Anthropogenic Noise Interact to Influence Bird Abundance Over a Continental Scale. Global Change Biology 27:3987-4004. https://doi.org/10.1111/gcb.15663
Wolf, MM & CD Francis. 2025. Eye Catching Light: Anthropogenic Light at Night and its Evolutionary Influence on the Avian Eye. iScience 28:112039. https://doi.org/10.1016/j.isci.2025.112039

