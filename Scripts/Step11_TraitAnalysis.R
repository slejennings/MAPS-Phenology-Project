###### MAPS Phenology #######
### Script name: Step11_TraitAnalyses.R
## Author(s): LP, CDF, SLJ

########## Objective/Description of Script #####################
# this script compares the long-term and decision window model for each species using the loo package
# we do this using two model comparison metrics: ELPD and waic
#################################################################

#### Setup ####

# load packages
library(here)
library(performance)
library(tidyverse)
library(ggeffects)
library(ape)
library(geiger)
library(nlme)
library(patchwork)
library(broom.mixed)

#### Import files

# decision window models
dw_models <- readRDS(here("Models/Model Outputs", "dw_model_summaries.rds"))

# long-term window models
lw_models <- readRDS(here("Models/Model Outputs", "lw_model_summaries.rds"))

# eye morphometrics
eye <- read.csv(here("Data", "species_eyes.csv"), header=T)

# species temperature indices
STI <- read.csv(here("Data", "species_STI.csv"), header=T)

# hand wing index (HWI), body mass, annual precip in range
traits <- read.csv(here("Data", "species_traits.csv"), header=T)

# preferred model (decision vs long-term) using ELPD
elpd <- readRDS(here("Outputs", "elpd_SE.rds")) %>% 
  select(SPEC, elpd_diff_model, SE_diff) %>%
  rename(elpd_SE_diff = SE_diff)

# preferred model (decision vs long-term) using WAIC # POSSIBLY DELETE THIS
waic <- readRDS(here("Outputs", "waic_SE.rds")) %>% 
  select(SPEC, waic_diff_model)

# phylogenetic tree for birds
tree <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

####################################################################################
#### Combine imported files to get data frame to use for trait models

# first, combine decision and long-term window models
# before doing this, we need to calculate t-statistics for the estimates generated in these models
# we will do this using the estimate divided by the standard error
# we also need to give the light t-statistic a unique name to differentiate the DW vs LW model

dw_models_wide <- dw_models %>%
  filter(Parameter %in% c(
    "scaletempanom_DW_tstat",
    "scaleprcp_DW_total_tstat",
    "scaleprcp_DW_cov_tstat",
    "scalelight_tstat"
  )) %>%
  rowwise() %>%
  mutate(tvalue = Estimate/Est.Error) %>%
  select(SPEC, Parameter, tvalue) %>%
  pivot_wider(
    names_from = Parameter,   # new column names from "Parameter"
    values_from = tvalue) %>%  # values to fill from "tvalue"
rename(scalelight_DW_tstat = scalelight_tstat) # give light t-statistic a unique name

lw_models_wide <- lw_models %>%
  filter(Parameter %in% c(
    "scaletempanom_LW_tstat",
    "scaleprcp_LW_total_tstat",
    "scaleprcp_LW_cov_tstat",
    "scalelight_tstat"
  )) %>%
  rowwise() %>%
  mutate(tvalue = Estimate/Est.Error) %>%
  select(SPEC, Parameter, tvalue) %>%
  pivot_wider(
    names_from = Parameter,   # new column names from "Parameter"
    values_from = tvalue) %>%   # values to fill from "tvalue"
  rename(scalelight_LW_tstat = scalelight_tstat) # give light t-statistic a unique name

# combine t-statistics from both long-term and decision window models
all_models_wide <- left_join(dw_models_wide, lw_models_wide)

# next, add traits and model preference scores
all_models_traits <- left_join(all_models_wide, eye, by="SPEC") %>% 
  select(-Species, -Source, -Museum.Collection.IDs) %>%
  left_join(., STI, by="SPEC") %>%
  select(-Species) %>%
  left_join(., traits, by="SPEC") %>%
  select(-Species_name) %>%
  rowwise() %>% mutate(AnnualPrecip_cm = AnnualPrecip/10) %>% # convert annual precip values to cm
  left_join(., elpd) %>%
  left_join(., waic) # POSSIBLY DELETE

# for most traits, we have values for all 20 species
# move the Tree_name column to rownames to facilitate pairing data with phylogenetic tree
traits_20spp <- all_models_traits %>%
  mutate(Species = str_replace(Tree_name, "_", " ")) %>% 
  column_to_rownames(., var = "Tree_name")

# for eye morphometrics, we only have values for 18 species
# move the Tree_name column to rownames to facilitate pairing data with phylogenetic tree
eye_18spp <- all_models_traits %>%
  filter(!is.na(C.T)) %>%
  mutate(Species = str_replace(Tree_name, "_", " ")) %>%
  column_to_rownames(., var = "Tree_name")


################################################################################
#### Pair data with phylogenetic tree

#### For all traits with 20 species ####
phydat_20spp <- geiger::treedata(tree, traits_20spp, sort=T) # join tree with data

birdtree_20spp <- phydat_20spp$phy # this is our trimmed tree for the 20 species

# these are the data associated with our trimmed tree
# NOTE: IF NOT USING WAIC, CHANGE THE NAME OF THE LAST COLUMN INSIDE MUTATE TO elpd_SE_diff
birddat_20spp <- as.data.frame(phydat_20spp$data) %>% # convert to df
  mutate(across(c(scaletempanom_DW_tstat:waic_diff_model), as.numeric)) # make sure columns that need to be are numeric
birddat_20spp$species2 <- rownames(birddat_20spp)

#### For 18 species with eye morphometrics ####
phydat_18spp <- geiger::treedata(tree, eye_18spp, sort=T) # join tree with data

birdtree_18spp <- phydat_18spp$phy # this is our trimmed tree for the 18 species

# these are the data associated with our trimmed tree
# NOTE: IF NOT USING WAIC, CHANGE THE NAME OF THE LAST COLUMN INSIDE MUTATE TO elpd_SE_diff
birddat_18spp <- as.data.frame(phydat_18spp$data) %>% # convert to df
  mutate(across(c(scaletempanom_DW_tstat:waic_diff_model), as.numeric)) # make sure columns that need to be are numeric

birddat_18spp$species2 <- rownames(birddat_18spp)

####################################################################################
##### Trait models with STI (species temperature index)
####################################################################################

# use birdtree_20spp as the tree and birddat_20spp for the data in these models
colnames(birddat_20spp)

######## Temperature anomaly in decision window ##############

# fixed lambda at zero because estimated at zero using phylolm
temp_STI_dw <- gls(scaletempanom_DW_tstat ~ STI + Body_mass_log, 
          data = birddat_20spp, 
          correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

# no relationships
summary(temp_STI_dw)
check_model(temp_STI_dw)

######## Temperature anomaly in long-term window ##############

# fixed lambda at zero because estimated at zero using phylolm
temp_STI_lw <- gls(scaletempanom_LW_tstat ~ STI + Body_mass_log, 
                   data = birddat_20spp, 
                   correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

# positive relationship between STI and response
summary(temp_STI_lw)
check_model(temp_STI_lw)

####################################################################################
##### Trait models with annual precipitation in range
####################################################################################

# use birdtree_20spp as the tree and birddat_20spp for the data in these models
colnames(birddat_20spp)

######## Total precipitation in decision window ##############

# # fixed lambda at zero because estimated at zero using phylolm
totalprcp_annualprcp_dw <- gls(scaleprcp_DW_total_tstat ~ AnnualPrecip_cm + Body_mass_log, 
                   data = birddat_20spp, 
                   correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

# no relationships
summary(totalprcp_annualprcp_dw)
check_model(totalprcp_annualprcp_dw)

######## Total precipitation in long-term window ##############

# # fixed lambda at zero because estimated at zero using phylolm
totalprcp_annualprcp_lw <- gls(scaleprcp_LW_total_tstat ~ AnnualPrecip_cm + Body_mass_log, 
                               data = birddat_20spp, 
                               correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form = ~species2), method = "ML")

summary(totalprcp_annualprcp_lw)
check_model(totalprcp_annualprcp_lw)

######## Precipitation variability in decision window ##############

# # fixed lambda at zero because estimated at zero using phylolm
covprcp_annualprcp_dw <- gls(scaleprcp_DW_cov_tstat ~ AnnualPrecip_cm + Body_mass_log, 
                               data = birddat_20spp, 
                               correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

#no relationships
summary(covprcp_annualprcp_dw)
check_model(covprcp_annualprcp_dw)

######## Precipitation variability in long-term window ##############

# # fixed lambda at zero because estimated at zero using phylolm
covprcp_annualprcp_lw <- gls(scaleprcp_LW_cov_tstat ~ AnnualPrecip_cm + Body_mass_log, 
                             data = birddat_20spp, 
                             correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

#no relationships
summary(covprcp_annualprcp_lw)
check_model(covprcp_annualprcp_lw)

####################################################################################
##### Trait models with HWI (hand wing index)
####################################################################################

# use birdtree_20spp as the tree and birddat_20spp for the data in these models
colnames(birddat_20spp)

######## Temperature anomaly in decision window ##############

# # fixed lambda at zero because estimated at zero using phylolm
temp_HWI_dw <- gls(scaletempanom_DW_tstat ~ HWI + Body_mass_log, 
                   data = birddat_20spp, 
                   correlation = corPagel(0, phy = birdtree_20spp, fixed=T,form =~species2), method = "ML")

#no relationships
summary(temp_HWI_dw)
check_model(temp_HWI_dw)

######## Temperature anomaly in long-term window ##############

# fixed lambda at zero because estimated at zero using phylolm
temp_HWI_lw <- gls(scaletempanom_LW_tstat ~ HWI + Body_mass_log, 
                   data = birddat_20spp, 
                   correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

#no relationships
summary(temp_HWI_lw)
check_model(temp_HWI_lw)

######## Preferred model (decision vs long-term) using ELPD ##############

# using the difference in SE for elpd as a weighting variable in this model

# fixed lambda at zero because estimated at zero using phylolm
elpddiff_HWI <- gls(elpd_diff_model ~ HWI, 
                    data = birddat_20spp, 
                    correlation = corPagel(0, phy = birdtree_20spp, fixed=T,form=~species2),
                    weights = varFixed(~1/sqrt(elpd_SE_diff)), method = "ML")

# no relationship
summary(elpddiff_HWI)
check_model(elpddiff_HWI)

######## Preferred model (decision vs long-term) using WAIC ##############

######## POSSIBLY DELETE THIS AND ONLY USE ELPD MODEL

# no weighting variable in this model as we can't get the difference in SE of waic (see notes in Step9 script)

# fixed lambda at zero because estimated at zero using phylolm
waicdiff_HWI <- gls(waic_diff_model ~ HWI, 
                   data = birddat_20spp, 
                   correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

#no relationship
summary(waicdiff_HWI)
check_model(waicdiff_HWI)

####################################################################################
##### Trait models with C.T (eye morphology)
####################################################################################

# use birdtree_18spp as the tree and birddat_18spp for the data in these models
colnames(birddat_18spp)


######## Light pollution in decision window ##############

# fixed lambda at zero because estimated at zero using phylolm
light_CT_dw <- gls(scalelight_DW_tstat ~ C.T + Body_mass_log, 
                   data = birddat_18spp, 
                   correlation = corPagel(0, phy = birdtree_18spp, fixed=T, form =~species2), method = "ML")

#no relationship
summary(light_CT_dw)
check_model(light_CT_dw)

######## Light pollution in long-term window ##############

# fixed lambda at zero because estimated at zero using phylolm
light_CT_lw <- gls(scalelight_LW_tstat ~ C.T + Body_mass_log, 
                   data = birddat_18spp, 
                   correlation = corPagel(0, phy = birdtree_18spp, fixed=T, form=~species2), method = "ML")

summary(light_CT_lw)
check_model(light_CT_lw)


####################################################################################
##### Fig 5: Panel of trait plots
####################################################################################

###### Panel A: using temp_STI_lw ##########

# get effect of STI on the model statistics for change in temperature anomalies
eff_temp_STI <- plot(ggeffects::predict_response(temp_STI_lw, terms =c("STI")), colors = "darkseagreen")

# add data, labels, nice formatting to plot
Panel_A <- eff_temp_STI +
  geom_point(data = birddat_20spp, aes(x = STI, y = scaletempanom_LW_tstat), color = "darkseagreen", size = 3.1, pch = 19)+
  labs(title= "", x = "Species temperature index (STI)", y = "Sensitivity to temperature anomalies")+
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), 
        axis.title.x = element_text(color = "black", size = 12), axis.title.y = element_text(color = "black", size = 12)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "gray", linewidth=.6)

Panel_A     

###### Panel B: using C.T ##########

# get effect of C.T on the model statistics for change in light pollution
eff_light_CT <- plot(ggeffects::predict_response(light_CT_lw, terms =c("C.T")), colors = "#DEB70D")

# add data, labels, nice formatting to plot
Panel_B <- eff_light_CT +
  geom_point(data = birddat_18spp, aes(x = C.T, y = scalelight_LW_tstat), color = "#DEB70D", size = 3.1, pch = 19)+
  labs(title= "", x = "Dim light vision", y = "Sensitivity to light pollution")+
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), 
        axis.title.x = element_text(color = "black", size = 12), axis.title.y = element_text(color = "black", size = 12)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "gray", linewidth=.6)

Panel_B

# combine panels and export
toprow <- (Panel_A + plot_spacer() + Panel_B) + plot_layout(widths = c(0.49, 0.01, 0.49)) # add a small space between the two plots


panelplot <- toprow +
             plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 14, face ="bold"))

ggsave(panelplot, filename = "Fig5_TraitPanel.pdf", path = here("Figures"), width=20, height=10, units = "cm")

####################################################################################
##### Fig 4: Model fit and dispersal ability (HWI)
####################################################################################

# get effect of HWI on the model fit (elpd)
eff_elpddiff_HWI <- plot(ggeffects::predict_response(elpddiff_HWI, terms =c("HWI")), colors = "#0B67A1")

# make plot with errorbars
Fig4 <- eff_elpddiff_HWI +
  geom_point(data = birddat_20spp, aes(x = HWI, y = elpd_diff_model), color = "#0B67A1", size = 3.1, pch = 19) +
  geom_errorbar(aes(ymin = (elpd_diff_model-SE_diff), ymax = (elpd_diff_model+SE_diff), x=STI),
                color="#0B67A1", data = combDAT,inherit.aes = FALSE, lwd=.9) +
  labs(title= "", x = "Hand-wing index", y = "Model strength")+
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), 
        axis.title.x = element_text(color = "black", size = 12), axis.title.y = element_text(color = "black", size = 12)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "gray", linewidth=.6)

# save plot
ggsave(HWI2, filename = "Fig4_HWI.pdf", path = here("Figures"), width=20, height=12, units = "cm")
ggsave(HWI2, filename = "Fig4_HWI.jpg", path = here("Figures"), width=20, height=12, units = "cm")


####################################################################################
##### Tidy models and save estimates
####################################################################################

# tidy models and save estimates

temp_STI_dw_tidy <- tidy(temp_STI_dw, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "STI",
         response_variable = "temperature anomaly",
         window = "decision",
         sample_size = temp_STI_dw$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

temp_STI_lw_tidy <- tidy(temp_STI_lw, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "STI",
         response_variable = "temperature anomaly",
         window = "long-term",
         sample_size = temp_STI_lw$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

totalprcp_annualprcp_dw_tidy <- tidy(totalprcp_annualprcp_dw, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "annual precip in range",
         response_variable = "total precipitation",
         window = "decision",
         sample_size = totalprcp_annualprcp_dw$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

totalprcp_annualprcp_lw_tidy <- tidy(totalprcp_annualprcp_lw, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "annual precip in range",
         response_variable = "total precipitation",
         window = "long-term",
         sample_size = totalprcp_annualprcp_lw$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

covprcp_annualprcp_dw_tidy <- tidy(covprcp_annualprcp_dw, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "annual precip in range",
         response_variable = "precipitation COV",
         window = "decision",
         sample_size = covprcp_annualprcp_dw$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

covprcp_annualprcp_lw_tidy <- tidy(covprcp_annualprcp_lw, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "annual precip in range",
         response_variable = "precipitation COV",
         window = "long-term",
         sample_size = covprcp_annualprcp_lw$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

temp_HWI_dw_tidy <- tidy(temp_HWI_dw, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "HWI",
         response_variable = "temperature anomaly",
         window = "decision",
         sample_size = temp_HWI_dw$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

temp_HWI_lw_tidy <- tidy(temp_HWI_lw, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "HWI",
         response_variable = "temperature anomaly",
         window = "long-term",
         sample_size = temp_HWI_lw$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

elpddiff_HWI_tidy <- tidy(elpddiff_HWI, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "HWI",
         response_variable = "model preference (elpd)",
         window = "NA",
         sample_size = elpddiff_HWI$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

waicdiff_HWI_tidy <- tidy(waicdiff_HWI, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "HWI",
         response_variable = "model preference (waic)",
         window = "NA",
         sample_size = waicdiff_HWI$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

light_CT_dw_tidy <- tidy(light_CT_dw, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "C.T",
         response_variable = "light pollution",
         window = "NA",
         sample_size = light_CT_dw$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

light_CT_lw_tidy <- tidy(light_CT_lw, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "C.T",
         response_variable = "light pollution",
         window = "NA",
         sample_size = light_CT_lw$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

# combine into one single table of results to export
trait_results <- bind_rows(temp_STI_dw_tidy, temp_STI_lw_tidy, totalprcp_annualprcp_dw_tidy, totalprcp_annualprcp_lw_tidy, 
                           covprcp_annualprcp_dw_tidy, covprcp_annualprcp_lw_tidy, temp_HWI_dw_tidy, temp_HWI_lw_tidy, 
                           elpddiff_HWI_tidy, waicdiff_HWI_tidy, light_CT_dw_tidy, light_CT_lw_tidy) %>%
  select(trait, response_variable, window, sample_size, term:conf.high, lambda)

write.csv(trait_results, here("Models/Model Outputs", "traitmodels.csv"))



