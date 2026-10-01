###### MAPS Phenology #######
### Script name: Step11_TraitAnalyses.R
## Author(s): XXX (removed for peer review)

########## Objective/Description of Script #####################
# run trait-based phylogenetic models
# create Figure 4
# export model outputs
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
dw_models <- readRDS(here("Models", "dw_model_summaries.rds"))

# long-term window models
lw_models <- readRDS(here("Models", "lw_model_summaries.rds"))

# eye morphometrics
eye <- read.csv(here("Data", "species_eyes.csv"), header=T)

# species temperature indices
STI <- read.csv(here("Data", "species_STI.csv"), header=T)

# hand wing index (HWI), body mass, annual precip in range
traits <- read.csv(here("Data", "species_traits.csv"), header=T)

# phylogenetic tree for birds
tree <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

####################################################################################
#### Combine imported files to get data frame to use for trait models

# first, combine decision and long-term window models
# before doing this, we need to calculate t-statistics for the estimates generated in these models
# we will do this using the estimate divided by the standard error
# we also need to give the light t-statistic a unique name to differentiate the DW vs LW model

dw_models_wide <- dw_models %>%
  filter(term %in% c(
    "scale(tempanom_DW_tstat)",
    "scale(prcp_DW_total_tstat)",
    "scale(prcp_DW_cov_tstat)",
    "scale(light_tstat)"
  )) %>%
  rowwise() %>%
  mutate(tvalue = estimate/std.error) %>%
  select(Species, term, tvalue) %>%
  pivot_wider(
    names_from = term,   # new column names from "Parameter"
    values_from = tvalue) %>%  # values to fill from "tvalue"
rename(scalelight_DW_tstat = 'scale(light_tstat)',
       scaletempanom_DW_tstat = 'scale(tempanom_DW_tstat)',
       scaleprcp_DW_total_tstat = 'scale(prcp_DW_total_tstat)',
       scaleprcp_DW_cov_tstat = 'scale(prcp_DW_cov_tstat)'
       ) # give light t-statistic a unique name

lw_models_wide <- lw_models %>%
  filter(term %in% c(
    "scale(tempanom_LW_tstat)",
    "scale(prcp_LW_total_tstat)",
    "scale(prcp_LW_cov_tstat)",
    "scale(light_tstat)"
  )) %>%
  rowwise() %>%
  mutate(tvalue = estimate/std.error) %>%
  select(Species, term, tvalue) %>%
  pivot_wider(
    names_from = term,   # new column names from "Parameter"
    values_from = tvalue) %>%  # values to fill from "tvalue"
  rename(scalelight_LW_tstat = 'scale(light_tstat)',
         scaletempanom_LW_tstat = 'scale(tempanom_LW_tstat)',
         scaleprcp_LW_total_tstat = 'scale(prcp_LW_total_tstat)',
         scaleprcp_LW_cov_tstat = 'scale(prcp_LW_cov_tstat)'
  )
# combine t-statistics from both long-term and decision window models
all_models_wide <- left_join(dw_models_wide, lw_models_wide) %>%
  rename(SPEC = Species)

# next, add traits and model preference scores
all_models_traits <- left_join(all_models_wide, eye, by="SPEC") %>% 
  select(-Species, -Source, -Museum.Collection.IDs) %>%
  left_join(., STI, by="SPEC") %>%
  select(-Species) %>%
  left_join(., traits, by="SPEC") %>%
  select(-Species_name) %>%
  rowwise() %>% 
  mutate(AnnualPrecip_cm = AnnualPrecip/10) # convert annual precip values to cm
      
# for most traits, we have values for all 20 species
# move the Tree_name column to rownames to facilitate pairing data with phylogenetic tree
traits_20spp <- all_models_traits %>%
  mutate(Species = str_replace(Tree_name, "_", " ")) %>% 
  column_to_rownames(., var = "Tree_name")

# for eye morphometrics, we only have values for 19 species
# move the Tree_name column to rownames to facilitate pairing data with phylogenetic tree
eye_19spp <- all_models_traits %>%
  filter(!is.na(C.T)) %>%
  mutate(Species = str_replace(Tree_name, "_", " ")) %>%
  column_to_rownames(., var = "Tree_name")

################################################################################
#### Pair data with phylogenetic tree

#### For all traits with 20 species ####
phydat_20spp <- geiger::treedata(tree, traits_20spp, sort=T) # join tree with data

birdtree_20spp <- phydat_20spp$phy # this is our trimmed tree for the 20 species

# these are the data associated with our trimmed tree
birddat_20spp <- as.data.frame(phydat_20spp$data) %>% # convert to df
  mutate(across(c(scalelight_DW_tstat:AnnualPrecip_cm), as.numeric)) # make sure columns that need to be are numeric
birddat_20spp$species2 <- rownames(birddat_20spp)

#### For 19 species with eye morphometrics ####
phydat_19spp <- geiger::treedata(tree, eye_19spp, sort=T) # join tree with data

birdtree_19spp <- phydat_19spp$phy # this is our trimmed tree for the 19 species

# these are the data associated with our trimmed tree
birddat_19spp <- as.data.frame(phydat_19spp$data) %>% # convert to df
  mutate(across(c(scalelight_DW_tstat:AnnualPrecip_cm), as.numeric)) # make sure columns that need to be are numeric

birddat_19spp$species2 <- rownames(birddat_19spp)

####################################################################################
##### Trait models with STI (species temperature index)
####################################################################################

# use birdtree_20spp as the tree and birddat_20spp for the data in these models
colnames(birddat_20spp)

######## Temperature anomaly in decision window ##############

# fixed lambda at zero because estimated at zero using phylolm
# run model with and without body mass and compare them
temp_STI_dw <- gls(scaletempanom_DW_tstat ~ STI + Body_mass_log, 
          data = birddat_20spp, 
          correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

temp_STI_dw2 <- gls(scaletempanom_DW_tstat ~ STI, 
                   data = birddat_20spp, 
                   correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

anova(temp_STI_dw, temp_STI_dw2) # no advantage of model with body mass

# no relationships
summary(temp_STI_dw2)
check_model(temp_STI_dw2)
confint(temp_STI_dw2)


######## Temperature anomaly in long-term window ##############

# fixed lambda at zero because estimated at zero using phylolm
# run model with and without body mass and compare them
temp_STI_lw <- gls(scaletempanom_LW_tstat ~ STI + Body_mass_log, 
                   data = birddat_20spp, 
                   correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

temp_STI_lw2 <- gls(scaletempanom_LW_tstat ~ STI, 
                   data = birddat_20spp, 
                   correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

anova(temp_STI_lw, temp_STI_lw2) # no advantage of model with body mass

# positive relationship between STI and response
summary(temp_STI_lw2)
confint(temp_STI_lw2)
check_model(temp_STI_lw2)

####################################################################################
##### Trait models with annual precipitation in range
####################################################################################

# use birdtree_20spp as the tree and birddat_20spp for the data in these models
colnames(birddat_20spp)

######## Total precipitation in decision window ##############

# fixed lambda at zero because estimated at zero using phylolm
# run model with and without body mass and compare them
totalprcp_annualprcp_dw <- gls(scaleprcp_DW_total_tstat ~ AnnualPrecip_cm + Body_mass_log, 
                               data = birddat_20spp, 
                               correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

totalprcp_annualprcp_dw2 <- gls(scaleprcp_DW_total_tstat ~ AnnualPrecip_cm, 
                               data = birddat_20spp, 
                               correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")


anova(totalprcp_annualprcp_dw, totalprcp_annualprcp_dw2) # no support for model with body mass

# no relationships
summary(totalprcp_annualprcp_dw2)
confint(totalprcp_annualprcp_dw2)
check_model(totalprcp_annualprcp_dw2)

######## Total precipitation in long-term window ##############

# # fixed lambda at zero because estimated at zero using phylolm
totalprcp_annualprcp_lw <- gls(scaleprcp_LW_total_tstat ~ AnnualPrecip_cm + Body_mass_log, 
                               data = birddat_20spp, 
                               correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form = ~species2), method = "ML")

totalprcp_annualprcp_lw2 <- gls(scaleprcp_LW_total_tstat ~ AnnualPrecip_cm, 
                               data = birddat_20spp, 
                               correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form = ~species2), method = "ML")

anova(totalprcp_annualprcp_lw, totalprcp_annualprcp_lw2) #  no support for model with body mass

summary(totalprcp_annualprcp_lw2)
confint(totalprcp_annualprcp_lw2)
check_model(totalprcp_annualprcp_lw2)

######## Precipitation variability in decision window ##############

# # fixed lambda at zero because estimated at zero using phylolm
covprcp_annualprcp_dw <- gls(scaleprcp_DW_cov_tstat ~ AnnualPrecip_cm + Body_mass_log, 
                               data = birddat_20spp, 
                               correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

covprcp_annualprcp_dw2 <- gls(scaleprcp_DW_cov_tstat ~ AnnualPrecip_cm, 
                             data = birddat_20spp, 
                             correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

anova(covprcp_annualprcp_dw, covprcp_annualprcp_dw2) # no support for model with body mass

#no relationships
summary(covprcp_annualprcp_dw2)
confint(covprcp_annualprcp_dw2)
check_model(covprcp_annualprcp_dw2)

######## Precipitation variability in long-term window ##############

# # fixed lambda at zero because estimated at zero using phylolm
# run model with and without body mass and compare them
covprcp_annualprcp_lw <- gls(scaleprcp_LW_cov_tstat ~ AnnualPrecip_cm + Body_mass_log, 
                             data = birddat_20spp, 
                             correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

covprcp_annualprcp_lw2 <- gls(scaleprcp_LW_cov_tstat ~ AnnualPrecip_cm, 
                             data = birddat_20spp, 
                             correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

anova(covprcp_annualprcp_lw, covprcp_annualprcp_lw2) # no support for model with body mass

#no relationships
summary(covprcp_annualprcp_lw2)
confint(covprcp_annualprcp_lw2)
check_model(covprcp_annualprcp_lw2)

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

temp_HWI_dw2 <- gls(scaletempanom_DW_tstat ~ HWI, 
                   data = birddat_20spp, 
                   correlation = corPagel(0, phy = birdtree_20spp, fixed=T,form =~species2), method = "ML")

anova(temp_HWI_dw, temp_HWI_dw2) # no support for model with body mass

#no relationships
summary(temp_HWI_dw2)
confint(temp_HWI_dw2)
check_model(temp_HWI_dw)

######## Temperature anomaly in long-term window ##############

# fixed lambda at zero because estimated at zero using phylolm
temp_HWI_lw <- gls(scaletempanom_LW_tstat ~ HWI + Body_mass_log, 
                   data = birddat_20spp, 
                   correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

temp_HWI_lw2 <- gls(scaletempanom_LW_tstat ~ HWI, 
                   data = birddat_20spp, 
                   correlation = corPagel(0, phy = birdtree_20spp, fixed=T, form =~species2), method = "ML")

anova(temp_HWI_lw, temp_HWI_lw2) # no support for model with body mass

#no relationships
summary(temp_HWI_lw2)
confint(temp_HWI_lw2)
check_model(temp_HWI_lw2)


####################################################################################
##### Trait models with C.T (eye morphology)
####################################################################################

# use birdtree_19spp as the tree and birddat_19spp for the data in these models
colnames(birddat_19spp)

######## Light pollution in decision window ##############

# fixed lambda at zero because estimated at zero using phylolm
light_CT_dw <- gls(scalelight_DW_tstat ~ C.T + Body_mass_log, 
                   data = birddat_19spp, 
                   correlation = corPagel(0, phy = birdtree_19spp, fixed=T, form =~species2), method = "ML")

light_CT_dw2 <- gls(scalelight_DW_tstat ~ C.T , 
                   data = birddat_19spp, 
                   correlation = corPagel(0, phy = birdtree_19spp, fixed=T, form =~species2), method = "ML")

anova(light_CT_dw, light_CT_dw2)  # no support for model with body mass

#no relationship
summary(light_CT_dw2)
confint(light_CT_dw2)
check_model(light_CT_dw2)


######## Light pollution in long-term window ##############

# fixed lambda at zero because estimated at zero using phylolm
light_CT_lw <- gls(scalelight_LW_tstat ~ C.T + Body_mass_log, 
                   data = birddat_19spp, 
                   correlation = corPagel(0, phy = birdtree_19spp, fixed=T, form=~species2), method = "ML")

light_CT_lw2 <- gls(scalelight_LW_tstat ~ C.T, 
                   data = birddat_19spp, 
                   correlation = corPagel(0, phy = birdtree_19spp, fixed=T, form=~species2), method = "ML")

anova(light_CT_lw, light_CT_lw2) # no support for model with body mass

summary(light_CT_lw2)
confint(light_CT_lw2)
check_model(light_CT_lw2)


####################################################################################
##### Fig 4:STI trait plot
####################################################################################

###### using temp_STI_lw ##########

# get effect of STI on the model statistics for change in temperature anomalies
eff_temp_STI <- plot(ggeffects::predict_response(temp_STI_lw, terms =c("STI")), colors = "#699864")

# add data, labels, nice formatting to plot
Fig4 <- eff_temp_STI +
  geom_point(data = birddat_20spp, aes(x = STI, y = scaletempanom_LW_tstat), color = "#699864", size = 3.1, pch = 19)+
  labs(title= "", x = "Species temperature index (STI)", y = "Phenological responsiveness \n to temperature anomalies")+
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), 
        axis.title.x = element_text(color = "black", size = 12, margin = margin(t=0.3, unit="cm")), # add space between axis title and axis labels
        axis.title.y = element_text(color = "black", size = 12, margin = margin(r=0.3, unit="cm"))) + # add space between axis title and axis labels
  geom_hline(yintercept=0, linetype="dashed", color = "gray", linewidth=.6)

Fig4     



ggsave(Fig4, filename = "Fig4_TraitPanel.pdf", path = here("Figures"), width=11, height=10, units = "cm", device=cairo_pdf)
#ggsave(Fig4, filename = "Fig4_TraitPanel.png", path = here("Figures"), width=11, height=10, units = "cm")

####################################################################################
##### Tidy models and save estimates
####################################################################################

# tidy models and save estimates

temp_STI_dw_tidy <- tidy(temp_STI_dw2, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "STI",
         response_variable = "temperature anomaly",
         window = "decision",
         sample_size = temp_STI_dw2$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

temp_STI_lw_tidy <- tidy(temp_STI_lw2, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "STI",
         response_variable = "temperature anomaly",
         window = "long-term",
         sample_size = temp_STI_lw2$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

totalprcp_annualprcp_dw_tidy <- tidy(totalprcp_annualprcp_dw2, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "annual precip in range",
         response_variable = "total precipitation",
         window = "decision",
         sample_size = totalprcp_annualprcp_dw2$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

totalprcp_annualprcp_lw_tidy <- tidy(totalprcp_annualprcp_lw2, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "annual precip in range",
         response_variable = "total precipitation",
         window = "long-term",
         sample_size = totalprcp_annualprcp_lw2$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

covprcp_annualprcp_dw_tidy <- tidy(covprcp_annualprcp_dw2, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "annual precip in range",
         response_variable = "precipitation COV",
         window = "decision",
         sample_size = covprcp_annualprcp_dw2$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

covprcp_annualprcp_lw_tidy <- tidy(covprcp_annualprcp_lw2, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "annual precip in range",
         response_variable = "precipitation COV",
         window = "long-term",
         sample_size = covprcp_annualprcp_lw2$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

temp_HWI_dw_tidy <- tidy(temp_HWI_dw2, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "HWI",
         response_variable = "temperature anomaly",
         window = "decision",
         sample_size = temp_HWI_dw2$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

temp_HWI_lw_tidy <- tidy(temp_HWI_lw2, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "HWI",
         response_variable = "temperature anomaly",
         window = "long-term",
         sample_size = temp_HWI_lw2$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

light_CT_dw_tidy <- tidy(light_CT_dw2, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "C.T",
         response_variable = "light pollution",
         window = "decision",
         sample_size = light_CT_dw2$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

light_CT_lw_tidy <- tidy(light_CT_lw2, conf.int = T, conf.level = 0.95) %>%
  mutate(trait = "C.T",
         response_variable = "light pollution",
         window = "long-term",
         sample_size = light_CT_lw2$dims$N,
         lambda = 0) %>% # if lambda was fixed in the model, need to manually enter it here
  select(-p.value) %>%
  mutate_if(is.numeric, round, 3)

# combine into one single table of results to export
trait_results <- bind_rows(temp_STI_dw_tidy, temp_STI_lw_tidy, totalprcp_annualprcp_dw_tidy, totalprcp_annualprcp_lw_tidy, 
                           covprcp_annualprcp_dw_tidy, covprcp_annualprcp_lw_tidy, temp_HWI_dw_tidy, temp_HWI_lw_tidy, light_CT_dw_tidy, light_CT_lw_tidy) %>%
  select(trait, response_variable, window, sample_size, term:conf.high, lambda)

write.csv(trait_results, here("Models", "traitmodels.csv"))
