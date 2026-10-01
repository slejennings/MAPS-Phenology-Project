###### MAPS Phenology #######
### Script name: Step9_CompareDW&LWModels.R
## Author(s): XXX (removed for peer review)

########## Objective/Description of Script #####################
# compare the long-term and decision window model for each species using AICc
# perform a post-hoc analysis to examine the relationship between log body mass and delta AICc across species
# create figure S3
#################################################################

#### Setup ####
# load packages
library(tidyverse)
library(here)
library(spmodel)

# import decision window models 

dw_sppmodels <- readRDS(here("Models", "dw_models.rds")) %>%
  select(SPEC, model) %>% 
  rename(dw_model = model)

# import long-term window models

lw_sppmodels <- readRDS(here("Models", "lw_models.rds")) %>%
  select(SPEC, model) %>%
  rename(lw_model = model)

# put the two models together into a single df
models_df <- left_join(dw_sppmodels, lw_sppmodels)

# get AICc for each model
models_AICc <- models_df %>%
  mutate(AICc_dw = purrr::map(dw_model, ~spmodel::AICc(.)),
         AICc_lw = purrr::map(lw_model, ~spmodel::AICc(.))) %>%
  select(-dw_model, -lw_model)


# compare models
models_AICc_diff <- models_AICc %>%
  rowwise() %>%
  mutate(AICc_diff = round(abs(AICc_dw-AICc_lw),1)) %>%
  filter(AICc_diff >=4) %>%
  mutate(preferredmod = ifelse(AICc_dw < AICc_lw, "decision", "long-term"))

# conclusions of model selection:
# requiring difference of at least 4 to identify a difference in model performance
# ORJU, WIWA, HOWR, DOWO, LISP - long-term
# SOSP, AMRO, GRCA - decision


##################################################################################
### post-hoc analysis of delta AICc and body mass, create figure S3

# load a few additional packages
library(performance)
library(ggeffects)
library(ape)
library(geiger)
library(nlme)

# import species traits file that contain body mass
traits <- read.csv(here("Data", "species_traits.csv"), header=T)
eye <- read.csv(here("Data", "species_eyes.csv"), header=T) # has correct scientific names for joining data to the tree

# phylogenetic tree for birds
tree <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

# get delta AICc for all species
# join traits (including log body mass)
# move the Tree_name column to rownames to facilitate pairing data with phylogenetic tree
models_AICc_diff_all <- models_AICc %>%
  rowwise() %>%
  mutate(AICc_diff = round(AICc_dw-AICc_lw,1)) %>%
  left_join(., traits) %>%
  left_join(., eye) %>%
  select(SPEC, AICc_diff, Body_mass_log, Tree_name) %>%
  column_to_rownames(., var = "Tree_name")

# pair data with phylogenetic tree
phydat_AICc <- geiger::treedata(tree, models_AICc_diff_all, sort=T) # join tree with data

birdtree_AICc <- phydat_AICc$phy # this is our trimmed tree for the 20 species

# these are the data associated with our trimmed tree
dat_AICc <- as.data.frame(phydat_AICc$data) %>% # convert to df
  mutate(across(c(AICc_diff, Body_mass_log), as.numeric)) # make sure columns that need to be are numeric

# run PGLS model
AICc_mass <- gls(AICc_diff ~ Body_mass_log, 
                   data = dat_AICc, 
                   correlation = corPagel(0, phy = birdtree_AICc, fixed=T), method = "ML")
summary(AICc_mass)
confint(AICc_mass)

# plot the results

# get effect of body mass on delta AICc
eff_AICc_mass <- plot(ggeffects::predict_response(AICc_mass, terms =c("Body_mass_log")), colors = "#1a5988")

# add data, labels, nice formatting to plot
FigS3 <- eff_AICc_mass  +
  geom_point(data = dat_AICc, aes(x = Body_mass_log, y = AICc_diff), color = "#1a5988", size = 3.1, pch = 19)+
  labs(title= "", x = "Natural log of body mass (g)", y = "AICc difference")+
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), 
        axis.title.x = element_text(color = "black", size = 12, margin = margin(t=0.3, unit="cm")), # add space between axis title and axis labels
        axis.title.y = element_text(color = "black", size = 12, margin = margin(r=0.3, unit="cm"))) + # add space between axis title and axis labels
  geom_hline(yintercept=0, linetype="dashed", color = "gray", linewidth=.6)

FigS3    

ggsave(FigS3, filename = "FigS3_AICcBodyMass.pdf", path = here("Figures"), width=11, height=10, units = "cm", device=cairo_pdf)
#ggsave(FigS3, filename = "FigS3_AICcBodyMass.png", path = here("Figures"), width=11, height=10, units = "cm")
