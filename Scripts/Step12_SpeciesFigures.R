###### MAPS Phenology #######
### Script name: Step12_SpeciesFigures.R
## Author(s): CDF, LP, SLJ

########## Objective/Description of Script #####################
# this script creates multiple figures for the manuscript that accompanies this analysis
# Fig 3: heatmap depicting the results of the species' models from Step8_BayesianSpeciesModels
# Fig 2: a histogram showing changes in breeding phenology across all species and stations
# supplementary figures showing change in breeding phenology, change in light, change in climate for each species
# supplementary figure (Fig S1) showing average change in breeding phenology, light, and climate with all species/populations pooled
#################################################################

#### Setup ####

# load packages
library(tidyverse)
library(here)
library(patchwork)
library(colorspace)
library(ape)
library(geiger)
library(ggtree)

# import files 

# phylogenetic tree for birds
tree <- read.tree(here("Data", "Jetz_ConsensusPhy.tre"))

# decision window model results
dw_models <- readRDS(here("Models/Model Outputs", "dw_model_summaries.rds"))

# long-term window model results
lw_models <- readRDS(here("Models/Model Outputs", "lw_model_summaries.rds")) 

# Sheet to convert 4-letter bird codes to commonnames
codetospec <- read_csv(here("Data", "BirdCodetoSpecies.csv")) %>%
  select(SPEC, COMMONNAME)

# t-statistics for change in bird breeding phenology and environmental variables
tstats <- readRDS(here("Outputs", "combined_t_stats.rds"))

# species eye morphometrics (we will be using this here to get scientific names for each species)
eye <- read.csv(here("Data", "species_eyes.csv"))

########################################################################################################
######### Figure 3: heatmap 
########################################################################################################


# we need to combine decision and long-term window models
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

# join all_models_wide with eye df that contains scientific names for species
# also join to codetospec which shows how four letter bird codes translate to common names
allmods_spp <- full_join(all_models_wide, eye, by="SPEC") %>% 
  select(-C.T, -Source, -Museum.Collection.IDs) %>%
  left_join(., codetospec) %>%
  column_to_rownames(., var = "Tree_name")

# join tree with the data
phydat <- treedata(tree, allmods_spp, sort=T) # join tree with data

birdtree <- phydat$phy # this is our trimmed tree for the 20 species

# these are the data associated with our trimmed tree
birddat <- as.data.frame(phydat$data) %>% # convert to df
  mutate(across(c(2:8), as.numeric)) # make sure columns containing t-stats are numeric

# get names for labeling the tree. 
# we need a df with scientific names using the underscore (Tree_name) and the common names we want to use as labels
treenames <- allmods_spp %>% rownames_to_column(., var="Tree_name") %>% select(Tree_name, COMMONNAME)

# plot the phylogenetic tree
treeplot <- ggtree(birdtree) %<+% 
  treenames + geom_tiplab(aes(label=COMMONNAME), # name tree tips
                           size = 3, # font size for tip labels
                           color="black", # font color for tip labels
                           offset=1) + # amount of space between tips and labels
  hexpand(0.9) # specify a fraction of x range to expand the x-axis limit by. direction = 1 is for right side

treeplot

# get data for the heatmap
# NOTE: T-STATISTICS FOR LIGHT FROM BOTH MODELS SHOWN IN THIS VERSION
heatmap_dat <- allmods_spp %>% 
  select(scaletempanom_DW_tstat, scaletempanom_LW_tstat, scaleprcp_DW_total_tstat, scaleprcp_LW_total_tstat,
         scaleprcp_DW_cov_tstat, scaleprcp_LW_cov_tstat, scalelight_DW_tstat, scalelight_LW_tstat)
  
# now create heatmap to go with the tree
# plot species effects as heatmap with phylogenetic tree on left 

# set up divergent palette
colorspace::diverging_hcl(n = 13, h = c(245, 125), c = c(30, 55), l = c(15, 95), power = c(0.75, 0.85), 
                          register = "custom_bluegreen" )
swatchplot(custom = diverging_hcl(n = 13, "custom_bluegreen" )) # examine palette

treeplusheatmap <- gheatmap(
  treeplot, heatmap_dat, offset=50, width=1.2,
  colnames_angle = 75, 
  colnames_offset_y = -0.2, 
  colnames_position="top",
  font.size=3,
  hjust=0, 
  custom_column_labels = c("Temperature Anomaly - decision", "Temperature Anomaly - long", "Total Precipitation - decision", "Total Precipitation - long", "Precipitation Variability - decision", "Precipitation Variability - long", "Light Pollution - decision", "Light Pollution - long")) +
  scale_fill_continuous_diverging(palette = "custom_bluegreen", n_interp=13, mid = 0, limits=c(-4.25,4.25), name="t-statistic", rev=F)+
  vexpand(.2,1) +
  theme(axis.title.x = element_text(size = 8, face = "bold"),
        legend.title = element_text(size= 8),
        legend.text = element_text(size= 8))

treeplusheatmap

# save
ggsave(treeplusheatmap, filename = "Fig3_TreePlusHeatMap.pdf", path = here("Figures"), width=20, height=21, units = "cm", device=cairo_pdf)
ggsave(treeplusheatmap, filename = "Fig3_TreePlusHeatMap.png", path = here("Figures"), width=20, height=21, units = "cm")


########################################################################################################
######### Figure 2: Histogram of Change in Bird Breeding Phenology
########################################################################################################

range(tstats$FY_tstat)

FYchangestat <- tstats %>%
  ggplot(aes(x = FY_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-5,5), name="t-statistic", rev=F) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in breeding phenology", 
       y = "Number of models") +
  coord_cartesian(xlim=c(-5.5,5.5)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=3, y=50, xend=5, yend=50, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-3, y=50, xend=-5, yend= 50,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 4, y=54, label = "Delaying", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -4, y=54, label = "Advancing", color = "#00294D", size=4)


FYchangestat

ggsave(FYchangestat, filename = "Fig2_ChangeStat.pdf", path = here("Figures"), width=20, height=12, units = "cm", device=cairo_pdf)
#ggsave(FYchangestat, filename = "Fig2_ChangeStat.png", path = here("Figures"), width=20, height=12, units = "cm")

#########################################################################################################
##### Supplementary figures: Figure S1
########################################################################################################

# this figure shows the change in breeding phenology, change in light pollution, and change in climate variables with all species and populations pooled

# set up divergent palette
colorspace::diverging_hcl(n = 13, h = c(245, 125), c = c(30, 55), l = c(15, 95), power = c(0.75, 0.85), 
                          register = "custom_bluegreen" )
swatchplot(custom = diverging_hcl(n = 13, "custom_bluegreen" )) # examine palette

##### Change in Breeding Phenology

range(tstats$FY_tstat) # use range to set limits for color palette

FYchangestat <- tstats %>%
  ggplot(aes(x = FY_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-5,5), name="t-statistic", rev=F, breaks=c(-4,-2,0,2,4)) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in breeding phenology", 
       y = "Number of models",
       title = "Bird Breeding Phenology") +
  coord_cartesian(xlim=c(-6,6)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=3, y=50, xend=5, yend=50, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-3, y=50, xend=-5, yend= 50,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 4, y=60, label = "Delaying", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -4, y=60, label = "Advancing", color = "#00294D", size=4) +
  theme(legend.key.size = unit(.4, 'cm'))


FYchangestat


##### Change in Temp Anomalies - DW

range(tstats$tempanom_DW_tstat) # use range to set limits for color palette

TAchangestat_dw <- tstats %>%
  ggplot(aes(x = tempanom_DW_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-6,6), name="t-statistic", rev=F, breaks=c(-6,-3,0,3,6)) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in temperature anomalies", 
       y = "Number of models",
       title = "Temperature Anomaly: Decision Window") +
  coord_cartesian(xlim=c(-6,6)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=3, y=50, xend=5, yend=50, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-3, y=50, xend=-5, yend= 50,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 4, y=60, label = "Warming", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -4, y=60, label = "Cooling", color = "#00294D", size=4) +
  theme(legend.key.size = unit(.4, 'cm'))

TAchangestat_dw

##### Change in Temp Anomalies - LW

range(tstats$tempanom_LW_tstat) # use range to set limits for color palette

TAchangestat_lw <- tstats %>%
  ggplot(aes(x = tempanom_LW_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-8,8), name="t-statistic", rev=F, breaks=c(-8,-4,0,4,8)) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in temperature anomalies", 
       y = "Number of models",
       title = "Temperature Anomaly: Long-term Window") +
  coord_cartesian(xlim=c(-6,8)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=3, y=50, xend=5, yend=50, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-3, y=50, xend=-5, yend= 50,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 4, y=60, label = "Warming", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -4, y=60, label = "Cooling", color = "#00294D", size=4) +
  theme(legend.key.size = unit(.4, 'cm'))

TAchangestat_lw

##### Change in Cumulative (total) Precipitation - DW

range(tstats$prcp_DW_total_tstat) # use range to set limits for color palette

TPchangestat_dw <- tstats %>%
  ggplot(aes(x = prcp_DW_total_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-6,6), name="t-statistic", rev=F, breaks=c(-6,-3,0,3,6)) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in cumulative precipitation", 
       y = "Number of models",
       title = "Cumulative Precipitation: Decision Window") +
  coord_cartesian(xlim=c(-6,6)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=3, y=50, xend=5, yend=50, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-3, y=50, xend=-5, yend= 50,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 4, y=60, label = "Wetter", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -4, y=60, label = "Drier", color = "#00294D", size=4) +
  theme(legend.key.size = unit(.4, 'cm'))

TPchangestat_dw

##### Change in Cumulative (total) Precipitation - LW

range(tstats$prcp_LW_total_tstat) # use range to set limits for color palette

TPchangestat_lw <- tstats %>%
  ggplot(aes(x = prcp_LW_total_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-5,5), name="t-statistic", rev=F, breaks=c(-4,-2,0,2,4)) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in cumulative precipitation", 
       y = "Number of models",
       title = "Cumulative Precipitation: Long-term Window") +
  coord_cartesian(xlim=c(-6,6)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=3, y=60, xend=5, yend=60, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-3, y=60, xend=-5, yend= 60,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 4, y=72, label = "Wetter", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -4, y=72, label = "Drier", color = "#00294D", size=4) +
  theme(legend.key.size = unit(.4, 'cm'))

TPchangestat_lw


##### Change in Precipitation Variability - DW

range(tstats$prcp_DW_cov_tstat) # use range to set limits for color palette


PVchangestat_dw <- tstats %>%
  ggplot(aes(x = prcp_DW_cov_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-6,6), name="t-statistic", rev=F, breaks=c(-6,-3,0,3,6)) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in coefficient of variation \n for precipitation", 
       y = "Number of models",
       title = "Precipitation Variability: Decision Window") +
  coord_cartesian(xlim=c(-6,6), ylim=c(0,85)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=4, y=50, xend=5, yend=50, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-4, y=50, xend=-5, yend= 50,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 4.5, y=60, label = "More Variable", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -4.5, y=60, label = "Less Variable", color = "#00294D", size=4) +
  theme(legend.key.size = unit(.4, 'cm'))

PVchangestat_dw

##### Change in Precipitation Variability - LW

range(tstats$prcp_LW_cov_tstat) # use range to set limits for color palette

PVchangestat_lw <- tstats %>%
  ggplot(aes(x = prcp_LW_cov_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-8,8), name="t-statistic", rev=F, breaks=c(-8,-4,0,4,8)) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in coefficient of variation \n for precipitation", 
       y = "Number of models",
       title = "Precipitation Variability: Long-term window") +
  coord_cartesian(xlim=c(-8,8), ylim=c(0,85)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=4, y=50, xend=6, yend=50, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-4, y=50, xend=-6, yend= 50,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 5, y=60, label = "More Variable", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -5, y=60, label = "Less Variable", color = "#00294D", size=4) +
  theme(legend.key.size = unit(.4, 'cm'))

PVchangestat_lw


##### Change in Light

range(tstats$light_tstat) # use range to set limits for color palette

NLchangestat <- tstats %>%
  ggplot(aes(x = light_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-18,18), name="t-statistic", rev=F, breaks=c(-14,-7,0,7,14)) + 
  geom_vline(aes(xintercept = -1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = - 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.44), 
             color = "#ABABAB", 
             linetype = "dashed", 
             linewidth = 1) +
  geom_vline(aes(xintercept = 1.96), 
             color = "#717171", 
             linetype = "dashed", 
             linewidth = 1) +
  labs(x = "Change in nighttime light", 
       y = "Number of models",
       title = "Light Pollution") +
  coord_cartesian(xlim=c(-18,18), ylim=c(0,200)) +
  theme_classic() +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.title = element_text(size= 9),
        legend.text = element_text(size= 9)) +
  annotate("segment", x=5, y=140, xend=10, yend=140, # add arrows at top of plot
           col="#002F00", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("segment", x=-5, y=140, xend=-10, yend= 140,
           col="#00294D", arrow=arrow(length=unit(0.3, "cm"))) +
  annotate("text", x = 7, y=167, label = "Brighter", color =  "#002F00", size=4) + # add text annotations above arrows
  annotate("text", x = -7, y=167, label = "Darker", color = "#00294D", size=4) +
  theme(legend.key.size = unit(.4, 'cm'))

NLchangestat

# combine plots
(row1<-wrap_elements((FYchangestat + plot_spacer() + NLchangestat) + plot_layout(widths = c(0.48, 0.04, 0.48)) + 
                       plot_annotation(tag_levels = list(c('(a)', '(b)'))) & theme(plot.tag = element_text(size = 14, face ="bold"))))
(row2<-wrap_elements((TAchangestat_dw + plot_spacer() + TAchangestat_lw) + plot_layout(widths = c(0.48, 0.04, 0.48)) + 
                       plot_annotation(tag_levels = list(c('(c)', '(d)'))) & theme(plot.tag = element_text(size = 14, face ="bold"))))
(row3<-wrap_elements((TPchangestat_dw + plot_spacer() + TPchangestat_lw) + plot_layout(widths = c(0.48, 0.04, 0.48)) + 
                       plot_annotation(tag_levels = list(c('(e)', '(f)'))) & theme(plot.tag = element_text(size = 14, face ="bold"))))
(row4<-wrap_elements((PVchangestat_dw + plot_spacer() + PVchangestat_lw) + plot_layout(widths = c(0.48, 0.04, 0.48)) + 
                       plot_annotation(tag_levels = list(c('(g)', '(h)'))) & theme(plot.tag = element_text(size = 14, face ="bold"))))

(tstat_hist <- row1 / row2 / row3 / row4) 

ggsave(tstat_hist, filename = "FigS1_TStatHistograms.pdf", path = here("Figures"), width=26, height=32, units = "cm", device=cairo_pdf)
#ggsave(tstat_hist, filename = "FigS1_TStatHistograms.png", path = here("Figures"), width=26, height=32, units = "cm")

#########################################################################################################
##### Supplementary figures: Figure S2
########################################################################################################

# set up divergent palette
colorspace::diverging_hcl(n = 13, h = c(245, 125), c = c(30, 55), l = c(15, 95), power = c(0.75, 0.85), 
                          register = "custom_bluegreen" )
swatchplot(custom = diverging_hcl(n = 13, "custom_bluegreen" )) # examine palette


tstats_plot <- left_join(tstats, codetospec)

for(i in unique(tstats_plot$SPEC)){
  
  tstats_sub <- tstats_plot %>%
    filter(SPEC == i)
  
  name <- unique(tstats_sub$COMMONNAME)
  
##### Change in Breeding Phenology
  FY <- tstats_sub %>%
    ggplot(aes(x = FY_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-5,5), name="t-statistic", rev=F, breaks=c(-4,-2,0,2,4)) + 
    geom_vline(aes(xintercept = 0), 
               color = "#717171", 
               linetype = "dashed", 
               linewidth = 1) +
    labs(x = "Change in breeding phenology", 
         y = "Number of models",
         title = "Bird Breeding Phenology") +
    coord_cartesian(xlim=c(-6,6)) +
    theme_classic() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size= 9),
          legend.text = element_text(size= 9)) +
    theme(legend.key.size = unit(.4, 'cm'))
  
  ##### Change in Temp Anomalies - DW
  
  TA_dw <- tstats_sub %>%
    ggplot(aes(x = tempanom_DW_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-6,6), name="t-statistic", rev=F, breaks=c(-6,-3,0,3,6)) + 
    geom_vline(aes(xintercept = 0), 
               color = "#717171", 
               linetype = "dashed", 
               linewidth = 1) +
    labs(x = "Change in temperature anomalies", 
         y = "Number of models",
         title = "Temperature Anomaly: Decision Window") +
    coord_cartesian(xlim=c(-6,6)) +
    theme_classic() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size= 9),
          legend.text = element_text(size= 9)) +
    theme(legend.key.size = unit(.4, 'cm'))
  
  
  ##### Change in Temp Anomalies - LW
  
  TA_lw <- tstats_sub %>%
    ggplot(aes(x = tempanom_LW_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-8,8), name="t-statistic", rev=F, breaks=c(-8,-4,0,4,8)) + 
    geom_vline(aes(xintercept = 0), 
               color = "#717171", 
               linetype = "dashed", 
               linewidth = 1) +
    labs(x = "Change in temperature anomalies", 
         y = "Number of models",
         title = "Temperature Anomaly: Long-term Window") +
    coord_cartesian(xlim=c(-6,8)) +
    theme_classic() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size= 9),
          legend.text = element_text(size= 9)) +
    theme(legend.key.size = unit(.4, 'cm'))
  
  ##### Change in Cumulative (total) Precipitation - DW
  
  TP_dw <- tstats_sub %>%
    ggplot(aes(x = prcp_DW_total_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-6,6), name="t-statistic", rev=F, breaks=c(-6,-3,0,3,6)) + 
    geom_vline(aes(xintercept = 0), 
               color = "#717171", 
               linetype = "dashed", 
               linewidth = 1) +
    labs(x = "Change in cumulative precipitation", 
         y = "Number of models",
         title = "Cumulative Precipitation: Decision Window") +
    coord_cartesian(xlim=c(-6,6)) +
    theme_classic() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size= 9),
          legend.text = element_text(size= 9)) +
    theme(legend.key.size = unit(.4, 'cm'))
  
  ##### Change in Cumulative (total) Precipitation - LW
  
  TP_lw <- tstats_sub %>%
    ggplot(aes(x = prcp_LW_total_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-5,5), name="t-statistic", rev=F, breaks=c(-4,-2,0,2,4)) + 
    geom_vline(aes(xintercept = 0), 
               color = "#717171", 
               linetype = "dashed", 
               linewidth = 1) +
    labs(x = "Change in cumulative precipitation", 
         y = "Number of models",
         title = "Cumulative Precipitation: Long-term Window") +
    coord_cartesian(xlim=c(-6,6)) +
    theme_classic() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size= 9),
          legend.text = element_text(size= 9)) +
    theme(legend.key.size = unit(.4, 'cm'))
  
  ##### Change in Precipitation Variability - DW
  PV_dw <- tstats_sub %>%
    ggplot(aes(x = prcp_DW_cov_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-6,6), name="t-statistic", rev=F, breaks=c(-6,-3,0,3,6)) + 
    geom_vline(aes(xintercept = 0), 
               color = "#717171", 
               linetype = "dashed", 
               linewidth = 1) +
    labs(x = "Change in coefficient of variation \n for precipitation", 
         y = "Number of models",
         title = "Precipitation Variability: Decision Window") +
    coord_cartesian(xlim=c(-6,6)) +
    theme_classic() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size= 9),
          legend.text = element_text(size= 9)) +
    theme(legend.key.size = unit(.4, 'cm'))
  
  ##### Change in Precipitation Variability - LW
  
  PV_lw <- tstats_sub %>%
    ggplot(aes(x = prcp_LW_cov_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-8,8), name="t-statistic", rev=F, breaks=c(-8,-4,0,4,8)) + 
    geom_vline(aes(xintercept = 0), 
               color = "#717171", 
               linetype = "dashed", 
               linewidth = 1) +
    labs(x = "Change in coefficient of variation \n for precipitation", 
         y = "Number of models",
         title = "Precipitation Variability: Long-term window") +
    coord_cartesian(xlim=c(-8,8)) +
    theme_classic() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size= 9),
          legend.text = element_text(size= 9)) +
    theme(legend.key.size = unit(.4, 'cm'))
  
  ##### Change in Light
  
  NL <- tstats_sub %>%
    ggplot(aes(x = light_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "custom_bluegreen", mid = 0, limits=c(-18,18), name="t-statistic", rev=F, breaks=c(-14,-7,0,7,14)) + 
    geom_vline(aes(xintercept = 0), 
               color = "#717171", 
               linetype = "dashed", 
               linewidth = 1) +
    labs(x = "Change in nighttime light", 
         y = "Number of models",
         title = "Light Pollution") +
    coord_cartesian(xlim=c(-18,18), ylim=c(0,8)) +
    theme_classic() +
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.title = element_text(size= 9),
          legend.text = element_text(size= 9))
  theme(legend.key.size = unit(.4, 'cm'))
  
  # combine plots
  # 1st row
  R1 <-wrap_elements((FY + plot_spacer() + NL) + plot_layout(widths = c(0.48, 0.04, 0.48)) + 
                       plot_annotation(tag_levels = list(c('(a)', '(b)'))) & theme(plot.tag = element_text(size = 14, face ="bold")))
  # 2nd row
  R2 <-wrap_elements((TA_dw + plot_spacer() + TA_lw) + plot_layout(widths = c(0.48, 0.04, 0.48)) + 
                       plot_annotation(tag_levels = list(c('(c)', '(d)'))) & theme(plot.tag = element_text(size = 14, face ="bold")))
  # third row
  R3 <-wrap_elements((TP_dw + plot_spacer() + TP_lw) + plot_layout(widths = c(0.48, 0.04, 0.48)) + 
                       plot_annotation(tag_levels = list(c('(e)', '(f)'))) & theme(plot.tag = element_text(size = 14, face ="bold")))
  # fourth row
  R4<-wrap_elements((PV_dw + plot_spacer() + PV_lw) + plot_layout(widths = c(0.48, 0.04, 0.48)) + 
                      plot_annotation(tag_levels = list(c('(g)', '(h)'))) & theme(plot.tag = element_text(size = 14, face ="bold")))
  
  spp_tstat_hist <- R1 / R2 / R3 / R4 +
    plot_annotation(title = name,
                    theme = theme(plot.title = element_text(size = 14, face="bold", hjust=0.5))) 
  
  ggsave(spp_tstat_hist, filename = paste(i, "_FigS2.pdf", sep = "_"), path = here("Figures/Figure S2"), width=26, height=32, units = "cm", device=cairo_pdf)
  
}

