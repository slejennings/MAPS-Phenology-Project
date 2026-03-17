###### MAPS Phenology #######
### Script name: Step12_SpeciesFigures.R
## Author(s): CDF, LP, SLJ

########## Objective/Description of Script #####################
# this script creates multiple figures for the manuscript that accompanies this analysis
# Fig 3: heatmap depicting the results of the species' models from Step8_BayesianSpeciesModels
# Fig 2: a histogram showing changes in breeding phenology across all species and stations
# supplementary figures showing change in breeding phenology, change in light, change in climate for each species
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
# NOTE: MAY NEED TO ADJUST THE LIMITS OF THE DIVERGING PALETTE IF T-STATS EXCEED -3.5 OR 3.5
treeplusheatmap <- gheatmap(
  treeplot, heatmap_dat, offset=50, width=1.2,
  colnames_angle = 75, 
  colnames_offset_y = -0.2, 
  colnames_position="top",
  font.size=3,
  hjust=0, 
  custom_column_labels = c("Temp. Anomaly - decision", "Temp. Anomaly - long", "Precipitation - decision", "Precipitation - long", "Precipitation COV - decision", "Precipitation COV - long", "Light Pollution - decision", "Light Pollution - long")) +
  scale_fill_continuous_diverging(palette = "cork", n_interp=11, mid = 0, limits=c(-3.5,3.5), name="t-statistic", rev=F,
                                  p1 = 1, p2=1)+
  vexpand(.2,1) +
  theme(axis.title.x = element_text(size = 8, face = "bold"),
        legend.title = element_text(size= 8),
        legend.text = element_text(size= 8))

treeplusheatmap

# save
ggsave(treeplusheatmap, filename = "Fig3_TreePlusHeatMap.pdf", path = here("Figures"), width=20, height=20, units = "cm")
ggsave(treeplusheatmap, filename = "Fig3_TreePlusHeatMap.png", path = here("Figures"), width=20, height=20, units = "cm")


########################################################################################################
######### Figure 2: Histogram of Change in Bird Breeding Phenology
########################################################################################################

FYchangestat <- tstats %>%
  ggplot(aes(x = FY_tstat)) + 
  geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
  scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-5,5), name="t-statistic", rev=F) + 
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

ggsave(FYchangestat, filename = "Fig2_ChangeStat.pdf", path = here("Figures"), width=20, height=12, units = "cm")
ggsave(FYchangestat, filename = "Fig2_ChangeStat.png", path = here("Figures"), width=20, height=12, units = "cm")


#########################################################################################################
##### Supplementary figures: change in breeding and environment for each species
########################################################################################################

tstats_plot <- left_join(tstats, codetospec)

# make vector of species where long-term model was preferred
longspp <- c("ORJU", "WIWA", "AUWA", "SWTH", "DOWO",
"LISP", "YEWA", "HOWR", "COYE", "BEWR", "NOCA", "CARW",
"TUTI") 

# specify species where long-term window was preferred
lw_tstats_plot <- tstats_plot %>%
  filter(SPEC %in% longspp)

# specify species where decision window was preferred or species with no difference between decision and long-term window
dw_tstats_plot <- tstats_plot %>%
  filter(!SPEC %in% longspp) # keep all the species not in longspp

### make plots for decision window species

for(i in unique(dw_tstats_plot$SPEC)){
  
  tstats_sub <- dw_tstats_plot %>%
    filter(SPEC == i)
  
  name <- unique(tstats_sub$COMMONNAME) 
  
  # breeding phenology
  specstat <- tstats_sub %>%
    ggplot(aes(x = FY_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-5,5), name="t-statistic", rev=F) + 
    geom_vline(aes(xintercept = 0), 
               color = "#ABABAB", 
               linetype = "dashed", 
               linewidth = 1) +
    coord_cartesian(xlim=c(-5, 5)) +
    labs(x = "Change in breeding phenology", 
         y = "Number of models") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size= 7),
          legend.text = element_text(size= 6)) +
    theme(legend.key.size = unit(.4, 'cm'))
  # temp anomaly
  anomstat <- tstats_sub %>%
    ggplot(aes(x = tempanom_DW_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-8,8), name="t-statistic", rev=F) + 
    geom_vline(aes(xintercept = 0), 
               color = "#ABABAB", 
               linetype = "dashed", 
               linewidth = 1) +
    coord_cartesian(xlim=c(-8, 8)) +
    labs(x = "Change in temperature anomaly - decision window", 
         y = "Number of models") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size= 7),
          legend.text = element_text(size= 6)) +
    theme(legend.key.size = unit(.4, 'cm')) 
  # total precip 
  prcpsumstat <- tstats_sub %>%
    ggplot(aes(x = prcp_DW_total_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-8,8), name="t-statistic", rev=F) + 
    geom_vline(aes(xintercept = 0), 
               color = "#ABABAB", 
               linetype = "dashed", 
               linewidth = 1) +
    coord_cartesian(xlim=c(-8, 8)) +
    labs(x = "Change in total precipitation - decision window", 
         y = "Number of models") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size= 7),
          legend.text = element_text(size= 6)) +
    theme(legend.key.size = unit(.4, 'cm'))
  # cov for precip
  prcpcovstat <- tstats_sub %>%
    ggplot(aes(x = prcp_DW_cov_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-8,8), name="t-statistic", rev=F) + 
    geom_vline(aes(xintercept = 0), 
               color = "#ABABAB", 
               linetype = "dashed", 
               linewidth = 1) +
    coord_cartesian(xlim=c(-8, 8)) +
    labs(x = "Change in precipitation variability - decision window", 
         y = "Number of models") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size= 7),
          legend.text = element_text(size= 6)) +
    theme(legend.key.size = unit(.4, 'cm'))
  # light
  lightstat <- tstats_sub %>%
    ggplot(aes(x = light_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-10,10), name="t-statistic", rev=F) + 
    geom_vline(aes(xintercept = 0), 
               color = "#ABABAB", 
               linetype = "dashed", 
               linewidth = 1) +
    coord_cartesian(xlim=c(-10, 10)) +
    labs(x = "Change in light pollution", 
         y = "Number of models") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size= 7),
          legend.text = element_text(size= 6)) +
    theme(legend.key.size = unit(.4, 'cm'))
  
  # specify plot layout
  a <- wrap_elements((plot_spacer() +specstat + plot_spacer()) + plot_layout(widths = c(0.22, 0.56, 0.22)))
  b <- wrap_elements(anomstat + plot_spacer() + lightstat + plot_layout(widths = c(0.48, 0.04, 0.48)))
  c <- wrap_elements(prcpcovstat + plot_spacer() + prcpsumstat + plot_layout(widths = c(0.48, 0.04, 0.48)))
  
  
  sppplot <- a/b/c +
    plot_annotation(title = name,
                    theme = theme(plot.title = element_text(size = 14, face="bold", hjust=0.5)))
  
  
  ggsave(sppplot, filename = paste(i, "dw_changestat.png", sep = "_"), path = here("Figures/Spp change stats/Decision window"), width=23, height=16, units = "cm")
  
}


### make plots for long window species

for(i in unique(lw_tstats_plot$SPEC)){
  
  tstats_sub <- lw_tstats_plot %>%
    filter(SPEC == i)
  
  name <- unique(tstats_sub$COMMONNAME)
  
  # breeding phenology
  specstat <- tstats_sub %>%
    ggplot(aes(x = FY_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-5,5), name="t-statistic", rev=F) + 
    geom_vline(aes(xintercept = 0), 
               color = "#ABABAB", 
               linetype = "dashed", 
               linewidth = 1) +
    coord_cartesian(xlim=c(-5, 5)) +
    labs(x = "Change in breeding phenology", 
         y = "Number of models") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size= 7),
          legend.text = element_text(size= 6)) +
    theme(legend.key.size = unit(.4, 'cm'))
  # temp anomaly
  anomstat <- tstats_sub %>%
    ggplot(aes(x = tempanom_LW_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-8,8), name="t-statistic", rev=F) + 
    geom_vline(aes(xintercept = 0), 
               color = "#ABABAB", 
               linetype = "dashed", 
               linewidth = 1) +
    coord_cartesian(xlim=c(-8, 8)) +
    labs(x = "Change in temperature anomaly - long window", 
         y = "Number of models") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size= 7),
          legend.text = element_text(size= 6)) +
    theme(legend.key.size = unit(.4, 'cm')) 
  # total precip 
  prcpsumstat <- tstats_sub %>%
    ggplot(aes(x = prcp_LW_total_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-8,8), name="t-statistic", rev=F) + 
    geom_vline(aes(xintercept = 0), 
               color = "#ABABAB", 
               linetype = "dashed", 
               linewidth = 1) +
    coord_cartesian(xlim=c(-8, 8)) +
    labs(x = "Change in total precipitation - long window", 
         y = "Number of models") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size= 7),
          legend.text = element_text(size= 6)) +
    theme(legend.key.size = unit(.4, 'cm'))
  # cov for precip
  prcpcovstat <- tstats_sub %>%
    ggplot(aes(x = prcp_LW_cov_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-8,8), name="t-statistic", rev=F) + 
    geom_vline(aes(xintercept = 0), 
               color = "#ABABAB", 
               linetype = "dashed", 
               linewidth = 1) +
    coord_cartesian(xlim=c(-8, 8)) +
    labs(x = "Change in precipitation variability - long window", 
         y = "Number of models") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size= 7),
          legend.text = element_text(size= 6)) +
    theme(legend.key.size = unit(.4, 'cm'))
  # light
  lightstat <- tstats_sub %>%
    ggplot(aes(x = light_tstat)) + 
    geom_histogram(aes(fill = after_stat(x)), binwidth = 0.2, position = "identity") + 
    scale_fill_continuous_diverging(palette = "cork", mid = 0, limits=c(-10,10), name="t-statistic", rev=F) + 
    geom_vline(aes(xintercept = 0), 
               color = "#ABABAB", 
               linetype = "dashed", 
               linewidth = 1) +
    coord_cartesian(xlim=c(-10, 10)) +
    labs(x = "Change in light pollution", 
         y = "Number of models") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size= 7),
          legend.text = element_text(size= 6)) +
    theme(legend.key.size = unit(.4, 'cm'))
  
  # specify plot layout
  a <- wrap_elements((plot_spacer() +specstat + plot_spacer()) + plot_layout(widths = c(0.22, 0.56, 0.22)))
  b <- wrap_elements(anomstat + plot_spacer() + lightstat + plot_layout(widths = c(0.48, 0.04, 0.48)))
  c <- wrap_elements(prcpcovstat + plot_spacer() + prcpsumstat + plot_layout(widths = c(0.48, 0.04, 0.48)))
  
  
  sppplot <- a/b/c +
    plot_annotation(title = name,
                    theme = theme(plot.title = element_text(size = 14, face="bold", hjust=0.5)))
  
  
  ggsave(sppplot, filename = paste(i, "lw_changestat.png", sep = "_"), path = here("Figures/Spp change stats/Long window"), width=23, height=16, units = "cm")
  
}
