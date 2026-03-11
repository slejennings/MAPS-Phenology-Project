###### MAPS Phenology #######
### Script name: Step8A_BayesianSpeciesModels.R
## Author(s): LP, CDF, SLJ

########## Objective/Description of Script #####################
# the goal of this script is to run models for each species examining the relationship between change in breeding phenology and change in environmental variables (light, climate)
# we want to obtain one model for each species for the decision window and another for the long-term window
# we implement models using a Bayesian framework and the brms package
# because MAPS stations are distributed across North America, we may need to account for the spatial relationships in the data (i.e., the distance between MAPS stations) 
# to do this, we can use a Gaussian Process term based on the latitude and longitude of each station
# We compare models with and without a Gaussian Process term
# Specifically, we run models where an exact GP term is calculated, models with an approximate GP term and a third batch without the GP term
# in the subsequent script, we will compare the models to find the model specifications that work best for species across the long-term and decision window
#################################################################

#################################################################
#### Setup ####
#################################################################

# load packages
library(tidyverse)
library(here)
library(gridExtra)
library(patchwork)
library(colorspace)
library(bayesplot)
library(brms)
library(posterior)
library(rstan)

# import data 

# t-statistics for bird phenology and environmental change
tstats <- readRDS(here("Outputs", "combined_t_stats.rds"))

# MAPS stations and their latitude/longitude
stations <- readRDS(here("Outputs", "STA_finallist.rds")) %>%
  mutate(LatR = round(DECLAT, 3), # round coordinates to 3 digits
         LongR = round(DECLNG, 3),
         STA = as.factor(STA)) %>%
  select(STA, LatR, LongR)

# join t-stats and station info
# set STA and SPEC as factors
sppmodels_dat <- left_join(tstats, stations) %>%
  mutate(STA = as.factor(STA),
         SPEC = as.factor(SPEC))

summary(sppmodels_dat)

# examine samples sizes for each species
# every species should have data from at least 30 stations
sample_check <- sppmodels_dat %>%
  group_by(SPEC) %>%
  count()

#################################################################
##### Set up priors and specifications for models
#################################################################

# examine default priors
get_prior(data = sppmodels_dat, 
          family = gaussian(),
          bf(FY_tstat ~ scale(tempanom_DW_tstat) + scale(prcp_DW_total_tstat) + 
               scale(prcp_DW_cov_tstat) + scale(light_tstat) + gp(LongR, LatR)))


# specify weakly informative priors
custompriors <- c(
  prior(normal(0, 5), class = "b"),       # Fixed effects
  prior(normal(0, 3), class = "Intercept"), # Intercept
  prior(normal(0, 1), class = "sdgp")   # GP std dev  
)

# for the lscale prior, you can specify a global prior using class = "lscale", coef = ""
# or you can use class = "lscale", coef = "gpLongRLatR" to set a prior for the individual coefficient
# however, only the coefficient-specific prior will be applied
# there will be a warning message produced about the global prior not being applied when running the model
# for the coefficient-specific prior, we are using the default of inv_gamma(1.49, 0.06)

# model specifications
CHAINS <- 4
ITER <- 8000
WARMUP <- 2000
THIN <- 2
CORES <- 4
BAYES_SEED <- 579

#################################################################
##### Run the models
#################################################################

############### Decision window models #####################

# note: default setting for brms is to scale the predictors inside the gp() function
# see: https://paulbuerkner.com/brms/reference/gp.html and look at scale argument
scale(sppmodels_dat$LongR)
scale(sppmodels_dat$LatR)

# Decision Window Model 1: calculate an exact Gaussian Process for each model
# these are the settings that seem to work best in terms of reducing the number of models/species with divergent transitions
# this output has been saved and can be imported below to save re-running these models

dw_sppmodels <- sppmodels_dat %>%
  group_by(SPEC) %>% # group each species into its own data frame
  nest() %>%
  mutate(dat1 = purrr::map(data, ~as_tibble(droplevels(.x))), # make a version of each species data frame where unused factor levels for STA are dropped
         dw_model = purrr::map(dat1, # use this df in the model
                        ~ brms::brm(
                          bf(FY_tstat ~ scale(tempanom_DW_tstat) + scale(prcp_DW_total_tstat) + 
                               scale(prcp_DW_cov_tstat) + scale(light_tstat) + gp(LongR, LatR)),
                          data = .,
                          family = gaussian(),
                          prior = custompriors,
                          save_pars = save_pars(all = TRUE),
                          control = list(adapt_delta = 0.999),
                          chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED, cores = CORES)))

# Some species still have divergent transitions:
# BEWR - 1, HOWR - 3, ORJU - 1
# if we increase adapt_delta to 0.9999 then BEWR and HOWR have no divergent transitions
# this does not help ORJU though, which has fewer divergent transitions using 0.999

# save models
saveRDS(dw_sppmodels, here("Models/Model Outputs", "dw_sppmodels.rds"))

# import saved models
dw_sppmodels <- readRDS(here("Models/Model Outputs", "dw_sppmodels.rds"))

# pull out a single model and examine it
# change species code in filter() and change the last function in the pipe to look at different types of outputs

dw_sppmodels  %>% filter(SPEC == "HOWR") %>% pull(dw_model) %>% pluck(1) %>% summary() # model summary
dw_sppmodels  %>% filter(SPEC == "ORJU") %>% pull(dw_model) %>% pluck(1) %>% summary() 
dw_sppmodels  %>% filter(SPEC == "HOWR") %>% pull(dw_model) %>% pluck(1) %>% plot() # trace and density plots
dw_sppmodels  %>% filter(SPEC == "ORJU") %>% pull(dw_model) %>% pluck(1) %>% plot() 
dw_sppmodels  %>% filter(SPEC == "HOWR") %>% pull(dw_model) %>% pluck(1) %>% pp_check(., ndraws=20) # ppc density plot
dw_sppmodels  %>% filter(SPEC == "ORJU") %>% pull(dw_model) %>% pluck(1) %>% pp_check(., ndraws=20) 

###############################

# the trace and density plots for the model with the exact GP often show some extreme values for lscale
# This doesn't seem ideal
# the prior for lscale (which is using the default prior) is already very narrow, so I did not restrict it anymore
# here is a plot of that prior:
library(invgamma)
# Define parameters
shape_val_invg <- 1.49
rate_val_invg <- 0.06

ggplot(data.frame(x = c(0, 5)), aes(x = x)) +
  stat_function(fun = dinvgamma, args = list(shape = shape_val_invg, rate = rate_val_invg),
                color = "darkred", linewidth = 1) +
  labs(title = paste("Inv Gamma Prior (Shape =", shape_val_invg, ", Rate =", rate_val_invg, ")"),
       y = "Density", x = "Parameter") + theme_minimal()

# see the page about default priors in brms: https://paulbuerkner.com/brms/reference/set_prior.html
# for GP priors, it notes that the default lscale prior can be too restrictive
# So, we will examine what happens if we try making the lscale prior a bit wider

# plot a new slightly wider prior
shape_val_invg <- 2
rate_val_invg <- 0.5

ggplot(data.frame(x = c(0, 5)), aes(x = x)) +
  stat_function(fun = dinvgamma, args = list(shape = shape_val_invg, rate = rate_val_invg),
                color = "darkred", linewidth = 1) +
  labs(title = paste("Inv Gamma Prior (Shape =", shape_val_invg, ", Rate =", rate_val_invg, ")"),
       y = "Density", x = "Parameter") + theme_minimal()

# specify new priors
custompriors_2 <- c(
  prior(normal(0, 5), class = "b"),       # Fixed effects
  prior(normal(0, 3), class = "Intercept"), # Intercept - I tightened this somewhat after looking at the default settings
  prior(normal(0, 1), class = "sdgp"),   # GP std dev
  prior(inv_gamma(2, 0.5), class = "lscale", coef ="gpLongRLatR") # GP length-scale
)

# do a test run with just a few species
dw_sppmodels_cp2 <- sppmodels_dat %>%
  filter(SPEC %in% c("HOWR", "BEWR", "ORJU")) %>%
  group_by(SPEC) %>% # group each species into its own data frame
  nest() %>%
  mutate(dat1 = purrr::map(data, ~as_tibble(droplevels(.x))), # make a version of each species data frame where unused factor levels for STA are dropped
         dw_model = purrr::map(dat1, # use this df in the model
                               ~ brms::brm(
                                 bf(FY_tstat ~ scale(tempanom_DW_tstat) + scale(prcp_DW_total_tstat) + 
                                      scale(prcp_DW_cov_tstat) + scale(light_tstat) + gp(LongR, LatR)), 
                                 data = .,
                                 family = gaussian(),
                                 prior = custompriors_2,
                                 save_pars = save_pars(all = TRUE),
                                 init = "random",
                                 control = list(adapt_delta = 0.999),
                                 chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED, cores = CORES)))


# This did not help. In fact there are now more extreme values in the lscale plot
dw_sppmodels_cp2  %>% filter(SPEC == "HOWR") %>% pull(dw_model) %>% pluck(1) %>% summary() # 0 dts
dw_sppmodels_cp2  %>% filter(SPEC == "BEWR") %>% pull(dw_model) %>% pluck(1) %>% summary() # 0 dts
dw_sppmodels_cp2  %>% filter(SPEC == "ORJU") %>% pull(dw_model) %>% pluck(1) %>% summary() # 6 dts
dw_sppmodels_cp2  %>% filter(SPEC == "HOWR") %>% pull(dw_model) %>% pluck(1) %>% plot() 
dw_sppmodels_cp2  %>% filter(SPEC == "BEWR") %>% pull(dw_model) %>% pluck(1) %>% plot()
dw_sppmodels_cp2  %>% filter(SPEC == "ORJU") %>% pull(dw_model) %>% pluck(1) %>% plot() 

###########################################
# Decision Window Model 2: use an approximate Gaussian Process for each model

# Another thing we can try is approximating the GP rather than letting the model calculate an exact GP
# to do this, we need to specify k = and c = inside the gp() function
# k is the number of basis functions to use when approximating the GP
# c defines the multiplicative constant of the predictors' range over which the predictions should be computed

# run all decision window models with k = 20 and c = 5/4 (this is a Hilbert space approximation)
dw_sppmodels_k20 <- sppmodels_dat %>%
  group_by(SPEC) %>% # group each species into its own data frame
  nest() %>%
  mutate(dat1 = purrr::map(data, ~as_tibble(droplevels(.x))), # make a version of each species data frame where unused factor levels for STA are dropped
         dw_model = purrr::map(dat1, # use this df in the model
                               ~ brms::brm(
                                 bf(FY_tstat ~ scale(tempanom_DW_tstat) + scale(prcp_DW_total_tstat) + 
                                      scale(prcp_DW_cov_tstat) + scale(light_tstat) + gp(LongR, LatR, k=20, c=5/4)), 
                                 data = .,
                                 family = gaussian(),
                                 prior = custompriors,
                                 save_pars = save_pars(all = TRUE),
                                 init = "random",
                                 control = list(adapt_delta = 0.999),
                                 chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED, cores = CORES)))


# save and export models
saveRDS(dw_sppmodels_k20, here("Models/Model Outputs", "dw_sppmodels_k20.rds"))

##################
# run all decision window models without GP

# specify prior for non-GP model (drop GP-specific parameters)
custompriors_noGP <- c(
    prior(normal(0, 5), class = "b"),       # Fixed effects
    prior(normal(0, 3), class = "Intercept") # Intercept
  )


dw_sppmodels_noGP <- sppmodels_dat %>%
  group_by(SPEC) %>% # group each species into its own data frame
  nest() %>%
  mutate(dat1 = purrr::map(data, ~as_tibble(droplevels(.x))), # make a version of each species data frame where unused factor levels for STA are dropped
         dw_model = purrr::map(dat1, # use this df in the model
                               ~ brms::brm(
                                 bf(FY_tstat ~ scale(tempanom_DW_tstat) + scale(prcp_DW_total_tstat) + 
                                      scale(prcp_DW_cov_tstat) + scale(light_tstat)), 
                                 data = .,
                                 family = gaussian(),
                                 prior = custompriors_noGP,
                                 save_pars = save_pars(all = TRUE),
                                 init = "random",
                                 control = list(adapt_delta = 0.999),
                                 chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED, cores = CORES)))


# save and export models
saveRDS(dw_sppmodels_noGP, here("Models/Model Outputs", "dw_sppmodels_noGP.rds"))

######################################################################################
############### Long-term window models #####################

# long-term window model for all species
# these are the settings that seem to work best in terms of reducing the number of models/species with divergent transitions
# this output has been saved and can be imported below to save re-running these models

lw_sppmodels <- sppmodels_dat %>%
  ungroup() %>%
  group_by(SPEC) %>% # group each species into its own data frame
  nest() %>%
  mutate(dat1 = purrr::map(data, ~as_tibble(droplevels(.x))), # make a version of each species data frame where unused factor levels for STA are dropped
         lw_model = purrr::map(dat1, # use this df in the model
                               ~ brms::brm(
                                 bf(FY_tstat ~ scale(tempanom_LW_tstat) + scale(prcp_LW_total_tstat) + 
                                      scale(prcp_LW_cov_tstat) + scale(light_tstat) + gp(LongR, LatR)),
                                 data = .,
                                 family = gaussian(),
                                 prior = custompriors_lsc,
                                 init = "random",
                                 save_pars = save_pars(all = TRUE),
                                 control = list(adapt_delta = 0.999),
                                 chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED, cores = CORES)))


# NOTES: getting warning message - The global prior 'gamma(2, 1)' of class 'lscale' will not be used in the model as all related coefficients have individual priors already. 
# Some species still have divergent transitions:
# HOWR - 1, ORJU - 3
# if we increase adapt_delta to 0.9999 then HOWR has no divergent transitions

# save models
saveRDS(lw_sppmodels, here("Models/Model Outputs", "lw_sppmodels.rds"))

# import saved models
lw_sppmodels <- readRDS(here("Models/Model Outputs", "lw_sppmodels.rds"))

# pull a specific model and examine the outputs
lw_sppmodels  %>% filter(SPEC == "BCCH") %>% pull(lw_model) %>% pluck(1) %>% summary()
lw_sppmodels  %>% filter(SPEC == "BEWR") %>% pull(lw_model) %>% pluck(1) %>% summary()
lw_sppmodels  %>% filter(SPEC == "HOWR") %>% pull(lw_model) %>% pluck(1) %>% summary()
lw_sppmodels %>% filter(SPEC == "ORJU") %>% pull(lw_model) %>% pluck(1) %>% summary() 

lw_sppmodels  %>% filter(SPEC == "BCCH") %>% pull(lw_model) %>% pluck(1) %>% plot()
lw_sppmodels  %>% filter(SPEC == "BEWR") %>% pull(lw_model) %>% pluck(1) %>% plot()
lw_sppmodels  %>% filter(SPEC == "HOWR") %>% pull(lw_model) %>% pluck(1) %>% plot()
lw_sppmodels %>% filter(SPEC == "ORJU") %>% pull(lw_model) %>% pluck(1) %>% plot() 

lw_sppmodels  %>% filter(SPEC == "BCCH") %>% pull(lw_model) %>% pluck(1) %>% pp_check()
lw_sppmodels  %>% filter(SPEC == "BEWR") %>% pull(lw_model) %>% pluck(1) %>% pp_check()
lw_sppmodels  %>% filter(SPEC == "HOWR") %>% pull(lw_model) %>% pluck(1) %>% pp_check()
lw_sppmodels %>% filter(SPEC == "ORJU") %>% pull(lw_model) %>% pluck(1) %>% pp_check()

##################
# run all long-term window models with k = 20 and c = 5/4

lw_sppmodels_k20 <- sppmodels_dat %>%
  group_by(SPEC) %>% # group each species into its own data frame
  nest() %>%
  mutate(dat1 = purrr::map(data, ~as_tibble(droplevels(.x))), # make a version of each species data frame where unused factor levels for STA are dropped
         lw_model = purrr::map(dat1, # use this df in the model
                               ~ brms::brm(
                                 bf(FY_tstat ~ scale(tempanom_LW_tstat) + scale(prcp_LW_total_tstat) + 
                                      scale(prcp_LW_cov_tstat) + scale(light_tstat) + gp(LongR, LatR, k=20, c=5/4)), 
                                 data = .,
                                 family = gaussian(),
                                 prior = custompriors,
                                 save_pars = save_pars(all = TRUE),
                                 init = "random",
                                 control = list(adapt_delta = 0.999),
                                 chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED, cores = CORES)))


# save and export models
saveRDS(lw_sppmodels_k20, here("Models/Model Outputs", "lw_sppmodels_k20.rds"))

##################
# run all decision window models without GP

# specify prior for non-GP model (drop GP-specific parameters)
custompriors_noGP <- c(
  prior(normal(0, 5), class = "b"),       # Fixed effects
  prior(normal(0, 3), class = "Intercept") # Intercept
)


lw_sppmodels_noGP <- sppmodels_dat %>%
  group_by(SPEC) %>% # group each species into its own data frame
  nest() %>%
  mutate(dat1 = purrr::map(data, ~as_tibble(droplevels(.x))), # make a version of each species data frame where unused factor levels for STA are dropped
         lw_model = purrr::map(dat1, # use this df in the model
                               ~ brms::brm(
                                 bf(FY_tstat ~ scale(tempanom_LW_tstat) + scale(prcp_LW_total_tstat) + 
                                      scale(prcp_LW_cov_tstat) + scale(light_tstat)), 
                                 data = .,
                                 family = gaussian(),
                                 prior = custompriors_noGP,
                                 save_pars = save_pars(all = TRUE),
                                 init = "random",
                                 control = list(adapt_delta = 0.999),
                                 chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED, cores = CORES)))


# save and export models
saveRDS(lw_sppmodels_noGP, here("Models/Model Outputs", "lw_sppmodels_noGP.rds"))
