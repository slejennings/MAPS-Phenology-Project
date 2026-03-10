###### MAPS Phenology #######
### Script name: Step8_BayesianSpeciesModels.R
## Author(s): LP, CDF, SLJ

########## Objective/Description of Script #####################
# run models for each species examining the relationship between change in breeding phenology and change in environmental variables (light, climate)
# we run one model for each species for the decision window and another for the long-term window
# we implement models using a Bayesian framework and the brms package
# all models will account for the spatial relationships in the data (i.e., the distance between MAPS stations)
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
  prior(normal(0, 3), class = "Intercept"), # Intercept - I tightened this somewhat after looking at the default settings
  prior(normal(0, 1), class = "sdgp"),   # GP std dev
  prior(gamma(2, 1), class = "lscale", coef ="") # GP length-scale
)

# for the lscale prior, you can specify a global prior as shown above
# or you can use class = "lscale", coef = "gpLongRLatR" to set a prior for the individual coefficient
# however, only the coefficient-specific prior will be applied. There will be a warning message produced about this when running the model
# the old batch of models used the default prior of inv_gamma(1.49, 0.057) for  class = "lscale", coef = "gpLongRLatR"
# so even though we specified a global lscale prior of gamma(2,1), this did not get implemented
# therefore, we can also use the priors below to achieve the same thing

custompriors <- c(
  prior(normal(0, 5), class = "b"),       # Fixed effects
  prior(normal(0, 3), class = "Intercept"), # Intercept
  prior(normal(0, 1), class = "sdgp")   # GP std dev  
)


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

# decision window model for all species
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

# NOTES: getting warning message - The global prior 'gamma(2, 1)' of class 'lscale' will not be used in the model as all related coefficients have individual priors already. 
# Some species still have divergent transitions:
# BEWR - 1, HOWR - 3, ORJU - 1
# if we increase adapt_delta to 0.9999 then BEWR and HOWR have no divergent transitions. This does not help ORJU though, which has fewer dts using 0.999


# save models
saveRDS(dw_sppmodels, here("Models/Model Outputs", "dw_sppmodels.rds"))

# import saved models
dw_sppmodels <- readRDS(here("Models/Model Outputs", "dw_sppmodels.rds"))

# pull out a single model and examine it
# change species code in filter() and change last function to look at different outputs

dw_sppmodels  %>% filter(SPEC == "HOWR") %>% pull(dw_model) %>% pluck(1) %>% summary() # model summary
dw_sppmodels  %>% filter(SPEC == "ORJU") %>% pull(dw_model) %>% pluck(1) %>% summary() 
dw_sppmodels  %>% filter(SPEC == "HOWR") %>% pull(dw_model) %>% pluck(1) %>% plot() # trace and density plots
dw_sppmodels  %>% filter(SPEC == "ORJU") %>% pull(dw_model) %>% pluck(1) %>% plot() 
dw_sppmodels  %>% filter(SPEC == "HOWR") %>% pull(dw_model) %>% pluck(1) %>% pp_check(., ndraws=20) # ppc density plot
dw_sppmodels  %>% filter(SPEC == "ORJU") %>% pull(dw_model) %>% pluck(1) %>% pp_check(., ndraws=20) 

###############################

# trace and density plots often show some extreme values for lscale
# I am not sure if this is a major problem but it doesn't seem ideal
# I tried experimenting with changing the priors
# changing intercept, slope and sdgp doesn't do anything helpful (not too surprisingly but it was worth trying)
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

# Out of curioisty...what if we try making the lscale prior a bit wider?
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


# test run with just a few species
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
# Another thing we can try is approximating the GP rather than letting the model calculate an exact GP
# to do this, we need to specify k = and c = inside the gp() function
# k is the number of basis functions to use when approximating the GP
# c defines the multiplicative constant of the predictors' range over which the predictions should be computed
# the reference doc suggests 5/4 could be a good default for c but also notes they are still working on providing better recommendations???

# I tried k = 20, 15, 10, 7 and 5 with c = 5/4 for a subset of species. 
# Approximating the GP does help eliminate the extreme values for lscale seen in the trace and density plots above
# However, we get a lot of divergent transitions (10s to 100s)
# There are fewer divergent transitions when using lower values of k but they are still not eliminated


# test run with just a few species
# modify k and c values
# I've tried k =5 with c = 5/4, 1.5 and 1.6. Notes below on how changing c alters the number of DTs
dw_sppmodels_k5 <- sppmodels_dat %>%
  filter(SPEC %in% c("HOWR", "BEWR", "ORJU")) %>%
  group_by(SPEC) %>% # group each species into its own data frame
  nest() %>%
  mutate(dat1 = purrr::map(data, ~as_tibble(droplevels(.x))), # make a version of each species data frame where unused factor levels for STA are dropped
        dw_model = purrr::map(dat1, # use this df in the model
                     ~ brms::brm(
                       bf(FY_tstat ~ scale(tempanom_DW_tstat) + scale(prcp_DW_total_tstat) + 
                            scale(prcp_DW_cov_tstat) + scale(light_tstat) + gp(LongR, LatR, k=5, c=1.4)), 
                       data = .,
                       family = gaussian(),
                       prior = custompriors,
                       save_pars = save_pars(all = TRUE),
                       init = "random",
                       control = list(adapt_delta = 0.999),
                       chains = CHAINS, iter = ITER, warmup = WARMUP, thin = THIN, seed = BAYES_SEED, cores = CORES)))


dw_sppmodels_k5  %>% filter(SPEC == "HOWR") %>% pull(dw_model) %>% pluck(1) %>% summary()
dw_sppmodels_k5  %>% filter(SPEC == "ORJU") %>% pull(dw_model) %>% pluck(1) %>% summary() 
dw_sppmodels_k5  %>% filter(SPEC == "BEWR") %>% pull(dw_model) %>% pluck(1) %>% summary() 

# with k = 5 and c = 5/4, the number of divergent transitions are HOWR - 41, ORJU - 32, BEWR - 31
# with k = 5 and c = 1.5, the number of divergent transitions are HOWR - 25, ORJU - 17, BEWR - 25
# with k = 5 and c = 1.6, the number of divergent transitions are HOWR - 21, ORJU - 23, BEWR - 25

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


##################################################################
##### Get the fixed effects, gp and sigma from each model summary 
###################################################################

# this section extracts the information from the summary for each model and organizes it into a data frame for export and use in subsequent steps
# we do this separately for the decision and long-term window models

# for the decision window models
dw_modelresults <- dw_sppmodels %>%
  select(SPEC, dw_model) %>%
  mutate(
    summary = purrr::map(dw_model, ~summary), # add model summary
    fixed = purrr::map(dw_model, ~ { summary(.x)$fixed %>% rownames_to_column(., "Parameter")}), # pull the fixed effects for each model
    gp = purrr::map(dw_model, ~ {summary(.x)$gp %>% bind_rows(., .id ="Parameter") %>% # pull the random effects for each model
        rownames_to_column(., "Intercept") %>% mutate(Parameter = paste0(Parameter, "_", str_sub(Intercept, end = -5))) %>%
        dplyr::select(-Intercept)}),
    sigma = purrr::map(dw_model, ~ {summary(.x)$spec_pars %>% bind_rows(., .id="Parameter")}) # pull sigma for each model
  )

# extract fixed effects, gp, and sigma from nested data frame for the decision window model
dw_fixed <- dw_modelresults %>% unnest(fixed) %>% select(- dw_model, - summary, -gp, -sigma)
dw_gp <- dw_modelresults %>% unnest(gp) %>% select(- dw_model, - summary, -fixed, -sigma)
dw_sigma <- dw_modelresults %>% unnest(sigma) %>% select( - dw_model, - summary, -gp, -fixed) %>%
  mutate(Parameter = paste0("sigma"))

# combine into new data frame and add rounding criteria
dw_model_summaries <- bind_rows(dw_fixed, dw_gp, dw_sigma) %>% 
  arrange(SPEC) %>%
  mutate(across(Estimate:Rhat, ~round(. ,4))) %>% # round estimate, SE, 95% CI boundaries and Rhat to 4 decimal places
  mutate(across(Bulk_ESS:Tail_ESS, ~ round(. ,0))) # round bulk and tail ESS to nearest integer

saveRDS(dw_model_summaries, here("Models/Model Outputs", "dw_model_summaries.rds"))
write.csv(dw_model_summaries, here("Models/Model Outputs", "dw_model_summaries.csv"))

################################
# for the long-term window models
lw_modelresults <- lw_sppmodels %>%
  select(SPEC, lw_model) %>%
  mutate(
    summary = purrr::map(lw_model, ~summary), # add model summary
    fixed = purrr::map(lw_model, ~ { summary(.x)$fixed %>% rownames_to_column(., "Parameter")}), # pull the fixed effects for each model
    gp = purrr::map(lw_model, ~ {summary(.x)$gp %>% bind_rows(., .id ="Parameter") %>% # pull the random effects for each model
        rownames_to_column(., "Intercept") %>% mutate(Parameter = paste0(Parameter, "_", str_sub(Intercept, end = -5))) %>%
        dplyr::select(-Intercept)}),
    sigma = purrr::map(lw_model, ~ {summary(.x)$spec_pars %>% bind_rows(., .id="Parameter")}) # pull sigma for each model
  )
# extract fixed effects, gp, and sigma from nested data frame for the long-term window model
lw_fixed <- lw_modelresults %>% unnest(fixed) %>% select(- lw_model, - summary, -gp, -sigma)
lw_gp <- lw_modelresults %>% unnest(gp) %>% select(- lw_model, - summary, -fixed, -sigma)
lw_sigma <- lw_modelresults %>% unnest(sigma) %>% select( - lw_model, - summary, -gp, -fixed) %>%
  mutate(Parameter = paste0("sigma"))

lw_model_summaries <- bind_rows(lw_fixed, lw_gp, lw_sigma) %>% 
  arrange(SPEC) %>%
  mutate(across(Estimate:Rhat, ~round(. ,4))) %>% # round estimate, SE, 95% CI boundaries and Rhat to 4 decimal places
  mutate(across(Bulk_ESS:Tail_ESS, ~ round(. ,0))) # round bulk and tail ESS to nearest integer

saveRDS(lw_model_summaries, here("Models/Model Outputs", "lw_model_summaries.rds"))
write.csv(lw_model_summaries, here("Models/Model Outputs", "lw_model_summaries.csv"))

#################################################################
##### Trace and density plots 
#################################################################

# this section makes trace and density plots for each model and exports them as pdf
# we do this separately for the decision and long-term window models

color_scheme_set("brightblue")

# decision window models
dw_tracedensity <- dw_sppmodels %>%
  mutate(tracedensity = purrr::map(dw_model, ~ plot(.x, theme=theme_minimal(),
                                       nvariables = 5, # put 5 variables/plots per page
                                       plot = T, newpage = T, ask= F)), .keep="none") %>% # .keep = "none" drops everything but the newly created plots
  deframe() # this converts the plots into a list

# export all plots as a single pdf and label each page with model/species
pdf(here("Models/Model Plots", "dw_tracedensityplots.pdf"), onefile = TRUE)
for (i in seq(length(dw_tracedensity))) {
  for(j in 1:2) { # there are two pages of plots for each species
    do.call("grid.arrange", c(dw_tracedensity[[i]][[j]], top = names(dw_tracedensity)[[i]]))   
  }
}
dev.off()   

###################################
# long-term window models
lw_tracedensity <- lw_sppmodels %>%
  mutate(tracedensity = purrr::map(lw_model, ~ plot(.x, theme=theme_minimal(),
                                             nvariables = 5, # put 5 variables/plots per page
                                             plot = T, newpage = T, ask= F)), .keep="none") %>% # .keep = "none" drops everything but the newly created plots
  deframe() # this converts the plots into a list

# export all plots as a single pdf and label each page with model/species
pdf(here("Models/Model Plots", "lw_tracedensityplots.pdf"), onefile = TRUE)
for (i in seq(length(lw_tracedensity))) {
  for(j in 1:2) { # there are two pages of plots for each species
    do.call("grid.arrange", c(lw_tracedensity[[i]][[j]], top = names(lw_tracedensity)[[i]]))   
  }
}
dev.off() 

#################################################################
##### Posterior predictive checks: density plots 
#################################################################

# this section of code gets posterior predictive check density plots for each model and exports them as pdf
# we do this separately for the decision and long-term window models

# for decision window model
dw_ppc_density <- dw_sppmodels %>%
  mutate(density = purrr::map(dw_model, ~ pp_check(.x, ndraws = 100) + xlim(-5, 5)), .keep="none") %>%
  deframe()

# export all plots as a single pdf and label each page with model/species
pdf(here("Models/Model Plots", "dw_ppc_densityplots.pdf"), onefile = TRUE)
for (i in seq(length(dw_ppc_density))) {
  do.call("grid.arrange", c(dw_ppc_density[i], top = names(dw_ppc_density)[i]))   
}
dev.off()

# for long-term window model
lw_ppc_density <- lw_sppmodels %>%
  mutate(density = purrr::map(lw_model, ~ pp_check(.x, ndraws = 100) + xlim(-5, 5)), .keep="none") %>%
  deframe()

# export all plots as a single pdf and label each page with model/species
pdf(here("Models/Model Plots", "lw_ppc_densityplots.pdf"), onefile = TRUE)
for (i in seq(length(lw_ppc_density))) {
  do.call("grid.arrange", c(lw_ppc_density[i], top = names(lw_ppc_density)[i]))   
}
dev.off()

