###### MAPS Phenology #######
### Script name: Step10_ModelValidation&Estimates.R
## Author(s): XXX (removed for peer review)

########## Objective/Description of Script #####################
# Extract parameter estimates for decision and long-term models
# Examine model convergence and validate that the models fit the data
#################################################################

#### Setup ####
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


# clear the environment to free up memory
rm(list = ls(all.names = TRUE))


# import the final models for the long-term and decision windows
# we are using the modes with no Gaussian Process 
dw_sppmodels <- readRDS(here("Models/Model Outputs", "dw_sppmodels_noGP.rds")) # models with no gaussian process

lw_sppmodels <- readRDS(here("Models/Model Outputs", "lw_sppmodels_noGP.rds")) # models with no gaussian process


##################################################################
##### Get the fixed effects and sigma from each model summary 
###################################################################

# this section extracts the information from the summary for each model and organizes it into a data frame for export and use in subsequent steps
# we do this separately for the decision and long-term window models

# for the decision window models
dw_modelresults <- dw_sppmodels %>%
  select(SPEC, dw_model) %>%
  mutate(
    summary = purrr::map(dw_model, ~summary), # add model summary
    fixed = purrr::map(dw_model, ~ { summary(.x)$fixed %>% rownames_to_column(., "Parameter")}), # pull the fixed effects for each model
    sigma = purrr::map(dw_model, ~ {summary(.x)$spec_pars %>% bind_rows(., .id="Parameter")}) # pull sigma for each model
  )

# extract fixed effects, gp, and sigma from nested data frame for the decision window model
dw_fixed <- dw_modelresults %>% unnest(fixed) %>% select(- dw_model, - summary, -sigma)
dw_sigma <- dw_modelresults %>% unnest(sigma) %>% select( - dw_model, - summary, -fixed) %>%
  mutate(Parameter = paste0("sigma"))

# combine into new data frame and add rounding criteria
dw_model_summaries <- bind_rows(dw_fixed, dw_sigma) %>% 
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
    sigma = purrr::map(lw_model, ~ {summary(.x)$spec_pars %>% bind_rows(., .id="Parameter")}) # pull sigma for each model
  )
# extract fixed effects, gp, and sigma from nested data frame for the long-term window model
lw_fixed <- lw_modelresults %>% unnest(fixed) %>% select(- lw_model, - summary, -sigma)
lw_sigma <- lw_modelresults %>% unnest(sigma) %>% select( - lw_model, - summary, -fixed) %>%
  mutate(Parameter = paste0("sigma"))

lw_model_summaries <- bind_rows(lw_fixed, lw_sigma) %>% 
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
                                                    nvariables = 3, # put 3 variables/plots per page
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
                                                    nvariables = 3, # put 3 variables/plots per page
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
