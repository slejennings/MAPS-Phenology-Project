
#### Make plots to visualize custom and default priors for Bayesian models

# load packages
library(invgamma)
library(ggplot2)
library(dplyr)
library(patchwork)


########## lscale ############

# plot custom prior for lscale
# gamma prior with (2, 0)
# define parameters for the distribution
shape_val <- 2
rate_val <- 1

gamma <- ggplot(data.frame(x = c(0, 7)), aes(x = x)) +
  stat_function(fun = dgamma, args = list(shape = shape_val, rate = rate_val),
                color = "darkblue", linewidth = 1) +
  labs(title = paste("Custom lscale prior: gamma (shape =", shape_val, ", rate =", rate_val, ")"),
       y = "Density", x = "Parameter") +
  theme_minimal()

# plot default prior for lscale
# inverse gamma prior

# define parameters for the distribution
shape_val_invg <- 1.49
rate_val_invg <- 0.06

invgamma <- ggplot(data.frame(x = c(0, 7)), aes(x = x)) +
  stat_function(fun = dinvgamma, args = list(shape = shape_val_invg, rate = rate_val_invg),
                color = "darkred", linewidth = 1) +
  labs(title = paste("Default lscale prior: inv gamma (shape =", shape_val_invg, ", rate =", rate_val_invg, ")"),
       y = "Density", x = "Parameter") +
  theme_minimal()

gamma/invgamma

########## sdgp (standard deviation of GP) ############

# plot default sdgp prior
# student-t with (3, 0, 2.5)
df_sdgp <- ggplot(data.frame(x = c(0, 7)), aes(x = x)) +
  stat_function(fun = dstudent_t, args = list(3, 0, 2.5),
                color = "gray40", linewidth = 1) +
  labs(title = paste("Default sdgp prior: student-t, df=3, mean = 0, scale = 2.5"),
       y = "Density", x = "Parameter") +
  theme_minimal()

# plot custom sdgp prior
# normal(0,1)
custom_sdgp <- ggplot(data.frame(x = c(0, 7)), aes(x = x)) +
  stat_function(fun = dnorm, args = list(0, 1),
                color = "darkseagreen", linewidth = 1) +
  labs(title = paste("Custom sdgp prior: normal, mean = 0, sd = 1"),
       y = "Density", x = "Parameter") +
  theme_minimal()

# combine two plots and compare
df_sdgp/custom_sdgp
