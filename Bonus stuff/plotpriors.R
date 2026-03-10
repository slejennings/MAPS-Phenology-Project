library(invgamma)
library(ggplot2)
library(patchwork)

# Plot Gamma Prior
# Define parameters
shape_val <- 2
rate_val <- 1

gamma <- ggplot(data.frame(x = c(0, 7)), aes(x = x)) +
  stat_function(fun = dgamma, args = list(shape = shape_val, rate = rate_val),
                color = "darkblue", linewidth = 1) +
  labs(title = paste("Gamma Prior (Shape =", shape_val, ", Rate =", rate_val, ")"),
       y = "Density", x = "Parameter") +
  theme_minimal()

# Plot Inverse Gamma Prior
# Define parameters
shape_val_invg <- 1.49
rate_val_invg <- 0.06

invgamma <- ggplot(data.frame(x = c(0, 7)), aes(x = x)) +
  stat_function(fun = dinvgamma, args = list(shape = shape_val_invg, rate = rate_val_invg),
                color = "darkred", linewidth = 1) +
  labs(title = paste("Inv Gamma Prior (Shape =", shape_val_invg, ", Rate =", rate_val_invg, ")"),
       y = "Density", x = "Parameter") +
  theme_minimal()

gamma/invgamma
