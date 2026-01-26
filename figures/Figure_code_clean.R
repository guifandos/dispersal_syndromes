# ==============================================================================
# CLEAN CODE FOR GENERATING FIGURES 1-3
# 
# Simple mechanistic traits outperform complex syndromes in predicting 
# avian dispersal distances
#
# Guillermo Fandos, Robert A. Robinson, Damaris Zurell
#
# Corresponding author: Guillermo Fandos (gfandos@ucm.es)
# ==============================================================================

library(tidyverse)
library(ggplot2)
library(ggdist)

# ==============================================================================
# FIGURE 1: Forest Plot - Dispersal Syndromes Coefficients
# ==============================================================================

# Load source data
fig1_data <- read.csv("Figure1_source_data.csv")

# Prepare data
fig1_data <- fig1_data %>%
  mutate(
    variable = factor(variable, levels = c(
      "HWI", "Habitat", "Diet", "Body mass", "Latitude", "Life history"
    )),
    dispersal_type = factor(dispersal_type, levels = c("Average Median", "Average Long", "Average long"))
  )

# Color palette
my_colors <- c(
  "Average Median" = "#ff7f0e",
  "Average Long" = "#1f77b4",
  "Average long" = "#1f77b4"
)

# Forest plot
figure1 <- ggplot(fig1_data, aes(y = variable, x = estimate, 
                                  xmin = CI_lower_95, xmax = CI_upper_95,
                                  color = dispersal_type)) +
  geom_pointrange(position = position_dodge(width = 0.5), size = 0.8) +
  geom_vline(xintercept = 0, lty = 2, color = "grey50") +
  scale_color_manual(
    values = my_colors,
    labels = c("Median dispersal", "Long-distance dispersal", "Long-distance dispersal"),
    name = NULL
  ) +
  labs(x = "Standardized coefficient", y = NULL) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 12)
  )

ggsave("Figure1_forest_plot.png", figure1, width = 8, height = 5, dpi = 600)


# ==============================================================================
# FIGURE 2: Scatterplots - Trait-Dispersal Relationships
# ==============================================================================

# Load source data
fig2_data <- read.csv("Figure2_source_data.csv")

# Common theme
theme_fig2 <- theme_classic() +
  theme(
    axis.text = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position = "bottom",
    legend.title = element_blank()
  )

# Color palette
colors_dispersal <- c("Median dispersal" = "#ff7f0e", "Long-distance dispersal" = "#1f77b4")

# Panel A: Life history (PC1)
fig2a <- ggplot(fig2_data, aes(x = life_history_PC1)) +
  geom_point(aes(y = median_dispersal_log, color = "Median dispersal"), 
             size = 2, alpha = 0.6) +
  geom_point(aes(y = long_dispersal_log, color = "Long-distance dispersal"), 
             size = 2, alpha = 0.6) +
  geom_smooth(aes(y = median_dispersal_log, color = "Median dispersal"),
              method = "lm", se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = long_dispersal_log, color = "Long-distance dispersal"),
              method = "lm", se = TRUE, alpha = 0.2, linetype = "dashed") +
  scale_color_manual(values = colors_dispersal) +
  labs(x = "Life history (slow-fast)", y = "Dispersal distance (log km)") +
  theme_fig2

# Panel B: Body mass
fig2b <- ggplot(fig2_data, aes(x = body_mass_log)) +
  geom_point(aes(y = median_dispersal_log, color = "Median dispersal"), 
             size = 2, alpha = 0.6) +
  geom_point(aes(y = long_dispersal_log, color = "Long-distance dispersal"), 
             size = 2, alpha = 0.6) +
  geom_smooth(aes(y = median_dispersal_log, color = "Median dispersal"),
              method = "lm", se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = long_dispersal_log, color = "Long-distance dispersal"),
              method = "lm", se = TRUE, alpha = 0.2, linetype = "dashed") +
  scale_color_manual(values = colors_dispersal) +
  labs(x = "Body mass (log g)", y = "Dispersal distance (log km)") +
  theme_fig2

# Panel C: HWI
fig2c <- ggplot(fig2_data, aes(x = HWI_log)) +
  geom_point(aes(y = median_dispersal_log, color = "Median dispersal"), 
             size = 2, alpha = 0.6) +
  geom_point(aes(y = long_dispersal_log, color = "Long-distance dispersal"), 
             size = 2, alpha = 0.6) +
  geom_smooth(aes(y = median_dispersal_log, color = "Median dispersal"),
              method = "lm", se = TRUE, alpha = 0.2) +
  geom_smooth(aes(y = long_dispersal_log, color = "Long-distance dispersal"),
              method = "lm", se = TRUE, alpha = 0.2, linetype = "dashed") +
  scale_color_manual(values = colors_dispersal) +
  labs(x = "Hand-Wing Index (log)", y = "Dispersal distance (log km)") +
  theme_fig2

# Combine panels
library(patchwork)
figure2 <- (fig2a | fig2b | fig2c) + 
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("Figure2_trait_dispersal.png", figure2, width = 14, height = 5, dpi = 600)


# ==============================================================================
# FIGURE 3: Predictive Accuracy - Within and Between Order
# ==============================================================================

# Load source data
fig3_data <- read.csv("Figure3_source_data.csv")

# Calculate summary statistics (mean and SD per model type and validation type)
fig3_summary <- fig3_data %>%
  group_by(model_type, validation_type) %>%
  summarise(
    mean_R2 = mean(R2, na.rm = TRUE),
    sd_R2 = sd(R2, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# Order model types
model_order <- c("life_history", "latitude", "body_mass", "HWI", 
                 "migration_distance", "diet", "habitat", "dispersal_syndrome")
fig3_summary$model_type <- factor(fig3_summary$model_type, levels = model_order)

# Clean labels
model_labels <- c(
  "life_history" = "Life history",
  "latitude" = "Latitude", 
  "body_mass" = "Body mass",
  "HWI" = "HWI",
  "migration_distance" = "Migration distance",
  "diet" = "Diet",
  "habitat" = "Habitat",
  "dispersal_syndrome" = "Dispersal syndrome"
)

# Colors
colors_validation <- c("within" = "#ff7f0e", "between" = "#1f77b4")

# Plot
figure3 <- ggplot(fig3_summary, aes(x = model_type, y = mean_R2, 
                                     fill = validation_type, color = validation_type)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  geom_errorbar(aes(ymin = mean_R2 - sd_R2, ymax = mean_R2 + sd_R2),
                width = 0.3, position = position_dodge(width = 0.7)) +
  scale_fill_manual(values = colors_validation, 
                    labels = c("Between-order", "Within-order"),
                    name = NULL) +
  scale_color_manual(values = colors_validation,
                     labels = c("Between-order", "Within-order"),
                     name = NULL) +
  scale_x_discrete(labels = model_labels) +
  scale_y_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = expression(R^2)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 11, color = "black", angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title = element_text(size = 14, color = "black"),
    legend.position = "bottom",
    legend.text = element_text(size = 12)
  )

ggsave("Figure3_predictive_accuracy.png", figure3, width = 10, height = 6, dpi = 600)


# ==============================================================================
# SUMMARY
# ==============================================================================

cat("\n=== FIGURES GENERATED ===\n")
cat("Figure 1: Forest plot of dispersal syndrome coefficients\n")
cat("Figure 2: Trait-dispersal relationships (3 panels)\n")
cat("Figure 3: Predictive accuracy within and between orders\n")
cat("\nSource data files required:\n")
cat("  - Figure1_source_data.csv\n")
cat("  - Figure2_source_data.csv\n")
cat("  - Figure3_source_data.csv\n")
