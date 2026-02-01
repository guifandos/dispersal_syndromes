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
library(Rmisc)
library(rcartocolor)

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
fig3_data <- read_csv("Figure3_source_data_b.csv")

# Fix typo in type column
fig3_data$type <- gsub("betweeen", "between", fig3_data$type)

# Calculate mean R2 per order, type, and complexity
summary_order <- fig3_data %>%
  dplyr::group_by(type, complexity, order) %>%
  dplyr::summarise(Mean_R2 = mean(R2, na.rm = TRUE), .groups = "drop")

# Rename complexity levels for display
table_prediction <- summary_order %>%
  mutate(complexity = as.factor(complexity))

table_prediction$complexity <- recode_factor(
  table_prediction$complexity,
  "only_PC1"          = "Life history",
  "only_latitude"     = "Latitude",
  "only_Latitude"     = "Latitude",
  "only_body_mass"    = "Body mass",
  "only_HWI"          = "HWI",
  "only_distance_mig" = "Distance migration",
  "only_diet"         = "Diet",
  "only_habitat"      = "Habitat",
  "only_habitat_for"  = "Habitat",
  "only_habita_for"   = "Habitat",
  "total"             = "Dispersal syndrome"
)

# Calculate summary statistics (mean and SE across orders)
sumrepdat2 <- summarySE(table_prediction,
                        measurevar = "Mean_R2",
                        groupvars = c("type", "complexity"))

# Define order of x-axis
level_order <- c("Life history", "Latitude", "Body mass", "HWI",
                 "Distance migration", "Diet", "Habitat", "Dispersal syndrome")

sumrepdat2$complexity <- factor(sumrepdat2$complexity, levels = level_order)
table_prediction$complexity <- factor(table_prediction$complexity, levels = level_order)

# Create plot
figure3 <- ggplot(sumrepdat2, aes(x = complexity, y = Mean_R2, fill = type, group = type)) +
  # Mean points
  geom_point(aes(color = type),
             shape = 16,
             size = 4,
             position = position_dodge(0.5)) +
  # Error bars (SE)
  geom_errorbar(aes(ymin = Mean_R2 - se, ymax = Mean_R2 + se, color = type),
                width = 0.2,
                linewidth = 0.8,
                position = position_dodge(0.5)) +
  # Individual order points
  geom_point(data = table_prediction,
             aes(x = complexity, y = Mean_R2, group = type, colour = type, shape = order),
             size = 3,
             alpha = 0.5,
             position = position_dodge(0.5)) +
  # Scales
  scale_colour_brewer(palette = "Dark2") +
  scale_fill_brewer(palette = "Dark2") +
  scale_shape_manual(values = c(15, 17, 4, 8)) +
  # Axis
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = NULL, y = expression(R^2)) +
  # Theme
  theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

# Save as PNG
ggsave("Figure3_predictive_accuracy.png", figure3, width = 10, height = 6, dpi = 600)

# Save as PDF (vectorial - recommended for publication)
ggsave("Figure3_predictive_accuracy.pdf", figure3, width = 8, height = 5, device = cairo_pdf)

# Save as TIFF (alternative)
ggsave("Figure3_predictive_accuracy.tiff", figure3, width = 8, height = 5, device = "tiff", dpi = 300)


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
cat("  - Figure3_source_data_b.csv\n")
