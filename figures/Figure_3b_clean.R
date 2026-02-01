#' ####################################################################### #
#' PROJECT: Dispersal syndromes on European birds
#' CONTENTS: Figure 3b - Predictive accuracy of single-trait and multi-trait models
#' AUTHOR: Guillermo Fandos
#' ####################################################################### #

# =============================================================================
# 1. SETUP
# =============================================================================

# Clear environment
rm(list = ls())

# Load packages
library(tidyverse)
library(Rmisc)
library(rcartocolor)

# Set theme
theme_set(theme_classic())

# =============================================================================
# 2. LOAD DATA
# =============================================================================

# Load predictions data
validation_long <- read_csv("Figure3_source_data_b.csv")

# =============================================================================
# 3. DATA PROCESSING
# =============================================================================

# Fix typo in type column
validation_long$type <- gsub("betweeen", "between", validation_long$type)

# Calculate mean R2 per order, type, and complexity
summary_order <- validation_long %>% 
  dplyr::group_by(type, complexity, order) %>% 
  dplyr::summarise(Mean_R2 = mean(R2, na.rm = TRUE), .groups = "drop")

summary_order

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

# =============================================================================
# 4. CREATE FIGURE
# =============================================================================

# Define order of x-axis
level_order <- c("Life history", "Latitude", "Body mass", "HWI", 
                 "Distance migration", "Diet", "Habitat", "Dispersal syndrome")

sumrepdat2$complexity <- factor(sumrepdat2$complexity, levels = level_order)
table_prediction$complexity <- factor(table_prediction$complexity, levels = level_order)

# Create plot
fig3b <- ggplot(sumrepdat2, aes(x = complexity, y = Mean_R2, fill = type, group = type)) +
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
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    legend.position = "right"
  )

# Display plot
print(fig3b)

# =============================================================================
# 5. SAVE FIGURE
# =============================================================================

# Save as PDF (vectorial - recommended for publication)
ggsave(filename = "Figure_3b.pdf", 
       plot = fig3b,
       width = 8, 
       height = 5, 
       device = cairo_pdf)

# Save as high-resolution TIFF (alternative)
ggsave(filename = "Figure_3b.tiff", 
       plot = fig3b,
       width = 8, 
       height = 5, 
       device = "tiff", 
       dpi = 300)

cat("\n✓ Figures saved:\n")
cat("  - Figure_3b.pdf (vectorial)\n")
cat("  - Figure_3b.tiff (300 DPI)\n")
