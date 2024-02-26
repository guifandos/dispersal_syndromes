#' ####################################################################### #
#' PROJECT: Dispersal syndromes on European birds
#' CONTENTS: 
#'  - This code is to analyse the correlation between the different traits
#'  DEPENDENCIES:
#'  - Run 1_Join_dispersal_traits
#'  - Data files: species dispersal distances and trait databases
#'  
#' AUTHOR: Guillermo Fandos
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list=ls())
## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "readr","dplyr","tidyr", "purrr", "reshape2", "caret", "readxl", "HH", "ggplot2", "tidymodels" # names of the packages required placed here as character objects
)

sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

tidymodels_prefer(quiet = TRUE) #uses the `conflicted` package to handle common conflicts with tidymodels and other packages.

# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}


# DATA ====================================================================

## Loading ----------------------------------------------------------------

dispersal_traits_total <- read.csv("data/data_process/dispersal_traits_avonet.csv")
names(dispersal_traits_total)

# Select traits to be analysed. 

traits_analysed <- c("median","upper_distance","HWI", "Body.mass..log." , "Range.Size", "Migration.2", "Diet","Habitat", "AnnualTemp", "TempRange", "PrecipRange" , "AnnualPrecip", "Clutch_MEAN", "Broods.per.year", "Age.of.independence", "Age.of.first.breeding", "Life.span")

traits_analysed <- c("median","upper_distance","LengthU_MEAN", "WeightU_MEAN" , "TarsusU_MEAN", "Clutch_MEAN", "Broods.per.year","Mortality.of.adults", "Post.fledging.mortality", "Hand.Wing.Index", "Mass" , "Age.of.independence", "Age.of.first.breeding", "Life.span")

# Transform the median and long distance dispersal to logarithm.

analysis_table <- dispersal_traits_total %>%
  dplyr::select(traits_analysed) %>% 
  dplyr::mutate(dispersal_distance= log(median +1),
         lon_dispersal_distance= log(upper_distance +1) )

# ANALYSIS ================================================================

### Correlation of the quantitative traits ####

# Select quantitative variables

analysis_table_quant <- analysis_table %>% 
  dplyr::select("median","upper_distance","LengthU_MEAN", "WeightU_MEAN" , "TarsusU_MEAN", "Clutch_MEAN", "Broods.per.year","Mortality.of.adults", "Post.fledging.mortality", "Hand.Wing.Index", "Mass" , "Age.of.independence", "Age.of.first.breeding", "Life.span") %>% 
  modify_if(., is.character, as.numeric) %>% 
  scale(.) %>% 
  as.data.frame(.) %>% 
  drop_na(.)

str(analysis_table_quant)

cormat <- round(cor(analysis_table_quant),2)
melted_cormat <- melt(cormat)
head(melted_cormat)
ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()


upper_tri <- get_upper_tri(cormat)
upper_tri

# Melt the correlation matrix

melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Heatmap
ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()



# Reorder the correlation matrix
cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Print the heatmap
print(ggheatmap)

############################
# Vif analysis####

resultado.vif<- HH::vif(analysis_table_quant)
resultado.vif


analysis_table_quant <- analysis_table %>% 
  dplyr::select("median","upper_distance", "TarsusU_MEAN", "Clutch_MEAN", "Broods.per.year","Mortality.of.adults", "Post.fledging.mortality", "Hand.Wing.Index", "Mass" , "Age.of.independence", "Age.of.first.breeding", "Life.span") %>% 
  modify_if(., is.character, as.numeric) %>% 
  scale(.) %>% 
  as.data.frame(.) %>% 
  drop_na(.)

resultado.vif<- HH::vif(analysis_table_quant)
resultado.vif
