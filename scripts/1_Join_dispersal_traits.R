#' ####################################################################### #
#' PROJECT: Dispersal syndromes on European birds
#' CONTENTS: 
#'  - This code is to join the dispersal distance dataset with the different traits databases 
#'  DEPENDENCIES:
#'  - Code documents needed to execute this document
#'  - Data files: species dispersal distances and trait databases
#'  
#' AUTHOR: Guillermo Fandos
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list=ls())

## Directories ------------------------------------------------------------
### Define dicrectories in relation to project directory
Dir.Base <- getwd()
Dir.Data <- file.path(Dir.Base, "data")
Dir.Exports <- file.path(Dir.Base, "Exports")
### Create directories which aren't present yet
Dirs <- c(Dir.Data, Dir.Exports)
CreateDir <- sapply(Dirs, function(x) if(!dir.exists(x)) dir.create(x))

## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "readr","dplyr","tidyr", "taxize", "taxadb", "magrittr", "readxl" # names of the packages required placed here as character objects
)

sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function


# DATA ====================================================================


## Loading ----------------------------------------------------------------

# Load dispersal data####

### Total ####
distance_total <- read_csv("data/dispersal_distance/species_dispersal_distances.csv") %>% 
  group_by(species, type) %>% 
  filter(function_comparison== max(function_comparison)) %>% 
  mutate_if(is.character,as.factor) %>% 
  filter(type== "total")
### Breeding ####
distance_breeding <- read_csv("data/dispersal_distance/species_dispersal_distances.csv") %>% 
  group_by(species, type) %>% 
  filter(function_comparison== max(function_comparison)) %>% 
  mutate_if(is.character,as.factor) %>% 
  filter(type== "breeding")
### Natal ####
distance_natal <- read_csv("data/dispersal_distance/species_dispersal_distances.csv") %>% 
  group_by(species, type) %>% 
  filter(function_comparison== max(function_comparison)) %>% 
  mutate_if(is.character,as.factor) %>% 
  filter(type== "natal")

#  all the tables

total_dispersal_all <- read_csv("data/dispersal_distance/species_dispersal_distances.csv") %>% 
  group_by(species, type) %>% 
  filter(function_comparison== max(function_comparison)) %>% 
  mutate_if(is.character,as.factor)

# Load trait data ####
# We are going to use taxadb package to join the different datasets and avoid species name issues
#Sys.setenv(TAXADB_HOME="D:/Github_clon/dispersal_traits/")
#taxadb:::taxadb_dir()

td_create("col", overwrite=FALSE)
td_create("gbif", overwrite = TRUE)
td_create("itis", overwrite=FALSE)

# Different traits databases ####
# Traits for European birds (Storchová, L., & Hořák, D. (2018). Life‐history characteristics of European birds. Global Ecology and Biogeography, 27(4), 400-406.)
bird_traits <- read.delim("data/traits/Life-history characteristics of European birds.txt")
# Life history traits already calculated (Hanzelka, J., Horká, P., & Reif, J. (2019). Spatial gradients in country‐level population trends of European birds. Diversity and Distributions, 25(10), 1527-1536.)
bird_traits_devel <- read_excel("bibliography/Hanzelka_et_al_2019_D_D/ddi12945-sup-0002-tables3.xls")
names(bird_traits_devel)
# Population trends European birds 
trends <- read.csv("bibliography/Brlík_2021/trends2017.csv")
names(trends)
head(trends)
# Eltonian traits
BirdFuncDat <- read.delim("bibliography/Willmann_et al_2014/BirdFuncDat.txt")
# Sheard_2020_traits
hwi_index <- read_excel("bibliography/Sheard_2020/Dataset HWI 2020-04-10.xlsx")

#avonet
avonet <- read_excel("data/traits/avonet/Supplementary dataset 1.xlsx", 
                     sheet = "AVONET3_BirdTree")

## Manipulation -----------------------------------------------------------

############################################
### Join Trait Databases with Dispersal ####
############################################

# 1. Recoding name of species 
dispersal_distance <- total_dispersal_all

dispersal_traits <- left_join(dispersal_distance, bird_traits, by= c("species"="Species"))

species_name_change <- dispersal_traits %>% filter(is.na(Family))
# Check what species we need to recode for a perfect match
unique(species_name_change$species) #Names to recode 23
dispersal_distance$species_name_recoded <- dispersal_distance$species


dispersal_distance <-  dispersal_distance %>%
  mutate(species_name_recoded = recode(species_name_recoded, 
                                       "Egretta alba" = "Ardea alba",
                                       "Anas penelope" =  "Mareca penelope"  ,          
                                       "Anas strepera" =  "Mareca strepera",           
                                       "Anas clypeata" = "Spatula clypeata",
                                       "Mergus albellus" = "Mergellus albellus",
                                       "Aquila pomarina" = "Clanga pomarina",
                                       "Tetrao tetrix" = "Lyrurus tetrix",             
                                       "Philomachus pugnax"= "Calidris pugnax", 
                                       "Phalaropus fulicaria" = "Phalaropus fulicarius",
                                       "Stercorarius skua" = "Catharacta skua",
                                       "Larus minutus" = "Hydrocoloeus minutus",           
                                       "Larus argentatus cachinnans" = "Larus argentatus",
                                       "Sterna caspia" = "Hydroprogne caspia", 
                                       "Sterna sandvicensis" = "Thalasseus sandvicensis",
                                       "Sterna albifrons" = "Sternula albifrons",          
                                       "Chlidonias hybridus" = "Chlidonias hybrida", 
                                       "Apus melba" = "Tachymarptis melba",      
                                       "Dendrocopos medius" = "Leiopicus medius",
                                       "Dendrocopos minor"  = "Dryobates minor",        
                                       "Dendrocopos tridactylus" = "Picoides tridactylus",
                                       "Calandrella rufescens" = "Alaudala rufescens",
                                       "Delichon urbica"  = "Delichon urbicum" ,
                                       "Anthus spinoletta petrosus" = "Anthus petrosus",
                                       "Saxicola torquata"= "Saxicola torquatus",   
                                       "Hippolais pallida" = "Iduna pallida",      
                                       "Hippolais pallida opacus" = "Iduna pallida",   
                                       "Regulus ignicapillus" = "Regulus ignicapilla",
                                       "Parus palustris" = "Poecile palustris",     
                                       "Parus montanus"  = "Poecile montanus",     
                                       "Parus cinctus" = "Poecile cinctus",             
                                       "Parus cristatus" = "Lophophanes cristatus",
                                       "Parus ater" = "Periparus ater",         
                                       "Parus caeruleus" = "Cyanistes caeruleus",      
                                       "Cyanopica cyana"   = "Cyanopica cyanus",         
                                       "Corvus corone cornix" = "Corvus corone",   
                                       "Serinus citrinella" = "Carduelis citrinella", 
                                       "Carduelis chloris" = "Chloris chloris",
                                       "Carduelis spinus" = "Spinus spinus",          
                                       "Carduelis cannabina"  = "Linaria cannabina",
                                       "Carduelis flavirostris" = "Linaria cannabina",
                                       "Carduelis flammea cabaret"  = "Linaria cannabina",
                                       "Carduelis flammea"  = "Linaria cannabina",
                                       "Miliaria calandra" = "Emberiza calandra",
                                       "Luscinia svecica" = "Luscinia luscinia"
  ))

dispersal_traits <- left_join(dispersal_distance, bird_traits, by= c("species_name_recoded"="Species"))
species_name_change <- dispersal_traits %>% filter(is.na(Family))
unique(species_name_change$species) #Names to recode 3
# Corvus cornix is very similar to corvus corone. 
# Casmerodius albus is Ardea alba
# Larus cachinnans is very similar to Larus michahellis

dispersal_distance <-  dispersal_distance %>%
  mutate(species_name_recoded = recode(species_name_recoded, 
                                       "Casmerodius albus" = "Ardea alba",
                                       "Corvus cornix" =  "Corvus corone",
                                       "Larus cachinnans" =  "Larus michahellis"))
dispersal_traits <- left_join(dispersal_distance, bird_traits, by= c("species_name_recoded"="Species"))
species_name_change <- dispersal_traits %>% filter(is.na(Family))
unique(species_name_change$species) #Names to recode 3

##################
# 2. Using taxadb and getting the code for each specie

species_name <- unique(dispersal_traits$species_name_recoded)
species_identification <- as.data.frame(species_name) %>% 
  mutate(taxon = get_ids(species_name, "col"))

species_wo_code <- species_identification %>% 
  filter(is.na(taxon))
ff <- filter_name(species_wo_code[5,]) %>%
  mutate(acceptedNameUsage = get_names(acceptedNameUsageID)) %>% 
  select(scientificName, taxonomicStatus, acceptedNameUsage, acceptedNameUsageID) %>% 
  drop_na(scientificName)
rr <- ff %>% 
  mutate(taxon = get_ids(acceptedNameUsage, "col"))


# Note: Still fix the name of several species
############################################################

traits_devel <- bird_traits_devel %>% mutate(taxon = get_ids(Species, "col"))
dispersal <- dispersal_traits %>% left_join(.,species_identification, by= c("species"= "species_name"))
trends_id <- trends %>% mutate(taxon = get_ids(species, "col"))
hwi_index_id <- hwi_index %>% 
  rename(species = "Species name") %>% 
  mutate(taxon = get_ids(species, "col"))
bird_func_trait <- BirdFuncDat %>% select(Scientific, BodyMass.Value, Diet.5Cat) %>%  
  mutate(taxon = get_ids(Scientific, "col"))

#save(dispersal, traits_devel, hwi_index_id, bird_func_trait, trends_id, file = "data/traits/trait_data_id.RData")
#load("data/traits/trait_data_id.RData")


#joined_traits <- left_join(dispersal[,c(1:14, 97 )], traits_devel, by = "taxon") %>% 
joined_traits <- left_join(dispersal[,c(1:56, 97 )], traits_devel, by = "taxon") %>% 
  left_join(., trends_id, by = "taxon") %>% 
  left_join(., bird_func_trait, by = "taxon") %>% 
  drop_na(ID)  %>% 
  distinct(species.x, .keep_all = TRUE) %>% 
  left_join(., hwi_index_id, by = "taxon") %>% 
  drop_na(ID)  %>% 
  distinct(species.x, .keep_all = TRUE) %>% 
  dplyr::select(-"Order.y", -"Family.y", -"note", -"Notes")

dispersal_traits_total <- joined_traits

# EXPORT ==================================================================

## Data -------------------------------------------------------------------

write.csv(dispersal_traits_total, "data/data_process/dispersal_traits_total_preliminar.csv")

