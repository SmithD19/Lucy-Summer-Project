##################################################
## Project:
## Script purpose:
## Date: 2020-08-24
## Author: Daniel Smith
## Contact: dansmi@ceh.ac.uk
## Licence: MIT
##################################################

library(tidyverse)
library(read.so)

# This package reads MD formatted tables like the ones Access outputs 
# to .txt files
oviposition_data <- read.so::read_md("Oviposition.txt")

# New code to change artificial and ground pool to better names
oviposition_data <- oviposition_data %>% mutate(Microsite = str_replace(Microsite, "Artificial", "ArtificialCon"),
                                                Microsite = str_replace(Microsite, "Ground Pool", "GroundPool"))

# Are any values NA in the ID column? (Blank lines)
if (any(is.na(oviposition_data$Oviposition_ID))) {
  stop(print("NA lines in the data file"))
}

# Function to run whole process of:
  # Splitting rows
  # Calculating frequency
  # Generating DF to work with this

cleanmossie <- function(data, var, group) {
  # Make them names so they can be passed to dplyr
  group = as.name(group)
  var = as.name(var)
  # Seperate rows of DF
  sep = separate_rows(data, {{var}})
  # Count pairings
  cnt = count(sep, {{group}}, {{var}})
  # Relative frequencies
  freq = cnt %>%
    # Gruop by
    group_by({{group}}) %>%
    # Now mutate works along the groups
    mutate(
      # Get the total number of records
      total = sum(n),
      # Use that to calculate new frequency column
      freq = round(n / total * 100, digits = 2),
      # you can check this worked by summing and seeing if it equals 100:
      check = sum(freq)
    )
  return(freq)
}

# Summarise Each mosquito variable we looked at using a for loop:

# These are the variables we looked at:
vars <- c("Broad_Habitat_Level_1", "Broad_Habitat_Level_2", "Microsite",
          "Water_Source", "Water_Permanence")

# Initialise a list to store these in:
mosquito_summary <- list()

# Now run the loop
for (i in seq_along(vars)) {
  mosquito_summary[[i]] <- cleanmossie(oviposition_data, vars[[i]], "Mosquito_Species")
}

# Is by country useful? If so it's here:
country_summary <- list()

for (i in seq_along(vars)) {
  country_summary[[i]] <- cleanmossie(oviposition_data, vars[[i]], "Country")
  names(country_summary[i]) <- vars[[i]]
}

##################################################
## Project:
## Script purpose:
## Date: 2020-08-24
## Author: Daniel Smith
## Contact: dansmi@ceh.ac.uk
## Licence: MIT
##################################################

library(tidyverse)
library(janitor)
library(ape)

# Summary Data ------------------------------------------------------------

# Make the data one line for each mosquito speceis, turn these values 
# into 1/0 binary values for a logistic regression

# Broad habitat level 1
Broad_Habitat_Level_1 <- mosquito_summary[[1]] %>%
  na.omit() %>%
  mutate(binary = ifelse(freq > 0, 1, 0)) %>%
  select(Mosquito_Species, Broad_Habitat_Level_1, binary) %>%
  pivot_wider(
    names_from = Broad_Habitat_Level_1,
    values_from = binary,
    values_fill = 0,
    names_prefix = "hab_1_"
  ) %>%
  janitor::clean_names()

# Level 2 habitat
Broad_Habitat_Level_2 <- mosquito_summary[[2]] %>%
  na.omit() %>%
  mutate(binary = ifelse(freq > 0, 1, 0)) %>%
  select(Mosquito_Species, Broad_Habitat_Level_2, binary) %>%
  pivot_wider(
    names_from = Broad_Habitat_Level_2,
    values_from = binary,
    values_fill = 0,
    names_prefix = "hab_2_"
  ) %>%
  janitor::clean_names()

# Microsites
Microsite <- mosquito_summary[[3]] %>%
  na.omit() %>%
  mutate(binary = ifelse(freq > 0, 1, 0)) %>%
  select(Mosquito_Species, Microsite, binary) %>%
  pivot_wider(names_from = Microsite,
              values_from = binary,
              values_fill = 0) %>%
  janitor::clean_names()

# Water sources
Water_Source <- mosquito_summary[[4]] %>%
  mutate(binary = ifelse(freq > 0, 1, 0)) %>%
  select(Mosquito_Species, Water_Source, binary) %>%
  pivot_wider(names_from = Water_Source,
              values_from = binary,
              values_fill = 0) %>%
  janitor::clean_names()

# Water permanences
Water_Permanence <- mosquito_summary[[5]] %>%
  mutate(binary = ifelse(freq > 0, 1, 0)) %>%
  select(Mosquito_Species, Water_Permanence, binary) %>%
  pivot_wider(names_from = Water_Permanence,
              values_from = binary,
              values_fill = 0) %>%
  janitor::clean_names()

# Phylogenetic Data -------------------------------------------------------

# Read in the generated phylogenetic tree from COI genes
phylo <- read.tree("phylotree_COI.new")

# Here are the tip labels for the tree
mosquitoes <-
  phylo$tip.label %>% str_replace_all("_", " ") %>% word(1, 2)

# We need to change Ochlerotatus to Aedes
mosquitoes <- mosquitoes %>% str_replace("Ochlerotatus", "Aedes")

# Replace the tiplabels with the new changed names
phylo$tip.label <- mosquitoes

# Virus interaction -------------------------------------------------------

# These are flavivirus interactions from EID2 database
int <-
  read.csv(
    "https://raw.githubusercontent.com/SmithD19/MiniProject1.1/master/data/species-master-interaction.csv"
  )

# Mosquito names are in the variable scientficname.y
interactions <- int %>%
  # Lets get only the ones that are present in our database and phylogeny
  filter(scientificname.y %in% mosquitoes) %>%
  # Any NA values wave goodbye
  na.omit() %>%
  # The names of mosquito and virus names
  select(mosquito = scientificname.y, virus = virus_name) %>%
  # Replace Ochlerotatus with Aedes
  mutate(mosquito = str_replace(mosquito, "Ochlerotatus", "Aedes"))

# Which mosquitoes aren't in the data? negate the %in% function
# this makes it do the opposite: so... %not-in% 
`%nin%` <- negate(`%in%`)

# Use %nin% to get mosquitoes that aren't in any virus interactions
virdat1 <-
  data.frame(mosquito_species = mosquitoes[mosquitoes %nin% interactions$mosquito])

# Get the ones that are in virus interactions
virdat2 <-
  interactions %>% select(mosquito_species = mosquito, virus = virus)

# Bind these together and make a matrix that shows 1 or 0 depending on whether
# it transmits
virus_mat <- bind_rows(virdat1, virdat2) %>%
  mutate(vector = ifelse(is.na(virus), 0, 1)) %>%
  pivot_wider(names_from = virus,
              values_from = vector,
              values_fill = 0) %>%
  janitor::clean_names() %>%
  select(-na)

# these are the virus colnames
colnames(virus_mat)[-1]

# Which ones are vectors?
vector <- virus_mat %>%
  group_by(mosquito_species) %>%
  summarise(# If the sum of all the vector counts is above 1, just give it a 1, else a 0
    vector = if_else(sum(across(
      colnames(virus_mat)[-1]
    )) > 1, 1, 0))

# Geodata -----------------------------------------------------------------

# Here is the data that Francis gave me from his paper about occurrence of
# some mosquito species across Europe
geo <-
  read.csv(
    "https://raw.githubusercontent.com/SmithD19/MiniProject1.1/master/data/occurance.csv"
  )

# Make it so the species names are only two words to get rid of any subspecies
geo <- geo %>% mutate(species = word(species, 1, 2))

# Lets only take the records for species that we have in the database
geo <- geo %>% filter(species %in% mosquitoes)

# Some data cleaning
geodata <- geo %>%
  # Change it from wide > long format
  pivot_longer(-species) %>%
  # Mutate when the numbers for this value = this thing
  mutate(
    value = case_when(
      value == "0" ~ "Absent",
      value == "1" ~ "Extinct",
      value == "2" ~ "Uncertain",
      value == "3" ~ "Introduced",
      value == "4" ~ "Native",
    )
  ) %>%
  # Give me a count of all the combinations (How many places is a mozzie absent? etc)
  count(species, value) %>%
  # Only take these Introduced and Native values. Others won't tell us much
  filter(value == "Introduced" | value == "Native") %>%
  # Now change the data back to a wide format ready for modeling
  pivot_wider(
    species,
    names_from = value,
    values_from = n,
    values_fill = 0
  ) %>%
  # Rename species to mosquito_species
  rename(mosquito_species = species)

# Habitat/country breadth and binary invasiveness in Europe
geodata <- geodata %>%
  mutate(
    geo_breadth = Introduced + Native,
    Introduced = if_else(Introduced > 1, 1, 0),
    Native = if_else(Native > 1, 1, 0)
  )

# Lat/Lon Data ------------------------------------------------------------

# # Max min lat of countries
# minmaxlat <- read.csv("qgis/MinMaxLat.csv")
# library(countrycode)
#
# # I've realized that the  geographical units in frnacis' data arent countries. 
# # This makes sense, Francis paper split up the big Greek islands and such.
# # obviously the country boundaries are useless for mosquitoes.
# # Unsure how to get this working because island max/min lats are hard to get
# # But if you found the added islands in this dataset 
# # and then addded them to the data you could use these
#
# countries1 <- minmaxlat$CNTRY_NAME
#
# countries2 <- geo %>% pivot_longer(-species) %>% pull(name) %>% unique()
#
# code1 <- countrycode::countrycode(countries1, origin = "country.name", destination = "country.name")
#
# code2 <- countrycode::countrycode(countries2, origin = "country.name", destination = "country.name")

# Assemble all data -------------------------------------------------------

# Join all the data frames together (except whole vrius mat)
df <- Water_Permanence %>%
  full_join(Water_Source) %>%
  full_join(Broad_Habitat_Level_1) %>%
  full_join(Microsite) %>% 
  full_join(vector) %>%
  full_join(virus_mat) %>% 
  full_join(geodata) %>%
  # Last three have NA values in them because not all species in database
  na.omit() %>%
  # Also ungroup or leads to nastiness
  ungroup()

# This one uses broader habitat characteristics 
write_csv(df, "Oviposition_Summary_Hab1.csv")

# Join all the data frames together (except whole vrius mat)
df <- Water_Permanence %>%
  full_join(Water_Source) %>%
  full_join(Broad_Habitat_Level_2) %>%
  full_join(Microsite) %>% 
  full_join(vector) %>%
  full_join(virus_mat) %>% 
  full_join(geodata) %>%
  # Last three have NA values in them because not all species in database
  na.omit() %>%
  # Also ungroup or leads to nastiness
  ungroup()

# This one uses less broad (?) habitat characteristics 
write_csv(df, "Oviposition_Summary_Hab2.csv")





















