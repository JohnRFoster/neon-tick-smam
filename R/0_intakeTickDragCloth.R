# this script is for processing neon tick data
# download with neonstore
# clean data into usable csv


library(tidyverse)
library(neonstore)
library(lubridate)
library(uuid)

# NEONSTORE_HOME, NEONSTORE_DB, and NEON_TOKEN defined in .Renviron


#### Last check 2021-10-14
# product <- "DP1.10093.001"
# neon_download(product = product)

tick.field.raw <- neon_read("tck_fielddata")
tick.taxon.raw <- neon_read("tck_taxonomyProcessed")

# there are lots of reasons why sampling didn't occur (logistics, too wet, too cold, etc.)
# so, keep records when sampling occurred
tick.field <- tick.field.raw %>% 
  filter(totalSampledArea > 0)


tick.taxon.wide <- tick.taxon.raw %>% 
  filter(sampleCondition == "OK") %>% # remove taxonomy samples with quality issues
  mutate(sexOrAge = if_else(sexOrAge == "Female" | sexOrAge == "Male", 
                            "Adult",     # convert to Adult
                            sexOrAge)) %>% 
  pivot_wider(id_cols = sampleID, # make wide by species and life stage
              names_from = c(acceptedTaxonID, sexOrAge),
              values_from = individualCount, 
              names_sep = "_",
              values_fn = {sum}, # duplicates occur because of Adults that where F/M - add them 
              values_fill = 0)

tick.joined <- left_join(tick.taxon.wide, tick.field, by = "sampleID") %>% 
  select(-NA_NA, -geodeticDatum, -samplingImpractical, -targetTaxaPresent,
         -adultCount, -nymphCount, -larvaCount, -samplingProtocolVersion, 
         -measuredBy, -sampleCode)

# all the species column names
spp.cols <- tick.joined %>% 
  select(contains("Larva"), contains("Nymph"), contains("Adult")) %>% 
  colnames()

# get matching taxon ids
taxon.ids <- tick.taxon.raw %>%
  filter(!is.na(acceptedTaxonID)) %>% 
  select(acceptedTaxonID, scientificName, taxonRank) %>% 
  distinct() 

# make longer
tick.long <- tick.joined %>% 
  pivot_longer(cols = all_of(spp.cols), 
               names_to = "taxonAge",
               values_to = "processedCount") %>% 
  separate(col = taxonAge, into = c("acceptedTaxonID", "lifeStage"), sep = "_")

# add taxon ids
tick.long <- left_join(tick.long, taxon.ids, by = "acceptedTaxonID")

write_csv(tick.joined, "Data/tickWide.csv")
write_csv(tick.long, "Data/tickLong.csv")
write_csv(taxon.ids, "Data/tickTaxonID.csv")


