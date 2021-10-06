# this script is for processing neon small mammal data
# download with neonstore
# clean data into usable csv 
# most identifiers kept, removing body condition stuff
# keeping records of animals found in traps, dead or alive
# removes individuals only identified to Order or Class
# marked individuals have same tag if tag was replaced


library(tidyverse)
library(neonstore)
library(lubridate)
library(uuid)

dir.neonstore <- "/projectnb/dietzelab/fosterj/Data/neonstore/"
Sys.setenv("NEONSTORE_HOME" = dir.neonstore)
Sys.setenv("NEONSTORE_DB" = dir.neonstore)

smam.product <- "DP1.10072.001"

#### Last check 2021-10-06
# neon_download(product = smam.product,
#               .token = Sys.getenv("NEON_TOKEN"))

# tables - can join by "nightuid"
# plot.night <- neon_read("mam_perplotnight") # per trapping grid data
trap.night <- neon_read("mam_pertrapnight") # per trap night data

smam.data <- trap.night %>% 
  mutate(collectYear = year(collectDate),
         collectMonth = month(collectDate),
         collectYearMonth = paste(collectYear, collectMonth, sep = "-")) %>% 
  select(all_of(c("uid", # unique id
                  "nightuid", # night id, maps to plot.night
                  "collectDate",
                  "collectYear",
                  "collectMonth",
                  "collectYearMonth",
                  "domainID",
                  "siteID",
                  "plotID",
                  "nlcdClass", # land cover classification
                  "decimalLatitude",
                  "decimalLongitude",
                  "trapStatus", 
                  "trapCoordinate", # trap location on grid - A1
                  "tagID",
                  "taxonID", # four letter species/ID code??
                  "scientificName",
                  "taxonRank",
                  "identificationQualifier", # not sure
                  "recapture",
                  "fate",
                  "replacedTag",
                  "lifeStage",
                  "tickNumber",
                  "larvalTicksAttached",
                  "nymphalTicksAttached",
                  "adultTicksAttached"))) %>% 
  filter((!taxonRank %in% c("order", "class"))) %>% #, # get rid of high taxon classes
         # trapStatus %in% c("4 - more than 1 capture in one trap", 
                           # "5 - capture")) # keep only captured individuals
  mutate(animalStatus = if_else(trapStatus == "4 - more than 1 capture in one trap" | 
                                trapStatus == "5 - capture",
                                1, 0))

# want the number of unique animals each day
# will go by unique tag
# some tags have been replaced, deal with those first
# the majority of tags replaced are just left or right
# and the tag number is the same
# Some have a new tag number, need to demarcate those
# the tagID in the replacedTag column is the old tag

# arrange by collect date to make replacing easier
smam.data <- smam.data %>% 
  arrange(collectDate) 

tag.pattern <- "[[:upper:]]\\d{4}" # end of tagID - what is replaced; old tag
new.tag.rows <- grep(tag.pattern, smam.data$replacedTag) # rows when new tags used
old.tag.id <- smam.data$replacedTag[new.tag.rows] # old tag Ids
new.tag.id <- smam.data$tagID[new.tag.rows] # the tags that were replaced

for(i in seq_along(new.tag.id)){
  
  smam.subset <- smam.data %>%
    filter(tagID == new.tag.id[i])
  
  if(nrow(smam.subset) > 0){ # revert back to old tag
    index <- which(smam.data$tagID == new.tag.id[i])
    smam.data$tagID[index] <- gsub(tag.pattern, old.tag.id[i], new.tag.id[i]) 
  }
}

# TO_DO
# check to make sure at least one animal was captured each trap night
smam.data %>% 
  group_by(plotID, collectDate) %>% 
  summarise()



write_csv(smam.data, "Data/allSmallMammals.csv")

