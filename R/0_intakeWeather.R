# ====================================================== #
#              NEON driver data download                 #
# ====================================================== #

library(neonstore)
library(tidyverse)
library(lubridate)

download.neon <- TRUE

jobs <- tibble(
  products = c(
    "DP1.00098.001", # relative humidity
    "DP1.00003.001", # triple aspirated air temperature
    "DP1.00005.001", # IR biological temperature (ground temp)
    "DP1.00094.001", # soil water content and water salinity
    "DP1.00041.001"),  # soil temperature
  table = c("RH_30min-basic","TAAT_30min-basic","IRBT_30_minute-basic","ST_30_minute-basic","SWS_30_minute-basic"),
  QF.col = c("RHFinalQF", "finalQF", "finalQF","finalQF","VSWCFinalQF"),
  min.col = c("RHMinimum", "tempTripleMinimum", "bioTempMinimum","soilTempMinimum","VSWCMinimum"),
  max.col = c("RHMaximum", "tempTripleMaximum", "bioTempMaximum","soilTempMaximum","VSWCMaximum"),
  var.col = c("RHVariance", "tempTripleVariance", "bioTempVariance","soilTempVariance","VSWCVariance"),
  exp.unc.col = c("RHExpUncert", "tempTripleExpUncert", "bioTempExpUncert","soilTempExpUncert","VSWCExpUncert"),
  new.min.var = c("RHMinimumVariance", "airTempMinimumVariance", "bioTempMinimumVariance","soilTempMinimumVariance","VSWCMinimumVariance"),
  new.min.exp.unc = c("RHMinimumExpUncert", "airTempMinimumExpUncert", "bioTempMinimumExpUncert","soilTempMinimumExpUncert","VSWCMinimumExpUncert"),
  new.max.var = c("RHMaximumVariance", "airTempMaximumVariance", "bioTempMaximumVariance","soilTempMaximumVariance","VSWCMaximumVariance"),
  new.max.exp.unc = c("RHMaximumExpUncert", "airTempMaximumExpUncert", "bioTempMaximumExpUncert","soilTempMaximumExpUncert","VSWCMaximumExpUncert"),
  csv.name = c("RelativeHumidityDaily.csv", "airTempDaily.csv", "bioTempDaily.csv", "soilTempDaily.csv", "soilWaterContentDaily.csv")
)

array.num <- as.numeric(Sys.getenv("SGE_TASK_ID")) # read array job number
variable <- jobs %>% slice(array.num)
product <- variable %>% pull(products)

site.coord <- readr::read_csv("Data/siteLatLon.csv") %>% suppressMessages()
site.vec <- site.coord %>% 
  pull(siteID) %>% 
  unique()

if(download.neon){
  message("Updating NEON Data")
  print(product)
  for(s in seq_along(site.vec)){
  message(paste("Downloading", site.vec[s], "..."))
    neon_download(product = product, 
                  # type = "expanded",
                  site = site.vec[s],
                  dir = Sys.getenv("NEONSTORE_HOME"),
                  .token = Sys.getenv("NEON_TOKEN"))
    print(warnings())  
    message(paste("Download Complete", site.vec[s]))
  }
  
  message("---- NEON Download complete ----")  
}


make_dataset <- function(table, QF.col, min.col, max.col, var.col, exp.unc.col, 
                         new.min.var, new.min.exp.unc, new.max.var, new.max.exp.unc, QF.flag = 0){
  
  cols.select <- c("Date", "Year", "domainID", "siteID",
                   min.col, max.col, var.col, exp.unc.col)
  
  cols.select <- cols.select[!is.na(cols.select)] # remove NAs
  
  df <- table %>% 
    filter(.data[[QF.col]] == QF.flag) %>% # remove observations that fail QF
    mutate(Date = date(endDateTime), # get date
           Year = year(Date)) %>%    # add year column
    select(all_of(cols.select))
  
  daily.min <- df %>% 
    group_by(siteID, Date) %>%
    slice(which.min(.data[[min.col]])) %>% # use slice to retain variance column
    select(-c(.data[[max.col]])) %>% 
    rename_at(vars(var.col), ~ new.min.var) %>% 
    rename_at(vars(exp.unc.col), ~ new.min.exp.unc)
  
  daily.max <- df %>% 
    group_by(siteID, Date) %>%
    slice(which.max(.data[[max.col]])) %>% # use slice to retain variance column
    select(-c(.data[[min.col]])) %>% 
    rename_at(vars(var.col), ~ new.max.var) %>% 
    rename_at(vars(exp.unc.col), ~ new.max.exp.unc)
  
  daily.data <- left_join(daily.max, daily.min)
  
  return(daily.data)
  
}


message(paste("Extracting", pull(variable, table)))
data <- neon_read(pull(variable, table))
df <- make_dataset(
  table = data,
  QF.col = pull(variable, QF.col),
  min.col = pull(variable, min.col),
  max.col = pull(variable, max.col),
  var.col = pull(variable, var.col),
  exp.unc.col = pull(variable, exp.unc.col),
  new.min.var = pull(variable, new.min.var),
  new.min.exp.unc = pull(variable, new.min.exp.unc),
  new.max.var = pull(variable, new.max.var),
  new.max.exp.unc = pull(variable, new.max.exp.unc)
)

message("Writing CSV")
write_csv(df, file.path("Data", pull(variable, csv.name)))

message("--- DONE ---")