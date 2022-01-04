library(tidyverse)
library(lubridate)
library(daymetr)


site.coord <- readr::read_csv("Data/siteLatLon.csv") %>% suppressMessages()

mouse.data <- read_csv("Data/allSmallMammals.csv")
mouse.plots <- mouse.data %>% 
  group_by(plotID) %>% 
  select(plotID, decimalLatitude, decimalLongitude) %>% 
  distinct() %>% 
  mutate(plotID = paste0("smam", plotID))

tick.data <- read_csv("Data/tickLong.csv")
tick.plots <- tick.data %>% 
  group_by(plotID) %>% 
  select(plotID, decimalLatitude, decimalLongitude) %>% 
  distinct() %>% 
  mutate(plotID = paste0("tick", plotID))

plot.df <- bind_rows(mouse.plots, tick.plots)
plot.df <- plot.df %>% 
  rename_with( ~ c("site", "latitude", "longitude"), everything())
write_csv(plot.df, path = "Data/plotLatLon.csv")

s <- 1
dm <- download_daymet(
  site = plot.df$site[s],
  lat = plot.df$latitude[s],
  lon = plot.df$longitude[s],
  start = 2014,
  end = 2020,
  internal = TRUE,
  silent = FALSE,
  force = FALSE,
  simplify = TRUE
)
# 
# dm <- download_daymet_batch(
#   file_location = 'Data/plotLatLon.csv',
#   start = 2014,
#   end = 2021,
#   internal = TRUE
# )

write_csv(dm, path = "Data/daymet.csv")
