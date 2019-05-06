library(ggplot2)  # FYI you need v2.0
library(dplyr)    # yes, i could have not done this and just used 'subset' instead of 'filter'
library(ggalt)    # devtools::install_github("hrbrmstr/ggalt")
library(ggthemes) # theme_map and tableau colors
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(googleway)
library(ggspatial)
library(rgeos)
library(hrbrthemes)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

# source base theme and colors for plots
source("Scripts/Figure_Themes.R")

# source scripts 
source("Scripts/PD_Talk_Functions.R")


world <- ne_countries(scale='medium',returnclass = 'sf')

wrld_pt <- ggplot(data = world) +
  base_theme +
  geom_sf(fill = background_color, inherit.aes = T) +
  coord_sf(ylim = c(-50, 90)) +
  theme(panel.background = element_rect(fill = background_color, colour = NA),
        panel.border = element_rect(fill = NA),
        plot.background = element_rect(fill = background_color))

europe <- world[world$region_un=="Europe"&world$name!='Russia',]

# A bounding box for continental Europe.
europe.bbox <- st_polygon(list(
  matrix(c(-25,29,45,29,45,75,-25,75,-25,29),byrow = T,ncol = 2)))
europe.clipped <- suppressWarnings(st_intersection(europe, st_sfc(europe.bbox, crs=st_crs(europe))))

euro_pt <- ggplot(europe.clipped) +
  geom_sf(fill = "antiquewhite", inherit.aes = T) +
  base_theme +
  coord_sf(crs="+proj=aea +lat_1=36.333333333333336 +lat_2=65.66666666666667 +lon_0=14") +
  theme(panel.background = element_rect(fill = background_color, colour = NA),
        panel.border = element_rect(fill = NA),
        plot.background = element_rect(fill = background_color))

(hawaii  <- ggplot(data = world) +
    geom_sf(fill = "antiquewhite", inherit.aes = T) +
    coord_sf(crs = st_crs(4135), xlim = c(-161, -154), ylim = c(18, 23), expand = FALSE, datum = NA)+
    base_theme +
    theme(panel.background = element_rect(fill = background_color),
          panel.border = element_rect(fill = NA),
          plot.background = element_rect(fill = background_color)))

# download data from CeNDR
isolation_info <- readr::read_tsv("https://elegansvariation.org/strain/strain_data.tsv")

strains_330 <- isolation_info%>%
  dplyr::filter(reference_strain == "True")%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude) %>%
  dplyr::filter(long != "None")

sites <- st_as_sf(data.frame(longitude = as.numeric(strains_330$long), latitude = as.numeric(strains_330$lat)), coords = c("longitude", "latitude"), crs = 4326, 
                  agr = "constant")

ggplot(data = world) +
  base_theme +
  geom_sf(fill = "antiquewhite", inherit.aes = T) +
  geom_sf(data = sites, size = 2, shape = 25, fill = highlight_color) +
  # coord_sf(crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs ") +
  theme(panel.background = element_rect(fill = "aliceblue"))

ggsave(filename = "Plots/World_Plot.pdf", height = 10, width = 14)

ggplot(data = world) +
  base_theme +
  geom_sf(fill = "antiquewhite", inherit.aes = T) +
  geom_sf(data = sites, size = 2, shape = 25, fill = highlight_color) +
  coord_sf(crs = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs ") +
  theme(panel.background = element_rect(fill = "aliceblue"))

ggsave(filename = "Plots/World_Plot_ROUND.pdf", height = 10, width = 10)

hawaii + geom_sf(data = sites, size = 2, shape = 25, fill = highlight_color) +
  coord_sf(crs = st_crs(4135), xlim = c(-161, -154), ylim = c(18, 23), expand = FALSE, datum = NA)+
  theme(panel.background = element_rect(fill = "aliceblue"))

ggsave(filename = "Plots/Hawaii.pdf", height = 10, width = 10)

strains_330 <- isolation_info%>%
  dplyr::filter(reference_strain == "True")%>%
  dplyr::select(strain = isotype, long = longitude, lat = latitude) %>%
  dplyr::filter(long != "None")

euro_pt + 
  geom_sf(data = sites, size = 2, shape = 25, fill = highlight_color) + 
  coord_sf(crs = st_crs(4135), xlim = c(-15, 40), ylim = c(30, 70), expand = FALSE, datum = NA) +
  theme(panel.background = element_rect(fill = "aliceblue"))

ggsave(filename = "Plots/Europe.pdf", height = 10, width = 10)
