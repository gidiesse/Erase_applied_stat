# world_plots: This is the file we used to make the plots where we plotted our 
# observations on top of maps from google maps

# loading the required packages
library(ggplot2)
library(ggmap)
setwd("~/Desktop/ERASE_project2023_desgrp")


# import data
individuals.deconvolution <- read.csv("corr_convolution.csv")
# drop NAs
na.lat_lon <- is.na(individuals.deconvolution$lat_num)
individuals.deconvolution <- individuals.deconvolution[!na.lat_lon,]
na.macroarea <- is.na(individuals.deconvolution$Macroarea)
individuals.deconvolution <- individuals.deconvolution[!na.macroarea,]
na.macroperiod <- is.na(individuals.deconvolution$Macro_Period)
individuals.deconvolution <- individuals.deconvolution[!na.macroperiod,]

# df for mapping
df <- data.frame(lon=individuals.deconvolution$long_num, 
                 lat=individuals.deconvolution$lat_num)

# getting the map
mapgilbert <- get_map(location = c(lon = mean(df$lon), lat = mean(df$lat)), 
                      zoom = 2, maptype = "satellite", scale = 2)

# plotting the map with some points on it
ggmap(mapgilbert) +
  geom_point(data = df, aes(x = lon, y = lat, 
  color = as.factor(individuals.deconvolution$Macroarea),
  alpha = 1.0), 
  size = 1.5, shape = 21, show.legend=T) +
  scale_x_continuous(limits = c(-38, 115), expand = c(0, 0)) +
  scale_y_continuous(limits = c(25, 72), expand = c(0, 0)) +
  scale_color_manual(values = c('orange', 'chartreuse', 'forestgreen', 'cyan', 
                                'yellow', 'steelblue1', 'red', 'mediumblue', 
                                'purple', 'white', 'tan', 'magenta')) +
  guides(fill=FALSE, alpha=FALSE, size=FALSE) +
  labs(colour="Macroarea")

# Now we want to plot only the individuals with an ancestry label
our_ancestries <- c('CHG', 'EHG', 'WHG', 'Anatolia_N', 'North_Africa')
idxs <- c()
for (i in our_ancestries) {
  idxs <- c(idxs, which(individuals.deconvolution$Main_Ancestry == i))
}
df_labeled <- data.frame(lon=individuals.deconvolution$long_num[idxs], 
                         lat=individuals.deconvolution$lat_num[idxs])

# getting the map
mapgilbert <- get_map(location = c(lon = mean(df_labeled$long), lat = mean(df_labeled$lat)), 
                      zoom = 2, maptype = "satellite", scale = 2)
# plotting the map with some points on it
ggmap(mapgilbert) +
  geom_point(data = df_labeled, aes(x = lon, y = lat, 
                            color = as.factor(individuals.deconvolution$Main_Ancestry[idxs]),
                            fill = as.factor(individuals.deconvolution$Main_Ancestry[idxs]),
                            alpha = 1.0), 
             size = 1.5, shape = 21, show.legend=T) +
  scale_x_continuous(limits = c(-38, 115), expand = c(0, 0)) +
  scale_y_continuous(limits = c(25, 72), expand = c(0, 0)) +
  scale_color_manual(values = c('red', 'orange', 'yellow', 'chartreuse', 
                                'cyan')) +
  scale_fill_manual(values = c('red', 'orange', 'yellow', 'chartreuse', 
                               'cyan')) + 
  guides(fill=FALSE, alpha=FALSE, size=FALSE) +
  labs(colour="Ancestry")


# getting the map
mapgilbert <- get_map(location = c(lon = 12.49, lat = 41.9), 
                      zoom = 3, maptype = "toner", scale = 2)
ggmap(mapgilbert) + scale_x_continuous(limits = c(-38, 60), expand = c(0, 0)) +
  scale_y_continuous(limits = c(25, 70), expand = c(0, 0)) 

### Map between NAs ###

metadata <- read.csv('data/erase_metadata.csv')
metadata.anc <- metadata[metadata$project=='ancient', ]
metadata.anc <- metadata.anc[!is.na(metadata.anc$long_num),]

# df for mapping
df <- data.frame(lon=metadata.anc$long_num, 
                 lat=metadata.anc$lat_num)

# getting the map
mapgilbert <- get_map(location = c(lon = mean(df$lon), lat = mean(df$lat)), 
                      zoom = 2, maptype = "satellite", scale = 2)

# plotting the map with some points on it
ggmap(mapgilbert) +
  geom_point(data = df, aes(x = lon, y = lat, 
                            color = as.factor(metadata.anc$Macroarea),
                            alpha = 1.0), 
             size = 1.5, shape = 21, show.legend=T) +
  scale_x_continuous(limits = c(-38, 115), expand = c(0, 0)) +
  scale_y_continuous(limits = c(25, 72), expand = c(0, 0)) +
  scale_color_manual(values = c('chartreuse', 'forestgreen', 'cyan', 
                                'yellow', 'steelblue1', 'red', 'mediumblue', 
                                'purple', 'white', 'tan', 'magenta'),
                     na.value = 'black') +
  guides(fill=FALSE, alpha=FALSE, size=FALSE) +
  labs(colour="Macroarea")



