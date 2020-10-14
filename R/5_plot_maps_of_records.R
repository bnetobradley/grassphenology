library(readr)
library(tidyverse)
library(maps)

## read in data
herb <- read_csv("data/merged_herbarium_climate_na.csv")

## clean up taxonomy and remove non flowering records
herb$binomial <- sub("Deschampsia_caespitosa", "Deschampsia_cespitosa", herb$binomial)
herb <- herb %>% filter(anthers == "anthers")

## check lat and lon
summary(herb$latitude)
summary(herb$longitude)

library(viridis)
library(gridExtra)

## plot map of all herbarium records
world <- map_data("world") 
nam <- world %>% filter(region == "Canada" | region == "USA")

dev.off()
png(file = "all_species_phenology_map.png", units="in", width=7, height=5, res=300)
ggplot() + geom_polygon(data = nam, aes(x=long, y = lat, group = group), fill = "#DDDDDD") + geom_jitter(data = herb, aes(x = longitude, y = latitude, color = doy), size = 0.1) + theme_minimal() + scale_color_viridis(begin = 1, end = 0.2, limits = c(72, 334)) +  xlim(-170, -50) + ylim(25, 85) + labs(x = "Longitude", y = "Latitude")+ scale_y_continuous(limits=c(25, 85), breaks=seq(40,80, 20)) + scale_x_continuous(limits=c(-170, -50), breaks=seq(-140, -60, 40)) + theme_minimal()
dev.off()

## species by species maps
sp <- sort(unique(herb$binomial))
p <- list()
for(i in 1:length(sp)) {
 p[[i]] <- ggplot() + geom_polygon(data = nam, aes(x=long, y = lat, group = group), fill = "#DDDDDD") + geom_jitter(data = herb %>% filter(binomial == sp[i]), aes(x = longitude, y = latitude, color = doy), size = 0.25) + scale_color_viridis(begin = 1, end = 0.2, limits = c(72, 334)) + labs(subtitle = paste(gsub("_", " ", sp[i])), x = element_blank(), y = element_blank()) + scale_y_continuous(limits=c(25, 85), breaks=seq(40,80, 20)) + scale_x_continuous(limits=c(-170, -50), breaks=seq(-140, -60, 40)) + theme_minimal() + theme(legend.position = "none") 
}
dev.off()
png(file = "species_phenology_map_1of6.png", units="in", width=7, height=9, res=300)
do.call(grid.arrange, c(ncol = 3, p[1:12], left = "Latitude", bottom = "Longitude"))
dev.off()
png(file = "species_phenology_map_2of6.png", units="in", width=7, height=9, res=300)
do.call(grid.arrange, c(ncol = 3, p[13:24], left = "Latitude", bottom = "Longitude"))
dev.off()
png(file = "species_phenology_map_3of6.png", units="in", width=7, height=9, res=300)
do.call(grid.arrange, c(ncol = 3, p[25:36], left = "Latitude", bottom = "Longitude"))
dev.off()
png(file = "species_phenology_map_4of6.png", units="in", width=7, height=9, res=300)
do.call(grid.arrange, c(ncol = 3, p[37:48], left = "Latitude", bottom = "Longitude"))
dev.off()
png(file = "species_phenology_map_5of6.png", units="in", width=7, height=9, res=300)
do.call(grid.arrange, c(ncol = 3, p[49:60], left = "Latitude", bottom = "Longitude"))
dev.off()
png(file = "species_phenology_map_6of6.png", units="in", width=7, height=9, res=300)
do.call(grid.arrange, c(ncol = 3, p[61:72], left = "Latitude", bottom = "Longitude"))
dev.off()