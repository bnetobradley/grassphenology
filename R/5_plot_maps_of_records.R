library(readr)
library(tidyverse)
library(maps)
library(mapr)
library(spocc)
library(viridis)
library(gridExtra)

## read in data
herb <- read_csv("data/merged_herbarium_climate_na.csv")

## clean up taxonomy and remove non flowering records
herb$binomial <- sub("Deschampsia_caespitosa", "Deschampsia_cespitosa", herb$binomial)
herb <- herb %>% filter(anthers == "anthers")

## check lat and lon
summary(herb$latitude)
summary(herb$longitude)

## plot map of all herbarium records
world <- map_data("world") 
nam <- world %>% filter(region == "Canada" | region == "USA")

dev.off()
png(file = "all_species_phenology_map.png", units="in", width=7, height=5, res=300)
ggplot() + geom_polygon(data = nam, aes(x=long, y = lat, group = group), fill = "#DDDDDD") + geom_jitter(data = herb, aes(x = longitude, y = latitude, color = doy), size = 0.1) + theme_minimal() + scale_color_viridis(begin = 1, end = 0.2, limits = c(72, 334)) +  xlim(-170, -50) + ylim(25, 85) + labs(x = "Longitude", y = "Latitude")+ scale_y_continuous(limits=c(25, 85), breaks=seq(40,80, 20)) + scale_x_continuous(limits=c(-170, -50), breaks=seq(-140, -60, 40)) + theme_minimal()
dev.off()

## create species maps that show known occurrence range using R package 'spocc'

# create lists/vectors to search and dump data into
sp <- sort(unique(herb$binomial))
p <- list()
occurence_dat<- c()
clean_occ <- c()

## create dataframe with north american gbif species occurrences (lat- and longitude) including the species name by which the record was originally searched for (noted in the column 'binomial')
for(i in 1:length(sp)) {
occurence_dat <- occ(query = sp[i], from = c('gbif'), geometry = c(-170, 25, -50, 85))
clean <- occ2df(occurence_dat)
clean$binomial <- paste(sp[i])
clean_occ <- rbind(clean_occ, clean)
}

## plot herbarium records in flower, with 2d density plots of those species' records from gbif in the background (as a proxy for range)
for(i in 1:length(sp)) {
  num <- herb %>% filter(binomial == sp[i]) %>% count()
 p[[i]] <- ggplot() + geom_polygon(data = nam, aes(x=long, y = lat, group = group), fill = "#DDDDDD") + geom_jitter(data = herb %>% filter(binomial == sp[i]), aes(x = Longitude, y = Latitude, color = doy), size = 0.25) + scale_color_viridis(begin = 1, end = 0.2, limits = c(72, 334)) + labs(title = paste(gsub("_", " ", sp[i])), subtitle = paste0(num, " specimens"), x = element_blank(), y = element_blank()) +
 scale_y_continuous(limits=c(25, 85), breaks=seq(40,80, 20)) + scale_x_continuous(limits=c(-170, -50), breaks=seq(-140, -60, 40)) + theme_minimal() + theme(legend.position = "none") + theme(plot.title = element_text(face = 'italic', size = 10), plot.subtitle = element_text(size = 10)) + stat_density_2d(data = clean_occ %>% filter(binomial == sp[i]), aes(x = longitude, y= latitude), geom = "polygon",alpha = 0.1, fill = "#403891") 
 print(sp[i])
 }

## save plots
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
