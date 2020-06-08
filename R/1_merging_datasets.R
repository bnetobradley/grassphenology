## In this script, we merge files returned from ClimateNA with herbarium specimen records. We also:
# i) clean up binomials
# ii) standardize dates
# iii) remove records without dates/climate data
# iv) remove species with fewer than 10 flowering specimens

#### load packages needed ####
library(tidyverse)
library(data.table)
library(ggplot2)

#### read in year by year climateNA files and merge into one file ####
setwd("data/files_from_climate_na/")
clim_dat <- rbindlist(lapply(list.files(), fread))
setwd("../..")

#### read in herbarium database ####
herb_dat <- read_csv("data/herbarium_climate_data/20180612_ubcgrassdatabase.csv")
## clean up column names and remove obvious duplicates
herb_dat <- herb_dat %>% rename("source" = "X9", "anthers" = "X10")
herb_dat <- herb_dat[- which(duplicated(herb_dat) ==  TRUE), ]

#### merge herbarium data with climate data ####
flr_dat <- inner_join(herb_dat, clim_dat, by = c("accession" = "ID1"))
flr_dat <- mutate(flr_dat, binomial=paste(genus, species, sep = " "))

#### standardize dates ####
aux <- as.Date(flr_dat$date, format="%Y %b %d")
flr_dat$doy <- as.numeric(strftime(aux, format = "%j"))

#### remove unusable data
flr_dat <- flr_dat[-which(is.na(flr_dat$doy) == TRUE),]
flr_dat <- flr_dat[-which(flr_dat$Tmin02 == -9999.0),]

#### plot proportion of anthers from each species' specimens
flr_prop <- flr_dat %>% group_by(binomial) %>% count(anthers == "anthers")
flr_prop <- pivot_wider(flr_prop, names_from = `anthers == "anthers"`, values_from = n )
flr_prop$`FALSE` <- flr_prop$`FALSE` %>% replace_na(0)
flr_prop$`TRUE` <- flr_prop$`TRUE` %>% replace_na(0)
flr_prop$total <- flr_prop$`FALSE` + flr_prop$`TRUE`
flr_prop <- flr_prop[which(flr_prop$`TRUE` > 9),]
flr_prop <- pivot_longer(flr_prop, cols = `FALSE`:`TRUE`, names_to = "anthers")
#png("figures/s_fig_1.png", units="in", width=8, height=10, res=300)
ggplot(flr_prop, aes(x = reorder(binomial, desc(binomial)) , y = value, alpha = anthers)) + geom_col(stat = "identity", position = 'stack', fill = "#3a4660") + coord_flip() + scale_alpha_manual(values = c(0.5,1), breaks = c("FALSE", "TRUE"), name = "", labels = c("Total", "Anthers present")) + theme_classic() + theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = c(0.8, 0.9))+ labs(x = "", y = "Number of Herbarium Specimens") + theme(axis.text.y = element_text(margin = margin(r = 0))) 
#dev.off()

## save cleaned version of merged data in data folder
flr_dat <- flr_dat[flr_dat$binomial %in% flr_prop$binomial,]
flr_dat <- mutate(flr_dat, binomial=paste(genus, species, sep = "_"))
write.csv(flr_dat, file = "data/merged_herbarium_climate_na.csv", row.names=FALSE) 
