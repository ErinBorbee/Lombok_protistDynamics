#Load packages
library(maptools)
library(raster)
library(maps)
library(GISTools)
library(prettymapr)

#Set working directory
setwd("~/Desktop/Dissertation/Chapter3/")

#Load metadata and subset to Lombok only
meta <- read.csv("~/Desktop/DataAnalysis/Master_V9_MappingFile.csv")
LBK_points <- subset(meta, region == "Lombok")

#Load country shapefile data
Indo <- getData('GADM', country ='IDN', level = 0)

#Plot Lombok map with site points and record plot
plot(Indo, col = "grey57", border = NA, bg = "white", axes=T, xlim = c(115.5,117), ylim = c(-9, -8))
points(LBK_points$long, LBK_points$lat, pch = 21, col = "white", bg = "black", cex = 3, lwd = 1)
LBK_map <- recordPlot()

#Plot Indonesia map for inset
#Load country shapfiles (level zero is country level, subsequent levels are states, provinces, etc. depending on the level of detail you need)
Indo <- getData('GADM', country ='IDN', level = 0)
Malay <- getData('GADM', country = 'MYS', level = 0)
Timor <- getData('GADM', country = 'TLS', level = 0)
Papua <- getData('GADM', country = 'PNG', level = 0)
Solomon <- getData('GADM', country = 'SLB', level = 0)
China <- getData('GADM', country = 'CHN', level = 0)
Taiwan <- getData('GADM', country = 'TWN', level = 0)
Phil <- getData('GADM', country = 'PHL', level = 0)
Vietnam <- getData('GADM', country = 'VNM', level = 0)
Camb <- getData('GADM', country = 'KHM', level = 0)
Laos <- getData('GADM', country = 'LAO', level = 0)
Aus <- getData('GADM', country = 'AUS', level = 0)
Thai <- getData('GADM', country = 'THA', level = 0)
Sing <- getData('GADM', country = 'SGP', level = 0)
Myan <- getData('GADM', country = 'MMR', level = 0)
Brun <- getData('GADM', country = 'BRN', level = 0)
Bhut <- getData('GADM', country = 'BTN', level = 0)
Bang <- getData('GADM', country = 'BGD', level = 0)
India <- getData('GADM', country = 'IND', level = 0)
Nepal <- getData('GADM', country = 'NPL', level = 0)

#Plot the different countries in designated colors, using xlim and ylim to help plot the region of interest faster
plot(Indo, col = "grey25", border = NA, bg = "white", axes=T, xlim = c(95,145), ylim = c(-12, 8))
plot(Malay, col = "grey90", border = "white", lwd = 4, add = T)
plot(Timor, col = "grey90", border = "white", lwd = 4, add = T)
plot(Papua, col = "grey90", border = "white", lwd = 4, add = T)
plot(Solomon, col = "grey90", border = "white", lwd = 4, add = T)
plot(Phil, col = "grey90", border = "white", lwd = 4, add = T)
plot(China, col = "grey90", border = "white", lwd = 4, add = T)
plot(Taiwan, col = "grey90", border = "white", lwd = 4, add = T)
plot(Vietnam, col = "grey90", border = "white", lwd = 4, add = T)
plot(Camb, col = "grey90", border = "white", lwd = 4, add = T)
plot(Laos, col = "grey90", border = "white", lwd = 4, add = T)
plot(Aus, col = "grey90", border = "white", lwd = 4, add = T)
plot(Thai, col = "grey90", border = "white", lwd = 4, add = T)
plot(Sing, col = "grey90", border = "white", lwd = 4, add = T)
plot(Brun, col = "grey90", border = "white", lwd = 4, add = T)
plot(Myan, col = "grey90", border = "white", lwd = 4, add = T)
plot(India, col = "grey90", border = "white", lwd = 4, add = T)
plot(Bhut, col = "grey90", border = "white", lwd = 4, add = T)
plot(Bang, col = "grey90", border = "white", lwd = 4, add = T)
plot(Nepal, col = "grey90", border = "white", lwd = 4, add = T)
Indo_map <- recordPlot()

plot.new()
png("Indo_Map.png", height = 100, width = 100)
Indo_map
dev.off()

