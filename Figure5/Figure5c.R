#Load packages
library(tidyverse)
library(patchwork)
library(hrbrthemes)
library(circlize)
library(networkD3)

#Set working directory
setwd("~/Desktop/Dissertation/Chapter3")

#Import data and subset by geographic site groupings
data_long <- read.csv("CiliateCerc_table_long.csv")
data_long_NE <- subset(data_long, location == "Northeast")
data_long_NE <- data_long_NE[1:3]

data_long_NW <- subset(data_long, location == "Northwest")
data_long_NW <- data_long_NW[1:3]

data_long_SW <- subset(data_long, location == "Southwest")
data_long_SW <- data_long_SW[1:3]

# From these flows we need to create a node data frame: it lists every entities involved in the flow
nodes_SW <- data.frame(name=c(as.character(data_long_SW$CercozoaID), as.character(data_long_SW$CiliophoraID)) %>% unique())
nodes_NW <- data.frame(name=c(as.character(data_long_NW$CercozoaID), as.character(data_long_NW$CiliophoraID)) %>% unique())
nodes_NE <- data.frame(name=c(as.character(data_long_NE$CercozoaID), as.character(data_long_NE$CiliophoraID)) %>% unique())

# With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
data_long_SW$IDsource=match(data_long_SW$CercozoaID, nodes_SW$name)-1 
data_long_SW$IDtarget=match(data_long_SW$CiliophoraID, nodes_SW$name)-1

data_long_NW$IDsource=match(data_long_NW$CercozoaID, nodes_NW$name)-1 
data_long_NW$IDtarget=match(data_long_NW$CiliophoraID, nodes_NW$name)-1

data_long_NE$IDsource=match(data_long_NE$CercozoaID, nodes_NE$name)-1 
data_long_NE$IDtarget=match(data_long_NE$CiliophoraID, nodes_NE$name)-1

#define color scale
ColorScal ='d3.scaleOrdinal() .range(["#800026","fc6f35","#d23428","#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF","#440154FF"])'

# Make the Network
SW_Net <- sankeyNetwork(Links = data_long_SW, Nodes = nodes_SW,
              Source = "IDsource", Target = "IDtarget",
              Value = "value", NodeID = "name", 
              sinksRight=FALSE, colourScale=ColorScal, nodeWidth=60, fontSize=20, nodePadding=3)
SW_Net

NW_Net <- sankeyNetwork(Links = data_long_NW, Nodes = nodes_NW,
                        Source = "IDsource", Target = "IDtarget",
                        Value = "value", NodeID = "name", 
                        sinksRight=FALSE, colourScale=ColorScal, nodeWidth=60, fontSize=20, nodePadding=3)
NW_Net

NE_Net <- sankeyNetwork(Links = data_long_NE, Nodes = nodes_NE,
                        Source = "IDsource", Target = "IDtarget",
                        Value = "value", NodeID = "name", 
                        sinksRight=FALSE, colourScale=ColorScal, nodeWidth=60, fontSize=20, nodePadding=3)
NE_Net

