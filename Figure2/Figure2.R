#Load packages
library(ggplot2)

#Set working directory
setwd("~/Desktop/Dissertation/Chapter3/")

#Load OTU table
data <- read.csv("SAR_OTU_percent.csv")
#Convert OTU table to long format
data_long <- melt(data)
#Load metadata
meta <- read.csv("~/Desktop/DataAnalysis/Master_V9_MappingFile.csv")

#Rename column 2 as "sampleID" and merge OTU table and metadata
names(data_long)[2] <- "sampleID"
data_long <- merge(data_long, meta, by = "sampleID")
#Subset merged table to only Lombok sites
data_long <- subset(data_long, region == "Lombok")

#Reorder merged table rows by siteID and taxonomic group
data_long$siteID <- factor(data_long$siteID, levels = c("LBK_01","LBK_02","LBK_03",
                                                        "LBK_04","LBK_05","LBK_06",
                                                        "LBK_07","LBK_08","LBK_09",
                                                        "LBK_10","LBK_11","LBK_12",
                                                        "LBK_13","LBK_14","LBK_15",
                                                        "LBK_16","LBK_17","LBK_18"))
data_long$name <- factor(data_long$name, levels = c("Apicomplexa",
                                                    "Syndiniales",
                                                    "Dinophyceae",
                                                    "Bacillariophyta",
                                                    "Pelagophyceae",
                                                    "MAST",
                                                    "Cercozoa",
                                                    "Foraminifera",
                                                    "Radiolaria",
                                                    "Ciliophora",
                                                    "Other"))

#Define colors for taxa
colors <- c('Apicomplexa'="#fb4451",
            'Syndiniales'="#be0027",
            'Dinophyceae'="#a7c564",
            'Bacillariophyta'="#82a732",
            'Pelagophyceae'="#648746",
            'MAST'="#85cbda",
            'Cercozoa'="#0baeca",
            'Foraminifera'="#008bb8",
            'Radiolaria'="#005f84",
            'Ciliophora'="#003957",
            'Copepod'="#6473b4",
            'Other'="#999999")

#Plot community bar plots by sample type and filter size and save
plot <- ggplot(data_long, aes(siteID, value, fill = name)) + geom_bar(stat = "identity", position = "fill") +
  facet_grid(depth ~ filter) + theme(panel.background = element_rect(color = "black", fill = "white"), panel.grid = element_blank()) +
  theme(panel.border = element_rect(color = "black", fill = NA)) + theme(strip.background = element_rect(fill = "black", color = "black"), strip.text = element_text(color = "white", face = "bold")) +
  scale_fill_manual(values = colors) + theme(axis.text.x = element_blank()) +
  ylab("relative abundance") + xlab("\n") + scale_y_continuous(expand = c(0,0), limits = c(0,NA)) + theme(legend.position = "bottom")
plot
ggsave("LBK_community_barplot.pdf", plot, height = 6, width = 6)

#plot heatmap of rubble coverage per site and save
rubble_plot <- ggplot(data_long, aes(siteID, region, fill = rubble_percent)) +
  geom_point(pch = 22, size = 6, color = "black") + theme(panel.background = element_blank(), panel.grid = element_blank()) +
  theme(axis.title = element_blank(), axis.text = element_blank()) + theme(axis.ticks = element_blank()) +
  scale_fill_gradient(low = "white", high = "#002139")
rubble_plot
ggsave("LBK_rubble_heatmap.pdf", rubble_plot, height = 2.5, width = 5)
