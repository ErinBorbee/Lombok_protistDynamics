#load packages
library(ggplot2)

#Set working directory
setwd("~/Desktop/DataAnalysis/V9mock/")
#Import data
data <- read.csv("level-9.csv")
#Reorder fill column
data$fill <- factor(data$fill, levels = c("Bacillariophyta",
                                          "Dictyochophyceae",
                                          "Chrysophyceae",
                                          "Phaeothamnion",
                                          "Xanthophyceae",
                                          "Eustigmatophyceae",
                                          "Alveolata",
                                          "Archaeplastida",
                                          "Opisthokonta",
                                          "Eukaryota undefined"))

#Define colors
colors <- c('Bacillariophyta' = "#00909f",
            'Dictyochophyceae' = "#6c073a",
            'Chrysophyceae' = "#f4545b",
            'Phaeothamnion' = "#ffaf52",
            'Xanthophyceae' = "#6bae5e",
            'Eustigmatophyceae' = "#004e41",
            'Alveolata' = "grey80",
            'Archaeplastida' = "grey20",
            'Opisthokonta' = "grey40",
            'Eukaryota undefined' = "grey60")

plot <- ggplot(data, aes(name, percent, fill = fill)) + 
  geom_bar(stat = "identity", position = "stack", width = 0.25) + 
  theme_bw() + theme(panel.grid = element_blank()) + 
  theme(panel.border = element_blank()) +
  scale_fill_manual(values = colors) + 
  theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(),
        axis.text.x = element_text(face = "bold", size = 12)) +
  theme(axis.ticks.length.y = unit(0.1, "in"),
        axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12)) + 
  ylab("Relative Abundance \n") +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 12), 
        legend.position = "right") + 
  theme(axis.line.y = element_line(color = "black")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))
plot
ggsave("V9_MockCommunity_barplot.png", width = 6, height = 8)

