#Load packages
library(ggplot2)
library(reshape2)
library(dplyr)
library(ggpubr)

#Set working directory
setwd("~/Desktop/Dissertation/Chapter3/")

#Load metadata and subset to benthic data of interest (make sure sites are in numeric order 01-18)
meta <- read.csv("~/Desktop/Dissertation/Chapter3/benthic_data.csv")
meta <- subset(meta, substrate %in% c("Hard Coral",
                                      "Soft Coral",
                                      "Rubble","Algae"))

#Add column for coverage grouping to combine soft and hard coral cover
meta <- meta %>% mutate(coverage = case_when(substrate %in% c("Hard Coral","Soft Coral") ~ "Coral Cover",
                                     substrate == "Algae" ~ "Macroalgal Cover",
                                     substrate == "Rubble" ~ "Rubble Cover") )

#Subset individual coverage data by type
coral_meta <- subset(meta, substrate %in% c("Hard Coral","Soft Coral"))
algae_meta <- subset(meta, substrate == "Algae")
rubble_meta <- subset(meta, substrate == "Rubble")

#Define colors for plot
colors <- c('Hard Coral' = "#e7404b",
            'Soft Coral' = "#ff867c",
            'Rubble' = "#8270b2",
            'Algae' = "#a6c54b")

#plot each coverage type as separate bar plot (rubble plot includes horizontal line for 30% rubble threshold established in Sawall, etal 2013)
coral_plot <- ggplot(coral_meta, aes(site, average, fill = substrate)) +
  geom_bar(stat = "identity", color = "black") + theme_bw() +
  theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) + scale_fill_manual(values = colors) +
  theme(strip.background.x = element_blank(), strip.background.y = element_rect(fill = "black"), strip.text.x = element_blank(), strip.text.y = element_text(color = "white", face = "bold")) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + ylab("percent coverage") + theme(legend.position = "none") +
  facet_grid(coverage~region) + theme(axis.ticks.x = element_blank())
coral_plot

algae_plot <- ggplot(algae_meta, aes(site, average, fill = substrate)) +
  geom_bar(stat = "identity", color = "black", position = "stack", alpha = 0.7) + theme_bw() +
  theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) + scale_fill_manual(values = colors) +
  theme(strip.background.x = element_blank(), strip.background.y = element_rect(fill = "black"), strip.text.x = element_blank(), strip.text.y = element_text(color = "white", face = "bold")) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + ylab("percent coverage") + theme(legend.position = "none") +
  facet_grid(coverage~region) + theme(axis.ticks.x = element_blank())
algae_plot

rubble_plot <- ggplot(rubble_meta, aes(site, average, fill = substrate)) +
  geom_bar(stat = "identity", color = "black", position = "stack") + theme_bw() +
  theme(panel.grid = element_blank()) + theme(axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,100)) + scale_fill_manual(values = colors) +
  theme(strip.background.x = element_blank(), strip.background.y = element_rect(fill = "black"), strip.text.x = element_blank(), strip.text.y = element_text(color = "white", face = "bold")) +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank()) + ylab("percent coverage") + theme(legend.position = "none") +
  facet_grid(coverage~region) + geom_hline(yintercept = 30, color = "black", linetype = "dashed", size = 0.5)
rubble_plot

#arrange bar plots to stack on top of one another and save
cover_plot <- ggarrange(algae_plot, coral_plot, rubble_plot,
                        nrow = 3, ncol = 1, align = "hv") + theme_transparent()
cover_plot
ggsave("~/Desktop/NE_benthicCoverage_plot.png", cover_plot, width = 4, height = 4, bg = "transparent")
