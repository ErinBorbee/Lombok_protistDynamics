#Load packages
library(phyloseq)
library(ape)
library(ggplot2)
library(vegan)
library(ggpubr)

#set working directory
setwd("~/Desktop/Dissertation/Chapter3/")

#Import data
otu_table <- read.csv("~/Desktop/DataAnalysis/threshold95/SAR/SAR-percent-table-95.csv", row.names = 1)
otu_table <- otu_table[,-c(1:8)]
sample_data <- read.csv("~/Desktop/DataAnalysis/Master_V9_MappingFile.csv", row.names = 1)
tax_table <- read.csv("~/Desktop/DataAnalysis/pr2-taxonomy.csv", row.names = 1)
taxmat <- as.matrix(tax_table)
tree <- read.tree("~/Desktop/DataAnalysis/tree.nwk")

#Convert data into phyloseq format
otu <- otu_table(otu_table, taxa_are_rows = TRUE)
meta <- sample_data(sample_data)
tax <- tax_table(taxmat)

#Generate phyloseq object
physeq <- phyloseq(otu,meta,tax,tree)

#Subset phyloseq object
physeqM <- subset_samples(physeq, depth == "M")
physeqS <- subset_samples(physeq, depth == "S")

physeq04 <- subset_samples(physeq, filter == "0.4um")
physeq12 <- subset_samples(physeq, filter == "12um")

physeqM04 <- subset_samples(physeqM, filter == "0.4um")
physeqM12 <- subset_samples(physeqM, filter == "12um")
physeqS04 <- subset_samples(physeqS, filter =="0.4um")
physeqS12 <- subset_samples(physeqS, filter == "12um")

physeqM04_LBK <- subset_samples(physeqM04, region == "Lombok")
physeqM12_LBK <- subset_samples(physeqM12, region == "Lombok")
physeqS04_LBK <- subset_samples(physeqS04, region == "Lombok")
physeqS12_LBK <- subset_samples(physeqS12, region == "Lombok")

#Run CAP ordination and build CAP plot
CAP_M04_LBK <- ordinate(physeqM04_LBK, "CAP","bray", formula = ~herbAbund + coraliAbund + 
                          copepod_percent_v9 + sponge_percent + anthozoa_percent +
                          rubble_percent + npp_mean, na.action = na.exclude)
anova(CAP_M04_LBK)
anova.cca(CAP_M04_LBK, by = "axis", step = 1000)
summary(CAP_M04_LBK)

CAP_M04_LBK_plot <- plot_ordination(physeqM04_LBK, CAP_M04_LBK, axes = c(1,2)) +
  geom_point(color = "black", size = 4, alpha = 0.7, aes(fill = npp_mean, shape = location)) + scale_fill_viridis_c() +
  scale_shape_manual(values = c('Northeast Lombok' = 21, 'Northwest Lombok' = 22, 'South Lombok' = 24))

#Add environmental variables as arrows to make a biplot
arrowmat <- vegan::scores(CAP_M04_LBK, display = "bp")
#Add labels to matrix to create data frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
#Define arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.5 * CAP1, 
                 y = 1.5 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

#Create biplot grpahic
CAP_M04_LBK_biplot <- CAP_M04_LBK_plot + geom_segment(mapping = arrow_map, size = .5, data = arrowdf, color = "black", arrow = arrowhead) + 
  geom_text(mapping = label_map, size = 4, data = arrowdf, show.legend = FALSE, color = "black") + 
  theme_bw() + theme(panel.grid = element_blank()) + ylim(-1.5,2) + xlim(-1.5,1.5)
CAP_M04_LBK_biplot

ggsave("~/Desktop/CAP_M04_LBK_biplot.pdf", CAP_M04_LBK_biplot, width = 7, height = 6)

#Extract coordinates for points on CAP axes
CAP_M04_LBK_scores <- scores(CAP_M04_LBK, display = "sites")
CAP_M04_LBK_scores <- merge(CAP_M04_LBK_scores, sample_data, by = "row.names")

#Plot CAP1 vs. rubble coverage data
M04_plot <- ggplot(CAP_M04_LBK_scores, aes(CAP1, rubble_percent)) +
  geom_smooth(method = "lm", color = "grey25", linetype = 5, alpha = 0.15) + 
  geom_point(color = "black", alpha = 0.7, size = 4, aes(fill = npp_mean, shape = location)) + 
  theme_bw() + theme(panel.grid = element_blank()) + 
  scale_fill_viridis_c() + scale_shape_manual(values = c('Northeast Lombok' = 21, 'Northwest Lombok' = 22, 'South Lombok' = 24)) +
  stat_regline_equation() + stat_cor(label.y = 72) + ylab("mean NPP")
M04_plot
ggsave("~/Desktop/M04_CAP_rubbleNPP_plot_LBK.pdf", M04_plot, width = 7, height = 6)

#Run CAP ordination and build CAP plot
CAP_M12_LBK <- ordinate(physeqM12_LBK, "CAP","bray", formula = ~invertAbund + coraliAbund + 
                          anthozoa_percent +
                          npp_mean + rubble_percent, na.action = na.exclude)
anova(CAP_M12_LBK)
anova.cca(CAP_M12_LBK, by = "axis", step = 1000)
summary(CAP_M12_LBK)

CAP_M12_LBK_plot <- plot_ordination(physeqM12_LBK, CAP_M12_LBK, axes = c(1,2)) +
  geom_point(color = "black", size = 4, alpha = 0.7, aes(fill = npp_mean, shape = location)) + scale_fill_viridis_c() +
  scale_shape_manual(values = c('Northeast Lombok' = 21, 'Northwest Lombok' = 22, 'South Lombok' = 24))

#Add environmental variables as arrows to make a biplot
arrowmat <- vegan::scores(CAP_M12_LBK, display = "bp")
#Add labels to matrix to create data frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
#Define arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.5 * CAP1, 
                 y = 1.5 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

#Create biplot grpahic
CAP_M12_LBK_biplot <- CAP_M12_LBK_plot + geom_segment(mapping = arrow_map, size = .5, data = arrowdf, color = "black", arrow = arrowhead) + 
  geom_text(mapping = label_map, size = 4, data = arrowdf, show.legend = FALSE, color = "black") + 
  theme_bw() + theme(panel.grid = element_blank()) + ylim(-2,1.5) + xlim(-1.5, 1.5)
CAP_M12_LBK_biplot
ggsave("CAP_M12_LBK_biplot.pdf", CAP_M12_LBK_biplot, width = 7, height = 6)

#Extract coordinates for points on CAP axes
CAP_M12_LBK_scores <- scores(CAP_M12_LBK, display = "sites")
CAP_M12_LBK_scores <- merge(CAP_M12_LBK_scores, sample_data, by = "row.names")

#Plot CAP1 vs. copepod data
M12_plot <- ggplot(CAP_M12_LBK_scores, aes(CAP1, rubble_percent)) +
  geom_smooth(method = "lm", color = "grey25", linetype = 5, alpha = 0.15) + geom_point(color = "black", alpha = 0.7, size = 4, aes(fill = npp_mean, shape = location)) + 
  theme_bw() + theme(panel.grid = element_blank()) + 
  scale_fill_viridis_c() + scale_shape_manual(values = c('Northeast Lombok' = 21, 'Northwest Lombok' = 22, 'South Lombok' = 24)) +
  stat_regline_equation() + stat_cor(label.y = 72) + ylab("rubble percent coverage")
M12_plot
ggsave("M12_CAP_rubble_plot_LBK.pdf", M12_plot, width = 8, height = 7)

#Run CAP ordination and build CAP plot
CAP_S04_LBK <- ordinate(physeqS04_LBK, "CAP","bray", formula = ~herbAbund + coraliAbund + anthozoa_percent +
                          copepod_percent_v9 + rubble_percent + npp_mean, na.action = na.exclude)
anova(CAP_S04_LBK)
anova.cca(CAP_S04_LBK, by = "axis", step = 1000)
summary(CAP_S04_LBK)

CAP_S04_LBK_plot <- plot_ordination(physeqS04_LBK, CAP_S04_LBK, axes = c(1,2)) +
  geom_point(color = "black", size = 4, alpha = 0.7, aes(fill = npp_mean, shape = location)) + scale_fill_viridis_c() +
  scale_shape_manual(values = c('Northeast Lombok' = 21, 'Northwest Lombok' = 22, 'South Lombok' = 24))

#Add environmental variables as arrows to make a biplot
arrowmat <- vegan::scores(CAP_S04_LBK, display = "bp")
#Add labels to matrix to create data frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
#Define arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.5 * CAP1, 
                 y = 1.5 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

#Create biplot grpahic
CAP_S04_LBK_biplot <- CAP_S04_LBK_plot + geom_segment(mapping = arrow_map, size = .5, data = arrowdf, color = "black", arrow = arrowhead) + 
  geom_text(mapping = label_map, size = 4, data = arrowdf, show.legend = FALSE, color = "black") + 
  theme_bw() + theme(panel.grid = element_blank()) + ylim(-3,1.5) + xlim(-1.5, 1.5)
CAP_S04_LBK_biplot

ggsave("CAP_S04_LBK_biplot.pdf", CAP_S04_LBK_biplot, width = 7, height = 6)

#Extract coordinates for points on CAP axes
CAP_S04_LBK_scores <- scores(CAP_S04_LBK, display = "sites")
CAP_S04_LBK_scores <- merge(CAP_S04_LBK_scores, sample_data, by = "row.names")

#Plot CAP1 vs. copepod data
S04_plot <- ggplot(CAP_S04_LBK_scores, aes(CAP1, rubble_percent)) +
  geom_smooth(method = "lm", color = "grey25", linetype = 5, alpha = 0.15) + geom_point(color = "black", alpha = 0.7, size = 4, aes(fill = npp_mean, shape = location)) + 
  theme_bw() + theme(panel.grid = element_blank()) + 
  scale_fill_viridis_c() + scale_shape_manual(values = c('Northeast Lombok' = 21, 'Northwest Lombok' = 22, 'South Lombok' = 24)) +
  stat_regline_equation() + stat_cor(label.y = 72) + ylab("rubble percent coverage")
S04_plot
ggsave("S04_CAP_rubble_plot_LBK.pdf", S04_plot, width = 8, height = 7)

#Run CAP ordination and build CAP plot
CAP_S12_LBK <- ordinate(physeqS12_LBK, "CAP","bray", formula = ~detritAbund + 
                          copepod_percent_v9 + anthozoa_percent +
                          sponge_percent + rubble_percent + npp_mean, na.action = na.exclude)

anova(CAP_S12_LBK)
anova.cca(CAP_S12_LBK, by = "axis", step = 1000)
summary(CAP_S12_LBK)

CAP_S12_LBK_plot <- plot_ordination(physeqS12_LBK, CAP_S12_LBK, axes = c(1,2)) +
  geom_point(color = "black", size = 4, alpha = 0.7, aes(fill = npp_mean, shape = location)) + scale_fill_viridis_c() +
  scale_shape_manual(values = c('Northeast Lombok' = 21, 'Northwest Lombok' = 22, 'South Lombok' = 24))

#Add environmental variables as arrows to make a biplot
arrowmat <- vegan::scores(CAP_S12_LBK, display = "bp")
#Add labels to matrix to create data frame
arrowdf <- data.frame(labels = rownames(arrowmat), arrowmat)
#Define arrow aesthetic mapping
arrow_map <- aes(xend = CAP1, 
                 yend = CAP2, 
                 x = 0, 
                 y = 0, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

label_map <- aes(x = 1.5 * CAP1, 
                 y = 1.5 * CAP2, 
                 shape = NULL, 
                 color = NULL, 
                 label = labels)

arrowhead = arrow(length = unit(0.02, "npc"))

#Create biplot grpahic
CAP_S12_LBK_biplot <- CAP_S12_LBK_plot + geom_segment(mapping = arrow_map, size = .5, data = arrowdf, color = "black", arrow = arrowhead) + 
  geom_text(mapping = label_map, size = 4, data = arrowdf, show.legend = FALSE, color = "black") + 
  theme_bw() + theme(panel.grid = element_blank())
CAP_S12_LBK_biplot
ggsave("CAP_S12_LBK_biplot.pdf", CAP_S12_LBK_biplot, width = 7, height = 6)

#Extract coordinates for points on CAP axes
CAP_S12_LBK_scores <- scores(CAP_S12_LBK, display = "sites")
CAP_S12_LBK_scores <- merge(CAP_S12_LBK_scores, sample_data, by = "row.names")

#Plot CAP1 vs. copepod data
S12_plot <- ggplot(CAP_S12_LBK_scores, aes(CAP1, rubble_percent)) +
  geom_smooth(method = "lm", color = "grey25", linetype = 5, alpha = 0.15) + geom_point(color = "black", alpha = 0.7, size = 4, aes(fill = npp_mean, shape = location)) + 
  theme_bw() + theme(panel.grid = element_blank()) + 
  scale_fill_viridis_c() + scale_shape_manual(values = c('Northeast Lombok' = 21, 'Northwest Lombok' = 22, 'South Lombok' = 24)) +
  stat_regline_equation() + stat_cor(label.y = 72) + ylab("rubble percent coverage")
S12_plot
ggsave("S12_CAP_copepod_plot_LBK.pdf", S12_plot, width = 8, height = 7)

#Combine plots into one 2x4 grid and save
plot <- ggarrange(CAP_M04_LBK_biplot, M04_plot,
                  CAP_M12_LBK_biplot, M12_plot,
                  CAP_S04_LBK_biplot, S04_plot,
                  CAP_S12_LBK_biplot, S12_plot,
                  ncol = 2, nrow = 4, legend = "bottom",
                  align = "hv", common.legend = TRUE)
plot
ggsave("CAP_rubble_plots.pdf", plot, width = 8, height = 12)
