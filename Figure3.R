#load packages
library(vegan)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(grid)

setwd("~/Desktop/Dissertation/Chapter3/")

#import data and subset by taxonomic group
data <- read.csv("~/Desktop/DataAnalysis/threshold95/richness/LBK_richness_by_site.csv")

api_data <- subset(data, taxa == "Apicomplexa")
syn_data <- subset(data, taxa == "Syndiniales")
dino_data <- subset(data, taxa == "Dinophyceae")
diatom_data <- subset(data, taxa == "Bacillariophyta")
pelago_data <- subset(data, taxa == "Pelagophyceae")
MAST_data <- subset(data, taxa == "MAST")
cerc_data <- subset(data, taxa == "Cercozoa")
foram_data <- subset(data, taxa == "Foraminifera")
radio_data <- subset(data, taxa == "Radiolaria")
ciliate_data <- subset(data, taxa == "Ciliophora")

#Since not all variables were normal after transformation we used nonparametric stats for consistency throughout

#test for normality
hist(api_data$ASV_richness, main = "Richness", breaks = 10)
shapiro.test(api_data$ASV_richness)
#Normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(ASV_richness ~ location, data = api_data)
#p-value is significant (p = 0.01053)
pairwise.wilcox.test(api_data$ASV_richness, api_data$location, p.adjust.method = "BH")
#SW significantly different from NE and NW (p= 0.015 and 0.046 respectively)

#test for normality
hist(syn_data$ASV_richness, main = "Richness", breaks = 10)
shapiro.test(syn_data$ASV_richness)
#Normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(ASV_richness ~ location, data = syn_data)
#p-value not significant (p = 0.1477)
pairwise.wilcox.test(syn_data$ASV_richness, syn_data$location, p.adjust.method = "BH")
#regions not significantly different

#test for normality
hist(dino_data$ASV_richness, main = "Richness", breaks = 10)
shapiro.test(dino_data$ASV_richness)
#Not normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(ASV_richness ~ location, data = dino_data)
#p-value is significant (p = 0.04534)
pairwise.wilcox.test(dino_data$ASV_richness, dino_data$location, p.adjust.method = "BH")
#SW significantly different from NE (p= 0.026)

#test for normality
hist(diatom_data$ASV_richness, main = "Richness", breaks = 10)
shapiro.test(diatom_data$ASV_richness)
#Not normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(ASV_richness ~ location, data = diatom_data)
#p-value not significant (p = 0.2095)
pairwise.wilcox.test(diatom_data$ASV_richness, diatom_data$location, p.adjust.method = "BH")
#loactions not significantly different

#test for normality
hist(pelago_data$ASV_richness, main = "Richness", breaks = 10)
shapiro.test(pelago_data$ASV_richness)
#Not normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(ASV_richness ~ location, data = pelago_data)
#p-value not significant (p = 0.06428)
pairwise.wilcox.test(pelago_data$ASV_richness, pelago_data$location, p.adjust.method = "BH")
#SW significantly different from NE (p= 0.02)

#test for normality
hist(MAST_data$ASV_richness, main = "Richness", breaks = 10)
shapiro.test(MAST_data$ASV_richness)
#Normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(ASV_richness ~ location, data = MAST_data)
#p-value not significant (p = 0.3016)
pairwise.wilcox.test(MAST_data$ASV_richness, MAST_data$location, p.adjust.method = "BH")
#locations not significantly different

#test for normality
hist(cerc_data$ASV_richness, main = "Richness", breaks = 10)
shapiro.test(cerc_data$ASV_richness)
#Normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(ASV_richness ~ location, data = cerc_data)
#p-value not significant (p = 0.6072)
pairwise.wilcox.test(cerc_data$ASV_richness, cerc_data$location, p.adjust.method = "BH")
#locations not significantly different

#test for normality
hist(foram_data$ASV_richness, main = "Richness", breaks = 10)
shapiro.test(foram_data$ASV_richness)
#Normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(ASV_richness ~ location, data = foram_data)
#p-value not significant (p = 0.2732)
pairwise.wilcox.test(foram_data$ASV_richness, foram_data$location, p.adjust.method = "BH")
#locations not significantly different

#test for normality
hist(radio_data$ASV_richness, main = "Richness", breaks = 10)
shapiro.test(radio_data$ASV_richness)
#Normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(ASV_richness ~ location, data = radio_data)
#p-value not significant (p = 0.09928)
pairwise.wilcox.test(radio_data$ASV_richness, radio_data$location, p.adjust.method = "BH")
#locations not significantly different

#test for normality
hist(ciliate_data$ASV_richness, main = "Richness", breaks = 10)
shapiro.test(ciliate_data$ASV_richness)
#Normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(ASV_richness ~ location, data = ciliate_data)
#p-value not significant (p = 0.515)
pairwise.wilcox.test(ciliate_data$ASV_richness, ciliate_data$location, p.adjust.method = "BH")
#locations not significantly different

#test for normality
hist(api_data$OTU_richness, main = "Richness", breaks = 10)
shapiro.test(api_data$OTU_richness)
#Normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(OTU_richness ~ location, data = api_data)
#p-value not significant (p = 0.2831)
pairwise.wilcox.test(api_data$OTU_richness, api_data$location, p.adjust.method = "BH")
#locations not significantly different

#test for normality
hist(syn_data$OTU_richness, main = "Richness", breaks = 10)
shapiro.test(syn_data$OTU_richness)
#Normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(OTU_richness ~ location, data = syn_data)
#p-value not significant (p = 0.1913)
pairwise.wilcox.test(syn_data$OTU_richness, syn_data$location, p.adjust.method = "BH")
#regions not significantly different

#test for normality
hist(dino_data$OTU_richness, main = "Richness", breaks = 10)
shapiro.test(dino_data$OTU_richness)
#Not normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(OTU_richness ~ location, data = dino_data)
#p-value is significant (p = 0.1456)
pairwise.wilcox.test(dino_data$OTU_richness, dino_data$location, p.adjust.method = "BH")
#locations not significantly different

#test for normality
hist(diatom_data$OTU_richness, main = "Richness", breaks = 10)
shapiro.test(diatom_data$OTU_richness)
#Normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(OTU_richness ~ location, data = diatom_data)
#p-value not significant (p = 0.1159)
pairwise.wilcox.test(diatom_data$OTU_richness, diatom_data$location, p.adjust.method = "BH")
#loactions not significantly different

#test for normality
hist(pelago_data$OTU_richness, main = "Richness", breaks = 10)
shapiro.test(pelago_data$OTU_richness)
#Not normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(OTU_richness ~ location, data = pelago_data)
#p-value not significant (p = 0.156)
pairwise.wilcox.test(pelago_data$OTU_richness, pelago_data$location, p.adjust.method = "BH")
#locations not significantly different

#test for normality
hist(MAST_data$OTU_richness, main = "Richness", breaks = 10)
shapiro.test(MAST_data$OTU_richness)
#Not normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(OTU_richness ~ location, data = MAST_data)
#p-value not significant (p = 0.2825)
pairwise.wilcox.test(MAST_data$OTU_richness, MAST_data$location, p.adjust.method = "BH")
#locations not significantly different

#test for normality
hist(cerc_data$OTU_richness, main = "Richness", breaks = 10)
shapiro.test(cerc_data$OTU_richness)
#Normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(OTU_richness ~ location, data = cerc_data)
#p-value not significant (p = 0.4047)
pairwise.wilcox.test(cerc_data$OTU_richness, cerc_data$location, p.adjust.method = "BH")
#locations not significantly different

#test for normality
hist(foram_data$OTU_richness, main = "Richness", breaks = 10)
shapiro.test(foram_data$OTU_richness)
#Not normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(OTU_richness ~ location, data = foram_data)
#p-value not significant (p = 0.4815)
pairwise.wilcox.test(foram_data$OTU_richness, foram_data$location, p.adjust.method = "BH")
#locations not significantly different

#test for normality
hist(radio_data$OTU_richness, main = "Richness", breaks = 10)
shapiro.test(radio_data$OTU_richness)
#Normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(OTU_richness ~ location, data = radio_data)
#p-value not significant (p = 0.0583)
pairwise.wilcox.test(radio_data$OTU_richness, radio_data$location, p.adjust.method = "BH")
#locations not significantly different

#test for normality
hist(ciliate_data$OTU_richness, main = "Richness", breaks = 10)
shapiro.test(ciliate_data$OTU_richness)
#Normal without transformation

#For comparison of richness by designated metadata variable we use Kruskal-Wallis
#If Kruskal-Wallis is significant (p < 0.05) then we use Wilcoxon for pairwise comparisons
kruskal.test(OTU_richness ~ location, data = ciliate_data)
#p-value not significant (p = 0.6431)
pairwise.wilcox.test(ciliate_data$OTU_richness, ciliate_data$location, p.adjust.method = "BH")
#locations not significantly different

#Order table rows by taxonomic group
data$taxa <- factor(data$taxa, levels = c("Apicomplexa",
                                          "Syndiniales",
                                          "Dinophyceae",
                                          "Bacillariophyta",
                                          "Pelagophyceae",
                                          "MAST",
                                          "Cercozoa",
                                          "Foraminifera",
                                          "Radiolaria",
                                          "Ciliophora"))
#Define color of boxplots by significance (p<0.05)
sig <- c('high' = "#f73e42",
         'low' = "#007cb5",
         'none' = "white")

#Plot boxplots by taxonomic group and geographic site grouping (i.e. NE, NW, SW)
plot <- ggplot(data, aes(location, OTU_richness, fill = NULL)) +
  geom_boxplot(color = "black", alpha = 0.7) + geom_point(position = "jitter", pch = 21, color = "black", alpha = 0.7) + 
  facet_wrap(~taxa, scales = "free") + theme_bw() + theme(panel.grid = element_blank()) +
  stat_cor(size = 2) + scale_fill_manual(values = sig) + theme(strip.text = element_text(color = "white", face = "bold")) +
  ylab("OTU richness") + theme(axis.title.x = element_blank())
plot

#Define colors for strip backgrounds to make the strips match colors of taxonomic groups in Figure 2
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

#Plot new strips colored by taxonomic group and save
dummy <- ggplot(data, aes(location, OTU_richness))+ facet_wrap(~taxa) + 
  geom_rect(aes(fill=taxa), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf, color = "black") +
  theme_minimal() + scale_fill_manual(values = colors) + theme(strip.text = element_text(color = "white", face = "bold"))
dummy

g1 <- ggplotGrob(plot)
g2 <- ggplotGrob(dummy)

gtable_select <- function (x, ...) 
{
  matches <- c(...)
  x$layout <- x$layout[matches, , drop = FALSE]
  x$grobs <- x$grobs[matches]
  x
}

panels <- grepl(pattern="panel", g2$layout$name)
stript <- grepl(pattern="strip-t", g2$layout$name)
strips <- grepl(pattern="strip_t", g2$layout$name)
g2$layout$t[panels] <- g2$layout$t[panels] - 1
g2$layout$b[panels] <- g2$layout$b[panels] - 1

new_strips <- gtable_select(g2, panels | strips | stript)
grid.newpage()
grid.draw(new_strips)

gtable_stack <- function(g1, g2){
  g1$grobs <- c(g1$grobs, g2$grobs)
  g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
  g1$layout <- rbind(g1$layout, g2$layout)
  g1
}

## ideally you'd remove the old strips, for now they're just covered
new_plot <- gtable_stack(g1, new_strips)
grid.newpage()
grid.draw(new_plot)

ggsave("OTU_richness_byTaxa.png", new_plot, width = 10, height = 7)
