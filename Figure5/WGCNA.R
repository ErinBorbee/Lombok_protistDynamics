#load packages
library(WGCNA)

#set working directory
setwd("~/Desktop/DataAnalysis/threshold95/copepodData/")
#import data (rows = taxa; columns = samples, column 1 = taxa names)
options(stringsAsFactors = FALSE)

table <- read.csv("SAR-copepod-OTU-percent-table-bySite.csv", row.names = 1)
table <- as.data.frame(t(table))

meta <- read.csv("~/Desktop/DataAnalysis/Master_V9_MappingFile.csv")
LBK_meta <- subset(meta, region == "Lombok")
LBK_meta <- subset(LBK_meta, depth == "M")
LBK_M04_meta <- subset(LBK_meta, filter == "0.4um")
LBK_M12_meta <- subset(LBK_meta, filter == "12um")

LBK_M04_IDs <- LBK_M04_meta$sampleID
LBK_M12_IDs <- LBK_M12_meta$sampleID

LBK_M04_table <- as.data.frame(subset(table, row.names(table) %in% LBK_M04_IDs))
LBK_M12_table <- as.data.frame(subset(table, row.names(table) %in% LBK_M12_IDs))

NW_LBK_M04_meta <- subset(LBK_M04_meta, location == "Northwest Lombok")
NW_LBK_M12_meta <- subset(LBK_M12_meta, location == "Northwest Lombok")
NW_LBK_M04_IDs <- NW_LBK_M04_meta$siteID
NW_LBK_M12_IDs <- NW_LBK_M12_meta$siteID

NW_LBK_M04_table <- subset(LBK_M04_table, row.names(LBK_M04_table) %in% NW_LBK_M04_IDs)
NW_LBK_M12_table <- subset(LBK_M12_table, row.names(LBK_M12_table) %in% NW_LBK_M12_IDs)

NE_LBK_M04_meta <- subset(LBK_M04_meta, location == "Northeast Lombok")
NE_LBK_M12_meta <- subset(LBK_M12_meta, location == "Northeast Lombok")
NE_LBK_M04_IDs <- NE_LBK_M04_meta$sampleID
NE_LBK_M12_IDs <- NE_LBK_M12_meta$sampleID

NE_LBK_M04_table <- subset(LBK_M04_table, row.names(LBK_M04_table) %in% NE_LBK_M04_IDs)
NE_LBK_M12_table <- subset(LBK_M12_table, row.names(LBK_M12_table) %in% NE_LBK_M12_IDs)

SW_LBK_M04_meta <- subset(LBK_M04_meta, location == "South Lombok")
SW_LBK_M12_meta <- subset(LBK_M12_meta, location == "South Lombok")
SW_LBK_M04_IDs <- SW_LBK_M04_meta$sampleID
SW_LBK_M12_IDs <- SW_LBK_M12_meta$sampleID

SW_LBK_M04_table <- subset(LBK_M04_table, row.names(LBK_M04_table) %in% SW_LBK_M04_IDs)
SW_LBK_M12_table <- subset(LBK_M12_table, row.names(LBK_M12_table) %in% SW_LBK_M12_IDs)

LBK_IDs <- LBK_meta$siteID

LBK_table <- subset(table, row.names(table) %in% LBK_IDs)

NW_LBK_meta <- subset(LBK_meta, location == "Northwest Lombok")
NW_LBK_IDs <- NW_LBK_meta$siteID

NW_LBK_table <- subset(table, row.names(table) %in% NW_LBK_IDs)

NE_LBK_meta <- subset(LBK_meta, location == "Northeast Lombok")
NE_LBK_IDs <- NE_LBK_meta$siteID

NE_LBK_table <- subset(table, row.names(table) %in% NE_LBK_IDs)

SW_LBK_meta <- subset(LBK_meta, location == "South Lombok")
SW_LBK_IDs <- SW_LBK_meta$siteID

SW_LBK_table <- subset(table, row.names(table) %in% SW_LBK_IDs)

#Look at data
dim(NW_LBK_table)
names(NW_LBK_table)

#Check for samples with too many missing values
missing_taxa <- goodSamplesGenes(LBK_table, verbose = 3)
missing_taxa$allOK

if (!missing_taxa$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!missing_taxa$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(LBK_table)[!missing_taxa$goodGenes], collapse = ", ")));
  if (sum(!missing_taxa$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(LBK_table)[!missing_taxa$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  table = LBK_table[missing_taxa$goodSamples, missing_taxa$goodGenes]
}

#Construct dendrogram
sampleTree = hclust(dist(table), method = "average");

# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);

par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 15, col = "red");

# Determine cluster under the line
clust <- cutreeStatic(sampleTree, cutHeight = 1, minSize = 0)
table(clust)

# clust 1 contains the samples we want to keep.
keepSamples <- (clust==1)
table <- table[keepSamples, ]

nGenes <- ncol(table)
nSamples <- nrow(table)

metadata <- read.csv("~/Desktop/DataAnalysis/Master_V9_MappingFile.csv")
dim(LBK_meta)
names(LBK_meta)

# remove columns that hold information we do not need (remove categorical data columns)
filteredMeta <- LBK_meta[-c(2,4,6,12,14,16,18,20,22,24,26,
                            28,30,32,33), 
                         -c(1:3,5:11,13:14,16,30:31,33:34,43,47)];
dim(filteredMeta)
names(filteredMeta)

# Form a data frame analogous to expression data that will hold the clinical traits.
samples <- rownames(table);
metadataRows <- match(samples, filteredMeta$siteID);
metadata_var <- filteredMeta[metadataRows, -1];
rownames(metadata_var) <- filteredMeta[metadataRows, 1];
collectGarbage();

# Re-cluster samples
sampleTree2 <- hclust(dist(table), method = "average")

# Convert traits to a color representation: white means low, red means high, grey means missing entry
metaColors <- numbers2colors(metadata_var, signed = FALSE);

# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, metaColors,
                    groupLabels = names(metadata_var),
                    main = "Sample dendrogram and trait heatmap")

# Choose a set of soft-thresholding powers
powers <- c(c(1:10), seq(from = 12, to=20, by=2))


# Call the network topology analysis function
sft <- pickSoftThreshold(table, powerVector = powers, verbose = 5)

# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

# this line corresponds to using an R^2 cut-off of h
abline(h=0.15,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#One-step network construction and module detection
net <- blockwiseModules(table, power = 16,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "sampleTOM",
                       verbose = 3)

#Identify how many modules were identified and how large they are
table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)

# Convert labels to colors for plotting
mergedColors <- labels2colors(net$colors)

# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs;
geneTree <- net$dendrograms[[1]];

# Define numbers of genes and samples
nGenes <- ncol(table);
nSamples <- nrow(table);

# Recalculate MEs with color labels
MEs0 <- moduleEigengenes(table, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, metadata_var, use = "p");
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples);

sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) <- dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
pdf('LBK-OTU-heatmap.pdf')
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(metadata_var),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

# Define variable dist containing the dist column of metadata
npp_mean <- as.data.frame(metadata_var$npp_mean);
land_15km <- as.data.frame(metadata_var$land_area_15km);
reef_15km <- as.data.frame(metadata_var$reef_area_15km);
wave_mean <- as.data.frame(metadata_var$wave_mean);
pop_25km <- as.data.frame(metadata_var$pop2020_25km);
dist <- as.data.frame(metadata_var$dist_market);
rubble <- as.data.frame(metadata_var$rubble_percent);
hardCoral <- as.data.frame(metadata_var$hardCoral_percent);
softCoral <- as.data.frame(metadata_var$softCoral_percent);
abundance <- as.data.frame(metadata_var$abundance);
target_biomass <- as.data.frame(metadata_var$target_biomass);
nontarget_biomass <- as.data.frame(metadata_var$nontarget_biomass);
total_biomass <- as.data.frame(metadata_var$total_biomass);

names(npp_mean) <- "npp_mean"
names(land_15km) <- "land_15km"
names(reef_15km) <- "reef_15km"
names(wave_mean) <- "wave_mean"
names(pop_25km) <- "pop_25km"
names(dist) <- "dist"
names(rubble) <- "rubble"
names(hardCoral) <- "hardCoral"
names(softCoral) <- "softCoral"
names(abundance) <- "abundance"
names(target_biomass) <- "target_biomass"
names(nontarget_biomass) <- "nontarget_biomass"
names(total_biomass) <- "total_biomass"

merged <- cbind(npp_mean,land_15km,reef_15km,wave_mean,pop_25km,
                dist,rubble, hardCoral, softCoral, abundance,
                target_biomass, nontarget_biomass, total_biomass)

# names (colors) of the modules (MM = module membership, GS = gene significance)
modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(table, MEs, use = "p"));
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) <- paste("MM", modNames, sep="");
names(MMPvalue) <- paste("p.MM", modNames, sep="");

geneTraitSignificance <- as.data.frame(cor(table, merged, use = "p"));
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) <- paste("GS.", names(merged), sep="");
names(GSPvalue) <- paste("p.GS.", names(merged), sep="");

#Identifying taxa with high module membership and gene significance
module = "grey60"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Taxa significance for dist",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#Summary of network analysis results
names(table)
names(table)[moduleColors=="grey"]

#Assigning taxonomy to ASV identifiers
taxonomy <- read.csv("~/Desktop/DataAnalysis/WGCNA/pr2-taxonomy.csv")
dim(taxonomy)
names(taxonomy)
probes <- names(table)
probes2taxa <- match(probes, taxonomy$ID)

# The following is the number or probes without annotation (should return 0)
sum(is.na(probes2taxa))

# Create the starting data frame
asvInfo0 = data.frame(ID = probes,
                       geneSymbol = taxonomy$name[probes2taxa],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)

# Order modules by their significance for dist
modOrder = order(-abs(cor(MEs, dist, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(asvInfo0)
  asvInfo0 = data.frame(asvInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(asvInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder <- order(asvInfo0$moduleColor, -abs(asvInfo0$GS.dist));
asvInfo <- asvInfo0[geneOrder, ]

write.csv(asvInfo, file = "LBK-otuInfo.csv")

# Get the corresponding Locus Link IDs
allIDs = taxonomy$ID[probes2taxa];

# Choose interesting modules
intModules = c("turquoise")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modIDs = allIDs[modGenes];
  # Write them into a file
  fileName = paste("ID-", module, ".txt", sep="");
  write.table(as.data.frame(modIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}

# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(table, power = 16);

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^16;

# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;

# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all taxa")

# Recalculate module eigengenes
MEs2 <- moduleEigengenes(table, moduleColors)$eigengenes

# Isolate weight from the clinical traits
npp <- as.data.frame(metadata_var$npp_mean);
names(npp) <- "npp"

# Add the weight to existing module eigengenes
MET <- orderMEs(cbind(MEs, npp))

# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle = 90)

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)

# Recalculate topological overlap
TOM <- TOMsimilarityFromExpr(table, power = 16);

# Select module
modules <- c("turquoise","blue","brown",
             "yellow","green","red","black");

# Select module probes
probes <- names(table)
inModule <- (moduleColors==modules);
modProbes <- probes[inModule];
modGenes = taxonomy$name[match(modProbes, taxonomy$ID)];

# Select the corresponding Topological Overlap
modTOM <- TOM[inModule, inModule];
dimnames(modTOM) <- list(modProbes, modProbes)

# Export the network into edge and node list files Cytoscape can read
cyt <- exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-LBK-OTUs.txt"),
                               nodeFile = paste("CytoscapeInput-nodes-LBK-OTUs.txt"),
                               weighted = TRUE,
                               threshold = 0.02,
                               nodeNames = modProbes,
                               altNodeNames = modGenes,
                               nodeAttr = moduleColors[inModule]);

