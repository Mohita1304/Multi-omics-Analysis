setwd("/path/")

library(DESeq2)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(CorLevelPlot)
library(WGCNA)

d<-read.csv("CGs_RawData_1.csv", header = TRUE, sep=",")
#View(d)
duplicate_genes=duplicated(d$Gene)
d= d[!duplicate_genes, ]
rownames(d)<-d$Gene
d1<-d[, -1]
#View(d1)
samples.to.be.excluded <- c('TCGA.AZ.6599.01A', 'TCGA-A6-2680-11A')
d1.subset <- d1[,!(colnames(d1) %in% samples.to.be.excluded)]

colData <- read.csv('Sample_Annot_1.csv', row.names = 1) #
all(colnames(d1) %in% rownames(colData)) # to check if column name are matching with row name 
all(colnames(d1)==rownames(colData)) # in the same order

dds <- DESeqDataSetFromMatrix(countData = d1,
                              colData = colData,
                              design = ~ 1) # not spcifying model

dds75 <- dds[rowSums(counts(dds) >= 15) >= 24,]
nrow(dds75) 

dds_norm <- vst(dds75)

norm.counts <- assay(dds_norm) %>% 
  t()

norm.counts[] <- sapply(norm.counts, as.numeric)

power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)


sft.data <- sft$fitIndices
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Soft Threshold Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),  
        axis.text = element_text(size = 12),   
        plot.title = element_text(size = 16)) 

a1

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic() +
  theme(axis.title = element_text(size = 14),  
        axis.text = element_text(size = 12),   
        plot.title = element_text(size = 16)) 

install.packages("gridExtra")
library(gridExtra)
grid.arrange(a1, a2, nrow = 2)
ggsave("plot_output.png", arrangeGrob(a1, a2, nrow = 2))


soft_power <- 18  # selected the soft threshold value 
temp_cor <- cor
cor <- WGCNA::cor

# memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                          maxBlockSize = 14000,
                          TOMType = "signed",
                          corType = "pearson",
                          power = soft_power,
                          mergeCutHeight = 0.25,
                          numericLabels = FALSE,
                          randomSeed = 1234,
                          verbose = 3)


cor <- temp_cor
module_eigengenes <- bwnet$MEs
head(module_eigengenes)
table(bwnet$colors)

plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

traits <- colData %>%
  mutate(Annotation = ifelse(grepl('Normal', Annotation), 1, 0)) %>%
  select(Annotation)

colData$Satge_Types <- factor(colData$Satge_Types, levels = c("Normal","EarlyStage","LateStage"))

severity.out <- binarizeCategoricalColumns(colData$Satge_Types,
                                           includePairwise = FALSE,
                                           includeLevelVsAll = TRUE,
                                           minCount = 1)

traits <- cbind(traits, severity.out)

nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)


module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

# visualize module-trait association as a heatmap

heatmap.data <- merge(module_eigengenes, traits, by = 'row.names')

head(heatmap.data)

names(heatmap.data)

heatmap.data <- heatmap.data %>% 
  column_to_rownames(var = 'Row.names')

library(CorLevelPlot)

CorLevelPlot(heatmap.data,
             x = names(heatmap.data)[8:14],
             y = names(heatmap.data)[1:7],
             col = c("blue1", "skyblue", "white", "pink", "red"))

module.gene.mapping <- as.data.frame(bwnet$colors)


module.membership.measure <- cor(module_eigengenes$MEgreen, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)

gene.signf.corr <- cor(norm.counts, traits, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

g<-gene.signf.corr.pvals %>% 
  as.data.frame()

#names (colors) of the modules
modNames = substring(names(module_eigengenes), 3)

geneModuleMembership = as.data.frame(cor(norm.counts, module_eigengenes$MEturquoise, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

geneTraitSignificance = as.data.frame(cor(norm.counts, traits$data.EarlyStage.vs.all, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

new_col_name <- paste0(colnames(GSPvalue), "_Pvalue")

# Rename the column in MMPvalue dataframe
colnames(GSPvalue) <- new_col_name

# Merge the two data frames by row index
merged_df <- cbind(geneTraitSignificance, GSPvalue)

#write.csv(merged_df, file="GeneSignificance_Pvalue.csv")

module = "blue"#########################putting the color below the plot
column = match(module, modNames)
moduleGenes = bwnet$colors ==module

#sizeGrWindow(7, 7)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))  # Adjust the margin to leave space for axis labels

verboseScatterplot(abs(geneModuleMembership[moduleGenes, "V1"]),
                   abs(geneTraitSignificance[moduleGenes, "V1"]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Normal state",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
abline(h = 0.4, col = "red")
abline(v = 0.8, col = "red")
text(abs(geneModuleMembership[moduleGenes, "V1"]),
     abs(geneTraitSignificance[moduleGenes, "V1"]),
     labels = rownames(geneModuleMembership)[moduleGenes],
     pos = 1, cex = 0.8)


########## Merging these two file in one 
new_col_name <- paste0(colnames(MMPvalue), "_Pvalue")

# Rename the column in MMPvalue dataframe
colnames(MMPvalue) <- new_col_name

# Merge the two data frames by row index
merged_df <- cbind(geneModuleMembership, MMPvalue)

write.csv(merged_df, file="/home/sukanta-mondal/Documents/cBioPortal/Multi_Omics_Analysis/WGCNA/WGCNA/Final_WGCNA/ModulesMembership_Pvalue.csv")



# TO get the hub genes based on the Gene Signifcacne and modulemembership

# Define the thresholds
mm_threshold <- 0.8
gs_threshold <- 0.5

# Extract the necessary values
mm_values <- abs(geneModuleMembership[moduleGenes, "V1"])
gs_values <- abs(geneTraitSignificance[moduleGenes, "V1"])
gene_names <- rownames(geneModuleMembership)[moduleGenes]

# Filter for genes that meet the thresholds
filter_indices <- which(mm_values > mm_threshold & gs_values > gs_threshold)
filtered_mm_values <- mm_values[filter_indices]
filtered_gs_values <- gs_values[filter_indices]
filtered_gene_names <- gene_names[filter_indices]

# Print the names of the filtered genes
print(filtered_gene_names)

# Plot the scatterplot with all points
verboseScatterplot(mm_values, gs_values,
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Normal state",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
abline(h = gs_threshold, col = "red")
abline(v = mm_threshold, col = "red")

# Add text labels for the filtered genes
text(filtered_mm_values, filtered_gs_values,
     labels = filtered_gene_names,
     pos = 1, cex = 0.8)

