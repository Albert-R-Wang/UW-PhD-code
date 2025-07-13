# ------------------------------------------------------------------------------
# Title:       4-Quadrant Plot for Visualizing gene expressions Between Two DE results
# Author:      Albert Wang
# Last updated:        2025-07-12
# ------------------------------------------------------------------------------


### Summary
# This script generates a 4-quadrant scatter plot to compare gene expressions from two differential expression (DE) analyses (two comparisons) based on log2 fold change and adjusted p-values.
# It loads two expression datasets, applies user-defined significance thresholds, and categorizes genes into one of seven expression patterns: Both Upregulated, Both Downregulated, Opposite Fold Change, Upregulated in comp1 only, Downregulated in comp1 only, Upregulated in comp2 only, and Downregulated in comp2 only.
# Users can choose to work with differentially expressed genes (DEGs) from each dataset, their union, or a custom gene list. The script also supports optional subsetting by category, and it assigns labels and colors for clear visualization.
# Final outputs include a composite 4-quadrant scatter plot with density maps or volcano plots for each comparison, and an Excel file containing the categorized DEG table.


###=============================================================================
### Required packages

library(ggplot2)
library(ggpubr)
library(EnhancedVolcano)
library(openxlsx)
library(dplyr)

## References
# http://www.sthda.com/english/wiki/scatter-plots-r-base-graphs
# http://www.sthda.com/english/wiki/ggplot2-scatter-plots-quick-start-guide-r-software-and-data-visualization
# http://www.sthda.com/english/wiki/ggplot2-legend-easy-steps-to-change-the-position-and-the-appearance-of-a-graph-legend-in-r-software
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/

### ============================================================================
### 4-Quadrant Scatter Plot

# Sections:
# 1. Load Expression Data and Define Parameters
# 2. Define the gene list to use for downstream plotting
# 3. Assign Gene Expression Categories
# 4. Configure Plot Labels and Color Schemes
# 5. Subset Data for Visualization
# 6. Generate 4-quadrant Scatter Plot with Density Maps or Volcano Plots



### ============================================================================
### 1. Load Expression Data and Define Parameters

# select the first expression file
RNAseqFile1 <- file.choose()
RNAseqData1 <- read.table(RNAseqFile1,header=TRUE ,sep = "\t", quote="\"", stringsAsFactors = FALSE, check.names = FALSE)
#use the normalized count from DESeq2 (or other methods that account for between-sample comparison)

# select the second expression file
RNAseqFile2 <- file.choose()
RNAseqData2 <- read.table(RNAseqFile2,header=TRUE ,sep = "\t", quote="\"", stringsAsFactors = FALSE, check.names = FALSE)

# Define the names of each comparison group (modify as needed)
comp1 <- "LTvsLN"
comp2 <- "BrMvsLT"

# Set DEG threshold
pValue <- 0.05
log2FC <- 0.6

# Find DEGs
sigGene1 <- subset(RNAseqData1, abs(empiricalLFC)>=log2FC & padj<=pValue)
sigGene2 <- subset(RNAseqData2, abs(empiricalLFC)>=log2FC & padj<=pValue)


### ============================================================================
### 2. Define the gene list to use for downstream plotting

## Three options for defining the genelist. Choose one option to run and skip the others.

# Option 1: Use the union of significant genes from both comparisons
sigGene_all <- union(sigGene1$gene, sigGene2$gene) # obtain a compiled list of DEG from both sigGene1 and sigGene2 (no duplicates)

# Option 2: Use significant genes from only one comparison (e.g., sigGene2)
sigGene_all = sigGene2$gene

# Option 3: Load a custom gene list from file
geneFile = file.choose()
geneList = read.table(geneFile, header=F ,sep = "\t", quote="\"", stringsAsFactors = FALSE, check.names = FALSE)
sigGene_all = geneList$V1


## Subset expression data to retain only genes of interest

# Extract only relevant columns: gene name, adjusted p-value (padj), and log2 fold change (empiricalLFC)
Data1_subset <- RNAseqData1[, c("gene", "padj", "empiricalLFC")]
Data2_subset <- RNAseqData2[, c("gene", "padj", "empiricalLFC")]

# Filter to retain only the selected genes
Data1_sigGene <- Data1_subset[Data1_subset$gene %in% sigGene_all, ]
Data2_sigGene <- Data2_subset[Data2_subset$gene %in% sigGene_all, ]


# Merge into one dataframe
DataMerge <- merge(Data1_sigGene, Data2_sigGene, by= "gene")

### ============================================================================
### 3. Assign Gene Expression Categories

## SKIP: Assign expression categories based on fold change thresholds only (legacy method)
# This block classifies genes based on log2FC alone, without considering statistical significance.
# Not recommended for current use — retained for reference.
DataMerge$Category[DataMerge$empiricalLFC.x>log2FC & DataMerge$empiricalLFC.y>log2FC] <- "Both Upregulated"
DataMerge$Category[DataMerge$empiricalLFC.x<(-log2FC) & DataMerge$empiricalLFC.y<(-log2FC)] <- "Both Downregulated"
DataMerge$Category[DataMerge$empiricalLFC.x>log2FC & DataMerge$empiricalLFC.y<(-log2FC)] <- "Opposite effect"
DataMerge$Category[DataMerge$empiricalLFC.x<(-log2FC) & DataMerge$empiricalLFC.y>log2FC] <- "Opposite effect"
DataMerge$Category[DataMerge$empiricalLFC.x>log2FC & abs(DataMerge$empiricalLFC.y)<log2FC] <- paste("Upregulated", comp1, "only", sep=" ")
DataMerge$Category[DataMerge$empiricalLFC.x<(-log2FC) & abs(DataMerge$empiricalLFC.y)<log2FC] <- paste("Downregulated", comp1, "only", sep=" ")
DataMerge$Category[abs(DataMerge$empiricalLFC.x)<log2FC & DataMerge$empiricalLFC.y>log2FC] <- paste("Upregulated", comp2, "only", sep=" ")
DataMerge$Category[abs(DataMerge$empiricalLFC.x)<log2FC & DataMerge$empiricalLFC.y<(-log2FC)] <- paste("Downregulated", comp2, "only",sep=" ")


## Assign expression categories using both log2FC and adjusted p-value
# Classify DEG status for each dataset individually
DataMerge <- DataMerge %>%
  mutate(
    Category.x = case_when(
      empiricalLFC.x >  log2FC & padj.x < pValue ~ "Upregulated",
      empiricalLFC.x < -log2FC & padj.x < pValue ~ "Downregulated",
      TRUE ~ "Not DE"
      # TRUE acts as the default case: if none of the above conditions are met 
      # (e.g., genes that do not meet thresholds or have NA values), assign "Not DE"
    ),
    Category.y = case_when(
      empiricalLFC.y >  log2FC & padj.y < pValue ~ "Upregulated",
      empiricalLFC.y < -log2FC & padj.y < pValue ~ "Downregulated",
      TRUE ~ "Not DE"
    )
  )

# Assign final expression category based on combined DEG status
DataMerge <- DataMerge %>%
  mutate(
    Category = case_when(
      Category.x == "Upregulated"   & Category.y == "Upregulated"   ~ "Both Upregulated",
      Category.x == "Downregulated" & Category.y == "Downregulated" ~ "Both Downregulated",
      Category.x == "Upregulated"   & Category.y == "Downregulated" ~ "Opposite FC",
      Category.x == "Downregulated" & Category.y == "Upregulated"   ~ "Opposite FC",
      Category.x == "Upregulated"   & Category.y == "Not DE"        ~ paste("Upregulated", comp1, "only"),
      Category.x == "Downregulated" & Category.y == "Not DE"        ~ paste("Downregulated", comp1, "only"),
      Category.x == "Not DE"        & Category.y == "Upregulated"   ~ paste("Upregulated", comp2, "only"),
      Category.x == "Not DE"        & Category.y == "Downregulated" ~ paste("Downregulated", comp2, "only"),
      TRUE ~ NA_character_
    )
  )


## Calculate number of genes in each category
geneCount = as.data.frame(table(t(DataMerge$Category)))
geneCount$percentOverall = geneCount$Freq/sum(geneCount$Freq)


### ============================================================================
### 4. Configure Plot Labels and Color Schemes

## Define thresholds for gene labeling on the scatter plot
SelectFC <- 4 # Use SelectFC = 4 for full dataset; for subsets, consider lowering to 3 or lower. Modify if needed.
SelectP <- 0.001 # Use SelectP = 0.001 for full dataset; for subsets, consider 0.01. Modify if needed

# Assign a default padj value of 1 to genes with missing p-values (ensures they are excluded from labeling)
DataMerge$padj.x[is.na(DataMerge$padj.x)] <- 1
DataMerge$padj.y[is.na(DataMerge$padj.y)] <- 1

# Select genes to label based on fold change and significance thresholds
DataMerge$LabelSelect[abs(DataMerge$empiricalLFC.x)>SelectFC & DataMerge$padj.x<SelectP] <- 
  DataMerge$gene[abs(DataMerge$empiricalLFC.x)>SelectFC & DataMerge$padj.x<SelectP] # Labels for first comp
DataMerge$LabelSelect[abs(DataMerge$empiricalLFC.y)>SelectFC & DataMerge$padj.y<SelectP] <- 
  DataMerge$gene[abs(DataMerge$empiricalLFC.y)>SelectFC & DataMerge$padj.y<SelectP] # Labels for second comp

# Sort data by expression category to allow manually define color when plotting
DataMerge <- DataMerge[order(DataMerge$Category),]

# Store the original DataMerge (with all data) prior to subsetting in the later section
DataMerge.orig <- DataMerge


## Assign colors to expression categories for scatter plot
colorCategory <- data.frame(
  Category = unique(DataMerge$Category),
  stringsAsFactors = FALSE
)

colorCategory <- colorCategory %>%
  mutate(
    Color = case_when(
      Category == "Both Upregulated"                         ~ "red3",
      Category == "Both Downregulated"                       ~ "royalblue",
      Category == paste("Downregulated", comp1, "only")      ~ "cyan3",
      Category == paste("Downregulated", comp2, "only")      ~ "grey30",
      Category == "Opposite FC"                              ~ "green4",
      Category == paste("Upregulated", comp1, "only")        ~ "orange2",
      Category == paste("Upregulated", comp2, "only")        ~ "magenta2",
      TRUE                                                   ~ NA_character_  # catch unrecognized categories
    )
  )


# Warn if any categories are unmatched
if (any(is.na(colorCategory$Color))) {
  warning("The following categories were not assigned a color:\n",
          paste(colorCategory$Category[is.na(colorCategory$Color)], collapse = ", "))
}


colorCategory.orig = colorCategory # Store the original colorCategory (with all data) prior to subsetting in the later section


### ============================================================================
### 5. Subset Data for Visualization

## Optionally subset DataMerge to plot a specific group; otherwise, all data will be used for visualization.

# NOTE: Need to run the following 'readline()' by itself and enter group in the console. 'readline()' does not pause if run all the code together.
# See discussion: https://stackoverflow.com/questions/27112370/make-readline-wait-for-input-in-r
# If run all code together, 'readline()' will return an empty string (""), and the script will default to using all data.
subsetGroup = readline(paste("Enter group to be subset (options: All, ", comp1, " only, ", comp2, " only, ", "Same/Opposite): ", sep = ""))


if (subsetGroup == "All"){
  message("Using all data")
  DataMerge <- DataMerge.orig
  colorCategory <- colorCategory.orig
  labelName <- "all"
  
} else if (subsetGroup == paste(comp1, " only", sep = "")){
  message(paste("Using", comp1, "only data"))
  DataMerge <- subset(DataMerge.orig, Category %in% c(paste("Downregulated", comp1, "only"), paste("Upregulated", comp1, "only")))
  colorCategory <- subset(colorCategory.orig, Category %in% c(paste("Downregulated", comp1, "only"), paste("Upregulated", comp1, "only")))
  labelName <- paste(comp1, "only")
  
} else if (subsetGroup == paste(comp2, " only", sep = "")){
  message(paste("Using", comp2, "only data"))
  DataMerge <- subset(DataMerge.orig, Category %in% c(paste("Downregulated", comp2, "only"), paste("Upregulated", comp2, "only")))
  colorCategory <- subset(colorCategory.orig, Category %in% c(paste("Downregulated", comp2, "only"), paste("Upregulated", comp2, "only")))
  labelName <- paste(comp2, "only")
  
} else if (subsetGroup == "Same/Opposite"){
  essage("Using Same/Opposite FC data")
  DataMerge <- subset(DataMerge.orig, Category %in% c("Both Upregulated", "Both Downregulated", "Opposite FC"))
  colorCategory <- subset(colorCategory.orig, Category %in% c("Both Upregulated", "Both Downregulated", "Opposite FC"))
  labelName <- "Same-Opposite change"
  
} else {
  message("Invalid input. Using all data.")
  DataMerge <- DataMerge.orig
  colorCategory <- colorCategory.orig
  labelName <- "all"
}

### ============================================================================
### 6. Generate 4-quadrant Scatter Plot with Density Maps or Volcano Plots

### Scatter plot ###

SP <- ggplot(DataMerge, aes(x=empiricalLFC.x, y=empiricalLFC.y, color=Category)) + geom_point(size=5, alpha = 0.6) +
  scale_color_manual(values=colorCategory$Color)+ # color for "Both Downregulated", "Both Upregulated", "Downregulated comp1", "Downregulated comp2", "Opposite effect", "Upregulated comp1", "Upregulated comp2"
  geom_vline(xintercept = c(log2FC, -log2FC), size=1, color = 'black', linetype = "dashed") + 
  geom_hline(yintercept = c(log2FC, -log2FC), size=1, color = 'black', linetype = "dashed") +
  #geom_smooth(method="lm", aes(color=NULL), color = "grey30", size = 2)+ # add a linear fit line (aes color=NULL to prevent drawing line for individual color groups defined previously)
  #stat_cor(aes(color=NULL), size = 10)+ #add R and p-value from Pearson correlation
  theme_classic() +
  labs(x=paste(comp1,"log2FC",sep=" "), y= paste(comp2,"log2FC",sep=" ")) +
  theme(text = element_text(size=35, family = "sans"), # sans = Arial, serif = Times New Roman, mono = Courier New
        axis.text.x = element_text(size=35, family = "sans"), 
        axis.text.y = element_text(size=35, family = "sans"), 
        axis.line = element_line(size=3),
        legend.title = element_text(size=20),
        legend.text = element_text(size=17),
        plot.margin = unit(c(0,0,0,0),"cm"),
        legend.position = c(0.99,0.99),
        legend.justification=c("right","top"))+ 
  xlim(-8, 11) + ylim(-12, 11)+ #control x and y scale
  geom_label_repel(aes(label=LabelSelect), max.overlaps = 30, box.padding = 0.5, size=6) + # Control data point labels
  border(size=3) #+ geom_rug()
# https://stackoverflow.com/questions/58984274/how-do-i-label-a-point-with-ggplot-depending-on-two-different-conditions-being-t


# Figure with just legend
SP.legend <- ggplot(DataMerge, aes(x=empiricalLFC.x, y=empiricalLFC.y, color=Category)) + geom_point(alpha = 0.75) +
  scale_color_manual(values=colorCategory$Color)+ # color for "Both Downregulated", "Both Upregulated", "Downregulated comp1", "Downregulated comp2", "Opposite effect", "Upregulated comp1", "Upregulated comp2"
  theme_void() +
  lims(x = c(0,0), y = c(0,0))+
  theme(legend.key.size = unit(1,"cm"),
        legend.title = element_text(size=20),
        legend.text = element_text(size=16),
        plot.margin = unit(c(0,0,0,0),"cm"),
        legend.position = c(0.5,0.5)) +
  guides(colour = guide_legend(override.aes = list(size=8)))

# Figure without legend
SP.noLegend <- SP+theme(legend.position = "none")


### Option 1: Density plot ###
# Marginal density plot of x (top panel)
xdensity <- ggplot(DataMerge, aes(empiricalLFC.x, fill=Category.x, ..count..)) + # can substitute ..count.. with ..density.. or ..scaled..
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c("Downregulated"='blue3', "Upregulated" = 'red3',"Not DE" = 'grey30')) + theme_classic() +
  labs(x=NULL, y="Count")+
  theme(legend.position = "none",
        text = element_text(size=35, family = "sans"),
        plot.margin = unit(c(0,0,0,0),"cm")) + border(size=3)
xdensity
# Marginal density plot of y (right panel)
ydensity <- ggplot(DataMerge, aes(empiricalLFC.y, fill=Category.y, ..count..)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c("Downregulated"='blue3', "Upregulated" = 'red3',"Not DE" = 'grey30')) + theme_classic() +
  labs(x=NULL, y="Count")+ scale_y_continuous(breaks = seq(0, 2500, by = 1000))+ #to adjust y scale so that number won't be out of the edge of figure
  theme(legend.position = "none",
        plot.margin = unit(c(0,0,0,0),"cm"),
        text = element_text(size = 35, family = "sans")) + border(size=3) + rotate()
ydensity

## Export scatter plot with density plot
tiff(paste('LBM_4QuadrantScatter_',comp1,'_',comp2,'_densityMap_', labelName, '.tiff', sep=""), width = 10000, height = 9000, units = "px", res = 600, compression = 'lzw')
ggarrange(xdensity, SP.legend, SP.noLegend, ydensity, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2,1), heights = c(1, 2),
          common.legend = F)
dev.off()




### Option 2: Volcano plot ###
# create custom key-value pairs for 'Up-regulated', 'Down-regulated', 'NS' expression by fold-change
# set the base colour as 'black'
keyvals1 <- rep('black', nrow(sigGene1))
keyvals2 <- rep('black', nrow(sigGene2))

# set the base name/label as 'Mid'
names(keyvals1) <- rep('Not significant', nrow(sigGene1))
names(keyvals2) <- rep('Not significant', nrow(sigGene2))

# modify keyvals for transcripts with fold change > log2FC
keyvals1[which(sigGene1$empiricalLFC > log2FC)] <- 'red3' #color the genes red if its log2FC is greater than 0.6
names(keyvals1)[which(sigGene1$empiricalLFC > log2FC)] <- 'Up-regulated'
keyvals2[which(sigGene2$empiricalLFC > log2FC)] <- 'red3' #color the genes red if its log2FC is greater than 0.6
names(keyvals2)[which(sigGene2$empiricalLFC > log2FC)] <- 'Up-regulated'

# modify keyvals for transcripts with fold change < -log2FC
keyvals1[which(sigGene1$empiricalLFC < -log2FC)] <- 'royalblue4'
names(keyvals1)[which(sigGene1$empiricalLFC < -log2FC)] <- 'Down-regulated'
keyvals2[which(sigGene2$empiricalLFC < -log2FC)] <- 'royalblue4'
names(keyvals2)[which(sigGene2$empiricalLFC < -log2FC)] <- 'Down-regulated'

# modify keyvals for transcripts with p-value > pValue
keyvals1[which(sigGene1$padj > pValue)] <- 'black' #color the genes that have log2FC>0.6 but also adjP>0.05 back to black (not significant genes)
names(keyvals1)[which(sigGene1$padj > pValue)] <- 'Not significant'
keyvals2[which(sigGene2$padj > pValue)] <- 'black' #color the genes that have log2FC>0.6 but also adjP>0.05 back to black (not significant genes)
names(keyvals2)[which(sigGene2$padj > pValue)] <- 'Not significant'

unique(names(keyvals1))
unique(names(keyvals2))



# Input data format: rownames = gene name, x axis = log2FC, y axis = adjusted p-value
xVolcano <- EnhancedVolcano(sigGene1,
                lab = sigGene1$gene, #lab = label
                x = 'empiricalLFC', 
                y = 'padj',
                selectLab =  DataMerge$gene[abs(DataMerge$empiricalLFC.x)>SelectFC & DataMerge$padj.x<SelectP],
                #selectLab = c('Gene1'), # only label key genes, can type 'nothing' (or anything random) to not display any gene name
                #selectLab = sigGenes, #only label the genes defined by top significant genes
                labSize = 5, #gene label size
                labFace = "plain",
                #labhjust = 0.5,
                #labvjust = 1,
                boxedLabels = T,
                drawConnectors = T,
                widthConnectors = 1,
                typeConnectors = 'closed',
                arrowheads = T,
                xlim = c(min(sigGene1$empiricalLFC, na.rm = TRUE) - 0.5, max(sigGene1$empiricalLFC, na.rm = T) + 0.5), #control x axis range
                ylim = c(0, max(-log10(sigGene1$padj), na.rm = TRUE) + 0.5),
                borderWidth = 1.5,
                xlab = NULL, #x and y-axis labels
                ylab = bquote(~-Log[10]~adj~italic(P)),
                axisLabSize = 35,
                pCutoff = pValue, #p-value cutoff
                FCcutoff = log2FC, #fold change cutoff
                cutoffLineWidth = 1,
                pointSize = 4,
                #legend=c('NS','Log2 fold-change','Adjusted P value',
                #         'Adjusted P value & Log2 fold-change'), #legend for pre-set colors
                legendPosition = 'none',
                legendLabSize = 23,
                colAlpha = 0.75, #transparency (1=100% opaque, 0=100% transparent)
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                #colOverride = keyvals#custom color defined in keyvals above
                colCustom = keyvals1,
                title = NULL,
                subtitle = NULL,
                caption = NULL
                #col = c('black', 'black', 'black', 'red3') #color of the 4 regions (NS, >log2FC, >-log10P, >log2FC&-log10P)
                )

yVolcano <- EnhancedVolcano(sigGene2,
                lab = sigGene2$gene, #lab = label
                x = 'empiricalLFC', 
                y = 'padj',
                selectLab = DataMerge$gene[abs(DataMerge$empiricalLFC.y)>SelectFC & DataMerge$padj.y<SelectP],
                #selectLab = c('Gene1'), # only label key genes, can type 'nothing' (or anything random) to not display any gene name
                #selectLab = sigGenes, #only label the genes defined by top significant genes
                labSize = 5, #gene label size
                labFace = "plain",
                #labhjust = 0.5,
                #labvjust = 1,
                boxedLabels = T,
                drawConnectors = T,
                widthConnectors = 1,
                typeConnectors = 'closed',
                arrowheads = T,
                xlim = c(min(sigGene2$empiricalLFC, na.rm = TRUE) - 0.5, max(sigGene2$empiricalLFC, na.rm = T) + 0.5), #control x axis range
                ylim = c(0, max(-log10(sigGene2$padj), na.rm = TRUE) + 0.5),
                borderWidth = 1.5,
                xlab = NULL, #x and y-axis labels
                ylab = bquote(~-Log[10]~adj~italic(P)),
                axisLabSize = 35,
                pCutoff = pValue, #p-value cutoff
                FCcutoff =log2FC, #fold change cutoff
                cutoffLineWidth = 1,
                pointSize = 4,
                legendPosition = 'none',
                legendLabSize = 23,
                colAlpha = 0.75, #transparency (1=100% opaque, 0=100% transparent)
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                #colOverride = keyvals#custom color defined in keyvals above
                colCustom = keyvals2,
                title = NULL,
                subtitle = NULL,
                caption = NULL
                #col = c('black', 'black', 'black', 'red3') #color of the 4 regions (NS, >log2FC, >-log10P, >log2FC&-log10P)
                ) + rotate()

## Export scatter plot with volcano plot
tiff(paste('LBM_4QuadrantScatter_',comp1,'_',comp2,'_volcano_', labelName, '.tiff', sep=""), width = 10000, height = 10000, units = "px", res = 600, compression = 'lzw')
ggarrange(xVolcano, SP.legend, SP.noLegend, yVolcano, 
          ncol = 2, nrow = 2,  align = "hv", 
          widths = c(2,1), heights = c(1, 2))
dev.off()



### Export full data as Excel file ###
# Rename comparison-specific columns in DataMerge.orig
DataMerge.export <- DataMerge.orig
colnames(DataMerge.export)[colnames(DataMerge.export) == "empiricalLFC.x"] <- paste("log2FC", comp1, sep = ".")
colnames(DataMerge.export)[colnames(DataMerge.export) == "padj.x"]         <- paste("padj", comp1, sep = ".")
colnames(DataMerge.export)[colnames(DataMerge.export) == "Category.x"]     <- paste("Category", comp1, sep = ".")
colnames(DataMerge.export)[colnames(DataMerge.export) == "empiricalLFC.y"] <- paste("log2FC", comp2, sep = ".")
colnames(DataMerge.export)[colnames(DataMerge.export) == "padj.y"]         <- paste("padj", comp2, sep = ".")
colnames(DataMerge.export)[colnames(DataMerge.export) == "Category.y"]     <- paste("Category", comp2, sep = ".")

# Export as Excel
exportList = list('Data' = DataMerge.export, 'Count' = geneCount)
write.xlsx(exportList, file = paste("LBM_4QuadrantScatter_",comp1,"_",comp2,"_data.xlsx", sep=""))
