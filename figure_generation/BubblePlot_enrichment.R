# ------------------------------------------------------------------------------
# Title:       Bubble Plot for Visualizing Enrichment Results
# Author:      Albert Wang
# Last updated:        2025-07-12
# ------------------------------------------------------------------------------


### Summary
# This script generates a bubble plot to visualize enrichment analysis results (e.g. Gene Ontology and KEGG pathway analysis) from a tab-delimited input file.
# It creates a figure where pathways (or GO terms) are plotted along one axis and experimental conditions along the other.
# The size of each bubble represents the percentage of differentially expressed genes (%DEG), and the fill color reflects the FDR-adjusted p-value, allowing for intuitive comparison of enrichment significance and gene representation across conditions.


###=============================================================================
### Required packages

library(ggplot2)
library(stringr)
library(ggpubr)
library(reshape2)
library(scales)

## References
# https://jkzorz.github.io/2019/06/05/Bubble-plots.html


### ============================================================================

# Load enrichment result file (e.g., LBM_meta_pathway_FDR.txt)
dataFile <- file.choose()
df <- read.table(dataFile,
                     sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# Generate export filename based on input file name
fileName = basename(dataFile) #get file name
fileName = tools::file_path_sans_ext(fileName) #remove extension
exportName = paste(fileName,"_bubblePlot.tiff", sep="") #for naming export tiff file



# Save plot as high-resolution TIFF image
tiff(exportName, width = 2700, height = 2800, units = "px", res = 600, compression = 'lzw')
# Generate bubble plot
ggplot(df, aes(x = reorder(str_wrap(name,width=17),setOrder, decreasing = T), y = reorder(condition, setOrder, decreasing = F))) + 
  geom_point(aes(size = CountPercent, fill = FDR), alpha = 0.75, shape = 21) + 
  labs( x= "", y = "", size = "    %DEG", fill = "p-value")  + 
  theme(legend.key=element_blank(), 
        axis.text.x = element_text(colour = "black", size = 14, face = "bold", angle = 90, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black", face = "bold", size = 14),
        axis.ticks = element_blank(), # set no ticks on the axis
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.title = element_text(size = 12, face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2), 
        legend.position = "right") + 
  #scale_fill_manual(values = colours, guide = FALSE) + 
  scale_fill_gradientn(limits = c(1e-5,0.11), # Gradient fill scale for FDR values. Modify the scale if needed
    #breaks = c(1,0.05,0.01,1e-11),
    oob = squish, # handle out of bound values
    colours=c("red","orange", "black")) + # continuous values
  scale_size_continuous( # Scale for %DEG bubble size
    limits = c(round(min(df$CountPercent)/5), max(df$CountPercent)*1.2),
    range = c(1,20), 
    breaks = c(round(min(df$CountPercent)/5)*5, round(median(df$CountPercent)/5)*5,round(max(df$CountPercent)/5)*5) # break is rounded to the nearest 5
    ) + 
  scale_y_discrete(limits = rev(levels(df$CountPercent))) + coord_flip() # Flip coordinates for better readability

dev.off()
