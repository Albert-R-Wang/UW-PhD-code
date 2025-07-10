# ------------------------------------------------------------------------------
# Title:       Venn Diagram Analysis for DEG and GO Term
# Author:      Albert Wang
# Last updated:        2025-07-10
# ------------------------------------------------------------------------------


### Summary
# This R script processes differential gene expression (DEG) and gene ontology (GO) analysis results from up to four comparisons, generating Venn diagrams to visualize shared and unique features.
# It supports input of DEG or DAVID (https://davidbioinformatics.nih.gov/home.jsp) enrichment result files, filters significant results based on user-defined thresholds, and creates customized Venn diagrams using the ggvenn package.


###=============================================================================
### Required packages

library(ggplot2)
library(ggvenn)
#library(ggVennDiagram) # can plot color map in venn diagram based on number of count in each circle



## References
# https://github.com/yanlinlin82/ggvenn
# https://github.com/gaospecial/ggVennDiagram
# https://www.datanovia.com/en/blog/beautiful-ggplot-venn-diagram-with-r/
# https://www.datanovia.com/en/blog/venn-diagram-with-r-or-rstudio-a-million-ways/

### ============================================================================
### Venn Diagram

# Sections:
# 1. Setup Functions and Variables
# 2. Load Input Files
# 3. DEGs Venn Diagram: 2-Group
# 4. DEGs Venn Diagram: 3-Group
# 5. GO Terms from DAVID: 3-Group Venn Diagram
# 6. GO Terms from DAVID: 4-Group Venn Diagram
# 7. Generate and Save Venn Diagram Plot



### ============================================================================
### 1. Setup Functions and Variables

# Thresholds and column settings
pValue <- 0.05
log2FC <- 0.6
GO_columns <- c(2, 3, 6, 12)

# Define a helper function for loading tab-delimited text files
load_file <- function(prompt_message = "File selected"){
  filePath <- file.choose()
  message(prompt_message, ": ", basename(filePath))
  
  df <- read.table(filePath, header = TRUE, sep = "\t", quote = "\"", stringsAsFactors = FALSE, check.names = FALSE)
  
  return(df)
}

# Function to filter DEGs
filter_DEG <- function(df){
  subset(df, abs(empiricalLFC)>=log2FC & padj<=pValue)
}

# Function to filter GO terms
filter_GO <- function(df){
  # Filter for significant GO terms based on adjusted p-value
  DAVIDsubset <- subset(df, Benjamini <= 0.05)[, GO_columns]
  
  
  # Remove the GO ID prefix (e.g., "GO:0008150~") from each term using a regular expression.
  # The pattern "GO:[0-9]+~" matches:
  # - "GO:" — a constant prefix present in all terms,
  # - "[0-9]+" — one or more digits representing the GO identifier,
  # - "~" — a tilde that separates the ID from the actual term name.
  # The matched pattern is then replaced with an empty string (""), leaving only the descriptive term.
  DAVIDsubset$Term <- gsub("GO:+[0-9]+~", "", DAVIDsubset$Term)
  
  DAVIDsubset
}

### ============================================================================
### 2. Load Input Files

## For DEG result files: ensure they contain a column named "empiricalLFC" for log2 fold change and "padj" for adjusted p-values.
## For GO result files: directly load the text files downloaded from DAVID.


# Load DEG or GO result files (adjust the number of datasets as needed)
data1 <- load_file("Loading data1")
data2 <- load_file("Loading data2")
data3 <- load_file("Loading data3") # Optional
data4 <- load_file("Loading data4") # Optional


## NOTE: Execute the section corresponding to the number of comparisons in your analysis.


### ============================================================================
### 3. DEGs Venn Diagram: 2-Group


# Define the names of each comparison group (modify as needed)
comp1 <- "LTvsLN"
comp2 <- "BrMvsLT"


# Find DEGs
sigGene1 <- filter_DEG(data1)
sigGene2 <- filter_DEG(data2)

# Create a list of significant gene sets for each comparison group to use in the Venn diagram.
x <- list(
  A = sigGene1$gene,
  B = sigGene2$gene
)

names(x) <- c(comp1, comp2) # change the category names



### ============================================================================
### 4. DEGs Venn Diagram: 3-Group

# Define the names of each comparison group (modify as needed)
comp1 = "BrMvsLT (All)"
comp2 = "BrMvsLT (Meta)"
comp3 = "BrMvsLT (Sync)"


# Find DEGs
sigGene1 <- filter_DEG(data1)
sigGene2 <- filter_DEG(data2)
sigGene3 <- filter_DEG(data3)

# Create a list of significant gene sets for each comparison group to use in the Venn diagram.
x <- list(
  A = sigGene1$gene,
  B = sigGene2$gene,
  C = sigGene3$gene
)

names(x) = c(comp1, comp2, comp3)


### ============================================================================
### 5. GO Terms from DAVID: 3-Group Venn Diagram

# Define the names of each comparison group (modify as needed)
comp1 <- "LTvsLN upregulated"
comp2 <- "LTvsLN downregulated"
comp3 <- "BrMvsLT downregulated"

# Filter and clean GO terms
DAVIDsubset1 <- filter_GO(data1)

DAVIDsubset2 <- filter_GO(data2)

DAVIDsubset3 <- filter_GO(data3)




# Create a list of GO term sets for each comparison group to use in the Venn diagram.
x <- list(
  A = DAVIDsubset1$Term,
  B = DAVIDsubset2$Term,
  C = DAVIDsubset3$Term
)

names(x) <- c(comp1, comp2, comp3)

# Find any common terms
group1 <- intersect(DAVIDsubset1$Term, DAVIDsubset2$Term)
writeClipboard(group1)

group2 <- intersect(DAVIDsubset1$Term, DAVIDsubset3$Term)
writeClipboard(group2)

group3 <- intersect(DAVIDsubset2$Term, DAVIDsubset3$Term)
writeClipboard(group3)

### ============================================================================
### 6. GO Terms from DAVID: 4-Group Venn Diagram


# Define the names of each comparison group (modify as needed)
comp1 <- "LTvsLN upregulated"
comp2 <- "LTvsLN downregulated"
comp3 <- "BrMvsLT upregulated"
comp4 <- "BrMvsLT downregulated"

# Filter for significant GO terms based on adjusted p-value
DAVIDsubset1 <- filter_GO(data1)

DAVIDsubset2 <- filter_GO(data2)

DAVIDsubset3 <- filter_GO(data3)

DAVIDsubset4 <- filter_GO(data4)



# Create a list of GO term sets for each comparison group to use in the Venn diagram.
x <- list(
  A = DAVIDsubset1$Term,
  B = DAVIDsubset2$Term,
  C = DAVIDsubset3$Term,
  D = DAVIDsubset4$Term
)

names(x) <- c(comp1, comp2, comp3, comp4)


### ============================================================================
### 7. Generate and Save Venn Diagram Plot

tiff("LBM_VennDiagram_LTvsLN_BTvsLT_DEG.tiff", width = 4000, height = 3000, units = "px", res = 600, compression = 'lzw') # 4000px for 2 categories, 5000px for 3 categories
ggvenn(x, 
       fill_color = c("#E66100", "#5D3A9B", "blue", "red"),
       fill_alpha = 0.5,
       #show_elements = T,
       #label_sep = "\n",
       stroke_size = 2, # line width
       set_name_size = 7, # text size of category name
       text_size = 7,
       show_percentage = T,
       digits = 0
       )
dev.off()