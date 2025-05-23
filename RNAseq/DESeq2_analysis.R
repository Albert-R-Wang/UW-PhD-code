# https://bioconductor.org/packages/release/bioc/html/DESeq2.html

#### Data import ####
## import count data (NSCLC_data_RNA_counts_nameCorrected_Adeno.txt or NSCLC_data_RNA_counts_nameCorrected_Adeno_BTonly.txt)
countFile <- file.choose()
cts <- read.table(countFile,header=TRUE,row.names=1 ,sep = "\t", quote="\"", stringsAsFactors = FALSE, check.names = FALSE)

## import sample information
annFile <- file.choose()
coldata <- read.table(annFile,header=TRUE,row.names=1 ,sep = "\t", quote="\"", stringsAsFactors = FALSE, check.names = FALSE)

## The columns of the count matrix and the rows of the column data (information about samples) are the same order
#rownames(coldata) <- sub("fb", "", rownames(coldata)) #if need to delete certain text pattern in rowname
all(rownames(coldata) %in% colnames(cts)) # check if naming is consistent

cts <- cts[, rownames(coldata)] #reorder the count matrix to match the order of sample information file
all(rownames(coldata) == colnames(cts)) # check if order is the same

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ condition) #specify a column called "condition" in the sample information file. For example, "treated" and "untreated"
dds
# can add additional feature data as needed
#featureData <- data.frame(gene=rownames(cts))
#mcols(dds) <- DataFrame(mcols(dds), featureData)

#### Pre-filtering ####
# it's not necessary to pre-filter low count genes. Two useful reasons: 1) reduce memory size of the dds data object, thus increase speed within DESeq2, and 2) improve visualizations as low counts gene are not plotted
# filter by sum of counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
# alternative option: filter by ensuring at least X samples with a count of 10 or more
x = ncol(cts)*0.2 # 20% of the number of samples
keep <- rowSums(counts(dds) >= 10) >= X
dds <- dds[keep,]

#### factor levels ####
dds$condition <- relevel(dds$condition, ref = "High Immune Cluster") #specify reference level (ex. "untreated")
#dds$condition <- relevel(dds$condition, ref = "Lung Primary Tumor")

#### Differential expression analysis ####
dds <- DESeq(dds, betaPrior = F)
#Starting in version 1.16, betaPrior is F by default, this means that DESeq2 will not perform any shrinkage of log2FC. This is moved to the lfcShrink function so that newer shrinkage methods can be more easily apply as needed (https://support.bioconductor.org/p/95695/). Turn betaPrior on if want to match the result of earlier versions of DESEq2.
#The purpose of shrinkage of log2FC is to address the problem that lowly expressed genes tend to have high variability (https://www.biostars.org/p/340269/)
#It looks at the largest fold changes that are not due to low counts and uses these to inform a prior distribution. So the large fold changes from genes with lots of statistical information are not shrunk, while the imprecise fold changes are shrunk. This allows you to compare all estimated LFC across experiments, for example, which is not really feasible without the use of a prior. (https://support.bioconductor.org/p/77461/)
res <- results(dds, pAdjustMethod = "bonferroni") #default is Benjamini and Hochberg method (BH or fdr)
res

## LFC (log fold change) shrinkage
resultsNames(dds) #Show all possible comparisons that can be used in lfcShrink as coef. Can input name directly, or use the order number
resLFC = lfcShrink(dds, coef=2, type="apeglm", res = res) #The newer shrinkage methods/estimators (apeglm and ashr) outperform the Normal prior (or betaprior) in most cases (https://support.bioconductor.org/p/115900/). The default is "apeglm"
#Unlike betaprio, all estimators in lfcShrink does not change padj.betaPrior will give the p-value for the shrunken LFC, while lfcShrink is only giving the shruken LFC, and keeping the original p-value
resLFC



#### Export results ####
fileName = basename(countFile) #get file name
fileName = tools::file_path_sans_ext(fileName) #remove extension
write.csv(as.data.frame(resLFC), 
          file=paste("DESeq2Result_", fileName, ".csv", sep=""))


nc <- counts(dds, normalized = T)
write.csv(as.data.frame(nc), 
          file=paste("DESeq2NormCount_", fileName, ".csv", sep=""))

#### Differential expressed gene (DEG) threshold ####
FC = 2 # log2 fold change
adjp = 0.01 # adjusted p-values

# Determine significant DEGs based on fold-change and adjusted p-value cut-off
sigGenes <- rownames(subset(resLFC, (abs(log2FoldChange)>=FC & padj<=adjp )))

# Extract RNAseq data for the significant genes
sig.res.DE = subset(resLFC, rownames(resLFC) %in% sigGenes)


#### Variance Stablizing Transformation (VST) ####
# Adjusts for sequencing depth and stabilizes variance across the range of mean values.
# This transformation is important for visualization and clustering, as it reduces the influence of highly variable genes and makes expression values more comparable across samples.
vsd = vst(dds, blind = F) #transform counts into log2 scale (modeled to stablize variance) for visualization

# Extract the matrix of VST-transformed values
vst_mat = assay(vsd)
sig.vst_mat = subset(vst_mat, rownames(vst_mat) %in% sigGenes)


#### Volcano plot ####
vol = resLFC %>% filter(!is.na(padj)) # Exclude genes with NA padj
# Determine which genes are upregulated, downregulated, or not significant (NS) based on the DEG cutoff
vol$group = with(vol, ifelse(padj < adjp & log2FoldChange > FC, "Upregulated",
                             ifelse(padj < adjp & log2FoldChange < -FC, "Downregulated", "NS")))
numDEG = table(vol$group) # Count the number of genes that are upregulated, downregulated, or NS

ggplot(vol, aes(x=log2FoldChange, y=-log10(padj), color = group)) + geom_point(size = 2, alpha = 0.75) +
  scale_color_manual(values = c("Upregulated" = "red3", "Downregulated" = "royalblue3", "NS" = "gray20")) +
  annotate("text", x = max(vol$log2FoldChange), y = max(-log10(vol$padj), na.rm = TRUE), 
           label = paste("Upregulated:", numDEG["Upregulated"]), color = "red3", hjust = 1) +
  annotate("text", x = min(vol$log2FoldChange), y = max(-log10(vol$padj), na.rm = TRUE), 
           label = paste("Downregulated:", numDEG["Downregulated"]), color = "royalblue3", hjust = 0) +
  theme_minimal()
