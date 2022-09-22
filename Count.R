# Workflow for DE analysis using DESeq2 #
# Import quantification via salmon tool

## Author : Lucie RIGOLOT
## Hôpitaux Universitaires de Strasbourg

# Sources :
# http://master.bioconductor.org/packages/release/workflows/vignettes/rnaseqGene/inst/doc/rnaseqGene.html#the-deseqdataset-object-sample-information-and-the-design-formula
# http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#differential-expression-analysis
# https://hbctraining.github.io/DGE_workshop_salmon_online/schedule/links-to-lessons.html
# http://genomicsclass.github.io/book/pages/rnaseq_gene_level.html

# 15/09/2022 : add argparse options and refactor some code - TL

# Load libraries
library(tximport)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(RColorBrewer)
library(pheatmap)
library(DEGreport)
library(heatmaply)
library(EnhancedVolcano)
library(genefilter)
library(ReportingTools)
library(regionReport)
library(pcaExplorer)

######Parsing input options and setting defaults########
option_list<-list(
    make_option('--data',default='data', help='Folder were the data are stored',dest='datafolder'),
    make_option('--folder',default='results', help='Folder were to save the results',dest='resultfolder')
)

opt<-parse_args(OptionParser(option_list=option_list))

FolderOutput=opt$resultfolder
DataInput=opt$datafolder


### Step 1 ### Setup
# Check the working directory
getwd()
# To change the working directory, use : setwd("/PathToYourWorkingdirectory/")

#Create output dir
dir.create(FolderOutput)

##Create a dataframe with transcript ID & gene ID with ENSEMBL ref
library("EnsDb.Hsapiens.v75")
edb<-EnsDb.Hsapiens.v75
k <- keys(edb, keytype = "TXNAME")
tx2gene <- select(edb, k, "GENENAME", "TXNAME")
tx2gene<-tx2gene[,1:2]
head(tx2gene)
write.table(tx2gene, "./tx2gene_ENSEMBL_v75.txt", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

# Load data (Salmon quantification files = quant.sf files) using tximport

## 1 # Create a "data" folder in your working directory to import Salmon data (ENSEMBL reference)
# Salmon data format = quant.sf file in a folder, with folder name = sample name

## 2 # List all directories containing data
samples <- list.files(path = "./data", full.names = T) #For all directories in "data" folder
# Use option "pattern" if necessary. For example, to have all directories with "salmon" 
# in the name, use : salmon_samples <- list.files(path = "./data", full.names = T, pattern="salmon?") 

## 3 # Obtain a vector of all filenames including the path
files <- file.path(samples, "quant.sf")

## 4 # Manipulate the paths to create names for each element matching the sample names
names(files) <- str_replace(samples, "./data/", "") %>% 
  str_replace("/quant.sf", "")
files # Check the names of the files, they must match the sample names

## 5 # Import Annotation file (data frame linking transcript ID to gene ID = "tx2gene") from your working directory
tx2gene<-read.table("./tx2gene_ENSEMBL_v75.txt", header=T, stringsAsFactors = FALSE, sep="\t")
head(tx2gene) # Check column names : transcript name must be the first column

## 6 # Run tximport for salmon count
txi <- tximport(files, type = "salmon", tx2gene=tx2gene, ignoreTxVersion = TRUE, countsFromAbundance="lengthScaledTPM")
# Use option countsFromAbundance="lengthScaledTPM" to obtain "raw" counts values 
# (non-normalized) from the TPM (recommended for DE analysis)
# Option countsFromAbundance="scaledTPM" to take the TPM scaled up to library size as "raw" counts
# Default option : takes TPM as scaled values (Abundance) and NumReads as "raw" counts (Counts)

## 7 # View and export counts data
# Display the counts
txi$counts %>% View()

# Round the counts and write it into an object
data <- txi$counts %>% round(digits=0) 
data<-cbind(rownames(data),data)
colnames(data)[1]<-"Gene"
head(data)
# Save counts ("raw" counts) as text file
write.table(data, "./results/tximport_counts_per_gene.txt", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

## 8 # Create a sampletable/metadata 
# Import from the working directory a file matching the samples to the corresponding groups that are investigated
# TODO create the sample table with a sample row, and each condition's rows you need
sample_table <- read.csv("./sample_table.csv", sep=";")
meta <- data.frame(sample_table, row.names = colnames(txi$counts))
# Check that the row names of the metadata are the same as the column names of the counts data 
# and that the column names of the counts data are the same as the sample names of meta
all(colnames(txi$counts) == rownames(meta)) # Value must be "TRUE"
all(colnames(txi$counts) == meta$sample) # Value must be "TRUE"

### Step 2 ### Create the dds object
dds <- DESeqDataSetFromTximport(txi, 
                                colData = meta, 
                                design = ~ condition) # design is the name(s) of the column(s) of the sample_table with the info on the groups/conditions you want to compare

## 9 # *Optional* ### Pre-filtering the dataset (dds object) before DE analysis
## Number of rows in the dataset before filtering
nrow(dds)
## Delete row if 0 counts for all samples
keep <- rowSums(counts(dds)) >= 1  # Keep if sum of row values >= 1 
dds <- dds[keep,]
nrow(dds)

## Other filters (select only genes with significant amount of counts in a significant number of samples)
## Ex : Delete row if less than 3 samples have 3 counts or less for this gene
keep <- rowSums(counts(dds) >= 3) >= 3  # Keep if >= 3 counts in >= 3 samples
dds <- dds[keep,]
nrow(dds)


### Step 3 ### Exploratory data analysis (PCA & hierarchical clustering) - QC control
# Identifying outliers and sources of variation in the data
# Hierarchical clustering places similar samples together, represented by the 
# tree structure. High correlations across the board (> 0.999) suggest no outlying sample(s).
# The samples are supposed to be clustering together by sample group

# Transform counts for data visualization (rlog transformation)
rld <- rlog(dds, blind=TRUE)

# PCA uses the top 500 most variable genes to determine the similarity of the samples
pdf(file = "./results/PCA.pdf")
plotPCA(rld, intgroup="condition")
dev.off()

# Extract the rlog matrix from the object and compute pairwise correlation values
rld_mat <- assay(rld)
rld_cor <- cor(rld_mat)

# Plot heatmap for unsupervised hierarchical clustering (sample to sample distance)
pdf(file = "./results/HeatMap.pdf")
pheatmap(rld_cor, annotation = meta)
dev.off()

### Step 4 ### Differential expression analysis with DESeq2
# **Optional step** - Re-create DESeq2 dataset if the design formula has changed after QC analysis
# For example, include other sources of variation (if the info is in the sample_table) using
# "dds <- DESeqDataSetFromTximport(txi, colData = meta, design = ~ covariable + condition)"
# Run DESeq2 differential expression analysis 
dds <- DESeq(dds)
resultsNames(dds)  # Check the default group comparison
normalized_counts <- counts(dds, normalized=TRUE)
normalized_counts <- cbind(rownames(normalized_counts),normalized_counts) # Add a first column with gene names
colnames(normalized_counts)[1]<-"gene"
write.table(normalized_counts, file="./results/DESeq_normalized_counts.txt", sep="\t", quote=F, col.names=TRUE, row.names=FALSE)

### Step 5 ### Check the fit of the dispersion estimates (Quality check)
# The dispersion estimates reflect the variance in gene expression for a given
# mean value. With few replicates per group, the estimates of variation for each gene
# are often unreliable. DESeq2 shares information across genes to generate more accurate
# estimates of variation using a "shrinkage" method (DESeq2 assumes that genes with
# similar expression levels should have similar dispersion).
# This is a good plot to examine to ensure the data is a good fit for the DESeq2 model. 
# Black dots = dispersion estimates for each gene, Blue dots = shrunken dispersion values.
# We expect the data to generally scatter around the curve, with the dispersion decreasing 
# with increasing mean expression levels.
# With more replicates per condition, less shrinkage is applied to the dispersion estimates.

# Plot dispersion estimates
pdf(file = "./results/DESeq2_DispersionEstimates.pdf")
plotDispEsts(dds)
dev.off()

### Step 6 ### Create contrasts to perform Wald testing on the shrunken log2 foldchanges
# DESeq2 uses Wald test to identify genes that are differentially expressed 
# between two sample classes. Results for different comparisons can be extracted 
# depending on how many factors are used in the design formula and how many factor 
# levels are present.
# By default, DESeq compares the last condition against the first condition.
# Example of typical use : contrast <- c("condition", "level_to_compare", "base_level")
# Here, for the factor "condition", we want to compare "level_to_compare" versus "base_level"
# The fold change will be calculated with "level_to_compare" as numerator and "base_level"
# as denominator. Fold-change = level_to_compare/base_level.
# By default, DESeq2 allows for the shrinkage of the LFC estimates toward zero 
# when the information for a gene is low (low counts or high dispersion values). 
# LFC shrinkage uses information from all genes to generate more accurate estimates.
# The coef will be dependent on what your contrast was and should be identical to what 
# is stored in resultsNames()
# Use "resultsNames(dds)" to see the values that you can provide to the "coef" argument

contrast <- c("condition", "level_to_compare", "base_level") # Specify contrast for comparison of interest
res <- results(dds, contrast = contrast, alpha = 0.05) # Output results of Wald test for contrast of interest
summary(res)
res_tableOE_unshrunken <- res # Save unshrunken results
res <- lfcShrink(dds, coef = "condition_base_level_vs_level_to_compare", type = "apeglm") # Shrink the log2 fold changes (LFC) to be more accurate
res_tableOE_shrunken <- res #Save shrunken results

# Compare unshrunken vs shrunken results with MA plot
# MA plot shows the mean of the normalized counts versus the log2foldchanges
# for all genes tested. The genes that are significantly DE are colored.
pdf(file = "./results/MAplot_unshrunken_foldchanges.pdf")
plotMA(res_tableOE_unshrunken, ylim=c(-2,2)) # MA plot using unshrunken fold changes
dev.off()
pdf(file = "./results/MAplot_shrunken_foldchanges.pdf")
plotMA(res_tableOE_shrunken, ylim=c(-2,2)) # MA plot using shrunken fold changes
dev.off()


### Step 7 ### Output significant results
# This function reports the number of genes up- and down-regulated at the selected
# threshold (padj/alpha), the number of genes that were tested (genes with non-zero 
# total read count), and the number of genes not included in multiple test correction 
# due to a low mean count
# summary(res, alpha = 0.05)

# Set thresholds (padj cutoff = FDR)
padj.cutoff <- 0.01

# Turn the results object into a tibble for use with tidyverse functions
res_tbl <- res %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# Subset the significant results
sig_res_01 <- dplyr::filter(res_tbl, padj < padj.cutoff)
summary(sig_res_01)
sig_res_LFC1<- dplyr::filter(sig_res_05, log2FoldChange > 1)
summary(sig_res_LFC1)
sig_res_LFC0<- dplyr::filter(sig_res_05, log2FoldChange > 0)
sig_res_LFCmoins1 <- dplyr::filter(sig_res_05, log2FoldChange < -1)
summary(sig_res_LFCmoins1)
sig_res_LFCmoinsde0<- dplyr::filter(sig_res_05, log2FoldChange < 0)

# Save results
write.table(sig_res, "./results/DEgenes_sig_res05.txt", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
write.table(sig_res_LFC1, "./results/DEgenes_sig_res05_LFC1.txt", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)
write.table(sig_res_LFCmoins1, "./results/DEgenes_sig_res05_LFC-1.txt", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)


### Step 8 ### Visualize results: volcano plots, heatmaps, normalized counts plots of top genes...
## 1 # Create reports from DESeq2 results with special R packages
## 1.1 # ReportingTools = HTML report with counts plots for every DE gene 
## (use dds with default parameters)
des2Report <- HTMLReport(shortName = 'ReportingTools_DE', 
                         title = 'ReportingTools',
                         reportDirectory = "./results/DEreports")
publish(dds,des2Report, pvalueCutoff=0.05,
        annotation.db="EnsDb.Hsapiens.v75", 
        factor = colData(dds)$condition,
        reportDir="./results/DEreports")
finish(des2Report)

## 1.2 # RegionReport = HTML and PDF Report with visualization plots
regionreport <- DESeq2Report(dds, project = "DE report with RegionReport", res=NULL,
                             intgroup = "condition", outdir = "./results/DEreports",
                             output = "RegionReport", theme = theme_bw())

## *Info* : if res=NULL, then results will be used on dds with default parameters (last 
## condition vs first condition). Use the "res" object created at step 6 to generate 
## a RegionReport for the contrast of interest = replace res by is("res object", "class")
## Note : Error while executing regionReport with res object (incorrect name of color).

## 1.3 # pcaExplorer = Opens a webpage where the DESeq2 results are uploaded, 
# upload the annotation file and construct a report using the different options
pcaExplorer(dds)

## 2 # Plotting significant DE genes
## 2.1 # Plot expression for a single gene of interest (e.g. CD180)
# Use DESeq2 plotCounts function to plot expression for a single gene and save as PDF
# The specified gene needs to match with the original input to DESeq2 (GeneID, EnsemblID...)
pdf(file = "./results/CD180_expression.pdf")
plotCounts(dds, gene="CD180", intgroup="condition")
dev.off()
# Save plotcounts to a data frame object
d <- plotCounts(dds, gene="CD180", intgroup="condition", returnData=TRUE)
# Use ggplot2 to plot the normalized counts using the samplenames (rownames(d) as labels)
# TODO gene should be an arguement
pdf(file = "./results/CD180_expression_labeled.pdf")
ggplot(d, aes(x = condition, y = count, color = condition)) + 
  geom_point(position=position_jitter(w = 0.1,h = 0)) +
  geom_text_repel(aes(label = rownames(d))) + 
  theme_bw() +
  ggtitle("CD180") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


## 2.2 # Plot expression for multiple genes of interest (e.g. top 20 DE genes)

## Order results by padj values
top20_sigOE_genes <- res_tbl %>% 
  arrange(padj) %>% 	#Arrange rows by padj values
  pull(gene) %>% 		#Extract character vector of ordered genes
  head(n=20)		#Extract the first 20 genes

## Normalized counts for all significant DE genes
normalized_counts <- as_tibble(normalized_counts)
AllDEgenes<-sig_res$gene
AllDEgenes_norm <- normalized_counts %>%
  dplyr::filter(gene %in% AllDEgenes)
write.table(AllDEgenes_norm, "./results/AllDEgenes_norm_counts.txt", sep="\t", quote=FALSE, row.names = FALSE, col.names = TRUE)

## Normalized counts for top 20 significant genes
top20_sigOE_norm <- normalized_counts %>%
  dplyr::filter(gene %in% top20_sigOE_genes)

## Gathering the columns to have normalized counts to a single column
gathered_top20_sigOE <- top20_sigOE_norm %>%
  gather(colnames(top20_sigOE_norm)[2:7], key = "sample", value = "normalized_counts")
top20_sigOE_final <- gathered_top20_sigOE %>% dplyr::select(gene, sample) #keep gene and sample columns
top20_sigOE_final <- top20_sigOE_final %>% 
  cbind(as.numeric(gathered_top20_sigOE$normalized_counts)) #join normalized_counts column as numeric format
colnames(top20_sigOE_final)[3]<-"normalized_counts"

## **Optional checks**
# View(top20_sigOE_final)    #Check the column header in the final "gathered" data frame
# summary(top20_sigOE_final) #Check that the normalized_counts are numeric

## Merge metadata information with the normalized counts data 
top20_sigOE_final <- inner_join(meta, top20_sigOE_final)
#inner_join() will merge 2 data frames with respect to the "sample" column, 
#i.e. a column with the same column name in both data frames.

## Plot using ggplot2
pdf(file="./results/Top20DEgenes.pdf")
ggplot(top20_sigOE_final) +
  geom_point(aes(x = gene, y = normalized_counts, color = condition)) +
  scale_y_log10() +
  xlab("Genes") +
  ylab("log10 Normalized Counts") +
  ggtitle("Top 20 Significant DE Genes LLC LPS vs NS") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()


## 3 # Gene clustering (Heatmap) using Heatmaply

## Use transformed counts by rlog transformation ("rld") and select top 50 DE genes 
topVarGenes <- head(order(rowVars(assay(rld)), decreasing = TRUE), 50)
mat  <- assay(rld)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)

## Create interactive heatmap using heatmaply and save as html file 
heatmaply(mat, file = "./results/heatmaply_plot_defaultcolors.html") #Save heatmap

## Personalize colors (e.g. blue for low and red for high expression), save and open
heatmaply(mat, file = "./results/heatmaply_plot_50genes_redbluecolors.html", 
          scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(  
            low = "blue", 
            high = "red", 
            midpoint = 0,       #Change threshold value between low and high colors
            limits = c(-4, 4), #Change scale
             ))

## **Optional** - Open heatmap in browser with browseURL("path_to_html_file")
## e.g. browseURL("heatmaps/heatmaply_allexon05_red_blue.html")


## 4 # Volcano plot (EnhancedVolcano)

pdf(file="./results/DE_Volcanoplot_LBns_vs_LLCns.pdf")
EnhancedVolcano(sig_res_05,
                lab = sig_res_05$gene,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'LB LPS vs NS (LFC>1, FDR<5%)',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 4.0,
                legendLabSize = 9,
                col=c('grey','grey','grey', 'red2'),
                colAlpha = 1)
dev.off()

pdf(file="./results/DE_Volcanoplot_LBns_vs_LLCns_TLRlabeled_LFC1.pdf")
EnhancedVolcano(sig_res,
                lab = sig_res$gene,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'LB NS vs LLC NS (LFC > 1, FDR < 5 %)',
                selectLab = c('TLR1','TLR2','TLR3','TLR4', 'TLR5', 'TLR6', 'TLR7', 'TLR8', 'TLR9', 'TLR10', 'CD180', 'LY86', 'LY96', 'NFKB1', 'NFKBIA', 'NFKBID'),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                pointSize = 1.0,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                col=c('grey', 'grey', 'grey', 'red2'),
                colAlpha = 0.75,
                legendLabels=c('Not sig.','Log (base 2) FC','p-value',
                               'p-value & Log (base 2) FC'),
                legendPosition = 'right',
                legendLabSize = 10,
                legendIconSize = 5.0,
                drawConnectors = TRUE,
                widthConnectors = 1.0,
                colConnectors = 'black')
dev.off()


### Step 9 ### Functional analysis: GO over-representation, GSEA with clusterProfiler

## 1 # Gene Ontology (GO) over-representation analysis with clusterProfiler

# Load libraries
library(DOSE)
library(pathview)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggnewscale)

## Ensembl and EntrezID are needed for this functional analysis
annotations_ahb <- read.csv("annotations_ahb.csv") #Upload annotation file
# Keep the first identifier for these multiple mapping cases.
annotations_ahb$entrezid <- map(annotations_ahb$entrezid,1) %>%  unlist()
# Determine the indices for the non-duplicated genes
non_duplicates_idx <- which(duplicated(annotations_ahb$gene_name) == FALSE)
# Return only the non-duplicated genes using indices
annotations_ahb <- annotations_ahb[non_duplicates_idx, ]

## Merge the AnnotationHub dataframe with the results 
res_ids <- left_join(res_tbl, annotations_ahb, by=c("gene"="gene_name"))   

## Create background dataset for hypergeometric testing = all genes tested                 
allOE_genes <- as.character(res_ids$gene)

## Extract significant results
sigOE <- dplyr::filter(res_ids, padj < 0.05)
sigOE_genes <- as.character(sig_res_LFCmoins1$gene)

## Run GO enrichment analysis (significant OE genes within universe=all OE genes)
ego <- enrichGO(gene = sigOE_genes, 
                universe = allOE_genes,
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = F)

## Output results from GO analysis to a table
cluster_summary <- as.data.frame(ego)
write.csv(as.data.frame(cluster_summary), 
          file="./results/clusterProfiler_LB_lps_vs_ns_LFCunder1.csv")

## Visualize clusterProfiler results :

# Dotplot #

pdf(file = "./results/clusterprofiler_dotplot_top20_LBns_vs_LLCns.pdf")
dotplot(ego, showCategory=20)
dev.off()

# Enrichmap #

#Add similarity matrix to the termsim slot of enrichment result
ego <- enrichplot::pairwise_termsim(ego)
#Enrichmap clusters the 20 most significant (by padj) GO terms to visualize 
#relationships between terms
pdf(file = "./results/clusterprofiler_enrichmap_LBns_vs_LLCns.pdf")
emapplot(ego, showCategory = 20, label_size = 1)
dev.off()

# Category netplot #

#To color genes by log2 fold changes, we need to extract the log2 fold changes 
#from our results table creating a named vector
OE_foldchanges <- sigOE$log2FoldChange
names(OE_foldchanges) <- sigOE$gene
#Cnetplot details the genes associated with one or more terms - by default gives 
#the top 5 significant terms (by padj)
pdf(file = "./results/clusterprofiler_CategoryNetplot_LBns_vs_LLCns.pdf")
cnetplot(ego, 
         categorySize="pvalue", 
         showCategory = 10, 
         foldChange=OE_foldchanges, 
         vertex.label.font=2)
dev.off()

## 2 # GSEA analysis with clusterProfiler (work in progress)

## **Info**
# Functional class scoring (FCS) tools, such as GSEA, most often use the
# gene-level statistics or log2 fold changes for all genes from the differential
# expression results, then look to see whether gene sets for particular
# biological pathways are enriched among the large positive or negative fold changes.


### Step 10 ### Output the versions of all tools used in the DE analysis
sessionInfo()
