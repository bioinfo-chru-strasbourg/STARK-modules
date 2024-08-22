##########################################################################
# PlotCoverage script       Version: 1
# Description:         		R script to plot coverage 
##########################################################################

################## Context ##################
# Draw barplot and boxplot of coverage data in png format
########## Note ########################################################################################
# DEV v1 11/03/2024
# Changelog
#   - optparse script, bugfixes & update code
#   - change png to pdf to avoid some size limitation
########################################################################################################

suppressPackageStartupMessages({
	library(ggplot2)
	library(reshape)
	library(reshape2)
	library(dplyr)
	library(gtools)
	library(optparse)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="Input coverage file", dest='input'),
  make_option(c("-o", "--output1"), type="character", help="Output pdf boxplot", dest='output1'),
  make_option(c("-p", "--output2"), type="character", help="Output pdf barplot", dest='output2')
)

opt_parser <- OptionParser(option_list=option_list)
options <- parse_args(opt_parser)


BoxPlotCoverage <- function(SAMPLE, INPUT, OUTPUT){
    coverage <- read.csv(INPUT, header=TRUE, check.names=FALSE, row.names=1, sep="\t")
    coverageMeans <- colMeans(coverage)
    coverageSample <- coverageMeans[SAMPLE]
    facteur <- coverageSample / coverageMeans
    coveragefact <- t(t(coverage) * facteur)
    coveragefactdf <- data.frame(Sample=coveragefact)
    coveragefactdfexclude <- coveragefactdf[,-which(colnames(coveragefactdf)==paste("Sample.", SAMPLE, sep=""))]
    coveragefactdfexclude[coveragefactdfexclude == 0] <- NA
    coverageMed <- apply(coveragefactdfexclude, 1, median, na.rm = TRUE)
    coverageMedian <- coveragefact / coverageMed
    coverageNormal <- log2(coverageMedian)
    
    lab <- mixedsort(rownames(coverageNormal))
    is.na(coverageNormal) <- sapply(coverageNormal, is.infinite)
    coverageNa <- coverageNormal[is.na(coverageNormal[,SAMPLE]),,drop=FALSE]
    coverageOut <- which(rownames(coverageNormal[lab,]) %in% rownames(coverageNa))

    my_data <- data.frame(Genomic_Region=row.names(coverageNormal), Sample=coverageNormal[,SAMPLE])
    my_data_box <- data.frame(Genomic_Region=row.names(coverageNormal), Sample=coverageNormal)
    y_data.long <- melt(my_data_box, id.vars="Genomic_Region")

    sample_list <- paste(colnames(coverage), collapse=',')
    nrows <- min(dim(coverage)[1] * 10, 10000)  # Limit width to avoid excessive size

    p <- ggplot() + 
        geom_boxplot(data=y_data.long, aes(x=factor(Genomic_Region, levels=lab), y=value)) +
        geom_point(data=my_data, aes(x=factor(Genomic_Region, levels=lab), y=Sample), shape=17, color="green") +
        geom_vline(xintercept=coverageOut, linetype="solid", color = "red1", linewidth=4) +
        geom_hline(yintercept = c(1, -1, 2, -2), linetype = "dashed", color = "red2") +
        labs(title = "Boxplot coverage mean target sample vs all (normalised)", subtitle = paste("[INFO] Target Sample: ", SAMPLE, "\n[INFO] Processed samples: ", sample_list, "\n[INFO] Warning: Without sexual chromosomes, Normalisation: For each sample and for each region the coverage mean is divide by the total coverage mean of the sample and multiply by the total coverage mean of the target sample (values are divide by the median (Graph:median=0), in log2 )"), x="Genomic Region (chromosome:start-end)", y="Coverage (normalised) (log2)") +
        theme(axis.text=element_text(size=7), axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1), plot.title = element_text(size=9), plot.subtitle = element_text(size=6))

    tryCatch({
        ggsave(basename(OUTPUT), p, width = nrows, limitsize = FALSE, path = dirname(OUTPUT))
    }, error = function(e) {
        message("Error saving boxplot: ", e)
    })
}

BarPlotCoverage <- function(SAMPLE, INPUT, OUTPUT){
    coverage <- read.csv(INPUT, header=TRUE, check.names=FALSE, row.names=1, sep="\t")
    my_data <- data.frame(Genomic_Region=row.names(coverage), Coverage_Mean_Sample=coverage[,SAMPLE])
    sample_list <- paste(colnames(coverage), collapse=',')
    nrows <- min(dim(coverage)[1] * 10, 10000)  # Limit width to avoid excessive size
    lab <- mixedsort(rownames(coverage))

    p <- ggplot() + 
        geom_bar(data=my_data, stat="identity", fill="skyblue", aes(x=factor(Genomic_Region, levels=lab), y=Coverage_Mean_Sample)) +
        labs(title = "Barplot coverage mean per base per region", subtitle = paste("[INFO] Target Sample: ", SAMPLE), x="Genomic Region (chromosome:start-end)", y="Coverage mean s(bp)") +
        theme(axis.text=element_text(size=7), axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1), plot.title = element_text(size=9), plot.subtitle = element_text(size=6))

    tryCatch({
        ggsave(basename(OUTPUT), p, width = nrows, limitsize = FALSE, path = dirname(OUTPUT))
    }, error = function(e) {
        message("Error saving barplot: ", e)
    })
}

# Read the coverage file to extract sample names from the header
coverage <- read.csv(options$input, header=TRUE, check.names=FALSE, row.names=1, sep="\t")
sample.names <- colnames(coverage)

for (i in 1:length(sample.names)) {
    BoxPlotCoverage(sample.names[i], options$input, options$output1)
	BarPlotCoverage(sample.names[i], options$input, options$output2)
}

# Execute the functions with the samples derived from the header
for (SAMPLE in sample.names) {
    BoxPlotCoverage(SAMPLE, options$input, options$output1)
    BarPlotCoverage(SAMPLE, options$input, options$output2)
}
