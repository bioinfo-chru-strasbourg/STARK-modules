##########################################################################
# PlotCoverage script       Version: 1
# Description:         		R script to plot coverage 
##########################################################################

################## Context ##################
# Draw barplot and boxplot of coverage datas in png format
########## Note ########################################################################################
# DEV v1 11/03/2024
# Changelog
#   - optparse script /bugfixes
########################################################################################################

library(ggplot2)
library(reshape2)
library(dplyr)
library(gtools)
library(optparse)

# Define the option parser
option_list <- list(
  make_option(c("-s", "--sample"), type="character", help="Sample name"),
  make_option(c("-i", "--input"), type="character", help="Input file path"),
  make_option(c("-o", "--output"), type="character", help="Output directory path")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
options <- parse_args(opt_parser)

BoxPlotCoverage <- function(SAMPLE, INPUT, OUTPUT){
	# Boxplot coverage data from coverage statistique
	coverage <- read.csv(INPUT, header=TRUE, check.names=FALSE, row.names=1, sep="\t")
	coverageMeans <- colMeans(coverage)
	coverageSample <- coverageMeans[SAMPLE]
	facteur <- coverageSample/coverageMeans
	coveragefact <- t(t(coverage)*facteur)
	coveragefactdf <- data.frame(Sample=coveragefact)
	coveragefactdfexclude <- coveragefactdf[,-which(colnames(coveragefactdf)==paste("Sample.",SAMPLE, sep=""))]
	coveragefactdfexclude[coveragefactdfexclude == 0] <- NA
	coverageMed <- apply(coveragefactdfexclude, 1, median, na.rm = TRUE)
	coverageMedian <- coveragefact/coverageMed
	coverageNormal <- log2(coverageMedian)

	# get outliers
	lab <- mixedsort(rownames(coverageNormal))
	is.na(coverageNormal) <- sapply(coverageNormal, is.infinite)
	coverageNa <- coverageNormal[is.na(coverageNormal[,SAMPLE]) > 0,,drop=FALSE]
	covergeOut <- which(rownames(coverageNormal[lab,])  %in% rownames(coverageNa))

	my_data <- data.frame(Genomic_Region=row.names(coverageNormal), Sample=coverageNormal[,SAMPLE])
	my_data_box <- data.frame(Genomic_Region=row.names(coverageNormal), Sample=coverageNormal)
	y_data.long<-melt(select(my_data_box, -ends_with(SAMPLE)))
	sample_list <- paste(colnames(coverage), collapse=',')
	nrows <- dim(coverage)[1]*10

	# do boxplot
	p <- ggplot() + 
	geom_boxplot(data=y_data.long, aes(x=factor(Genomic_Region, levels=lab), y=value)) +
	geom_point(data=my_data, aes(x=factor(Genomic_Region, levels=lab), y=Sample), shape=17, color="green") +
	geom_vline(xintercept=covergeOut, linetype="solid", color = "red1", size=4) +
	geom_hline(yintercept = c(1, -1, 2, -2), linetype = "dashed", color = "red2") +
	labs(title = "Boxplot coverage mean taget sample vs all(normalised)", subtitle = paste("[INFO] Target Sample: ", SAMPLE, "\n[INFO] Processed samples: ", sample_list, "\n[INFO] Warning: Without sexual chromosomes, Normalisation: For each sample and for each region the coverage mean is divide by the total coverage mean of the sample and multply by the total coverage mean of the target sample (values are divide by the median(Graph:median=0), in log2 )"), x="Genomic Region(chromosome:start-end)", y="Coverage(normalised)(log2)") +
	theme(axis.text=element_text(linewidth=0.5), axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1), plot.title = element_text(linewidth=0.5), plot.subtitle = element_text(linewidth=0.5))
	#theme(axis.text=element_text(size=7),axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1), plot.title = element_text(size=9), plot.subtitle = element_text(size=6))
  
  	# save plot
  	png(OUTPUT, width=nrows, height=400)
  	print(p)
  	dev.off()

}

BarPlotCoverage <- function(SAMPLE, INPUT, OUTPUT){
	# Barplot coverage data from coverage statistique
	coverage <- read.csv(INPUT, header=TRUE, check.names=FALSE, row.names=1, sep="\t")
	my_data <- data.frame(Genomic_Region=row.names(coverage), Coverage_Mean_Sample=coverage[,SAMPLE])
	sample_list <- paste(colnames(coverage), collapse=',')
	nrows <- dim(coverage)[1]*10
	lab <- mixedsort(rownames(coverage))

	# do barplot
	p <- ggplot() + 
		geom_bar(data=my_data, stat="identity", fill="skyblue", aes(x=factor(Genomic_Region, levels=lab), y=Coverage_Mean_Sample)) +
		labs(title = "Barplot coverage mean per base per region", subtitle = paste("[INFO] Target Sample: ", SAMPLE), x="Genomic Region(chromosome:start-end)", y="Coverage mean(bp)") +
		theme(axis.text=element_text(linewidth=0.5) ,axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1), plot.title = element_text(linewidth=0.5), plot.subtitle = element_text(linewidth=0.5))
		#theme(axis.text=element_text(size=7) ,axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1), plot.title = element_text(size=9), plot.subtitle = element_text(size=6))
  
	# save plot
  	png(OUTPUT, width=nrows, height=400)
  	print(p)
  	dev.off()
}
# Execute the functions with command-line arguments
BoxPlotCoverage(options$sample, options$input, options$output)
BarPlotCoverage(options$sample, options$input, options$output)