##########################################################################
# DECON script          Version: 1
# Description:          R script to generate metrics from coverage data
##########################################################################

################## Context ##################
# Original R script from https://github.com/RahmanTeam/DECoN
# DECON is an ExomeDepth wrapper 
########## Note ########################################################################################
# PROD v1 21/11/2023
# Changelog
#   - Refactor code, remove install system, optparse script, update ExomeDepth 1.16
########################################################################################################

suppressPackageStartupMessages({
library(R.utils)
library(optparse)
library(ExomeDepth)
})

# Define options and parse command line arguments
option_list <- list(
    make_option("--rdata", help="Input summary RData file (required)", dest='data'),
    make_option("--mincorr", help='Minimum correlation to consider, default=0.98', default=0.98, dest='mincorr'),
    make_option("--mincov", help='Minimum coverage to consider, default=100', default=100, dest='mincov'),
    make_option("--tsv", default='./Metrics.tsv', help='Output metrics tsv file, default=./Metrics.tsv', dest='out')
)

opt_parser <- OptionParser(option_list=option_list)
options <- parse_args(opt_parser)

# Function to stop execution if a required value is missing
stop_if_missing <- function(val, message) {
    if (is.null(val) || length(val) == 0) {
        stop(message)
    }
}

# Function to load data from RData file
load_data <- function(rdata_file) {
    stop_if_missing(rdata_file, "ERROR: no Rdata summary file provided -- Execution halted")
    load(rdata_file)
    return(list(counts = counts, bed.file = bed.file, sample.names = sample.names))
}

# Function to prepare bed file (sorting by chromosome)
prepare_bed_file <- function(bed.file) {
    bed.file <- bed.file[order(as.numeric(gsub('chr', '', bed.file$chromosome))), ]
    return(bed.file)
}

# Function to calculate metrics from coverage data
calculate_metrics <- function(count, sample.names, min_corr, min_cov, bed.file) {
    
        # Check if sample.names are valid column names in count
    if (!all(sample.names %in% colnames(count))) {
        invalid_samples <- sample.names[!sample.names %in% colnames(count)]
        stop("Invalid sample names: ", paste(invalid_samples, collapse = ", "))
    }
    
    
    ReadDepths <- count[, sample.names]
    Corr <- cor(ReadDepths)
    MaxCorr <- apply(Corr, 1, function(x) max(x[x != 1]))
    SampleMedian <- apply(ReadDepths, 2, median)
    ExonMedian <- apply(ReadDepths, 1, median)
    
    # Initialize result vectors
    Sample <- Exon <- Details <- Types <- Gene <- character()
    
    # Detect low correlation samples
    low_corr_samples <- sample.names[MaxCorr < min_corr]
    if (length(low_corr_samples) > 0) {
        Sample <- c(Sample, low_corr_samples)
        Exon <- c(Exon, rep("All", length(low_corr_samples)))
        Types <- c(Types, rep("Whole sample", length(low_corr_samples)))
        Details <- c(Details, paste("Low correlation:", MaxCorr[MaxCorr < min_corr]))
        Gene <- c(Gene, rep("All", length(low_corr_samples)))
    }
    
    # Detect low median coverage samples
    for (i in seq_along(sample.names)) {
        if (SampleMedian[i] < min_cov) {
            Sample <- c(Sample, sample.names[i])
            Exon <- c(Exon, "All")
            Types <- c(Types, "Whole sample")
            Details <- c(Details, paste("Low median read depth (FPKM):", SampleMedian[i]))
            Gene <- c(Gene, "All")
        }
    }
    
    # Detect low median coverage exons
    low_cov_exons <- which(ExonMedian < min_cov)
    if (length(low_cov_exons) > 0) {
        Exon <- c(Exon, low_cov_exons)
        Sample <- c(Sample, rep("All", length(low_cov_exons)))
        Types <- c(Types, rep("Whole exon", length(low_cov_exons)))
        Details <- c(Details, paste("Low median read depth (FPKM):", ExonMedian[low_cov_exons]))
        Gene <- c(Gene, bed.file[low_cov_exons, 4])
    }
    
    return(data.frame(Sample, Exon, Types, Gene, Details, stringsAsFactors = FALSE))
}

# Function to write metrics to a file
write_metrics <- function(metrics, output_file, bed.file) {
    if ("exon" %in% colnames(bed.file) & any(metrics$Exon != "All")) {
        exons <- bed.file
        failed.calls <- bed.file[as.numeric(metrics$Exon[metrics$Exon != "All"]),]
        
        exonnumber <- sapply(failed.calls$chr, '==', exons$Chr) & 
                      sapply(failed.calls$start, '>=', exons$Start) & 
                      sapply(failed.calls$end, '<=', exons$End) |
                      sapply(failed.calls$chr, '==', exons$Chr) & 
                      sapply(failed.calls$start, '<=', exons$Start) & 
                      sapply(failed.calls$end, '>=', exons$Start) |
                      sapply(failed.calls$chr, '==', exons$Chr) & 
                      sapply(failed.calls$start, '<=', exons$End) & 
                      sapply(failed.calls$end, '>=', exons$End)
        
        exonlist <- which(colSums(exonnumber) != 0)
        a <- apply(exonnumber[, exonlist, drop=F], 2, which)
        Custom_Numbering <- rep("NA", nrow(metrics))
        Custom_Numbering[as.character(metrics$Exon) %in% row.names(failed.calls[exonlist,])] <- exons$Custom[a]
        Metrics <- data.frame(metrics, Custom_Numbering = Custom_Numbering)
        names(Metrics) <- c("Sample", "Exon", "Type", "Gene", "Custom.numbering", "Info")
        Metrics_custom <- Metrics[Metrics$Custom_Numbering != "NA" | Metrics$Types == "Whole sample",]
        write.table(Metrics_custom, file = output_file, quote = F, row.names = F, sep = "\t")    
    } else {
        write.table(metrics, file = output_file, quote = FALSE, row.names = FALSE, sep = "\t")
    }
}


# Main function to orchestrate the process
main <- function(rdata_file, min_corr, min_cov, output_file) {
    print("BEGIN IdentifyFailures script")
    
    # Load data from RData file
    data <- load_data(rdata_file)
    counts <- data$counts
    bed.file <- data$bed.file
    sample.names <- data$sample.names
    
    # Assign column names to each sample
    colnames(counts)[1:length(sample.names)+5] <- sample.names

    # Prepare bed file (sorting by chromosome)
    bed.file <- prepare_bed_file(bed.file)
    
    # Calculate metrics
    metrics <- calculate_metrics(counts, sample.names, min_corr, min_cov, bed.file)
    
    # Write metrics to output file
    write_metrics(metrics, output_file, bed.file)
    warnings()
    print("END IdentifyFailures script")
}

# Run the main function with parsed arguments
main(rdata_file = options$data, min_corr = as.numeric(options$mincorr), min_cov = as.numeric(options$mincov), output_file = options$out)