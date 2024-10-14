###############################################################################
# DECON script          Version: 1
# Description:          R script to extract coverage data from a list of BAMs
###############################################################################

################## Context ###############################################
# Original R script from https://github.com/RahmanTeam/DECoN
# DECON is an ExomeDepth wrapper 
########## Note ########################################################################################
# PROD v1 21/11/2023
# Changelog
#   - refactor code, remove install system, update ExomeDepth 1.16
#   - parallelisation to speed up the process
#   - detect bed file with a 5 column with exon numbers, remove custom exon options
########################################################################################################

suppressPackageStartupMessages({
  library(R.utils)
  library(optparse)
  library(ExomeDepth)
  library(doParallel)
  library(foreach)
  library(dplyr)
})

option_list <- list(
    make_option("--bams", help="Text file containing a list of bam files to process (required)", dest='bamfiles'),
    make_option("--refbams", help="Text file containing a list of reference bam files to process (optional)", dest='rbams'),
    make_option("--bed", help='Bed file with 4 or 5 columns (chr, start, stop, gene, +/- exon) used to generate coverage data (required)', dest='bed'),
    make_option("--fasta", help='Reference genome fasta file to use (required)', default=NULL, dest='fasta'),
    make_option("--rdata", default="./ReadInBams.Rdata", help="Output Rdata file, default: ./ReadInBams.Rdata", dest='data'),
    make_option("--maxcores", default=16, help="Maximum cores to use, default: 16", dest='mcore')
)

opt <- parse_args(OptionParser(option_list=option_list))

stop_if_missing <- function(val, message) {
    if (is.null(val) || length(val) == 0) {
        stop(message)
    }
}

process_bams <- function(bamfiles, rbams, bed, fasta, output, maxcores = 16) {
    # Read BAM files
    bams <- read_bam_files(bamfiles)
    
    # Optionally append reference BAMs
    if (!is.null(rbams) && file.exists(rbams)) {
        refbams <- read_reference_bams(rbams)
        bams <- append(refbams, bams)
    }
    
    # Get sample names
    sample.names <- get_sample_names(bams)
    
    # Read BED file without headers
    bed.file <- read_bed_file(bed)
    
    nfiles <- length(bams)
    message(paste('Parse', nfiles, 'BAM files'))
    
    # Set up parallel processing
    numCores <- min(detectCores(), maxcores)
    cl <- parallel::makeForkCluster(numCores)
    doParallel::registerDoParallel(cl)
    
    # Debugging: Check structure of bams and bed.file
    message("BAM files: ", paste(bams, collapse=", "))
    message("BED file structure: ", str(bed.file))
    head(bed.file)

    # Process each BAM file in parallel and get counts
    unfilteredcounts <- foreach(i = 1:nfiles, .combine = 'cbind') %dopar% {
        bam <- bams[[i]]
        message("Processing BAM: ", bam)
        
        getBamCounts(bed.frame = bed.file, bam.files = bam, include.chr = FALSE, referenceFasta = fasta)
    }
    
    parallel::stopCluster(cl)

    # Dynamically adjust the filtercount based on the number of columns in the BED file
    filtercount <- c("chromosome", "start", "end")
    if (ncol(bed.file) == 5) {
        filtercount <- append(filtercount, "exon_number")
    }

    # Extract the unique chromosome, start, end, gene (and exon if applicable) from the BED file
    bed_info <- bed.file[, c("V1", "V2", "V3", "V4")]
    if (ncol(bed.file) == 5) {
        bed_info <- cbind(bed_info, exon = bed.file$V5)
    }

    # Remove the metadata columns from the unfiltered counts
    counts_data <- unfilteredcounts[, !colnames(unfilteredcounts) %in% filtercount]

    # Combine bed_info (unique) with counts_data
    counts <- cbind(bed_info[!duplicated(bed_info), ], counts_data)
    
    # Rename the coverage columns to sample names
    colnames(counts)[-(1:(ncol(bed_info)))] <- sample.names
    
    # Append gene and exon_number (if applicable) at the start of counts
    if (ncol(bed.file) == 5) {
        final_counts <- cbind(gene = bed_info$V4, exon_number = bed_info$exon, counts)
    } else {
        final_counts <- cbind(gene = bed_info$V4, counts)
    }
    
    # Save the output
    save(final_counts, bams, bed.file, sample.names, fasta, file = output)
    
    warnings()
}





read_bam_files <- function(bam_file) {
    apply(read.table(bam_file), 1, toString)
}

read_reference_bams <- function(refbams_file) {
    refbams <- read.csv(refbams_file, header=TRUE, sep="\t")
    refbams$bam
}

get_sample_names <- function(bam_paths) {
    # Extract sample names for each BAM path
    sample_names <- sapply(bam_paths, function(bam_path) {
        # Extract the basename of the file
        base_name <- basename(bam_path)
        
        # Split the basename by '.' and take the first part
        sample_name <- strsplit(base_name, "\\.")[[1]][1]
        
        return(sample_name)
    })
    
    return(sample_names)
}


read_bed_file <- function(bedfile) {
    bed.file <- read.table(bedfile, sep="\t", header=FALSE, stringsAsFactors=FALSE)
    # Debugging output to check the structure
    message("Number of columns in bed file: ", ncol(bed.file))
    print(head(bed.file))
    
    # Remove any completely empty columns, just in case
    bed.file <- bed.file[, colSums(is.na(bed.file) | bed.file == "") != nrow(bed.file)]
    
    # Debugging output after cleanup
    message("Number of columns after cleanup: ", ncol(bed.file))
    
    # Set column names based on the number of columns
    if (ncol(bed.file) == 5) {
        colnames(bed.file) <- c("chromosome", "start", "end", "gene", "exon_number")
    } else if (ncol(bed.file) == 4) {
        colnames(bed.file) <- c("chromosome", "start", "end", "gene")
    } else {
        stop("Unexpected number of columns in BED file")
    }
    return(bed.file)  # Make sure to return the data frame
}

main <- function(bamfiles, rbams, bed, fasta, data, mcore) {
    print("BEGIN ReadInBams script")
    
    # Validate input
    stop_if_missing(bamfiles, "ERROR: BAM files must be provided -- Execution halted")
    stop_if_missing(bed, "ERROR: BED file must be provided -- Execution halted")
    stop_if_missing(fasta, "ERROR: No reference FASTA file detected -- Execution halted")
    
    process_bams(bamfiles, rbams, bed, fasta, data, mcore)
    warnings()
    print("END ReadInBams script")
}

# Run the main function with parsed arguments
main(opt$bamfiles, opt$rbams, opt$bed, opt$fasta, opt$data, opt$mcore)