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

# Function to sort a file by chromosome (any file with a chromosome column)
# usage sorted_df <- chr_sort_df(df, "chromosome")
chr_sort_df <- function(df, col_name, add_prefix = TRUE) {
    # Remove "chr" prefix temporarily
    chromosomes <- gsub("chr", "", df[[col_name]])

    # Separate numeric and non-numeric chromosomes
    numeric_chromosomes <- suppressWarnings(as.numeric(chromosomes))
    non_numeric_chromosomes <- chromosomes[is.na(numeric_chromosomes)]
    numeric_chromosomes <- numeric_chromosomes[!is.na(numeric_chromosomes)]
    
    # Sort numeric chromosomes and non-numeric chromosomes separately
    sorted_numeric <- sort(unique(numeric_chromosomes), na.last = TRUE)
    sorted_non_numeric <- sort(unique(non_numeric_chromosomes))
    
    # Combine sorted numeric and non-numeric chromosomes
    sorted_chromosomes <- c(sorted_numeric, sorted_non_numeric)
    
    # Reapply the "chr" prefix if needed
    if (add_prefix) {
        sorted_chromosomes <- paste0("chr", sorted_chromosomes)
    }

    # Keep only levels that are actually present in the data
    actual_levels <- intersect(sorted_chromosomes, unique(df[[col_name]]))
    df[[col_name]] <- factor(df[[col_name]], levels = actual_levels)
    
    # Order by the factor levels
    df <- df[order(df[[col_name]]), ]
    rownames(df) <- NULL  # Reset row indices if needed
    
    return(df)
}






process_bams <- function(bamfiles, rbams, bed, fasta, output, maxcores = 16) {
    bams <- read_bam_files(bamfiles)
    
    if (!is.null(rbams) && file.exists(rbams)) {
        refbams <- read_reference_bams(rbams)
        bams <- append(refbams, bams)
    }
    
    sample.names <- get_sample_names(bams)
    bed.file <- read_bed_file(bed)
    bed.file <- chr_sort_df(bed.file, "chromosome", add_prefix = TRUE)
    head(bed.file)
    nfiles <- length(bams)
    message(paste('Parse', nfiles, 'BAM files'))
    numCores <- min(detectCores(), maxcores)
    cl <- parallel::makeForkCluster(numCores)
    doParallel::registerDoParallel(cl)
    
    unfilteredcounts <- foreach(i = 1:nfiles, .combine='cbind') %dopar% {
        bam <- bams[i]
        getBamCounts(bed.frame = bed.file, bam.files = bam, include.chr = FALSE, referenceFasta = fasta)
    }
    
    parallel::stopCluster(cl)
    
    filtercount1 <- basename(bams)
    filtercount2 <- c("chromosome", "start", "end", "exon", "GC")
    filtercount <- append(filtercount2, filtercount1)
    
    counts <- dplyr::select(unfilteredcounts, all_of(filtercount))
    names(counts) <- unlist(lapply(strsplit(names(counts), "\\."), "[[", 1))
    colnames(counts)[colnames(counts) == "exon"] <- "gene"
   
    if (ncol(bed.file) == 5) {
        counts <- dplyr::bind_cols(counts, bed.file["exon_number"])
        counts <- counts %>%
            dplyr::relocate(exon_number, .after = gene)
        colnames(bed.file)[colnames(bed.file) == "exon_number"] <- "exon"
    }
    
    save(counts, bams, bed.file, sample.names, fasta, file=output)
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
    bed.file <- read.table(bedfile)
    if (ncol(bed.file) == 5) {
        colnames(bed.file) <- c("chromosome", "start", "end", "gene","exon_number")
    } else if (ncol(bed.file) == 4) {
        colnames(bed.file) <- c("chromosome", "start", "end", "gene")
    }
    bed.file
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