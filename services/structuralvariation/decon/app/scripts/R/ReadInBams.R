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
# usage sorted_df <- sort_chromosome_df(df)
sort_chromosome_df <- function(df) {
  # Check if 'chr' prefix is present in chromosome column
  chr_prefix <- grepl("^chr", df$chromosome[1])
  
  # Strip 'chr' prefix if present for numeric conversion
  df$chromosome <- gsub("chr", "", df$chromosome)
  
  # Convert chromosome column to numeric, handling X and Y as 23 and 24
  df$chromosome_numeric <- as.numeric(ifelse(df$chromosome == "X", 23,
                                             ifelse(df$chromosome == "Y", 24, df$chromosome)))
  
  # Sort by chromosome_numeric and start columns
  df <- df[order(df$chromosome_numeric, df$start), ]
  
  # Convert 23 and 24 back to "X" and "Y"
  df$chromosome <- ifelse(df$chromosome_numeric == 23, "X",
                          ifelse(df$chromosome_numeric == 24, "Y", df$chromosome))
  
  # Drop the helper column
  df$chromosome_numeric <- NULL
  
  # Add 'chr' prefix back if it was originally present
  if (chr_prefix) {
    df$chromosome <- paste0("chr", df$chromosome)
  }
  
  # Reset row indices
  row.names(df) <- NULL
  
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
  bed.file <- sort_chromosome_df(bed.file)

  message(paste('Parse', length(bams), 'BAM files'))

  # Run getBamCounts on all BAMs at once
  counts <- getBamCounts(
    bed.frame = bed.file,
    bam.files = bams,
    include.chr = FALSE,
    referenceFasta = fasta
  )

  # Clean column/sample names
  colnames(counts) <- sapply(basename(colnames(counts)), tools::file_path_sans_ext)

  # ---- Missing sample check ----
  non_sample_cols <- c("chromosome", "start", "end", "GC", "exon", "gene", "exon_number")
  present_samples <- setdiff(colnames(counts), non_sample_cols)
  missing_samples <- setdiff(sample.names, present_samples)

  if (length(missing_samples) > 0) {
    warning(sprintf("[WARN] %d sample(s) missing from counts table: %s",
                    length(missing_samples),
                    paste(missing_samples, collapse = ", ")))
  }

  # Rename "exon" to "gene" if needed
  colnames(counts)[colnames(counts) == "exon"] <- "gene"

  # Optionally add exon_number
  if ("exon_number" %in% colnames(bed.file)) {
    counts <- dplyr::bind_cols(counts, bed.file["exon_number"])
    counts <- counts %>%
      dplyr::relocate(exon_number, .after = gene)
    colnames(bed.file)[colnames(bed.file) == "exon_number"] <- "exon"
  }

  # Remove duplicate rows
  counts <- counts[!duplicated(counts[, c("chromosome", "start", "end", "gene")]), ]

  # Save results
  save(counts, bams, bed.file, sample.names, fasta, file = output)
  warnings()
}


read_bam_files <- function(bam_file) {
    # Log message before reading the file
    message("Starting to read BAM file: ", bam_file)
    
    # Check if file exists before reading
    if (!file.exists(bam_file)) {
        message("Error: File does not exist: ", bam_file)
        return(NULL)  # Return NULL if file doesn't exist
    }
    
    # Log message for successful file loading
    message("File exists, loading the file: ", bam_file)
    
    # Read the file and convert each row to a string
    result <- apply(read.table(bam_file), 1, toString)
    
    # Log message after processing the file
    message("Finished processing the BAM file: ", bam_file)
    
    # Return the result
    return(result)
}

read_reference_bams <- function(refbams_file) {
    refbams <- read.csv(refbams_file, header=TRUE, sep="\t")
    refbams$bam
}

get_sample_names <- function(bam_paths) {
    # Log the start of the function
    message("Starting to extract sample names for BAM files.")
    
    # Initialize a vector to store sample names
    sample_names <- sapply(bam_paths, function(bam_path) {
        # Log the current BAM file being processed
        message("Processing BAM file: ", bam_path)
        
        # Extract the basename of the file
        base_name <- basename(bam_path)
        
        # Split the basename by '.' and take the first part
        sample_name <- strsplit(base_name, "\\.")[[1]][1]
        
        # Log the extracted sample name
        message("Extracted sample name: ", sample_name)
        
        return(sample_name)
    })
    
    # Log the completion of the function
    message("Finished extracting sample names.")
    
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