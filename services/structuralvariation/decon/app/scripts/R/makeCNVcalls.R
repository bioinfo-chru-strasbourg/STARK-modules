##########################################################################
# DECON script          Version: 2
# Description:          R script to call CNVs from coverage data
##########################################################################

################## Context ###############################################
# Original R script from https://github.com/RahmanTeam/DECoN
# DECON is an ExomeDepth wrapper 
########## Note #########################################################
# PROD v2 21/11/2023
# Changelog
#   - Refactored code, removed install system, updated ExomeDepth 1.16
#   - optparse script, added option to use a list of ref bam files for comparison
#   - Detect bed file with a 5th column with exon numbers, removed custom exon options
#   - Mode to remove chrX or Y for calling
########################################################################################################

suppressPackageStartupMessages({
  library(R.utils)
  library(optparse)
  library(ExomeDepth)
})

###### Parsing input options and setting defaults ########
option_list <- list(
  make_option('--rdata', help='Input summary RData file containing coverage data, bed file, and GC content (required)', dest='data'),
  make_option("--chromosome", default="A", help='Perform calling for autosomes or chr XX or chr XY', dest='chromosome'),
  make_option('--transProb', default=0.01, help='Transition probability for the HMM statistical analysis, default=0.01', dest='transProb'),
  make_option("--refbams", default=NULL, help="Text file containing the list of reference bam files for calling (full path) (optional)", dest='refbams'),
  make_option("--samples", default=NULL, help="Text file containing the list of sample bams to analyse", dest='samples'),
  make_option("--tsv", default="./CNVcalls.tsv", help="Output tsv file, default: ./CNVcalls.tsv", dest='tsv'),
  make_option("--outrdata", default="./CNVcalls.Rdata", help="Output Rdata file, default: ./CNVcalls.Rdata", dest='outdata')
)

opt_parser <- OptionParser(option_list=option_list)
options <- parse_args(opt_parser)

# Function to stop execution if a required value is missing
stop_if_missing <- function(val, message) {
    if (is.null(val) || length(val) == 0) {
        stop(message)
    }
}

# Function to prepare bed file (sorting by chromosome)
prepare_bed_file <- function(bed.file) {
    bed.file <- bed.file[order(as.numeric(gsub('chr', '', bed.file$chromosome))), ]
    return(bed.file)
}

multi_strsplit <- function(x, splits, y) {
  for (i in seq_along(splits)) {
    x <- strsplit(x, splits[i], fixed = TRUE)[[1]][y[i]]
  }
  return(x)
}

process_refbams <- function(refbams.file, mode.chrom) {
  refsample.names <- vector()
  if (!is.null(refbams.file)) {
    raw.refbams <- read.csv(refbams.file, header=TRUE, sep="\t")
    refbams <- apply(raw.refbams, 1, toString)
    a <- length(strsplit(refbams[1], "/")[[1]])
    refsample.names <- sapply(refbams, multi_strsplit, c("/", "."), c(a, 1))
    names(refsample.names) <- NULL
    sample.names <- setdiff(sample.names, refsample.names)
    if ("gender" %in% colnames(raw.refbams)) {
      if (mode.chrom == "XX") {
        refbams <- subset(refbams, raw.refbams$gender == 'F')
      } else if (mode.chrom == "XY") {
        refbams <- subset(refbams, raw.refbams$gender == 'M')
      }
    } else if (mode.chrom %in% c("XX", "XY")) {
      stop('ERROR: No gender specified in the reference bam list, calling of chrX is not possible -- Execution halted')
    }
  }
  return(refsample.names)
}

process_samplebams <- function(samples.file) {
  if (!is.null(samples.file)) {
    samples.bams <- apply(read.table(samples.file), 1, toString)
    a <- length(strsplit(samples.bams[1], "/")[[1]])
    sample.names <- sapply(samples.bams, multi_strsplit, c("/", "."), c(a, 1))
    names(sample.names) <- NULL
  } else {
    stop('ERROR: No samples to analyse -- Execution halted')
  }
  return(sample.names)
}

filter_data_by_chromosome <- function(ExomeCount, bed.file, counts, mode.chrom) {
  if (mode.chrom == "A") {
    ExomeCount <- subset(ExomeCount, !chromosome %in% c("X", "Y"))
    bed.file <- subset(bed.file, !chromosome %in% c("chrX", "chrY"))
    counts <- subset(counts, !chromosome %in% c("chrX", "chrY"))
  } else if (mode.chrom %in% c("XX", "XY")) {
    ExomeCount <- subset(ExomeCount, chromosome == "X")
    bed.file <- subset(bed.file, chromosome == "chrX")
    counts <- subset(counts, chromosome == "chrX")
  }
  list(ExomeCount = ExomeCount, bed.file = bed.file, counts = counts)
}

perform_cnv_calling <- function(ExomeCount, sample.names, refsample.names, trans.prob, bed.file) {
  #cnv.calls <- NULL
  #refs <- list()
  #models <- list()
  
  for (i in seq_along(sample.names)) {
    print(paste("Processing sample:", sample.names[i], i, "/", length(sample.names)))
    my.test <- ExomeCount[, sample.names[i]]
    my.ref.samples <- if (length(refsample.names) == 0) sample.names[-i] else refsample.names
    my.reference.set <- as.matrix(ExomeCount[, my.ref.samples])
    
    my.choice <- select.reference.set(
      test.counts = my.test,
      reference.counts = my.reference.set,
      bin.length = (ExomeCount$end - ExomeCount$start) / 1000,
      n.bins.reduced = 10000
    )
    
    my.matrix <- as.matrix(ExomeCount[, my.choice$reference.choice, drop = FALSE])
    my.reference.selected <- rowSums(my.matrix)
    
    all.exons <- new('ExomeDepth', test = my.test, reference = my.reference.selected, formula = 'cbind(test, reference) ~ 1')
    all.exons <- CallCNVs(
      x = all.exons,
      transition.probability = trans.prob,
      chromosome = ExomeCount$chromosome,
      start = ExomeCount$start,
      end = ExomeCount$end,
      name = ExomeCount$gene
    )
    
    if (nrow(all.exons@CNV.calls) > 0) {
      Comparator.name <- paste(my.choice$reference.choice, collapse = ",")
      cnvs <- cbind(sample.names[i], cor(my.test, rowSums(my.matrix)), length(my.choice[[1]]), Comparator.name, all.exons@CNV.calls)
      cnv.calls <- rbind(cnv.calls, cnvs)
    }
    
    refs[[i]] <- my.choice$reference.choice
    models[[i]] <- c(all.exons@expected[1], all.exons@phi[1])
  }
  
  names(refs) <- sample.names
  names(models) <- sample.names

  if(!is.null(cnv.calls)){
    names(cnv.calls)[1] <- "Sample"
    cnv.calls$sample <- paste(cnv.calls$sample)
    names(cnv.calls)[2] <- "Correlation"
    names(cnv.calls)[3] <- "N.comp"
    names(cnv.calls)["start.p"] <- "Start.p"
    names(cnv.calls)["end.p"] <- "End.p"
    names(cnv.calls)["type"] <- "CNV.type"
    names(cnv.calls)["nexons"] <- "N.exons"
    names(cnv.calls)["start"] <- "Start"
    names(cnv.calls)["end"] <- "End"
    names(cnv.calls)["chromosome"] <- "Chromosome"
    names(cnv.calls)["id"] <- "Genomic.ID"
    names(cnv.calls)["reads.expected"] <- "Reads.expected"
    names(cnv.calls)["reads.observed"] <- "Reads.observed"
    names(cnv.calls)["start.b"] <- "Start.b"
    names(cnv.calls)["end.b"] <- "End.b"
 
        if(ncol(bed.file)>=4){
        genes <- apply(cnv.calls, 1, function(x) {paste(unique(bed.file[x['start.p']:x['end.p'], 4]), collapse = ", ")})
        cnv.calls<-cbind(cnv.calls,genes)
        names(cnv.calls)[ncol(cnv.calls)] <- "Gene"
        }
  
  
  }else{
    print('No CNV detected')
    cnv.calls=NULL
  }
  list(cnv.calls = cnv.calls, refs = refs, models = models)
}

calculate_confidence <- function(cnv.calls, bed.file) {
  Confidence <- rep("HIGH", nrow(cnv.calls))
  Confidence[cnv.calls$correlation < 0.985] <- "LOW"
  Confidence[cnv.calls$reads.ratio < 1.25 & cnv.calls$reads.ratio > 0.75] <- "LOW"
  Confidence[cnv.calls$N.comp <= 3] <- "LOW"
  Confidence[cnv.calls$genes == "PMS2"] <- "LOW"
  return(Confidence)
}

# Function to add custom exon numbers
add_custom_exon_numbers <- function(cnv.calls_ids, bed.file, counts) {
  if (colnames(counts)[5] == "exon") {
    Custom.first <- rep(NA, nrow(cnv.calls_ids))
    Custom.last <- rep(NA, nrow(cnv.calls_ids))
    exons <- bed.file
    exonnumber <- sapply(cnv.calls_ids$Gene, '==', exons$Gene)
    
    for (i in 1:nrow(exonnumber)) {
      for (j in 1:ncol(exonnumber)) {
        temp <- cnv.calls_ids$start[j] <= exons$End[i] & cnv.calls_ids$end[j] >= exons$Start[i]
        exonnumber[i, j] <- exonnumber[i, j] & temp
      }
    }
    
    exonlist <- which(colSums(exonnumber) != 0)
    
    if (length(exonlist) > 0) {
      a <- list(length = length(exonlist))
      for (i in 1:length(exonlist)) {
        a[[i]] <- which(exonnumber[, exonlist[i]])
      }
      
      # identifies the first and last Custom exon in the deletion/duplication.
      first_exon <- unlist(lapply(a, function(a, b) min(b[a,]$Custom), exons))
      last_exon <- unlist(lapply(a, function(a, b) max(b[a,]$Custom), exons))
      Custom.first[exonlist] <- first_exon
      Custom.last[exonlist] <- last_exon
    }
    
    cnv.calls_ids <- cbind(cnv.calls_ids, Custom.first, Custom.last)
  }
  
  return(cnv.calls_ids)
}

# Replaces single calls involving multiple genes with multiple calls with a single call ID/gene
split_multi_gene_calls <- function(cnv.calls, bed.file, counts) {
  
  # Add a new column named "CNV.ID" with increasing numbers from 1 to the number of rows in cnv.calls
  cnv.calls_ids <- cbind(CNV.ID = 1:nrow(cnv.calls), cnv.calls)
  
  trim <- function (x) gsub("^\\s+|\\s+$", "", x)
  
  for (i in 1:nrow(cnv.calls_ids)) {
    if (length(grep(',', cnv.calls_ids$Gene[i])) > 0) {
      temp <- cnv.calls_ids[i, ]
      genes <- unlist(strsplit(cnv.calls_ids$Gene[i], split = ','))
      temp <- temp[rep(seq_len(nrow(temp)), each = length(genes)), ]
      temp$Gene <- genes
      
      gene.index <- which(bed.file[, 4] %in% genes)
      whole.index <- which(bed.file[, 4] == genes[1])
      
      for (j in 1:length(genes)) {
        gene.index <- which(bed.file[, 4] == genes[j])
        overlap <- gene.index[gene.index %in% whole.index]
        temp[j, ]$start.p <- min(overlap)
        temp[j, ]$end.p <- max(overlap)
      }
      
      if (i == 1) {
        cnv.calls_ids <- rbind(temp, cnv.calls_ids[(i + 1):nrow(cnv.calls_ids), ])
      } else if (i == nrow(cnv.calls_ids)) {
        cnv.calls_ids <- rbind(cnv.calls_ids[1:(i - 1), ], temp)
      } else {
        cnv.calls_ids <- rbind(cnv.calls_ids[1:(i - 1), ], temp, cnv.calls_ids[(i + 1):nrow(cnv.calls_ids), ])
      }
    }
  }
  
  cnv.calls_ids <- add_custom_exon_numbers(cnv.calls_ids, bed.file, counts)
  
  cnv.calls_ids$Gene <- trim(cnv.calls_ids$Gene)
  #Gene.index <- vector()
  genes_unique <- unique(bed.file[, 4])
  
  for (i in 1:length(genes_unique)) {
    Gene.index <- c(Gene.index, 1:sum(bed.file[, 4] == genes_unique[i]))
  }
  
  start.b <- Gene.index[cnv.calls_ids$start.p]
  end.b <- Gene.index[cnv.calls_ids$end.p]
  cnv.calls_ids <- cbind(cnv.calls_ids, start.b, end.b)
  
  return(cnv.calls_ids)
}



save_results <- function(cnv.calls, cnv.calls_ids, ExomeCount, output, sample.names, bams, output.rdata,  bed.file, counts, refs, models, fasta) {
  if (!is.null(cnv.calls_ids)) {
   
    cnv.calls_ids$sample <- as.character(cnv.calls_ids$sample)
    cnv.calls_ids <- cbind(cnv.calls_ids, calculate_confidence(cnv.calls_ids, bed.file))
    colnames(cnv.calls_ids)[ncol(cnv.calls_ids)] <- "Confidence"
    
    # Handle exon numbers if present
    if (colnames(counts)[5] == "exon") {
      cnv.calls_ids <- cnv.calls_ids[, c(1:4, 14, 15, 5:13, 16)]
    }
    
   write.table(cnv.calls_ids, file = output, sep = "\t", row.names = FALSE, quote = FALSE)
  }
  
  save(ExomeCount,bed.file,counts,sample.names,bams,cnv.calls_ids,cnv.calls, refs, models, fasta, file=output.rdata)
}

main <- function(data_file, modechrom, samples, p_value, output_file, rdata_output = NULL, refbams_file = NULL) {
  print("BEGIN makeCNVCalls script")
    
  stop_if_missing(data_file, "ERROR: no Rdata summary file provided -- Execution halted")
  load(data_file)

  stop_if_missing(bed.file, "ERROR: bed.file object not found in the loaded RData -- Execution halted")
  stop_if_missing(counts, "ERROR: counts object not found in the loaded RData -- Execution halted")
  stop_if_missing(sample.names, "ERROR: sample.names object not found in the loaded RData -- Execution halted")

  if (!file.exists(dirname(output_file))) {
    dir.create(dirname(output_file))
  }

  bed.file <- prepare_bed_file(bed.file)
  ExomeCount <- as.data.frame(counts)
  colnames(ExomeCount)[1:length(sample.names) + 5] <- sample.names
  ExomeCount$chromosome <- gsub('chr', '', as.character(ExomeCount$chromosome))
  
  filtered.data <- filter_data_by_chromosome(ExomeCount, bed.file, counts, modechrom)
  ExomeCount <- filtered.data$ExomeCount
  bed.file <- filtered.data$bed.file
  counts <- filtered.data$counts
  
  refsample.names <- process_refbams(refbams_file, modechrom)
  sample.names <- process_samplebams(samples)
  
  result <- perform_cnv_calling(ExomeCount, sample.names, refsample.names, p_value, bed.file)
  cnv.calls <- result$cnv.calls
  refs <- result$refs
  models <- result$models
  
  # Split multi-gene calls
  cnv.calls_ids <- split_multi_gene_calls(cnv.calls, bed.file, counts)
  
  save_results(cnv.calls, cnv.calls_ids, ExomeCount, output_file, sample.names, bams, rdata_output, bed.file, counts, refs, models, fasta)

  warnings()
  print("END makeCNVCalls script")
}

# Parse command line arguments and call the main function
main(data_file = options$data, modechrom = options$chromosome, samples = options$samples, p_value = as.numeric(options$transProb), output_file = options$tsv, rdata_output = options$outdata, refbams_file = options$refbams)