###############################################################################
# DECON script          Version: 1
# Description:          R script to extract coverage datas from a list of bams
###############################################################################

################## Context ###############################################
# Original R script from https://github.com/RahmanTeam/DECoN
# DECON is an ExomeDepth wrapper 
########## Note ########################################################################################
# PROD v1 21/11/2023
# Changelog
#   - refactor code, remove install system, update for ExomeDepth 1.16
#   - parallelisation to speed up the process
#   - detect bed file with a 5 column with exon numbers, remove custom exon options
########################################################################################################


library(R.utils)
library(optparse)
library(ExomeDepth)
library(doParallel)
library(foreach)
library(dplyr)

print("BEGIN ReadInBams script")

option_list <- list(
    make_option("--bams", help="Text file containing a list of bam files to process (required)", dest='bamfiles'),
    make_option("--refbams", help="Text file containing a list of reference bam files to process (optional)", dest='rbams'),
    make_option("--bed", help='Bed file with 4 or 5 columns (chr, start, stop, gene, +/- exon) used to generate coverage datas (required)', dest='bed'),
    make_option("--fasta", help='Reference genome fasta file to use (required)', default=NULL, dest='fasta'),
    make_option("--rdata", default="./ReadInBams.Rdata", help="Output Rdata file, default: ./ReadInBams.Rdata", dest='data'),
    make_option("--maxcores", default=16, help="Maximum cores to use, default: 16", dest='mcore')
)

opt <- parse_args(OptionParser(option_list=option_list))
bam_file <- opt$bamfiles
bedfile <- opt$bed
fasta <- opt$fasta
output <- opt$data
refbams_file <- opt$rbams
maxcores <- opt$mcore

stop_if_missing <- function(val, message) {
    if (is.null(val) || length(val) == 0) {
        stop(message)
    }
}

stop_if_missing(bedfile, "ERROR : bed file must be provided. Execution halted")
stop_if_missing(fasta, "ERROR : no reference fasta files detected. Execution halted")
stop_if_missing(bam_file, "ERROR : bam files must be provided. Execution halted")

bams<-apply(read.table(paste(bam_file)),1,toString)

if (length(refbams_file) > 0) {
    refbams <- read.csv(paste(refbams_file), header=TRUE, sep="\t")
    bams <- append(refbams$bam, bams)
}

# to test
get_sample_names <- function(bam_paths) {
  sapply(strsplit(bam_paths, "/|\\."), function(x) x[length(x) - 1], simplify = TRUE)
}
#sample_names <- get_sample_names(bam)

# function which recursively splits x by an element of 'splits' then extracts the y element of the split vector
multi_strsplit<-function(x,splits,y){
	X<-x
	for(i in 1:length(splits)){X=strsplit(X,splits[i], fixed = TRUE)[[1]][y[i]]}
	return(X)
}

# get the sample names from the path of the bams
a<-length(strsplit(bams[1],"/")[[1]])
sample.names<-sapply(bams,multi_strsplit,c("/","."),c(a,1))
names(sample.names)<-NULL

bed.file <- read.table(paste(bedfile))
if (ncol(bed.file) == 4){
colnames(bed.file)<-c("chromosome","start","end","gene")
}
if (ncol(bed.file) == 5){
colnames(bed.file)<-c("chromosome","start","end","gene","exon_number")
}

nfiles <- length(bams)
message('Parse ', nfiles, ' BAM files')
numCores <- min(detectCores(), maxcores)
cl <- parallel::makeForkCluster(numCores)
doParallel::registerDoParallel(cl)

unfilteredcounts <- foreach(i = 1:nfiles, .combine='cbind') %dopar% {
    bam <- bams[i]
    getBamCounts(bed.frame = bed.file, bam.files = bam, include.chr = FALSE, referenceFasta = fasta)
}

parallel::stopCluster(cl)

filtercount1 <- basename(bams)
if (ncol(bed.file) == 4){
filtercount2 <- c("chromosome", "start", "end", "exon", "GC")
}
if (ncol(bed.file) == 5){
filtercount2 <- c("chromosome", "start", "end", "exon", "exon_number", "GC")
}
filtercount <- append(filtercount2, filtercount1)

counts <- dplyr::select(unfilteredcounts, all_of(filtercount))
names(counts) <- unlist(lapply(strsplit(names(counts), "\\."), "[[", 1))
colnames(counts)[colnames(counts) == "exon"] <- "gene"
colnames(counts)[colnames(counts) == "exon_number"] <- "exon"

if (ncol(bed.file) == 5){
    colnames(bed.file)[colnames(bed.file) == "exon_number"] <- "exon"
}

save(counts, bams, bed.file, sample.names, fasta, file=output)
warnings()
print("END ReadInBams script")