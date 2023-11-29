# Changelog 21/11/2023 TL
# Refactor some code, remove packrat ref & BRCA specific code
# add parallel to speed up the process

# TODO directory save and prefix should be different

print("BEGIN ReadInBams.R")

library(R.utils)
library(optparse)
library(ExomeDepth)
library(doParallel)
library(foreach)
library(dplyr)


######Parsing input options and setting defaults########
option_list<-list(
    make_option("--bams",help="txt file containing a list of bam files or a directory containing all the bam files",dest='bams'),
    make_option("--refbams",help="txt file containing a list of reference bam files, full path",dest='rbams'),
    make_option("--bed",help='bed file to be used',dest='bed'),
    make_option("--fasta",help='ref genome fasta file to be used, optional',default=NULL,dest='fasta'),
    make_option("--out",default="DECoN",help="output folder, default: DECoN",dest='out')
)
opt<-parse_args(OptionParser(option_list=option_list))
# location of bam files; can be a directory containing only bam files to be processed or the name of a file containing a list of bam files to be processed.
bam_file=opt$bams
bedfile=opt$bed
fasta=opt$fasta
output=opt$out
refbams_file=opt$rbams

if(bam_file=="NULL"){bam_file=NULL}
if(is.null(bam_file)){
print("ERROR bam files must be provided. Execution halted")
quit()
}

if(bedfile=="NULL"){bedfile=NULL}
if(is.null(bedfile)){
print("ERROR bed file must be provided. Execution halted")
quit()
}

# Read list of bams to process
bams<-apply(read.table(paste(bam_file)),1,toString)
if(length(refbam_file)==0){
bams<-apply(read.table(paste(refbams_file)),1,toString)
}

# Check if bams is not empty
if(length(bams)==0){
print("ERROR NO BAM FILES DETECTED")
quit()
}

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

# reads in the bedfile and gives each column a name - expects 4 columns: chr, start, stop, gene
bed.file<-read.table(paste(bedfile))
colnames(bed.file)<-c("chromosone","start","end","gene")

# parallelisation
# add some limitation to the numCores (12 to start)
nfiles <- length(bams)
message('Parse ', nfiles, ' BAM files')
numCores <- detectCores()
if (numCores > 12) {
    numCores <- 12
}
cl <- parallel::makeForkCluster(numCores)
doParallel::registerDoParallel(cl)

# reads in coverage info from each bam file in 'bams'; expects chromosomes to be given as numbers, eg 1, 2 etc, not chr1, chr2 etc.
unfilteredcounts <- foreach(i=1:nfiles, .combine='cbind') %dopar% {
	bam <- bams[i]
	unfilteredcounts <- getBamCounts(bed.frame = bed.file, bam.files = bam, include.chr = FALSE, referenceFasta = fasta)
}
parallel::stopCluster(cl)

# counts should be chromosome, start, end, exon, GC, sample1, sample2, etc., so we need to clean up the df
filtrecount1 <-basename(apply(read.table(paste(bam_file)),1,toString))
filtrecount2 <- c("chromosome", "start", "end", "exon", "GC")
filtrecount <-append(filtrecount2, filtrecount1)
counts<-dplyr::select(unfilteredcounts, all_of(filtrecount))
names(counts) <- unlist(lapply(strsplit(names(counts), "\\."), "[[", 1))
save(counts,bams,bed.file,sample.names,fasta,file=paste(output,"DECoN.ReadInBams.RData",sep=""))

print("END ReadInBams.R")