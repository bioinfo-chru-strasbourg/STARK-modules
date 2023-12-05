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

# TODO : add a function to validate the bed ?

print("BEGIN ReadInBams script")

suppressMessages(library(R.utils))
suppressMessages(library(optparse))
suppressMessages(library(ExomeDepth))
suppressMessages(library(doParallel))
suppressMessages(library(foreach))
suppressMessages(library(dplyr))


###### Parsing input options and setting defaults ########
option_list<-list(
    make_option("--bams",help="Text file containing a list of bam files to process (required)",dest='bamfiles'),
    make_option("--refbams",help="Text file containing a list of reference bam files to process (optional)",dest='rbams'),
    make_option("--bed",help='Bed file with 4 or 5 columns (chr, sart, stop, gene, +/- exon) used to generate coverage datas (required)',dest='bed'),
    make_option("--fasta",help='Reference genome fasta file to use (required)',default=NULL,dest='fasta'),
    make_option("--out",default="./coverage",help="Output folder for the ReadInBams.Rdata, default: ./coverage",dest='out'),
    make_option("--filename",default="ReadInBams.Rdata",help="File name of the .Rdata, default: ReadInBams.Rdata",dest='filename')
)
opt<-parse_args(OptionParser(option_list=option_list))
# location of bam files; can be a directory containing only bam files to be processed or the name of a file containing a list of bam files to be processed.
bam_file=opt$bamfiles
bedfile=opt$bed
fasta=opt$fasta
output=opt$out
file_name=opt$filename
refbams_file=opt$rbams

if(!file.exists(output)){dir.create(output)}

if(bedfile=="NULL"){bedfile=NULL}
if(is.null(bedfile)){
print("ERROR bed file must be provided. Execution halted")
quit()
}

if(length(fasta)==0){
print("ERROR No reference fasta files detected. Execution halted")
quit()
}

# Check if bams is not empty then read list of bams to process
if(bam_file=="NULL"){bam_file=NULL}
if(is.null(bam_file)){
print("ERROR bam files must be provided. Execution halted")
quit()
}

if(length(bam_file)==0){
print("ERROR No bam files found. Execution halted")
quit()
}else{
bams<-apply(read.table(paste(bam_file)),1,toString)
}

# Check if refbams is not empty then append refbams to bams
if(length(refbams_file)>0){
refbams<- read.csv(paste(refbams_file), header=TRUE, sep="\t")
bams <- append(refbams$bam,bams)
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

# reads in the bedfile and gives each column a name - expects 4 columns: chr, start, stop, gene, +/- exon
bed.file<-read.table(paste(bedfile))
if (ncol(bed.file) == 4){
colnames(bed.file)<-c("chromosome","start","end","gene")
}
if (ncol(bed.file) == 5){
colnames(bed.file)<-c("chromosome","start","end","gene","exon_number")
}

# parallelisation
nfiles <- length(bams)
message('Parse ', nfiles, ' BAM files')
numCores <- detectCores()
# add some limitation to the numCores (12 to start)
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

# counts should be chromosome, start, end, exon, +/- exon_number, GC, sample1, sample2, etc., so we need to clean up the df
filtercount1 <-basename(bams)
if (ncol(bed.file) == 4){
filtercount2 <- c("chromosome", "start", "end", "exon", "GC")
}
if (ncol(bed.file) == 5){
filtercount2 <- c("chromosome", "start", "end", "exon", "exon_number", "GC")
}
filtercount <-append(filtercount2, filtercount1)
counts<-dplyr::select(unfilteredcounts, all_of(filtercount))
names(counts) <- unlist(lapply(strsplit(names(counts), "\\."), "[[", 1))
colnames(counts)[colnames(counts) == "exon"] <- "gene"
colnames(counts)[colnames(counts) == "exon_number"] <- "exon"
if (ncol(bed.file) == 5){
    colnames(bed.file)[colnames(bed.file) == "exon_number"] <- "exon"
}
# Rdata counts table is chromosome, start, end, gene, +/- exon, GC, sample1, sample2 & bed.file table & sample names & path to genome ref fasta file
save(counts,bams,bed.file,sample.names,fasta,file=paste(output,file_name, sep=""))
warnings()
print("END ReadInBams script")