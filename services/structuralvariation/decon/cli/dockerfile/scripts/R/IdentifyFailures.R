##########################################################################
# DECON script          Version: 1
# Description:          R script to generate metrics from coverage datas
##########################################################################

################## Context ##################
# Original R script from https://github.com/RahmanTeam/DECoN
# DECON is an ExomeDepth wrapper 
########## Note ########################################################################################
# PROD v1 21/11/2023
# Changelog
#   - refactor code, remove install system, optparse script, update for ExomeDepth 1.16
########################################################################################################

library(R.utils)
library(optparse)
library(ExomeDepth)

print("BEGIN IdentifyFailures script")

option_list <- list(
    make_option("--rdata", help="Input summary RData file (required)", dest='data'),
    make_option("--mincorr", help='Minimum correlation to consider, default=0.98', default=0.98, dest='mincorr'),
    make_option("--mincov", help='Minimum coverage to consider, default=100', default=100, dest='mincov'),
    make_option("--tsv", default='./Metrics.tsv', help='Output metrics tsv file, default=./Metrics.tsv', dest='out')
)

opt <- parse_args(OptionParser(option_list=option_list))

count_data <- opt$data
corr_thresh <- as.numeric(opt$mincorr)
cov_thresh <- as.numeric(opt$mincov)
output <- opt$out

stop_if_missing <- function(val, message) {
    if (is.null(val) || length(val) == 0) {
        stop(message)
    }
}

stop_if_missing(count_data, "ERROR: no Rdata summary file provided -- Execution halted")

load(count_data)
ExomeCount <- as(counts, 'data.frame')
ExomeCount$chromosome <- gsub(as.character(ExomeCount$chromosome), pattern = 'chr', replacement = '')

# assigns the sample names to each column
colnames(ExomeCount)[1:length(sample.names)+5]=sample.names

Sample <- vector()
Exon <- vector()
Details <- vector()
Types <- vector()
Gene <- vector()

temp<-gsub('chr','',bed.file[,1])
temp1<-order(as.numeric(temp))
bed.file=bed.file[temp1,]

ReadDepths <- ExomeCount[, sample.names]
Corr <- cor(ReadDepths)
MaxCorr <- apply(Corr, 1, function(x) max(x[x != 1]))

for (i in 1:length(MaxCorr)) {
    if (MaxCorr[i] < corr_thresh) {
        Sample <- c(Sample, sample.names[i])
        Exon <- c(Exon, "All")
        Types <- c(Types, "Whole sample")
        Details <- c(Details, paste("Low correlation: ", MaxCorr[i], sep = ""))
        Gene <- c(Gene, "All")
    }
}

SampleMedian <- apply(ReadDepths, 2, median)

for (i in 1:length(SampleMedian)) {
    if (SampleMedian[i] < cov_thresh) {
        if (sample.names[i] %in% Sample) {
            k <- which(Sample == sample.names[i])
            Details[k] <- paste(Details[k], ", Low median read depth (FPKM): ", SampleMedian[i], sep = "")
        } else {
            Sample <- c(Sample, sample.names[i])
            Exon <- c(Exon, "All")
            Types <- c(Types, "Whole sample")
            Details <- c(Details, paste("Low median read depth (FPKM): ", SampleMedian[i], sep = ""))
            Gene <- c(Gene, "All")
        }
    }
}

ExonMedian <- apply(ReadDepths, 1, median)

for (i in 1:length(ExonMedian)) {
    if (ExonMedian[i] < cov_thresh) {
        Exon <- c(Exon, i)
        Sample <- c(Sample, "All")
        Types <- c(Types, "Whole exon")
        Details <- c(Details, paste("Low median read depth (FPKM): ", ExonMedian[i], sep = ""))
        Gene <- c(Gene, paste(bed.file[i, 4]))
    }
}

if("exon" %in% colnames(bed.file) & any(Exon!="All")){
	exons <- bed.file
	failed.calls=bed.file[Exon[Exon!="All"],]
	exonnumber<-sapply(failed.calls$chr, '==',exons$Chr) & sapply(failed.calls$start, '>=',exons$Start) & sapply(failed.calls$end, '<=',exons$End) | sapply(failed.calls$chr, '==',exons$Chr) & sapply(failed.calls$start, '<=',exons$Start) & sapply(failed.calls$end, '>=',exons$Start) | sapply(failed.calls$chr, '==',exons$Chr) & sapply(failed.calls$start, '<=',exons$End) & sapply(failed.calls$end, '>=',exons$End)
	exonlist<-which(colSums(exonnumber)!=0)
	a<-apply(exonnumber[,exonlist,drop=F],2,which)
	Custom_Numbering=rep("NA",length(Exon))
	Custom_Numbering[paste(Exon)%in%row.names(failed.calls[exonlist,])] = exons$Custom[a]
	Metrics<-data.frame(Sample,Exon,Types,Gene,Custom_Numbering,Details)
	names(Metrics)=c("Sample","Exon","Type","Gene","Custom.numbering","Info")
	Metrics_custom=Metrics[Metrics$Custom_Numbering!="NA" | Metrics$Types=="Whole sample",]
	write.table(Metrics_custom,file=output,quote=F,row.names=F,sep="\t")

}else{
	Metrics<-data.frame(Sample,Exon,Types,Gene,Details)
	names(Metrics)=c("Sample","Exon","Type","Gene","Info")
	write.table(Metrics,file=output,quote=F,row.names=F,sep="\t")
}
warnings()
print("END IdentifyFailures script")