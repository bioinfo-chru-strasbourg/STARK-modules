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

print("BEGIN IdentifyFailures script")

suppressMessages(library(R.utils))
suppressMessages(library(optparse))
suppressMessages(library(ExomeDepth))

###### Parsing input options and setting defaults ########
option_list<-list(
    make_option("--RData",help="Input summary RData file (required)",dest='Rdata'),
    make_option("--mincorr",help='Minimum correlation to consider, default=0.98',default=.98,dest='mincorr'),
    make_option("--mincov",help='Minimum coverage to consider, default=100',default=100,dest='mincov'),
    make_option("--out",default='./metrics',help='Output directory, default=./metrics',dest='out'),
	make_option("--filename",default="ReadInBams.Rdata",help="File name of the .Rdata, default: ReadInBams.Rdata",dest='filename')
)
opt<-parse_args(OptionParser(option_list=option_list))
count_data=opt$Rdata
if(count_data=="NULL"){count_data=NULL}
if(is.null(count_data)){
print("ERROR: no Rdata summary file provided -- Execution halted")
quit()
}
corr_thresh=as.numeric(opt$mincorr)
if(length(corr_thresh)==0){corr_thresh=0.98}
cov_thresh=as.numeric(opt$mincov)
if(length(cov_thresh)==0){cov_thresh=100}
output=opt$out
file_name=opt$filename
if(!file.exists(output)){dir.create(output)}


# R workspace with the coverage data, bedfile, GC content from ref genome saved in it.
load(count_data)
# converts counts, a ranged data object, to a data frame
ExomeCount<-as(counts, 'data.frame')
# remove any chr letters, and coerce to a string.
ExomeCount$chromosome <- gsub(as.character(ExomeCount$chromosome),pattern = 'chr',replacement = '') 
# assigns the sample names to each column
colnames(ExomeCount)[1:length(sample.names)+5]=sample.names

Sample<-vector()
Exon<-vector()
Details<-vector()
Types<-vector()
Gene<-vector()

####### Make sure bed file is in chromosome order ################
temp<-gsub('chr','',bed.file[,1])
temp1<-order(as.numeric(temp))
bed.file=bed.file[temp1,]


# extracts just the read depths
ReadDepths<-ExomeCount[,sample.names]

# calculates correlation matrix
Corr<-cor(ReadDepths)
# finds the maximum correlation for each sample
MaxCorr<-apply(Corr,1,function(x)max(x[x!=1]))

# tests correlation for each sample; if below corr_thresh, adds to list of fails
for(i in 1:length(MaxCorr)){
	if(MaxCorr[i]<corr_thresh){
			Sample<-c(Sample,sample.names[i])
			Exon<-c(Exon,"All")
			Types<-c(Types,"Whole sample")
			Details<-c(Details,paste("Low correlation: ", MaxCorr[i],sep=""))
			Gene<-c(Gene,"All")
	}
}

# calculates median coverage per sample
SampleMedian<-apply(ReadDepths,2,median)

# tests median coverage for each sample; if below cov_thresh, adds to list of fails
for(i in 1:length(SampleMedian)){
	if(SampleMedian[i]<cov_thresh){
		if(sample.names[i]%in%Sample){
			k=which(Sample==sample.names[i])
			Details[k] = paste(Details[k],", Low median read depth (FPKM): ", SampleMedian[i],sep="")
		}else{
			Sample<-c(Sample,sample.names[i])
			Exon<-c(Exon,"All")
			Types<-c(Types,"Whole sample")
			Details<-c(Details,paste("Low median read depth (FPKM): ", SampleMedian[i],sep=""))
			Gene<-c(Gene,"All")
		}
	}
}

# calculates median coverage per exon
ExonMedian<-apply(ReadDepths,1,median)

# tests median coverage for each exon; if below cov_thresh, adds to list of fails
for(i in 1:length(ExonMedian)){
	if(ExonMedian[i]<cov_thresh){
		Exon<-c(Exon,i)
		Sample<-c(Sample,"All")
		Types<-c(Types,"Whole exon")
		Details<-c(Details,paste("Low median read depth (FPKM): ",ExonMedian[i],sep=""))
		Gene<-c(Gene,paste(bed.file[i,4]))
	}
}

if((names(bed.file)[5]=="exon") & any(Exon!="All")){
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
	write.table(Metrics_custom,file=paste(output,"Custom.",file_name,sep=""),quote=F,row.names=F,sep="\t")

}else{
	Metrics<-data.frame(Sample,Exon,Types,Gene,Details)
	names(Metrics)=c("Sample","Exon","Type","Gene","Info")
	write.table(Metrics,file=paste(output,file_name,sep=""),quote=F,row.names=F,sep="\t")
}
warnings()
print("END IdentifyFailures script")