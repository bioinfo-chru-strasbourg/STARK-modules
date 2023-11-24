# Changelog 21/11/2023 TL
# Refactor some code, remove packrat ref & BRCA specific code

print("BEGIN IdentifyFailures script")

library(R.utils)
library(optparse)
library(ExomeDepth)

######Parsing input options and setting defaults########
option_list<-list(
    make_option("--RData",help="Summary RData file",dest='Rdata'),
    make_option("--mincorr",help='minimum correlation, default .98',default=.98,dest='mincorr'),
    make_option("--mincov",help='minimum coverage, default=100',default=100,dest='mincov'),
    make_option("--exons",default=NULL,help="File containing custom exon annotation",dest='exons'),
    make_option("--custom",default=FALSE,help='Use custom reporting, default=FALSE',dest='custom'),
    make_option("--out",default='./',help='output directory, default=./',dest='out')
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
exon_numbers=opt$exons
Custom=as.logical(opt$custom)
if(length(Custom)==0){Custom=FALSE}
Output=opt$out

# R workspace with the coverage data saved in - this workspace already has the bedfile, fasta file saved in it.
load(count_data)
# converts counts, a ranged data object, to a data frame
ExomeCount<-as(counts, 'data.frame')
# remove any chr letters, and coerce to a string.
ExomeCount$chromosome <- gsub(as.character(ExomeCount$chromosome),pattern = 'chr',replacement = '') 
# assigns the sample names to each column
if(names(ExomeCount)[5]=="GC"){ 
	colnames(ExomeCount)[1:length(sample.names)+5]=sample.names
	}else{
	colnames(ExomeCount)[1:length(sample.names)+4]=sample.names
	}

Sample<-vector()
Exon<-vector()
Details<-vector()
Types<-vector()
Gene<-vector()

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

if(!is.null(exon_numbers)&any(Exon!="All")){
	exons<-read.table(exon_numbers,sep="\t",header=T)
	failed.calls=bed.file[Exon[Exon!="All"],]
	exonnumber<-sapply(failed.calls$chr, '==',exons$Chr) & sapply(failed.calls$start, '>=',exons$Start) & sapply(failed.calls$end, '<=',exons$End) | sapply(failed.calls$chr, '==',exons$Chr) & sapply(failed.calls$start, '<=',exons$Start) & sapply(failed.calls$end, '>=',exons$Start) | sapply(failed.calls$chr, '==',exons$Chr) & sapply(failed.calls$start, '<=',exons$End) & sapply(failed.calls$end, '>=',exons$End)
	exonlist<-which(colSums(exonnumber)!=0)
	a<-apply(exonnumber[,exonlist,drop=F],2,which)
	Custom_Numbering=rep("NA",length(Exon))
	Custom_Numbering[paste(Exon)%in%row.names(failed.calls[exonlist,])] = exons$Custom[a]
	Metrics<-data.frame(Sample,Exon,Types,Gene,Custom_Numbering,Details)
	names(Metrics)=c("Sample","Exon","Type","Gene","Custom.numbering","Info")

if(Custom){
	Metrics_custom=Metrics[Metrics$Custom_Numbering!="NA" | Metrics$Types=="Whole sample",]
	if(nrow(Metrics_custom)>0){
	write.table(Metrics_custom,file=paste(Output,"Metrics.tsv",sep=""),quote=F,row.names=F,sep="\t")}
}
}else{
Metrics<-data.frame(Sample,Exon,Types,Gene,Details)
names(Metrics)=c("Sample","Exon","Type","Gene","Info")
}

if(nrow(Metrics)>0){
	write.table(Metrics,file=paste(Output,"Metrics.tsv",sep=""),quote=F,row.names=F,sep="\t")
}

print("END IdentifyFailures script")