##########################################################################
# DECON script          Version: 1
# Description:          R script to call CNVs from coverage datas
##########################################################################

################## Context ###############################################
# Original R script from https://github.com/RahmanTeam/DECoN
# DECON is an ExomeDepth wrapper 
########## Note ########################################################################################
# PROD v1 21/11/2023
# Changelog
#   - refactor code, remove install system, update for ExomeDepth 1.16
#   - optparse script, add a option to use a list of ref bam files for comparaison
#   - detect bed file with a 5 column with exon numbers, remove custom exon options
#   - mode to remove chrX or Y for calling
########################################################################################################

print("BEGIN makeCNVCalls script")

suppressMessages(library(R.utils))
suppressMessages(library(optparse))
suppressMessages(library(ExomeDepth))

###### Parsing input options and setting defaults ########
option_list<-list(
    make_option('--RData',help='Input summary RData file containing coverage datas, bed file and GC content (required)',dest='RData'),
    make_option("--chromosome",default="A",help='Performe calling for Autosomes, XX or XY only',dest='chromosome'),
    make_option('--transProb',default=.01,help='Transition probability for the HMM statistical analysis, default=0.01',dest='transProb'),
	make_option("--refbams",default=NULL,help="Text file containing a list of reference bam files for calling (full path) (optional)",dest='refbams'),
    make_option('--out',default='./results',help='Output folder, default=./results',dest='out')
)
opt<-parse_args(OptionParser(option_list=option_list))

# R workspace with the coverage data, the bedfile & ref genome fasta file saved in it
count_data=opt$RData
if(count_data=="NULL"){count_data=NULL}
if(is.null(count_data)){
    print("ERROR: no RData summary file provided -- Execution halted")
    quit()
}
load(count_data)

modechrom=opt$chromosome
refbams_file=opt$refbams
trans_prob=as.numeric(opt$transProb)
output=opt$out

if(!file.exists(output)){dir.create(output)}

# function which recursively splits x by an element of 'splits' then extracts the y element of the split vector
multi_strsplit<-function(x,splits,y){
	X<-x
	for(i in 1:length(splits)){X=strsplit(X,splits[i], fixed = TRUE)[[1]][y[i]]}
	return(X)
}

################# CNV CALLING #######################################
# To call autosomes you need to use a bed containing autosome only
# To call chrom X you need to use a bed containing chrom X only AND separate Male & Female patients 
ExomeCount<-as(counts, 'data.frame')
# remove any chr letters, and coerce to a string
ExomeCount$chromosome <- gsub(as.character(ExomeCount$chromosome),pattern = 'chr',replacement = '') 

# filter ExomeCount by modechrom
if (modechrom=="A"){
    ExomeCount<-subset(ExomeCount, chromosome!="X" & chromosome!="Y")
}

if (modechrom=="XX"){
    ExomeCount<-subset(ExomeCount, chromosome=="X")
}

if (modechrom=="XY"){
    ExomeCount<-subset(ExomeCount, chromosome=="X" | chromosome=="Y")
}

cnv.calls = NULL
refs<-list()
models<-list()
# for each sample :
# extracts the sample to be tested, uses all other samples as the reference set, places coverage info for all samples in the reference set into a matrix 
# from the reference set, selects correlated samples to be used.
# places selected, correlated samples into a matrix, sums the selected samples across each exon
# creates ExomeDepth object containing test data, reference data, and linear relationship between them. Automatically calculates likelihoods
# fits a HMM with 3 states to read depth data; transition.probability - transition probability for HMM from normal to del/dup. Returns ExomeDepth object with CNVcalls
# if refsample.names is not empty, the comparison will be made with those ref bams
# filter sample.names removing refsample.names

# Read list of refbams to process
if(length(refbams_file)>0){
refbams<-apply(read.table(paste(refbams_file)),1,toString)
# get the sample names (first part of the filename, separated by dot)
a<-length(strsplit(refbams[1],"/")[[1]])
refsample.names<-sapply(refbams,multi_strsplit,c("/","."),c(a,1))
sample.names <- sample.names[!sample.names %in% refsample.names]
}else{
refsample.names<-vector()
}

for(i in 1:length(sample.names)){
    print(paste("Processing sample: ",sample.names[i]," ", i,"/",length(sample.names),sep=""))
    my.test <- ExomeCount[,sample.names[i]]
	if(length(refsample.names)==0){
		my.ref.samples <- sample.names[-i]
    }else{
        my.ref.samples <- refsample.names
    }
    my.reference.set <- as.matrix(ExomeCount[,my.ref.samples])
    my.choice <- select.reference.set (test.counts = my.test,reference.counts = my.reference.set,bin.length = (ExomeCount$end - ExomeCount$start)/1000,n.bins.reduced = 10000)
    my.matrix <- as.matrix( ExomeCount[, my.choice$reference.choice, drop = FALSE])
    my.reference.selected <- apply(X = my.matrix,MAR = 1,FUN = sum)
    all.exons <- new('ExomeDepth', test = my.test, reference = my.reference.selected, formula = 'cbind(test, reference) ~ 1')
    all.exons <- CallCNVs(x = all.exons, transition.probability = trans_prob, chromosome = ExomeCount$chromosome, start = ExomeCount$start, end = ExomeCount$end, name = ExomeCount$gene)
    my.ref.counts <- apply(my.matrix, MAR = 1, FUN = sum)
    if(nrow(all.exons@CNV.calls)>0){
        cnvs= cbind(sample.names[i],cor(my.test, my.ref.counts),length( my.choice[[1]] ),all.exons@CNV.calls)
        cnv.calls = rbind(cnv.calls,cnvs)
    }
    refs[[i]] = my.choice$reference.choice
    models[[i]] = c(all.exons@expected[1],all.exons@phi[1])
}

names(refs) = sample.names
names(models) = sample.names

if(!is.null(cnv.calls)){
names(cnv.calls)[1] = "sample"
cnv.calls$sample = paste(cnv.calls$sample)
names(cnv.calls)[2] = "correlation"
names(cnv.calls)[3] = "N.comp"

##################### CONFIDENCE CALCULATION ##############################
Flag<-rep("HIGH",nrow(cnv.calls))

a<-which(cnv.calls$correlation<.985)
if(length(a)>0){Flag[a]="LOW"}

a<-which(cnv.calls$reads.ratio<1.25 & cnv.calls$reads.ratio>.75)
if(length(a)>0){Flag[a]="LOW"}

a<-which(cnv.calls$N.comp<=3)
if(length(a)>0){Flag[a]="LOW"}

if(ncol(bed.file)>=4){
    genes<-apply(cnv.calls,1,function(x)paste(unique(bed.file[x[4]:x[5],4]),collapse=", "))
    a<-which(genes=="PMS2")
    if(length(a)>0){Flag[a]="LOW"}
}

cnv.calls<-cbind(cnv.calls,genes)
names(cnv.calls)[ncol(cnv.calls)]="Gene"

cnv.calls<-cbind(cnv.calls,Flag)
names(cnv.calls)[ncol(cnv.calls)]="Flag"

# replaces single calls involving multiple genes with multiple calls with a single call ID/gene
cnv.calls_ids=cbind(1:nrow(cnv.calls),cnv.calls)
names(cnv.calls_ids)[1]="ID"
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

for(i in 1:nrow(cnv.calls_ids)){
    genes=strsplit(paste(cnv.calls_ids[i,]$Gene),",")[[1]]
    genes=trim(genes)
    whole.index=cnv.calls_ids[i,]$start.p:cnv.calls_ids[i,]$end.p
    if(length(genes)>1){
        temp=cnv.calls_ids[rep(i,length(genes)),]
        temp$Gene=genes
        for(j in 1:length(genes)){
            gene.index=which(bed.file[,4]==genes[j])
            overlap=gene.index[gene.index%in%whole.index]
            temp[j,]$start.p=min(overlap)
            temp[j,]$end.p=max(overlap)
        }
        if(i==1){
            cnv.calls_ids=rbind(temp,cnv.calls_ids[(i+1):nrow(cnv.calls_ids),])
        }else if(i==nrow(cnv.calls_ids)){
            cnv.calls_ids=rbind(cnv.calls_ids[1:(i-1),],temp)
        }else{
            cnv.calls_ids=rbind(cnv.calls_ids[1:(i-1),],temp,cnv.calls_ids[(i+1):nrow(cnv.calls_ids),])
        }
    }
}

cnv.calls_ids$Gene=trim(cnv.calls_ids$Gene)
Gene.index=vector()
genes_unique<-unique(bed.file[,4])

for(i in 1:length(genes_unique)){
    Gene.index=c(Gene.index,1:sum(bed.file[,4]==genes_unique[i]))
}

start.b<-Gene.index[cnv.calls_ids$start.p]
end.b<-Gene.index[cnv.calls_ids$end.p]
cnv.calls_ids<-cbind(cnv.calls_ids,start.b,end.b)
}

################# ADD Custom EXON NUMBERS #####################
if(colnames(counts)[5]=="exon"){
    Custom.first = rep(NA,nrow(cnv.calls_ids))
    Custom.last=rep(NA,nrow(cnv.calls_ids))
    # exons<-read.table(exon_numbers,sep="\t",header=T)
    exons <- bed.file
    exonnumber<-sapply(cnv.calls_ids$Gene, '==',exons$Gene)
    for(i in 1:nrow(exonnumber)){
        for(j in 1:ncol(exonnumber)){
            temp=cnv.calls_ids$start[j]<=exons$End[i] & cnv.calls_ids$end[j]>=exons$Start[i]
            exonnumber[i,j]=exonnumber[i,j] & temp
        }
    }
    exonlist<-which(colSums(exonnumber)!=0)
    if(length(exonlist)>0){
        a<-list(length=length(exonlist))
        for(i in 1:length(exonlist)){a[[i]]=which(exonnumber[,exonlist[i]])}
        first_exon<-unlist(lapply(a,function(a,b)min(b[a,]$Custom),exons))		#identifies the  first and last Custom exon in the deletion/duplication.
        last_exon<-unlist(lapply(a,function(a,b)max(b[a,]$Custom),exons))
        Custom.first[exonlist] = first_exon
        Custom.last[exonlist] = last_exon
    }
    cnv.calls_ids=cbind(cnv.calls_ids,Custom.first,Custom.last)
}

##################### Output #######################
if(colnames(counts)[5]=="exon"){
cnv.calls_ids_out<-data.frame(cnv.calls_ids$ID,cnv.calls_ids$sample,cnv.calls_ids$correlation,cnv.calls_ids$N.comp,cnv.calls_ids$start.p,cnv.calls_ids$end.p,cnv.calls_ids$type,cnv.calls_ids$nexons,cnv.calls_ids$start,cnv.calls_ids$end,cnv.calls_ids$chromosome,cnv.calls_ids$id,cnv.calls_ids$BF,cnv.calls_ids$reads.expected,cnv.calls_ids$reads.observed,cnv.calls_ids$reads.ratio,cnv.calls_ids$Gene,cnv.calls_ids$Custom.first,cnv.calls_ids$Custom.last)
names(cnv.calls_ids_out)<-c("CNV.ID","Sample","Correlation","N.comp","Start.b","End.b","CNV.type","N.exons","Start","End","Chromosome","Genomic.ID","BF","Reads.expected","Reads.observed","Reads.ratio","Gene","Custom.first","Custom.last")
save(ExomeCount,bed.file,counts,fasta,sample.names,bams,cnv.calls,cnv.calls_ids,refs,models,exon_numbers,exons,file=paste(output,"/CNVcall_custom.RData",sep=""))
custom_calls=cnv.calls_ids_out[!is.na(cnv.calls_ids$Custom.first),]
write.table(custom_calls,file=paste(output,"/CNVcalls_custom.tsv",sep=""),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)

}else{

cnv.calls_ids_out<-data.frame(cnv.calls_ids$ID,cnv.calls_ids$sample,cnv.calls_ids$correlation,cnv.calls_ids$N.comp,cnv.calls_ids$start.p,cnv.calls_ids$end.p,cnv.calls_ids$type,cnv.calls_ids$nexons,cnv.calls_ids$start,cnv.calls_ids$end,cnv.calls_ids$chromosome,cnv.calls_ids$id,cnv.calls_ids$BF,cnv.calls_ids$reads.expected,cnv.calls_ids$reads.observed,cnv.calls_ids$reads.ratio,cnv.calls_ids$Gene,cnv.calls_ids$Flag)
names(cnv.calls_ids_out)<-c("CNV.ID","Sample","Correlation","N.comp","Start.b","End.b","CNV.type","N.exons","Start","End","Chromosome","Genomic.ID","BF","Reads.expected","Reads.observed","Reads.ratio","Gene", "Flag")
save(ExomeCount,bed.file,counts,fasta,sample.names,bams,cnv.calls,cnv.calls_ids,refs,models,file=paste(output,"/CNVcall.RData",sep=""))
write.table(cnv.calls_ids_out,file=paste(output,"/CNVcalls.tsv",sep=""),sep="\t",quote=FALSE,col.names=TRUE,row.names=FALSE)
}

print("END makeCNVCalls script")