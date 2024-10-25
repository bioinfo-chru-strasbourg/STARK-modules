##########################################################################
# DECON script          Version: 1
# Description:          R script to generate plots
##########################################################################

################## Context ##################
# Original R script from https://github.com/RahmanTeam/DECoN
# DECON is an ExomeDepth wrapper 
########## Note ########################################################################################
# PROD v1 21/11/2023
# Changelog
#   - separate plot script from makeCNVcall, refactor code, remove install systemn, update ExomeDepth 1.16
#   - optparse script
########################################################################################################
 
print("BEGIN DECONPlot script")

suppressPackageStartupMessages({
    library(R.utils)
    library(ExomeDepth)
    library(optparse)
    library(grid)
    library(reshape)
    library(ggplot2)
    library(dplyr)
    library(GenomicRanges)
})

# Function to sanitize gene names for file naming
sanitize_filename <- function(name) {
    # Replace forbidden characters with an underscore
    gsub("[[:punct:] ]+", "_", name)
}

# Function to filter df based on overlapping intervals
filter_df <- function(input_df, filtering_df) {
  
  # Ensure the required columns exist in both data frames
  required_cols <- c("chromosome", "start", "end")

  # Check if the required columns exist (case-insensitive)
  input_cols <- tolower(colnames(input_df))
  filtering_cols <- tolower(colnames(filtering_df))
  
  if (!all(required_cols %in% input_cols)) {
    stop("CNV calls data frame must contain the following columns: chromosome, start, end")
  }
  
  if (!all(required_cols %in% filtering_cols)) {
    stop("BED file data frame must contain the following columns: chromosome, start, end")
  }
    #for debug
    #rdata_file <- "/app/res/input_and_filtering_data.RData"
    #save(input_df, filtering_df, file = rdata_file)

    # Detect the chromosome, start, and end columns and their original cases in input_df
    chrom_column <- if ("chromosome" %in% colnames(input_df)) "chromosome" else "Chromosome"
    start_column <- if ("Start" %in% colnames(input_df)) "Start" else "start"
    end_column <- if ("End" %in% colnames(input_df)) "End" else "end"

    # Save the original Chromosome, Start, and End values to restore later
    input_df$Original_Chromosome <- input_df[[chrom_column]]
    input_df$Original_Start <- input_df[[start_column]]
    input_df$Original_End <- input_df[[end_column]]

    # Standardize chromosome format in input_df
    input_df[[chrom_column]] <- ifelse(grepl("^chr", input_df[[chrom_column]], ignore.case = TRUE),
                                    input_df[[chrom_column]],
                                    paste0("chr", input_df[[chrom_column]]))

    # If `Start` and `End` columns were uppercase, temporarily rename them for processing
    if (start_column == "Start") colnames(input_df)[colnames(input_df) == "Start"] <- "start"
    if (end_column == "End") colnames(input_df)[colnames(input_df) == "End"] <- "end"

    # Standardize filtering_df column names and chromosome format
    colnames(filtering_df) <- tolower(colnames(filtering_df))  # Ensure lowercase column names
    filtering_df$chromosome <- ifelse(grepl("^chr", filtering_df$chromosome, ignore.case = TRUE),
                                    filtering_df$chromosome,
                                    paste0("chr", filtering_df$chromosome))

    # Convert input_df to GRanges using standardized column names
    input_gr <- GRanges(
    seqnames = input_df[[chrom_column]],
    ranges = IRanges(start = input_df$start, end = input_df$end),
    CNV.type = input_df$CNV.type
    )

    # Convert filtering_df to GRanges
    filter_gr <- GRanges(
    seqnames = filtering_df$chromosome,
    ranges = IRanges(start = filtering_df$start, end = filtering_df$end)
    )

    # Find overlaps
    overlaps <- findOverlaps(input_gr, filter_gr)

    # Filter input_df based on overlaps and remove duplicates
    input_df_filtered <- input_df[unique(queryHits(overlaps)), ]

    # Revert Chromosome, Start, and End columns to original values
    input_df_filtered[[chrom_column]] <- input_df_filtered$Original_Chromosome
    input_df_filtered[[start_column]] <- input_df_filtered$Original_Start
    input_df_filtered[[end_column]] <- input_df_filtered$Original_End

    # Drop temporary columns
    input_df_filtered$Original_Chromosome <- NULL
    input_df_filtered$Original_Start <- NULL
    input_df_filtered$Original_End <- NULL


  return(input_df_filtered)
}



###### Parsing input options and setting defaults ########
option_list <- list(
    make_option('--rdata', help = 'Input summary RData file (required)', dest = 'data'),
    make_option('--out', default = './plots', help = 'Output directory, default=./plots', dest = 'pfolder'),
    make_option('--prefix', default = '', help = 'Prefix for the files, default=None', dest = 'prefix'),
    make_option('--bedfile', default = NULL, help = 'Specify a BED file. Default is the one in RData', dest = 'bedfile'), # BED file option
    make_option("--outrdata", default="./DECONplot.Rdata", help="Output Rdata file, default: ./DECONplot.Rdata", dest='outdata')
)
opt <- parse_args(OptionParser(option_list = option_list))

rdata_output = opt$outdata
# Load R workspace with all the results saved : save(ExomeCount,bed.file,counts,fasta,sample.names,bams,cnv.calls,cnv.calls_ids,refs,models,exon_numbers,exons)
count_data=opt$data
# Check if count_data is NULL or an empty string
if(is.null(count_data) || count_data == "") {
    print("ERROR: no RData summary file provided -- Execution halted")
    quit()
}
message('Loading Rdata')
load(count_data)
message('Loading done')

if (length(cnv.calls_ids)==0){
    message('No CNVcalls in the Rdata file. No plot to draw')
    quit()
}

plotFolder=opt$pfolder
if(!file.exists(plotFolder)){dir.create(plotFolder)}

# Set the prefix for file names
prefixfile=opt$prefix
# Get current date and time
current_datetime <- Sys.time()
# Format date and time as "YYYYMMDD-HHMMSS"
formatted_datetime <- format(current_datetime, "%Y%m%d-%H%M%S")

# Check if prefixfile is NULL or an empty string
if (is.null(prefixfile) || prefixfile == "") {
    prefixfile <- formatted_datetime
}

###### Handling BED file option ########
if (!is.null(opt$bedfile)) {
    # Load the specified BED file
    message('Loading specified BED file: ', opt$bedfile)
    # Check and read the BED file, ensuring it has at least four columns
    bed.file <- read.table(opt$bedfile, header = FALSE, sep = "\t",
                        colClasses = c("character", "integer", "integer", "character"))

    # Ensure at least four columns are present
    if (ncol(bed.file) < 4) {
        stop("Error: The BED file must have at least four columns (chromosome, start, end, gene).")
    }

    # Keep only the first four columns
    bed.file <- bed.file[, 1:4]
    # Assign the required column names
    colnames(bed.file) <- c("chromosome", "start", "end", "gene")

    # Filtering
    counts <- filter_df(counts, bed.file)
    ExomeCount <- filter_df(ExomeCount, bed.file)
    cnv.calls <- filter_df(cnv.calls, bed.file)
    cnv.calls_ids <- filter_df(cnv.calls_ids, bed.file)
    save.image(file = rdata_output)
}


#################### Plots ####################################
message('Start generating plots')

if ("Custom.first" %in% colnames(cnv.calls_ids)){
    cnv.calls_plot=cnv.calls_ids[!is.na(cnv.calls_ids$Custom.first),] # filter NA
}else{
    cnv.calls_plot=cnv.calls_ids
}

cnv.calls_plot$chr=paste('chr',cnv.calls_plot$Chromosome,sep='')


Index=vector(length=nrow(bed.file))
Index[1]=1
for(i in 2:nrow(bed.file)){
    if(bed.file[i,4]==bed.file[i-1,4]){
        Index[i]=Index[i-1]+1
    }else{
        Index[i]=1
    }
}

if(colnames(counts)[5]=="exon_number"){
    exons <- bed.file[, c("chromosome", "start", "end", "exon")]
    for(i in 1:nrow(exons)){
        x=which(paste(bed.file[,4])==paste(exons[i,4]) & bed.file[,2]<=exons[i,3] & bed.file[,3]>=exons[i,2])
        Index[x]=exons[i,5]
    }
}

for(call_index in 1:nrow(cnv.calls_plot)){
    Sample<-cnv.calls_plot[call_index,]$Sample
    Gene<-unlist(strsplit(cnv.calls_plot[call_index,]$Gene,split=", "))
    exonRange<-which(bed.file[,4]%in%Gene)
    if((cnv.calls_plot[call_index,]$Start.p-5)<min(exonRange) & ((cnv.calls_plot[call_index,]$Start.p-5)>=1)){exonRange=(cnv.calls_plot[call_index,]$Start.p-5):max(exonRange)}
    if((cnv.calls_plot[call_index,]$Start.p-5)<min(exonRange) & ((cnv.calls_plot[call_index,]$Start.p-5)<=0)){exonRange=1:max(exonRange)}
    if((cnv.calls_plot[call_index,]$End.p+5)>max(exonRange) & ((cnv.calls_plot[call_index,]$End.p+5)<=nrow(bed.file))){exonRange=min(exonRange):(cnv.calls_plot[call_index,]$End.p+5)}
    if((cnv.calls_plot[call_index,]$End.p+5)>max(exonRange) & ((cnv.calls_plot[call_index,]$End.p+5)>nrow(bed.file))){exonRange=min(exonRange):nrow(bed.file)}
    singlechr=length(unique(bed.file[exonRange,1]))==1

    if(!singlechr){
        if(bed.file[exonRange[1],1]!=cnv.calls_plot[call_index,]$chr){
            prev=TRUE
            newchr=bed.file[exonRange[1],1]
        }else{
            prev=FALSE
            newchr=bed.file[exonRange[length(exonRange)],1]
        }
        exonRange=exonRange[bed.file[exonRange,1]==cnv.calls_plot[call_index,]$chr]
    }

    ###### Part of plot containing the coverage points ###############
    VariantExon<- unlist(mapply(function(x,y)x:y,cnv.calls[cnv.calls$Sample==Sample,]$Start.p,cnv.calls[cnv.calls$Sample==Sample,]$End.p))
    refs_sample<-refs[[Sample]]
    Data<-cbind(ExomeCount[exonRange,c(Sample,refs_sample)],exonRange)
    Data[,-ncol(Data)]=log(Data[,-ncol(Data)])
    Data1<-melt(Data,id=c("exonRange"))
    testref<-rep("gray",nrow(Data1))
    testref[Data1$variable==Sample]="blue"
    Data1<-data.frame(Data1,testref)
    levels(Data1$variable)=c(levels(Data1$variable),"VAR")
    Data1$testref=as.factor(Data1$testref)
    levels(Data1$testref)=c(levels(Data1$testref),"red")
    data_temp<-Data1[Data1$variable==Sample & Data1$exonRange%in%VariantExon,]

    if(nrow(data_temp)>0){
        data_temp$variable="VAR"
        data_temp$testref="red"
        Data1<-rbind(Data1,data_temp)
    }

    levels(Data1$testref)=c("Test Sample","Reference Sample","Affected exon")
    new_cols=c("blue","gray","red")
    A1<-ggplot(data=Data1,aes(x=exonRange,y=value,group=variable,colour=testref))
    A1<-A1 + geom_point(cex=2.5,size=1.5)                        #Have to set up points with scale to get correct legend, then re-plot in correct order etc.
    A1<-A1 + scale_colour_manual(values=new_cols)  
    A1<-A1 + geom_line(data=subset(Data1,testref=="Reference Sample"),lty="dashed",lwd=1.5,col="grey") 
    A1<-A1 + geom_point(data=subset(Data1,testref=="Reference Sample"),cex=2.5,col="grey")   
    A1<-A1+ geom_line(data=subset(Data1,testref=="Test Sample"),lty="dashed",lwd=1.5,col="blue")  
    A1<-A1 + geom_point(data=subset(Data1,testref=="Test Sample"),cex=2.5,col="blue") 
    A1<-A1 + geom_point(data=subset(Data1,testref=="Affected exon"),cex=3.5,col="red") 
    A1<-A1 + ylab("Log (Coverage)")  + theme_bw() + theme(legend.position="none")+ xlab(" ")
    Data2<-Data1[Data1$testref=="Affected exon",]
    if(nrow(Data2)>1){
        for(i in 1:(nrow(Data2)-1)){
            if((Data2$exonRange[i]+1)==Data2$exonRange[i+1]){A1<-A1 + geom_line(data=Data2[i:(i+1),],aes(x=exonRange,y=value,group=1),lwd=1.5,col="red")}
        }
    }

    if(!singlechr){
        if(prev){A1<-A1 + scale_x_continuous(breaks=(min(exonRange)-6):max(exonRange),labels=c(rep("",6),paste(Index[exonRange])),limits=c(min(exonRange)-6.75,max(exonRange)))
        }else{A1<-A1 + scale_x_continuous(breaks=min(exonRange):(max(exonRange)+6),labels=c(paste(Index[exonRange]),rep("",6)),limits=c(min(exonRange),max(exonRange)+6.75))}
    }else{A1<-A1 + scale_x_continuous(breaks=exonRange,labels=paste(Index[exonRange]))}


    ############## Part of plot containing the gene names ###########
    genes_sel = unique(bed.file[exonRange,4])
    temp<-cbind(1:nrow(bed.file),bed.file)[exonRange,]
    len<-table(temp$gene)
    mp<-tapply(exonRange,temp[,5],mean)
    mp<-mp[genes_sel]
    len<-len[genes_sel]
    Genes<-data.frame(c(genes_sel),c(mp),c(len-.5),1)
    names(Genes)=c("Gene","MP","Length","Ind")

    if(!singlechr){
        if(prev){
            fakegene=c(newchr,min(exonRange)-5,3.5,1)
        }else{
            fakegene=c(newchr,max(exonRange)+5,3.5,1)
        }
        Genes<-rbind(Genes,fakegene)
        levels(Genes$Gene)=c(levels(Genes$Gene),newchr)
        Genes$MP=as.numeric(Genes$MP)
        Genes$Length=as.numeric(Genes$Length)
    }

    GenesPlot<-ggplot(data=Genes, aes(x=MP,y=Ind,fill=Gene,width=Length,label=Gene)) +geom_tile() + geom_text() + theme_bw() + theme(legend.position="None",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),plot.margin=unit(c(.5,.5,.5,.55),"cm")) + ylab(" ") + xlab(" ")
    GenesPlot<-GenesPlot + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

    ####### Part of the plot containing the normalized ratio of the coverage #############
    Totals<-rowSums(ExomeCount[exonRange,c(Sample,refs_sample)])
    ratio = (ExomeCount[exonRange,Sample]/Totals)/models[[Sample]][1]
    mins <- vector(length=length(exonRange))
    maxs <-vector(length=length(exonRange))
    for(i in 1:length(exonRange)){
        temp = qbetabinom(p=0.025,Totals[i],models[[Sample]][2],models[[Sample]][1])
        mins[i] = (temp/Totals[i])/models[[Sample]][1]
        temp = qbetabinom(p=0.975,Totals[i],models[[Sample]][2],models[[Sample]][1])
        maxs[i] = (temp/Totals[i])/models[[Sample]][1]
    }       

    CIData<-data.frame(exonRange,ratio,mins,maxs)
    names(CIData)<-c("Exon","Ratio","Min","Max")
    CIPlot<-ggplot(CIData,aes(x=Exon,y=Ratio))+geom_ribbon(aes(ymin=Min,ymax=Max),fill="grey")+geom_point(col="blue",cex=3.5) + theme_bw() +xlab("")+ylab("Observed/Expected")
    temp = cnv.calls[cnv.calls$Sample==Sample,]
    if(sum(temp$Start.p%in%exonRange |temp$End.p%in%exonRange)>0){
        temp = temp[temp$Start.p%in%exonRange|temp$End.p%in%exonRange,]
        for(i in 1:nrow(temp)){
            start.temp = temp[i,]$Start.p
            end.temp = temp[i,]$End.p
            CIPlot<-CIPlot + geom_point(data=CIData[CIData$Exon%in%start.temp:end.temp,], aes(x=Exon,y=Ratio),color="red",cex=3.5)
        }
    }

    if(!singlechr){
        if(prev){CIPlot<- CIPlot + scale_x_continuous(breaks=(min(exonRange)-6):max(exonRange),labels=c(rep("",6),paste(Index[exonRange])),limits=c(min(exonRange)-6.75,max(exonRange)))
        }else{CIPlot<-CIPlot + scale_x_continuous(breaks=min(exonRange):(max(exonRange)+6),labels=c(paste(Index[exonRange]),rep("",6)),limits=c(min(exonRange),max(exonRange)+6.75))}
    }else{CIPlot<-CIPlot + scale_x_continuous(breaks=exonRange,labels=paste(Index[exonRange]))}

    ######### Save plot in pdf format ###########
    cnv_genes_sample=cnv.calls_plot[cnv.calls_plot$Sample==Sample,]$Gene
    cleaned_gene <- sanitize_filename(Gene)
    if (sum(cnv_genes_sample == cleaned_gene) == 1) {
        pdf(file = paste(plotFolder, "/DECON.", prefixfile, ".", Sample, ".", cleaned_gene, ".pdf", sep = ""), useDingbats = FALSE)
    } else {
        cnv_genes_sample_index = which(cnv.calls_plot$Sample == Sample & cnv.calls_plot$Gene == paste(cleaned_gene, collapse = ", "))
        pdf(file = paste(plotFolder, "/DECON.", prefixfile, ".", Sample, ".", paste(cleaned_gene, collapse = "_"), "_", which(cnv_genes_sample_index == call_index), ".pdf", sep = ""), useDingbats = FALSE)
    }

    if (sum(cnv_genes_sample == cleaned_gene) > 2) {
        print(paste("WARNING: more than 2 calls in ", cleaned_gene, ", could affect plotting", sep = ""))
    }

 
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(6, 1)))
    print(A1,vp=viewport(layout.pos.row=1:3,layout.pos.col=1))
    print(GenesPlot,vp=viewport(layout.pos.row=4,layout.pos.col=1))
    print(CIPlot,vp=viewport(layout.pos.row=5:6,layout.pos.col=1))
    dev.off()
}
warnings()
print("END DeconPlot script")