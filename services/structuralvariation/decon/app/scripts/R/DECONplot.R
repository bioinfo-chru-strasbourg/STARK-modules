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

# Function to sanitize gene names for file naming
sanitize_filename <- function(name) {
    # Replace forbidden characters with an underscore
    gsub("[[:punct:] ]+", "_", name)
}

# To exclude "X" and "Y" chromosomes: filtered_df <- filter_chromosomes(df, exclude.chrom = c("X", "Y"))
# To include only chromosome "X": filtered_df <- filter_chromosomes(df, include.chrom = c("X"))
# Work if chromosome is "chrX" or "X" in the df
filter_chromosomes <- function(df, include.chrom = NULL, exclude.chrom = NULL) {
  # Check if the dataframe has either "chromosome" or "Chromosome" column
  if ("chromosome" %in% colnames(df)) {
    chrom_col <- "chromosome"
  } else if ("Chromosome" %in% colnames(df)) {
    chrom_col <- "Chromosome"
  } else {
    stop("Data frame must contain a 'chromosome' or 'Chromosome' column.")
  }
  
  # Create versions of include/exclude lists with and without "chr" prefix
  if (!is.null(include.chrom)) {
    include.chrom <- unique(c(include.chrom, sub("^chr", "", include.chrom), paste0("chr", include.chrom)))
  }
  if (!is.null(exclude.chrom)) {
    exclude.chrom <- unique(c(exclude.chrom, sub("^chr", "", exclude.chrom), paste0("chr", exclude.chrom)))
  }
  
  # Apply inclusion and exclusion filters
  if (!is.null(include.chrom)) {
    df <- subset(df, df[[chrom_col]] %in% include.chrom)
  }
  if (!is.null(exclude.chrom)) {
    df <- subset(df, !df[[chrom_col]] %in% exclude.chrom)
  }
  
  return(df)
}


# Function to lowercase the specified columns name if necessary
lowercase_required_cols <- function(df, cols) {
  # Convert the column names and the specified columns to lowercase for case-insensitive matching
  colnames_lower <- tolower(colnames(df))
  cols_lower <- tolower(cols)
  
  # Find which columns should be lowercased by matching case-insensitively
  matching_cols <- colnames(df)[colnames_lower %in% cols_lower]
  
  # Lowercase the matching column names
  colnames(df)[colnames(df) %in% matching_cols] <- tolower(colnames(df)[colnames(df) %in% matching_cols])
  
  return(df)
}

# Function to standardize the chromosome format
standardize_chromosome <- function(df) {
  df$chromosome <- gsub("^chr(.*)$", "chr\\1", tolower(df$chromosome))
  df$chromosome[df$chromosome == "chrx"] <- "chr23"
  df$chromosome[df$chromosome == "chry"] <- "chr24"
  return(df)
}

remove_duplicate_rows <- function(df, subset_cols) {
  # Convert subset_cols to lower case
  subset_cols_lower <- tolower(subset_cols)
  
  # Find matching columns in the data frame (case-insensitive)
  matched_cols <- tolower(colnames(df)) %in% subset_cols_lower
  
  # If no matching columns are found, stop with an error
  if (all(!matched_cols)) {
    stop("None of the specified columns were found in the data frame.")
  }
  
  # Select the columns that match (case-insensitive)
  df_matched <- df[, matched_cols, drop = FALSE]
  
  # Remove duplicate rows based on only the matched columns
  df <- df[!duplicated(df_matched), ]
  
  return(df)
}

# Function to capitalize the first letter of specific columns
capitalize_first_letter <- function(df, cols) {
  colnames(df)[colnames(df) %in% cols] <- sub("^(.)", "\\U\\1", colnames(df)[colnames(df) %in% cols], perl = TRUE)
  return(df)
}


# Function to add or remove a prefix from a specific column in a dataframe
modify_prefix <- function(df, column_name, action = "remove", prefix = "chr") {
  # Check if the column exists in the dataframe
  if (!(column_name %in% colnames(df))) {
    stop("The specified column does not exist in the dataframe.")
  }
  
  # Perform the action based on the user's choice
  if (action == "remove") {
    df[[column_name]] <- sub(paste0("^", prefix), "", df[[column_name]])
  } else if (action == "add") {
    df[[column_name]] <- paste0(prefix, df[[column_name]])
  } else {
    stop("Invalid action. Choose 'remove' or 'add'.")
  }
  
  return(df)
}

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


###### Parsing input options and setting defaults ########
option_list <- list(
    make_option('--rdata', help = 'Input summary RData file (required)', dest = 'data'),
    make_option('--out', default = './plots', help = 'Output directory, default=./plots', dest = 'pfolder'),
    make_option("--chromosome", default="A", help='Perform plots for autosomes or chrX', dest='chromosome'),
    make_option('--prefix', default = '', help = 'Prefix for the files, default=None', dest = 'prefix'),
    make_option("--outrdata", default="./DECONplot.Rdata", help="Output Rdata file, default: ./DECONplot.Rdata", dest='outdata'),
	make_option('--debug', action="store_true", default=FALSE, help="Enable debug mode to save intermediate RData files", dest="debug")

)
opt <- parse_args(OptionParser(option_list = option_list))


# Get current date and time & format date and time as "YYYYMMDD-HHMMSS"
current_datetime <- Sys.time()
formatted_datetime <- format(current_datetime, "%Y%m%d-%H%M%S")

# Set the prefix for file names
prefixfile=opt$prefix
# Check if prefixfile is NULL or an empty string
if (is.null(prefixfile) || prefixfile == "") {
    prefixfile <- formatted_datetime
}

# Debug
if (opt$debug) {
debugFolder <- paste("/app/res/debug", prefixfile, sep = "/")
dir.create(debugFolder, recursive = TRUE, showWarnings = FALSE)
}

# set output
rdata_output = opt$outdata

# Check if count_data is NULL or an empty string
count_data=opt$data
if(is.null(count_data) || count_data == "") {
    print("[ERROR] no RData summary file provided -- Execution halted")
    quit()
}
message('[INFO] Loading Rdata')
load(count_data)
message('[INFO] Loading done')

if (length(cnv.calls_ids)==0){
    message('[ERROR] No CNVcalls in the Rdata file. No plot to draw')
    quit()
}

plotFolder=opt$pfolder
if(!file.exists(plotFolder)){dir.create(plotFolder)}

# Chromosome filtering
modechrom = opt$chromosome
if (modechrom == "XX" || modechrom == "XY") {
	bed.file <- filter_chromosomes(bed.file, include.chrom = c("chrX")) # we don't call Y
	ExomeCount <- filter_chromosomes(ExomeCount, include.chrom = c("chrX"))
	cnv.calls_ids <- filter_chromosomes(cnv.calls_ids, include.chrom = c("chrX"))
}
if (modechrom == "A") {
	bed.file <- filter_chromosomes(bed.file, exclude.chrom = c("chrX", "chrY"))
	ExomeCount <- filter_chromosomes(ExomeCount, exclude.chrom = c("chrX", "chrY"))
	cnv.calls_ids <- filter_chromosomes(cnv.calls_ids, exclude.chrom = c("chrX", "chrY"))
}

# Debug
if (opt$debug) {
rdata_file <- sprintf("%s/debug_filtering_data_%s.RData", debugFolder, modechrom)
save(bed.file , models, refs, ExomeCount, cnv.calls_ids, file = rdata_file)
}

# Creating cnv.calls_plot
cnv.calls_plot <- transform(cnv.calls_ids, Chromosome = paste('chr', Chromosome, sep = ''))

# Starting index
Index=vector(length=nrow(bed.file))
Index[1]=1
for(i in 2:nrow(bed.file)){
	if(bed.file[i,4]==bed.file[i-1,4]){
		Index[i]=Index[i-1]+1
	}else{
		Index[i]=1
	}   
}   

if (colnames(bed.file)[5] == "exon") {
    message('[INFO] Exon numbers detected')
    exons <- bed.file
    for (i in 1:nrow(exons)) {
        x = which(paste(bed.file[,5]) == paste(exons[i,5]) & 
                  bed.file[,2] <= exons[i,3] & 
                  bed.file[,3] >= exons[i,2])
        Index[x] = exons[i,5]
    }
}

if (opt$debug) {
rdata_file <- sprintf("%s/debug_start_plot_%s.RData", debugFolder, modechrom)
save(exons, cnv.calls_plot, bed.file , models, refs, ExomeCount, cnv.calls_ids, file = rdata_file)
}

for(call_index in 1:nrow(cnv.calls_plot)){
	
	Sample <- cnv.calls_plot[call_index,]$Sample
	Gene <- unlist(strsplit(cnv.calls_plot[call_index,]$Gene,split=", "))
	exonRange <- which(bed.file[,4]%in%Gene)

	if((cnv.calls_plot[call_index,]$Start.p-5) < min(exonRange) & ((cnv.calls_plot[call_index,]$Start.p-5) >= 1)){
		exonRange = (cnv.calls_plot[call_index,]$Start.p-5):max(exonRange)
		}
	
	if((cnv.calls_plot[call_index,]$Start.p-5) < min (exonRange) & ((cnv.calls_plot[call_index,]$Start.p-5) <= 0 )){
		exonRange = 1:max(exonRange)
		}
	
	if((cnv.calls_plot[call_index,]$End.p+5) > max(exonRange) & ((cnv.calls_plot[call_index,]$End.p+5) <= nrow (bed.file))){
		exonRange = min(exonRange):(cnv.calls_plot[call_index,]$End.p+5)
		}

	if((cnv.calls_plot[call_index,]$End.p+5) > max (exonRange) & ((cnv.calls_plot[call_index,]$End.p+5) > nrow (bed.file))){
		exonRange = min(exonRange):nrow(bed.file)
		}

	singlechr=length(unique(bed.file[exonRange,1]))==1

	if(!singlechr){
		if(bed.file[exonRange[1],1]!=cnv.calls_plot[call_index,]$Chromosome){
			prev=TRUE
			newchr=bed.file[exonRange[1],1]
		}else{
			prev=FALSE 
			newchr=bed.file[exonRange[length(exonRange)],1]
		}
		exonRange=exonRange[bed.file[exonRange,1]==cnv.calls_plot[call_index,]$Chromosome]
	}

	###### Part of plot containing the coverage points (upper part) ###############
	message('[INFO] Starting drawing coverage')
	VariantExon <- unlist(mapply(function(x,y)x:y,cnv.calls_plot[cnv.calls_plot$Sample==Sample,]$Start.p,cnv.calls_plot[cnv.calls_plot$Sample==Sample,]$End.p))
	refs_sample<-refs[[Sample]]
	Data <- cbind(ExomeCount[exonRange,c(Sample,refs_sample)],exonRange)
	Data[,-ncol(Data)] = log(Data[,-ncol(Data)])
	
	if (opt$debug) {
    rdata_file <- sprintf("%s/debug_before_melt_plot_%s_%s.RData", debugFolder, modechrom, call_index)
    save(Gene, VariantExon, Sample, refs_sample, Data, exons, cnv.calls_plot, bed.file , models, refs, ExomeCount, cnv.calls_ids, exonRange, file = rdata_file)
	}
	
	Data1<-melt(Data,id=c("exonRange"))
	testref<-rep("gray",nrow(Data1))
	testref[Data1$variable==Sample] = "blue"
	Data1<-data.frame(Data1,testref)
	levels(Data1$variable)=c(levels(Data1$variable),"VAR")
	Data1$testref=as.factor(Data1$testref)
	levels(Data1$testref)=c(levels(Data1$testref),"red")

	data_temp <- Data1[Data1$variable == Sample & Data1$exonRange%in%VariantExon,]

	if(nrow(data_temp)>0){
		Data1 <- rbind(Data1, transform(data_temp, variable = "VAR", testref = "red"))
	}
	# Include the sample name 
	sample_name_label <- paste("Test Sample (", Sample, ")", sep = "")
	levels(Data1$testref)=c(sample_name_label,"Reference Sample","Affected exon")
	new_cols=c("blue","gray","red")

	A1 <- ggplot(data=Data1,aes(x=exonRange,y=value,group=variable,colour=testref))
	A1 <- A1 + geom_point(cex=2.5,lwd=1.5)
	A1 <- A1 + scale_colour_manual(values=new_cols)  
	A1 <- A1 + geom_line(data=subset(Data1,testref=="Reference Sample"),lty="dashed",lwd=1.5,col="grey") 
	A1 <- A1 + geom_point(data=subset(Data1,testref=="Reference Sample"),cex=2.5,col="grey")   
	A1 <- A1 + geom_line(data=subset(Data1,testref==sample_name_label),lty="dashed",lwd=1.5,col="blue")  
	A1 <- A1 + geom_point(data=subset(Data1,testref==sample_name_label),cex=2.5,col="blue") 
	A1 <- A1 + geom_point(data=subset(Data1,testref=="Affected exon"),cex=3.5,col="red") 
    A1 <- A1 + ylab("Log (Coverage)")  + theme_bw() + theme(legend.position= "top") + xlab(" ") + guides(color = guide_legend(nrow = 1, title = NULL))

	Data2 <- Data1[Data1$testref=="Affected exon",]
	if(nrow(Data2)>1){
		for(i in 1:(nrow(Data2)-1)){
			if ((Data2$exonRange[i]+1)==Data2$exonRange[i+1]) {
				A1<-A1 + geom_line(data=Data2[i:(i+1),],aes(x=exonRange,y=value,group=1),lwd=1.5,col="red")
			}
		}
	}   

	if (!singlechr) {
		if (prev) {
			A1 <- A1 + scale_x_continuous(breaks=(min(exonRange)-6):max(exonRange), labels=c(rep("",6),paste(Index[exonRange])),limits=c(min(exonRange)-6.75,max(exonRange))) + theme(axis.text.x = element_text(angle = 45, size = 6, hjust = 1))
		} else {
			A1 <- A1 + scale_x_continuous(breaks=min(exonRange):(max(exonRange)+6),labels=c(paste(Index[exonRange]),rep("",6)),limits=c(min(exonRange),max(exonRange)+6.75)) + theme(axis.text.x = element_text(angle = 45, size = 6, hjust = 1))
			}
	} else {
	A1 <- A1 + scale_x_continuous(breaks=exonRange, labels = Index[exonRange]) + theme(axis.text.x = element_text(angle = 45, size = 6, hjust = 1))
	}

	############## Part of plot containing the gene names as legend ###########
	message('[INFO] Adding gene names')
	genes_sel = unique(bed.file[exonRange,4])
	temp <- cbind(1:nrow(bed.file),bed.file)[exonRange,]
	len <- table(temp$gene)
	mp <- tapply(exonRange,temp[,5],mean)
	mp <- mp[genes_sel]
	len <- len[genes_sel]
	Genes <- data.frame(c(genes_sel),c(mp),c(len-.5),1)
    names(Genes) = c("Gene","MP","Length","Ind")

	if(!singlechr){
		if(prev){
			fakegene=c(newchr,min(exonRange)-5,3.5,1)
		} else {
			fakegene=c(newchr,max(exonRange)+5,3.5,1)
		}
		Genes <- rbind(Genes,fakegene)
		Genes$MP = as.numeric(Genes$MP)
		Genes$Length = as.numeric(Genes$Length)
	}

	GenesPlot <- ggplot(data=Genes, aes(x=MP,y=Ind,fill=Gene,width=Length,label=Gene)) + geom_tile() + geom_text() + theme_bw() + theme(legend.position="None",panel.grid.major = element_blank(), panel.grid.minor = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),plot.margin=unit(c(.5,.5,.5,.5),"cm")) + ylab(" ") + xlab(" ")
	GenesPlot <- GenesPlot + theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())

    # Debug
	if (opt$debug) {
		rdata_file <- sprintf("%s/debug_after_adding_genes_names_%s_%s.RData", debugFolder, modechrom, call_index)
		save(exons, Index, singlechr, genes_sel, temp, len, mp, Genes, GenesPlot, A1, models, refs, bed.file, cnv.calls_plot, exonRange, Gene, Sample, VariantExon, refs_sample, ExomeCount, cnv.calls_ids, Data, file = rdata_file)
	}

	####### Part of the plot containing the normalized ratio of the coverage (lower part) #############
	message('[INFO] Adding normalized ratio')
	Totals <- rowSums(ExomeCount[exonRange,c(Sample,refs_sample)])
	ratio = (ExomeCount[exonRange,Sample]/Totals)/models[[Sample]][1]
	mins <- maxs <- rep(NA, length(exonRange))
	
	for(i in 1:length(exonRange)){
		temp = qbetabinom(p=0.025,Totals[i],models[[Sample]][2],models[[Sample]][1])
		mins[i] = (temp/Totals[i])/models[[Sample]][1]
		temp = qbetabinom(p=0.975,Totals[i],models[[Sample]][2],models[[Sample]][1])
		maxs[i] = (temp/Totals[i])/models[[Sample]][1]
	}       

	CIData <- setNames(data.frame(exonRange, ratio, mins, maxs), c("Exon", "Ratio", "Min", "Max"))
	CIPlot <- ggplot(CIData, aes(x=Exon,y=Ratio)) + geom_ribbon(aes(ymin=Min,ymax=Max),fill="grey") + geom_point(col="blue",cex=3.5) + theme_bw() + xlab("") + ylab("Observed/Expected")
	
	temp = cnv.calls_plot[cnv.calls_plot$Sample==Sample,]
	if (sum(temp$Start.p%in%exonRange | temp$End.p%in%exonRange) > 0){
		temp = temp[temp$Start.p%in%exonRange|temp$End.p%in%exonRange,]
		for(i in 1:nrow(temp)){
			start.temp = temp[i,]$Start.p
			end.temp = temp[i,]$End.p
			CIPlot <- CIPlot + geom_point(data=CIData[CIData$Exon%in%start.temp:end.temp,], aes(x=Exon,y=Ratio),color="red",cex=3.5)
		}
	}

	if (!singlechr) {
		if (prev) {
			CIPlot <- CIPlot + scale_x_continuous(breaks=(min(exonRange)-6):max(exonRange),labels=c(rep("",6),paste(Index[exonRange])),limits=c(min(exonRange)-6.75,max(exonRange))) + theme(axis.text.x = element_text(angle = 45, size = 6, hjust = 1))
		} else {
			CIPlot <- CIPlot + scale_x_continuous(breaks=min(exonRange):(max(exonRange)+6),labels=c(paste(Index[exonRange]),rep("",6)),limits=c(min(exonRange),max(exonRange)+6.75)) + theme(axis.text.x = element_text(angle = 45, size = 6, hjust = 1))
			}
	} else {
		CIPlot <- CIPlot + scale_x_continuous(breaks=exonRange, labels = Index[exonRange]) + theme(axis.text.x = element_text(angle = 45, size = 6, hjust = 1))
		}

	# Debug
	if (opt$debug) {
	rdata_file <- sprintf("%s/debug_end_plot_data_%s_%s.RData", debugFolder, modechrom, call_index)
    save(singlechr, A1, Index, CIPlot, CIData, GenesPlot, exonRange, Totals, ratio, mins, maxs, models, bed.file, ExomeCount, cnv.calls_ids, cnv.calls_plot, file = rdata_file)
	}

    ######### Save plot in pdf format ###########
    message('[INFO] Saving plots in pdf format')
    cnv_genes_sample=cnv.calls_plot[cnv.calls_plot$Sample==Sample,]$Gene
    cleaned_gene <- sanitize_filename(Gene) # for naming files

    if (sum(cnv_genes_sample == cleaned_gene) == 1) {
        pdf(file = paste(plotFolder, "/DECON.", prefixfile, ".", Sample, ".", cleaned_gene, ".pdf", sep = ""), useDingbats = FALSE)
    } else {
        cnv_genes_sample_index = which(cnv.calls_plot$Sample == Sample & cnv.calls_plot$Gene == paste(cleaned_gene, collapse = ", "))
        pdf(file = paste(plotFolder, "/DECON.", prefixfile, ".", Sample, ".", paste(cleaned_gene, collapse = "_"), "_", which(cnv_genes_sample_index == call_index), ".pdf", sep = ""), useDingbats = FALSE)
    }

if(sum(cnv_genes_sample==Gene)>2){
	print(paste("[WARNING] More than 2 calls in ",Gene,", could affect plotting",sep=""))
	}

	grid.newpage()
	pushViewport(viewport(layout = grid.layout(6, 1)))
	print(A1,vp=viewport(layout.pos.row=1:3,layout.pos.col=1))
	print(GenesPlot,vp=viewport(layout.pos.row=4,layout.pos.col=1))
	print(CIPlot,vp=viewport(layout.pos.row=5:6,layout.pos.col=1))
	dev.off()

}