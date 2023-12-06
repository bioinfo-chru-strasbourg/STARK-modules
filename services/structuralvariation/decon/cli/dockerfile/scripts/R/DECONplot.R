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
#   - separate plot script from makeCNVcall, refactor code, remove install system
#   - optparse script
########################################################################################################

print("BEGIN DECONPlot script")

suppressMessages(library(R.utils))
suppressMessages(library(ExomeDepth))
suppressMessages(library(optparse))
suppressMessages(library(grid))
suppressMessages(library(reshape))
suppressMessages(library(ggplot2))

###### Parsing input options and setting defaults ########
option_list<-list(
    make_option('--rdata',help='Input summary RData file (required)',dest='data'),
    make_option('--out',default='./plots',help='Output directory, default=./plots',dest='plotFolder')
)
opt<-parse_args(OptionParser(option_list=option_list))

# Load R workspace with all the results save : save(ExomeCount,bed.file,counts,fasta,sample.names,bams,cnv.calls,cnv.calls_ids,refs,models,exon_numbers,exons)
count_data=opt$data
if(count_data=="NULL"){count_data=NULL}
if(is.null(count_data)){
    print("ERROR: no RData summary file provided -- Execution halted")
    quit()
}
load(count_data)

plotFolder=opt$plotFolder
if(!file.exists(plotFolder)){dir.create(plotFolder)}

#################### Plots ####################################
message('Start generating plots')

if ("Custom.first" %in% colnames(cnv.calls_ids)){
    cnv.calls_plot=cnv.calls_ids[!is.na(cnv.calls_ids$Custom.first),]
}else{
    cnv.calls_plot=cnv.calls_ids
}

cnv.calls_plot$chr=paste('chr',cnv.calls_plot$chromosome,sep='')
Index=vector(length=nrow(bed.file))
Index[1]=1
for(i in 2:nrow(bed.file)){
    if(bed.file[i,4]==bed.file[i-1,4]){
        Index[i]=Index[i-1]+1
    }else{
        Index[i]=1
    }
}

if(colnames(counts)[5]=="exon"){
    for(i in 1:nrow(exons)){
        x=which(paste(bed.file[,4])==paste(exons[i,4]) & bed.file[,2]<=exons[i,3] & bed.file[,3]>=exons[i,2])
        Index[x]=exons[i,5]
    }
}

for(call_index in 1:nrow(cnv.calls_plot)){
    Sample<-cnv.calls_plot[call_index,]$sample
    Gene<-unlist(strsplit(cnv.calls_plot[call_index,]$Gene,split=", "))
    exonRange<-which(bed.file[,4]%in%Gene)
    if((cnv.calls_plot[call_index,]$start.p-5)<min(exonRange) & ((cnv.calls_plot[call_index,]$start.p-5)>=1)){exonRange=(cnv.calls_plot[call_index,]$start.p-5):max(exonRange)}
    if((cnv.calls_plot[call_index,]$start.p-5)<min(exonRange) & ((cnv.calls_plot[call_index,]$start.p-5)<=0)){exonRange=1:max(exonRange)}
    if((cnv.calls_plot[call_index,]$end.p+5)>max(exonRange) & ((cnv.calls_plot[call_index,]$end.p+5)<=nrow(bed.file))){exonRange=min(exonRange):(cnv.calls_plot[call_index,]$end.p+5)}
    if((cnv.calls_plot[call_index,]$end.p+5)>max(exonRange) & ((cnv.calls_plot[call_index,]$end.p+5)>nrow(bed.file))){exonRange=min(exonRange):nrow(bed.file)}
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
    VariantExon<- unlist(mapply(function(x,y)x:y,cnv.calls[cnv.calls$sample==Sample,]$start.p,cnv.calls[cnv.calls$sample==Sample,]$end.p))
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
temp = cnv.calls[cnv.calls$sample==Sample,]
    if(sum(temp$start.p%in%exonRange |temp$end.p%in%exonRange)>0){
        temp = temp[temp$start.p%in%exonRange|temp$end.p%in%exonRange,]
        for(i in 1:nrow(temp)){
            start.temp = temp[i,]$start.p
            end.temp = temp[i,]$end.p
            CIPlot<-CIPlot + geom_point(data=CIData[CIData$Exon%in%start.temp:end.temp,], aes(x=Exon,y=Ratio),color="red",cex=3.5)
        }
    }

    if(!singlechr){
        if(prev){CIPlot<- CIPlot + scale_x_continuous(breaks=(min(exonRange)-6):max(exonRange),labels=c(rep("",6),paste(Index[exonRange])),limits=c(min(exonRange)-6.75,max(exonRange)))
        }else{CIPlot<-CIPlot + scale_x_continuous(breaks=min(exonRange):(max(exonRange)+6),labels=c(paste(Index[exonRange]),rep("",6)),limits=c(min(exonRange),max(exonRange)+6.75))}
    }else{CIPlot<-CIPlot + scale_x_continuous(breaks=exonRange,labels=paste(Index[exonRange]))}

    ######### Save plot in pdf format ###########
    cnv_genes_sample=cnv.calls_plot[cnv.calls_plot$sample==Sample,]$Gene
    if(sum(cnv_genes_sample==Gene)==1){
        pdf(file=paste(plotFolder,"/",Sample,"_",Gene,".pdf",sep=""),useDingbats=FALSE)
    }else{
        cnv_genes_sample_index=which(cnv.calls_plot$sample==Sample & cnv.calls_plot$Gene==paste(Gene,collapse=", "))
        pdf(file=paste(plotFolder,"/",Sample,"_",paste(Gene,collapse="_"),"_",which(cnv_genes_sample_index==call_index),".pdf",sep=""),useDingbats=F)
    }

    if(sum(cnv_genes_sample==Gene)>2){print(paste("WARNING: more than 2 calls in ",Gene,", could affect plotting",sep=""))}

    grid.newpage()
    pushViewport(viewport(layout = grid.layout(6, 1)))
    print(A1,vp=viewport(layout.pos.row=1:3,layout.pos.col=1))
    print(GenesPlot,vp=viewport(layout.pos.row=4,layout.pos.col=1))
    print(CIPlot,vp=viewport(layout.pos.row=5:6,layout.pos.col=1))
    dev.off()
}
warnings()
print("END DeconPlot script")