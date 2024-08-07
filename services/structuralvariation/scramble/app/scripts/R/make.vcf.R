##########################################################################
# SCRAMBLE Rscript         Version: 2.0
# Description:             R script to make a vcf file from the scramble df
##########################################################################

################## Context ##################
# Original R script from https://github.com/GeneDx/scramble
# script is located in /opt/conda/share/scramble/bin/make.vcf.R
########## Note ########################################################################################
# DEV v1 11/07/2024
# Changelog
#   - refactor INFO fields, add all the values computed into the vcf (Number of supporting reads)
########################################################################################################

suppressPackageStartupMessages(library(stringr))

get_score = function(right_score, left_score){
  if(is.na(right_score)){
    return(left_score)
  }else if(is.na(left_score)){
    return(right_score)
  }else{
    return(mean(c(left_score, right_score)))
  }
}

get_refs = function(fa, chrom, start, end){
  if (missing(fa) | missing(chrom) | missing(start) | missing(end)) return('N')
  if (! chrom %in% names(fa)) return('N')
  fa = fa[chrom]
  seq = subseq(fa, start=start, end=end)
  return(as.vector(seq))
}

generate_contig_headers <- function(fa) {
  if (missing(fa) || length(names(fa)) == 0) return(NULL)
  
  contigs <- names(fa)
  contig_lengths <- width(fa)
  
  contig_headers <- sapply(seq_along(contigs), function(i) {
    paste('##contig=<ID=chr', contigs[i], ',length=', contig_lengths[i], '>', sep="")
  })
  
  return(contig_headers)
}

make.vcf.header = function(fa, blastRef=NULL){
  if (missing(fa)) return(NULL)
  contig_headers <- generate_contig_headers(fa)
	header <- c('##fileformat=VCFv4.3',
						 paste('##reference=', blastRef, sep=''),
						 contig_headers,
						 '##FILTER=<ID=PASS,Description="All filters passed">',
						 '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
						 '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
						 '##INFO=<ID=END,Number=1,Type=Integer,Description="End position for structural variants">',
						 '##INFO=<ID=MEINFO,Number=.,Type=String,Description="Mobile element info of the form CHR:NAME_START_DIRECTION">',
             '##INFO=<ID=STRAND,Number=.,Type=String,Description="Strand of the mobile element insertion; can be + or -">',
             '##INFO=<ID=COUNTS,Number=.,Type=Integer,Description="Number of supporting reads for the MEI">',
						 '##INFO=<ID=CLIPPED_READS_IN_CLUSTER,Number=.,Type=String,Description="Number of supporting reads in cluster">',
						 '##INFO=<ID=ALIGNMENT_PERCENT_LENGHT,Number=.,Type=Float,Description="Percent of clipped read consensus sequence involved in alignment to MEI reference sequence">',
						 '##INFO=<ID=ALIGNMENT_PERCENT_IDENTITY,Number=.,Type=Float,Description="Percent identify of alignment of clipped read consensus sequence with MEI reference sequence">',
						 '##INFO=<ID=CLIPPED_SEQUENCE,Number=.,Type=String,Description="Clipped cluster consensus sequences">',
						 '##INFO=<ID=CLIPPED_SIDE,Number=.,Type=String,Description="Left or right, side of read where soft-clipping ocurred">',
						 '##INFO=<ID=Start_In_MEI,Number=.,Type=Integer,Description="Left-most position of alignment to MEI reference sequence">',
						 '##INFO=<ID=Stop_In_MEI,Number=.,Type=Integer,Description="Right-most position of alignment to MEI reference sequence">',
						 '##INFO=<ID=polyA_Position,Number=.,Type=String,Description="Position of polyA clipped read cluster if found">',
						 '##INFO=<ID=polyA_Seq,Number=.,Type=String,Description="Clipped cluster consensus sequences of polyA clipped read cluster if found">',
						 '##INFO=<ID=polyA_SupportingReads,Number=.,Type=String,Description="Number of supporting reads in polyA clipped read cluster if found">',
						 '##INFO=<ID=TSD,Number=.,Type=String,Description="Target site duplication sequence if polyA clipped read cluster found">',
						 '##INFO=<ID=TSD_length,Number=.,Type=String,Description="Length of target site duplication if polyA clipped read cluster found">',
						 '##INFO=<ID=REF_ANCHOR_BASE,Number=.,Type=Integer,Description="Reference based at deletion start">',
						 '##INFO=<ID=RIGHT_CLUSTER,Number=.,Type=Integer,Description="Name of right cluster">',
						 '##INFO=<ID=RIGHT_CLUSTER_COUNTS,Number=.,Type=Integer,Description="Number of supporting reads in right cluster">',
						 '##INFO=<ID=LEFT_CLUSTER,Number=.,Type=Integer,Description="Name of left cluster">',
						 '##INFO=<ID=LEFT_CLUSTER_COUNTS,Number=.,Type=Integer,Description="Number of supporting reads in left cluster">',
						 '##INFO=<ID=LEN_RIGHT_ALIGNMENT,Number=.,Type=Integer,Description="Length of right-clipped consensus sequence involved in alignment">',
						 '##INFO=<ID=SCORE_RIGHT_ALIGNMENT,Number=.,Type=Float,Description="BLAST alignment bitscore for right-clipped consensus">',
						 '##INFO=<ID=PCT_COV_RIGHT_ALIGNMENT,Number=.,Type=Float,Description="Percent length of right-clipped consensus involved in alignment">',
						 '##INFO=<ID=PCT_IDENTITY_RIGHT_ALIGNMENT,Number=.,Type=Float,Description="Percent identity of right-clipped consensus in alignment">',
						 '##INFO=<ID=LEN_LEFT_ALIGNMENT,Number=.,Type=Integer,Description="Length of left-clipped consensus sequence involved in alignment">',
						 '##INFO=<ID=SCORE_LEFT_ALIGNMENT,Number=.,Type=Float,Description="BLAST alignment bitscore for left-clipped consensus">',
						 '##INFO=<ID=PCT_COV_LEFT_ALIGNMENT,Number=.,Type=Float,Description="Percent length of left-clipped consensus involved in alignment">',
						 '##INFO=<ID=PCT_IDENTITY_LEFT_ALIGNMENT,Number=.,Type=Float,Description="Percent identity of right-clipped consensus in alignment">',
						 '##INFO=<ID=INS_SIZE,Number=.,Type=Integer,Description="Length of insert within deleted sequence (for two-end deletions only)">',
						 '##INFO=<ID=INS_SEQ,Number=.,Type=String,Description="Inserted sequence (for two-end deletions only)">',
						 '##INFO=<ID=RIGHT_CLIPPED_SEQ,Number=.,Type=String,Description="Clipped consensus sequence for right-clipped cluster">',
						 '##INFO=<ID=LEFT_CLIPPED_SEQ,Number=.,Type=String,Description="Clipped consensus sequence for left-clipped cluster">',
             '##INFO=<ID=VAF,Number=.,Type=Float,Description="Variant Allele Frequency">',
						 '##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">',
						 '##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">',
						 '##ALT=<ID=INS:ME:SVA,Description="Insertion of SVA element">'
						 )
	return(header)
}

write.scramble.vcf = function(winners, fa, meis = FALSE) {
if (is.null(winners) || nrow(winners) == 0) {
    fixed = data.frame('#CHROM' = character(),
                       POS = character(),
                       ID = character(),
                       REF = character(),
                       ALT = character(),
                       QUAL = character(),
                       FILTER = character(),
                       INFO = character(),
                       stringsAsFactors = FALSE, check.names = FALSE)
    return(fixed)
  }

  if(deletions) {

    fixed = data.frame('#CHROM' = winners$CONTIG,
                       POS = winners$DEL.START,
                       ID = 'DEL',
                       QUAL = sapply(1:nrow(winners), function(i) get_score(winners$SCORE.RIGHT.ALIGNMENT[i], winners$SCORE.LEFT.ALIGNMENT[i])),
                       FILTER = 'PASS',
                       svtype = 'DEL',
                       stringsAsFactors = F, check.names = F)
    fixed$ALT = '<DEL>'
    fixed$svlen = -nchar(sapply(1:nrow(winners), function(i) {
            get_refs(fa, winners$CONTIG[i], winners$DEL.START[i], winners$DEL.END[i] + 1)
            }))
    fixed$end = fixed$POS + fixed$svlen
    fixed$INFO = paste0('SVTYPE=', fixed$svtype, ';', 'SVLEN=', -fixed$svlen, ';', 'END=', fixed$end, ';', 'REF_ANCHOR_BASE=', winners$REF.ANCHOR.BASE, ';', 'RIGHT_CLUSTER=', winners$RIGHT.CLUSTER, ';', 'RIGHT_CLUSTER_COUNTS=', winners$RIGHT.CLUSTER.COUNTS, ';', 'LEFT_CLUSTER=', winners$LEFT.CLUSTER, ';', 'LEFT_CLUSTER_COUNTS=', winners$LEFT.CLUSTER.COUNTS, ';', 'LEN_RIGHT_ALIGNMENT=', winners$LEN.RIGHT.ALIGNMENT, ';', 'SCORE_RIGHT_ALIGNMENT=', winners$SCORE.RIGHT.ALIGNMENT, ';', 'PCT_COV_RIGHT_ALIGNMENT=', winners$PCT.COV.RIGHT.ALIGNMENT, ';', 'PCT_IDENTITY_RIGHT_ALIGNMENT=', winners$PCT.IDENTITY.RIGHT.ALIGNMENT, ';', 'LEN_LEFT_ALIGNMENT=', winners$LEN.LEFT.ALIGNMENT, ';', 'SCORE_LEFT_ALIGNMENT=', winners$SCORE.LEFT.ALIGNMENT, ';', 'PCT_COV_LEFT_ALIGNMENT=', winners$PCT.COV.LEFT.ALIGNMENT, ';', 'PCT_IDENTITY_LEFT_ALIGNMENT=', winners$PCT.IDENTITY.LEFT.ALIGNMENT, ';', 'INS_SIZE=', winners$INS.SIZE, ';', 'RIGHT_CLIPPED_SEQ=', winners$RIGHT.CLIPPED.SEQ, ';', 'LEFT_CLIPPED_SEQ=', winners$LEFT.CLIPPED.SEQ)
    fixed$REF = sapply(1:nrow(fixed), function(i) get_refs(fa, fixed[i, '#CHROM'], fixed$POS[i], fixed$POS[i]))
  } 
  
  if (meis) {
    fixed = data.frame('#CHROM' =  gsub("(.*):(\\d*)$", "\\1", winners$Insertion),
                       POS = as.integer(gsub("(.*):(\\d*)$", "\\2", winners$Insertion)),
                       ID = 'INS:ME',
                       FILTER = 'PASS',
                       ALT = paste('<INS:ME:', toupper(winners$MEI_Family), '>', sep=''),
                       QUAL = winners$Alignment_Score,
                       name = paste(winners$Insertion, toupper(winners$MEI_Family), winners$Insertion_Direction, sep="_"),
                       polarity = ifelse(winners$Insertion_Direction == 'Plus', "+", "-"),
                       stringsAsFactors = F, check.names = F)
    fixed$start = fixed$POS
    fixed$INFO = paste0('MEINFO=', fixed$name, ';' , 'STRAND=', fixed$polarity, ';', 'COUNTS=', winners$Clipped_Reads_In_Cluster, ';','ALIGNMENT_PERCENT_LENGHT=' , winners$Alignment_Percent_Length , ';', 'ALIGNMENT_PERCENT_IDENTITY=', winners$Alignment_Percent_Identity, ';','CLIPPED_SEQUENCE=', winners$Clipped_Sequence, ';', 'CLIPPED_SIDE=', winners$Clipped_Side, ';', 'Start_In_MEI=', winners$Start_In_MEI, ';',  'Stop_In_MEI=', winners$Stop_In_MEI, ';', 'polyA_Position=',  winners$polyA_Position, ';', 'polyA_Seq=', winners$polyA_Seq, ';', 'polyA_SupportingReads=', winners$polyA_SupportingReads, ';', 'TSD=', winners$TSD, ';' , 'TSD_length=', winners$TSD_length)
    fixed$REF = sapply(1:nrow(fixed), function(i) get_refs(fa, fixed[i, '#CHROM'], fixed$POS[i], fixed$POS[i]))
  }   

  vcf.cols = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')
  return(fixed[,vcf.cols])

}