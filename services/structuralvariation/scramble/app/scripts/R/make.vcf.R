##############################
# convert to VCF
##############################

library(stringr)

get_score = function(right_score, left_score){
	if(is.na(right_score)){
		return(left_score)
	}else if(is.na(left_score)){
		return(right_score)
	}else{
		return(mean(c(left_score, right_score)))
	}
}
##############################
get_refs = function(fa, chrom, start, end){
	if (missing(fa) | missing(chrom) | missing(start) | missing(end)) return('N')
	if (! chrom %in% names(fa)) return('N')
	fa = fa[chrom]
	seq = subseq(fa, start=start, end=end)
	return(as.vector(seq))
}
##############################
make.vcf.header = function(fa, blastRef=None){
	if (missing(fa)) return(NULL)
	contigs = names(fa)
	header = c('##fileformat=VCFv4.3',
						 paste('##reference=', blastRef, sep=''),
						 paste('##contig=<ID=', contigs, '>', sep=""),
						 '##FILTER=<ID=PASS,Description="All filters passed">',
						 '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">',
						 '##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">',
						 '##INFO=<ID=END,Number=.,Type=Integer,Description="End position for structural variants">',
						 '##INFO=<ID=MEINFO,Number=4,Type=String,Description="Mobile element info of the form NAME,START,END,POLARITY">',
						 '##INFO=<ID=CLIPPED_READS_IN_CLUSTER,Number=.,Type=Integer,Description="Number of supporting reads in cluster">',
						 '##INFO=<ID=ALIGNMENT_PERCENT_LENGHT,Number=.,Type=Integer,Description="Percent of clipped read consensus sequence involved in alignment to MEI reference sequence">',
						 '##INFO=<ID=ALIGNMENT_PERCENT_IDENTITY,Number=.,Type=Integer,Description="Percent identify of alignment of clipped read consensus sequence with MEI reference sequence">',
						 '##INFO=<ID=CLIPPED_SEQUENCE,Number=.,Type=Integer,Description="Clipped cluster consensus sequences">',
						 '##INFO=<ID=CLIPPED_SIDE,Number=.,Type=Integer,Description="Left or right, side of read where soft-clipping ocurred">',
						 '##INFO=<ID=Start_In_MEI,Number=.,Type=Integer,Description="Left-most position of alignment to MEI reference sequence">',
						 '##INFO=<ID=Stop_In_MEI,Number=.,Type=Integer,Description="Right-most position of alignment to MEI reference sequence">',
						 '##INFO=<ID=polyA_Position,Number=.,Type=Integer,Description="Position of polyA clipped read cluster if found">',
						 '##INFO=<ID=polyA_Seq,Number=.,Type=Integer,Description="Clipped cluster consensus sequences of polyA clipped read cluster if found">',
						 '##INFO=<ID=polyA_SupportingReads,Number=.,Type=Integer,Description="Number of supporting reads in polyA clipped read cluster if found">',
						 '##INFO=<ID=TSD,Number=.,Type=Integer,Description="Target site duplication sequence if polyA clipped read cluster found">',
						 '##INFO=<ID=TSD_length,Number=.,Type=Integer,Description="Length of target site duplication if polyA clipped read cluster found">',
						 '##INFO=<ID=REF_ANCHOR_BASE,Number=.,Type=Integer,Description="Reference based at deletion start">',
						 '##INFO=<ID=DEL_LENGTH,Number=.,Type=Integer,Description="Deletion length">',
						 '##INFO=<ID=RIGHT_CLUSTER,Number=.,Type=Integer,Description="Name of right cluster">',
						 '##INFO=<ID=RIGHT_CLUSTER_COUNTS,Number=.,Type=Integer,Description="Number of supporting reads in right cluster">',
						 '##INFO=<ID=LEFT_CLUSTER,Number=.,Type=Integer,Description="Name of left cluster">',
						 '##INFO=<ID=LEFT_CLUSTER_COUNTS,Number=.,Type=Integer,Description="Number of supporting reads in left cluster">',
						 '##INFO=<ID=LEN_RIGHT_ALIGNMENT,Number=.,Type=Integer,Description="Length of right-clipped consensus sequence involved in alignment">',
						 '##INFO=<ID=SCORE_RIGHT_ALIGNMENT,Number=.,Type=Integer,Description="BLAST alignment bitscore for right-clipped consensus">',
						 '##INFO=<ID=PCT_COV_RIGHT_ALIGNMENT,Number=.,Type=Integer,Description="Percent length of right-clipped consensus involved in alignment">',
						 '##INFO=<ID=PCT_IDENTITY_RIGHT_ALIGNMENT,Number=.,Type=Integer,Description="Percent identity of right-clipped consensus in alignment">',
						 '##INFO=<ID=LEN_LEFT_ALIGNMENT,Number=.,Type=Integer,Description="Length of left-clipped consensus sequence involved in alignment">',
						 '##INFO=<ID=SCORE_LEFT_ALIGNMENT,Number=.,Type=Integer,Description="BLAST alignment bitscore for left-clipped consensus">',
						 '##INFO=<ID=PCT_COV_LEFT_ALIGNMENT,Number=.,Type=Integer,Description="Percent length of left-clipped consensus involved in alignment">',
						 '##INFO=<ID=PCT_IDENTITY_LEFT_ALIGNMENT,Number=.,Type=Integer,Description="Percent identity of right-clipped consensus in alignment">',
						 '##INFO=<ID=INS_SIZE,Number=.,Type=Integer,Description="Length of insert within deleted sequence (for two-end deletions only)">',
						 '##INFO=<ID=INS_SEQ,Number=.,Type=Integer,Description="Inserted sequence (for two-end deletions only)">',
						 '##INFO=<ID=RIGHT_CLIPPED_SEQ,Number=.,Type=Integer,Description="Clipped consensus sequence for right-clipped cluster">',
						 '##INFO=<ID=LEFT_CLIPPED_SEQ,Number=.,Type=Integer,Description="Clipped consensus sequence for left-clipped cluster">',
						 '##ALT=<ID=INS:ME:ALU,Description="Insertion of ALU element">',
						 '##ALT=<ID=INS:ME:L1,Description="Insertion of L1 element">',
						 '##ALT=<ID=INS:ME:SVA,Description="Insertion of SVA element">'
						 )
	return(header)
}

##############################
write.scramble.vcf = function(winners, fa, meis=F){
 
 # return empty fixed data when no variants found
		if(nrow(winners) == 0){
				fixed = data.frame('#CHROM' = character(),
							 POS = character(),
							 ID = character(),
							 REF = character(),
							 ALT = character(),
							 QUAL = character(),
							 FILTER = character(),
							 INFO = character(),
							 check.names = F)
				return(fixed)
		}
	#argument checks
	if (is.null(winners)) return(NULL)

fixed = data.frame('#CHROM' = ifelse(!meis, winners$CONTIG, gsub("(.*):(\\d*)$", "\\1", winners$Insertion)),
									 POS = ifelse(!meis, winners$DEL.START, as.integer(gsub("(.*):(\\d*)$", "\\2", winners$Insertion))),
									 ID = ifelse(!meis, 'DEL', 'INS:ME'),
									 QUAL = ifelse(!meis, sapply(1:nrow(winners), function(i) get_score(winners$SCORE.RIGHT.ALIGNMENT[i], winners$SCORE.LEFT.ALIGNMENT[i])), winners$Alignment_Score),
									 FILTER = 'PASS',
									 ALT = ifelse(!meis, '<DEL>' , paste('<INS:ME:', toupper(winners$MEI_Family), '>', sep='')),
									 stringsAsFactors = FALSE, check.names = FALSE)

if (!meis) { # for DEL
		fixed$REF = sapply(1:nrow(winners), function(i) get_refs(fa, winners$CONTIG[i], winners$DEL.START[i], winners$DEL.START[i] + 1))
		fixed$svtype = 'DEL'
		fixed$svlen = nchar(fixed$REF)
		fixed$end = fixed$POS + fixed$svlen
		fixed$INFO <- sprintf("SVTYPE=DEL;SVLEN=%d;END=%d;REF_ANCHOR_BASE=%s;DEL_LENGTH=%s;RIGHT_CLUSTER=%s;RIGHT_CLUSTER_COUNTS=%s;LEFT_CLUSTER=%s;LEFT_CLUSTER_COUNTS=%s;LEN_RIGHT_ALIGNMENT=%s;SCORE_RIGHT_ALIGNMENT=%s;PCT_COV_RIGHT_ALIGNMENT=%s;PCT_IDENTITY_RIGHT_ALIGNMENT=%s;LEN_LEFT_ALIGNMENT=%s;SCORE_LEFT_ALIGNMENT=%s;PCT_COV_LEFT_ALIGNMENT=%s;PCT_IDENTITY_LEFT_ALIGNMENT=%s;INS_SIZE=%d;RIGHT_CLIPPED_SEQ=%s;LEFT_CLIPPED_SEQ=%s",
													fixed$svlen, 
													fixed$end, 
													REF_ANCHOR_BASE$alignments, 
													DEL_LENGTH$alignments, 
													RIGHT_CLUSTER$alignments², 
													RIGHT_CLUSTER_COUNTS$alignments, 
													LEFT_CLUSTER$alignments, 
													LEFT_CLUSTER_COUNTS$alignments, 
													LEN_RIGHT_ALIGNMENT$alignments, 
													SCORE_RIGHT_ALIGNMENT$alignments, 
													PCT_COV_RIGHT_ALIGNMENT$alignments, 
													PCT_IDENTITY_RIGHT_ALIGNMENT$alignments, 
													LEN_LEFT_ALIGNMENT$alignments, 
													SCORE_LEFT_ALIGNMENT$alignments, 
													PCT_COV_LEFT_ALIGNMENT$alignments, 
													PCT_IDENTITY_LEFT_ALIGNMENT$alignments, 
													INS.SIZE, 
													RIGHT_CLIPPED_SEQ$alignments, 
													LEFT_CLIPPED_SEQ$alignments)
} else { # for MEI
		fixed$INFO <- sprintf("SVTYPE=INS:ME;MEINFO=%s,%d,%s;CLIPPED_READS_IN_CLUSTER=%s;ALIGNMENT_PERCENT_LENGHT=%s;ALIGNMENT_PERCENT_IDENTITY=%s;CLIPPED_SEQUENCE=%s;CLIPPED_SIDE=%s;Start_In_MEI=%s;Stop_In_MEI=%s;polyA_Position=%s;polyA_Seq=%s;polyA_SupportingReads=%s;TSD=%s;TSD_length=%s",
													paste(fixed$name, fixed$POS, fixed$polarity, sep = ","), 
													winners$counts, 
													winners$percent_clipped_read_aligned, 
													winners$percent_identity, 
													toupper(winners$clipped.consensus), 
													winners$clipped, 
													winners$starts, 
													winners$stops, 
													winners$polyA_Position, 
													winners$polyA_Seq, 
													winners$polyA_SupportingReads, 
													winners$TSD, 
													winners$TSD_length)
		fixed$REF = ifelse(meis, sapply(1:nrow(fixed), function(i) get_refs(fa, fixed[i, '#CHROM'], fixed$POS[i], fixed$POS[i])), "")
}

	vcf.cols = c('#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO')
	return(fixed[,vcf.cols])
	
}