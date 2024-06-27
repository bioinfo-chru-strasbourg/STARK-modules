##########################################################################
# pybedtools        	Version: 2.0
# Description:          Python script to extract coverage from a bam file or a list of bams, using a bed file
##########################################################################

########## Note ########################################################################################
# DEV v1 27/06/2024
# Changelog
#   - clean & simplify code, add option to input a tsv file with bam path

########################################################################################################

import argparse
import pandas as pd
import pysam
import os
from multiprocessing import Pool

def split_bed(bedfile):
	df = pd.read_csv(bedfile, sep='\t', header=None)
	return [group for _, group in df.groupby(df.iloc[:, 0])]

def multicov(bamfile, bedfile, mq):
	bam = pysam.AlignmentFile(bamfile, "rb")
	results = []

	for _, row in bedfile.iterrows():
		chrom, start, end = row[0], row[1], row[2]
		count_reads = sum(1 for read in bam.fetch(chrom, start, end) if read.mapping_quality > mq)
		bam_name = os.path.basename(bamfile).split('.')[0]
		results.append({'chrom': chrom, 'start': start, 'end': end, bam_name : count_reads})

	return pd.DataFrame(results)

def read_input_tsv(input_tsv):
	return pd.read_csv(input_tsv, sep='\t')['bam'].tolist()

def multicov_wrapper(bamfile, chr_split, mq):
	return pd.concat([multicov(bamfile, chr_df, mq) for chr_df in chr_split])

def main():
	args = parseargs()
	chr_split = split_bed(args.bed)
	bam_files = read_input_tsv(args.tsv) if args.tsv else [args.bam]
	list_cov = [(bam, chr_split, args.mq) for bam in bam_files]

	with Pool(args.threads) as pool:
		results = pool.starmap(multicov_wrapper, list_cov)

	df_final = pd.concat(results, ignore_index=True)
	df_final.to_csv(args.output, sep="\t", index=False, header=True)

def parseargs():
	parser = argparse.ArgumentParser()
	group = parser.add_mutually_exclusive_group(required=True)
	group.add_argument("-b", "--bam", type=str, help="Absolute path of bam alignment file")
	group.add_argument("-f", "--tsv", type=str, help="Absolute path of the input tsv file containing bam paths")
	
	parser.add_argument("-e", "--bed", type=str, required=True, help="Absolute path of the bed intervals file")
	parser.add_argument("-q", "--mq", type=int, required=True, help="Mapping Quality of read to consider (default None)")
	parser.add_argument("-o", "--output", type=str, required=True, help="Absolute path of the tsv file")
	parser.add_argument("-t", "--threads", type=int, default=6, help="Threads to use (default 6)")

	return parser.parse_args()

if __name__ == '__main__':
	main()