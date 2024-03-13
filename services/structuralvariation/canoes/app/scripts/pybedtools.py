from multiprocessing import Pool
import argparse
import pandas as pd
from inspect import getmembers, isfunction
import time
import pysam
import subprocess
import datetime
import os

def systemcall(command):
	"""
	Passing command to the shell, return list containing stdout lines
	"""
	p = subprocess.Popen([command], stdout=subprocess.PIPE, shell=True)
	return p.stdout.read().decode('utf8').strip().split('\n')

def split_bed(bedfile):
	list_df = []
	df = pd.read_csv(bedfile, sep='\t', header=None)
	grouped = df.groupby(df.iloc[:,0])
	for values in df.iloc[:,0].unique():
		list_df.append(grouped.get_group(values))
	return list_df

def multicov(bamfile, bedfile, mq):
	bam = pysam.AlignmentFile(bamfile, "rb")
	for j, rows in bedfile.iterrows():
		regions = rows.tolist()[0:3]
		reads_list = []
		for read in bam.fetch(str(regions[0]), regions[1], regions[2]):
			if read.mapping_quality > mq:
				reads_list.append(read)
		bedfile.loc[j, 6] = str(len(reads_list))

	return bedfile

def main():
	args = parseargs()
	t = time.process_time()
	chr_split = split_bed(args.bed)
	output = args.output+".tmp"
	list_cov = [(args.bam, chr, args.mq) for chr in chr_split]

	with Pool(args.threads) as pool:
		result = pool.starmap(multicov, list_cov)

	df_final = pd.concat(result)
	df_final.iloc[:, -1] = df_final.iloc[:, -1].apply(lambda x: str(x).split('.')[0])
	df_final.to_csv(output, sep="\t", index=False, header=False)
	print("#[INFO] Create unsorted multicov ")
	if os.path.exists(output):
		print("#[INFO] Sort output file ...")
		systemcall("sort -k1,1V -k2,2n "+output+" > "+args.output)
		os.remove(output)
	
def parseargs():
	'''
	parse arguments
	'''
	parser = argparse.ArgumentParser()
	parser.add_argument("-b", "--bam", type=str, help="Absolute path of alignment file, .bam")
	parser.add_argument("-e", "--bed", type=str, help="Absolute path of the intervals file, .bed")
	parser.add_argument("-q","--mq", type=int, help="Mapping Quality of read to consider")
	parser.add_argument("-o", "--output", type=str, help="Absolute path of the output file conf Snakefile")
	parser.add_argument("-t", "--threads", type=int, default=6, help="Threads to use, default 6")
	args = parser.parse_args()
	return args


if __name__ == '__main__':
	main()