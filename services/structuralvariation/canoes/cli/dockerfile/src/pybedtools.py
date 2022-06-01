#params.bedtools}/bedtools multicov -bams {input.bam} -bed {input.bed} -q 20 > {output} 2> {log}


from multiprocessing import Pool
import argparse
import pandas as pd
from inspect import getmembers, isfunction
import time
import pysam
import subprocess
import datetime
import os



#prepare args for multiprocessing
#print("#[INFO] Parse variants by chr")
#list_cov = [(osj(output.vcf+'.dict', covfiles), df, col) for covfiles in os.listdir(output.vcf+'.dict') if not covfiles.startswith('COVFILELocus')]
#print("#[INFO] chr nbr", len(list_cov))
#with Pool(8) as pool:
#	result = pool.starmap(parseSpeed, list_cov)
#	for locus in result:
#		list_var.extend(locus)
#bedfile = "SGT186075.bed"
#bamfile = "SGT186075.bwamem.bam"



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
		#print("#[INFO] Chromosomes "+values)
		list_df.append(grouped.get_group(values))
		#print("#[INFO] Number of regions ", len(grouped.get_group(values)))
	return list_df

def multicov(bamfile, bedfile, mq):
	bam = pysam.AlignmentFile(bamfile, "rb")
	for j, rows in bedfile.iterrows():
		#print("#[INFO] Regions", rows)
		regions = rows.tolist()[0:3]
		#print("#[INFO] Regions ", regions )
		reads_list = []
		for read in bam.fetch(str(regions[0]), regions[1], regions[2]):
			if read.mapping_quality > mq:
				#print(read)
				reads_list.append(read)
		#print("#[INFO] Number of reads ", 
		# len(reads_list))
		bedfile.loc[j, 6] = len(reads_list)
	
	bedfile.iloc[:, -1] = bedfile.iloc[:, -1].astype(int)
	return bedfile

def main():
	args = parseargs()
	t = time.process_time()
	chr_split = split_bed(args.bed)
	output = args.output+".tmp"
	#print(lis[6].head())
	list_cov = [(args.bam, chr, args.mq) for chr in chr_split]

	with Pool(args.threads) as pool:
		result = pool.starmap(multicov, list_cov)
	
	df_final = pd.concat(result)
	df_final.to_csv(output, sep="\t", index=False, header=False)
	print("#[INFO] Create unsorted multicov ")
	if os.path.exists(output):
		print("#[INFO] Sort output file ...")
		systemcall("sort -k1,1V -k2,2n "+output+" > "+args.output)
		os.remove(output)
	
	#time elapsed
	#elapsed = time.process_time() - t
	#hours, seconds = divmod(elapsed * 60, 3600)  # split to hours and seconds
	#minutes, seconds = divmod(seconds, 60)
	#timet = "{:02.0f}:{:02.0f}:{:02.0f}".format(hours, minutes, seconds)

	#print("#[INFO] Process "+os.path.abspath(args.bam)+" in ", timet)

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
	#test()
