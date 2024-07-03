##########################################################################
# merge_multicov.py     Version: 2.0
# Description:          Python script to merge coverage tsv
##########################################################################

########## Note ########################################################################################
# DEV v1 27/06/2024
# Changelog
#   - clean & simplify code

########################################################################################################

import argparse
import pandas as pd

def parse_arguments():
	parser = argparse.ArgumentParser()
	parser.add_argument("input", nargs="+", help="Input multicov files")
	parser.add_argument("-o", "--output", required=True, help="Output merged tsv file")
	return parser.parse_args()

def main():
	args = parse_arguments()
	data = pd.DataFrame(columns=["chrom", "start", "end"])
	
	for multicovfile in args.input:
		df = pd.read_csv(multicovfile, sep="\t")
		data = pd.merge(data, df, on=["chrom", "start", "end"], how="outer", validate="one_to_one")
	
	data.fillna("NA", inplace=True)
	data.to_csv(args.output, sep="\t", index=False)

if __name__ == "__main__":
	main()