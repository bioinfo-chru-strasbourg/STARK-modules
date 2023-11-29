#! /usr/bin/env python
# -*- coding: utf-8 -*-
###Author : RAUCH Mateusz###
###Maintainer : LAMOUCHE Jean-Baptiste###

# changelog add gene columns TL

from __future__ import division
from __future__ import print_function

import argparse
import os
import re
import time
import subprocess
import sys


def main(args):
	bedFile = args.input
	kmerSize = args.kmer
	if args.output_folder:
		if not os.path.exists(args.output):
			os.mkdir(args.output)
		elif os.path.isfile(args.output):
			print("ERROR "+args.output+" exists and is a file")
			exit()
		kmerBedFile = os.path.join(args.output, os.path.basename(bedFile)[:-4] + '_' + str(kmerSize) + '_kmer.bed')
	elif args.output:
		kmerBedFile = args.output
	else:
		kmerBedFile = bedFile[:-4] + '_' + str(kmerSize) + '_kmer.bed'
	check_bed(bedFile)

	count = 0
	print(str(kmerSize) + ' k-merisation of your bed ' + bedFile + ' in progress...')
	with open(bedFile, 'r') as readBed:
		with open(kmerBedFile, 'w+') as writeKbed:
			for line in readBed.readlines():
				count += 1
				if line.startswith('#'):
					continue
				else:
					line = line.split()
					diff = int(line[2])-int(line[1])
					start = line[1]
					end = line[2]
					gene = line[3]
					while diff >= kmerSize:
						newEnd = int(start) + int(kmerSize) - 1
						newLine = line[0] + '\t' + str(start) + '\t' + str(newEnd) + '\t' + gene
						diff = diff - kmerSize
						start = int(newEnd) + 1
						writeKbed.write(newLine + '\n')
					if diff > 0:
						newLine = line[0] + '\t' + str(start) + '\t' + str(end) + '\t' + gene
						writeKbed.write(newLine + '\n')

	sort_bed(kmerBedFile)
	print('Done !')


def sort_bed(kmerBedFile):
	print('Sorting your magnificent kbed...')
	tmpKbedSorted = kmerBedFile + '.tmp'

	with open(tmpKbedSorted, 'w') as writeSortedKbed:
		subprocess.call(['sort', kmerBedFile, '-k1,1V', '-k2,2n'], stdout=writeSortedKbed, stderr=subprocess.STDOUT, universal_newlines=True)
	os.remove(kmerBedFile)
	os.rename(tmpKbedSorted, kmerBedFile)


def check_bed(bedFile):
	count = 0

	with open(bedFile, 'r') as readBed:
		for line in readBed.readlines():
			count += 1
			if line.startswith('#'):
					continue
			else:
				line = line.split()
				if not re.match(r'chr[0-9XYxy]+', line[0]):
					print('[ERROR] The first column must be a valid chromosome format, exiting')
					sys.exit()
				elif not re.match(r'[0-9]+', line[1]):
					print('[ERROR] The second column must be the start region of your design, exiting')
					sys.exit()
				elif not re.match(r'[0-9]+', line[2]):
					print('[ERROR] The third column must be the end region of your design, exiting')
					sys.exit()
				elif int(line[2])-int(line[1]) < 0:
					print('[ERROR] The end and start value of your regions may be inverted at line ' + str(count) + ' please check your bedfile, exiting')


def check_args(args):
	if not os.path.isfile(args.input):
		print('[ERROR] The input file ' + args.input + ' doesn\'t exist, please check your input argument, exiting')
		sys.exit()

	if not args.input.endswith('.bed'):
		print('[ERROR] The input file ' + args.input + ' don\'t seems to be a .bed extension, please provide a file with correct extension, exiting')
		sys.exit()


if __name__ == '__main__':
	startTime = time.time()
	parser = argparse.ArgumentParser(description='kmerization script')
	parser.add_argument('-i', '--input', type=str, required=True, help='path to your original bed file')
	parser.add_argument('-k', '--kmer', type=int, required=True, help='length of kmer you want in your bed file')
	parser.add_argument('-f', '--output_folder', type=str, help='optionnal output folder, by default kmerized bed go to the same folder as original bed file')
	parser.add_argument('-o', '--output', type=str, help='Absolute path and name of the output bedfile .bed')
	args = parser.parse_args()
	check_args(args)
	main(args)
	print("--- %s seconds ---" % (time.time() - startTime))