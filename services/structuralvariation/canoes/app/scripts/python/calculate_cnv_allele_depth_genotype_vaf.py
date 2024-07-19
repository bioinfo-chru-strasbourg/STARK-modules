##########################################################################
# calculate_cnv_allele_depth_genotype_vaf.py       	Version: 1.0
# Description:          							Python script to add AD, VAF and GT to a vcf
##########################################################################

########## Note ########################################################################################
# DEV v1 19/07/2024
# need testing EXPERIMENTAL !
########################################################################################################

import pysam
import argparse

def calculate_cnv_allele_depth_genotype_vaf(bamfile, vcffile, output_vcf, hom_threshold=0.2, het_threshold=0.8):
	bam = pysam.AlignmentFile(bamfile, "rb")
	vcf = pysam.VariantFile(vcffile, "r")
	
	# Add VAF to the VCF header
	if 'VAF' not in vcf.header.formats:
		vcf.header.formats.add("VAF", "1", "Float", "Variant Allele Frequency")
	
	output = pysam.VariantFile(output_vcf, "w", header=vcf.header)
	
	for record in vcf:
		chrom = record.chrom
		start = record.pos - 1  # Convert to 0-based for pysam
		end = record.stop  # VCF stop position is inclusive, but pysam is 0-based exclusive
		
		# Initialize read counts
		ref_count = 0
		alt_count = 0
		
		# Iterate through pileup columns in the CNV region
		for pileupcolumn in bam.pileup(chrom, start, end, truncate=True):
			for pileupread in pileupcolumn.pileups:
				if not pileupread.is_del and not pileupread.is_refskip:
					base = pileupread.alignment.query_sequence[pileupread.query_position]
					if base == record.ref:
						ref_count += 1
					elif base in record.alts:
						alt_count += 1
		
		# Total read count
		total_count = ref_count + alt_count
		
		# Calculate average depth of coverage across the CNV region
		region_length = end - start
		avg_depth = total_count / region_length if region_length > 0 else 0
		
		# Calculate VAF
		vaf = alt_count / total_count if total_count > 0 else 0
		
		# Update AD field
		for sample in record.samples:
			record.samples[sample]['AD'] = [ref_count, alt_count]
			
			# Determine genotype based on depth thresholds
			if avg_depth <= hom_threshold:
				record.samples[sample]['GT'] = (0, 0)
			elif avg_depth >= het_threshold:
				record.samples[sample]['GT'] = (1, 1)
			else:
				record.samples[sample]['GT'] = (0, 1)
			
			# Add VAF field
			record.samples[sample]['VAF'] = vaf
		
		output.write(record)
	
	bam.close()
	vcf.close()
	output.close()

def main():
	parser = argparse.ArgumentParser(description='Calculate AD, GT, and VAF for CNVs in a VCF file.')
	parser.add_argument('-b', '--bamfile', required=True, help='Input BAM file')
	parser.add_argument('-v', '--vcffile', required=True, help='Input VCF file')
	parser.add_argument('-o', '--output_vcf', required=True, help='Output VCF file with updated AD, GT, and VAF fields')
	parser.add_argument('--hom_threshold', type=float, default=0.2, help='Threshold for homozygous reference (default: 0.2)')
	parser.add_argument('--het_threshold', type=float, default=0.8, help='Threshold for homozygous alternate (default: 0.8)')
	
	args = parser.parse_args()
	
	calculate_cnv_allele_depth_genotype_vaf(
		bamfile=args.bamfile,
		vcffile=args.vcffile,
		output_vcf=args.output_vcf,
		hom_threshold=args.hom_threshold,
		het_threshold=args.het_threshold
	)

if __name__ == "__main__":
	main()