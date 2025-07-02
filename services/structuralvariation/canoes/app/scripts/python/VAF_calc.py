##########################################################################
# calculate_cnv_allele_depth_genotype_vaf.py        Version: 1.0
# Description:                                       Python script to add AD, VAF, and GT to a VCF
##########################################################################

########## Note ########################################################################################
# DEV v1 19/07/2024
# Need testing EXPERIMENTAL!
########################################################################################################

import pysam
import argparse
import gzip
import os
import subprocess
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def calculate_cnv_allele_depth_genotype_vaf(bamfile, vcffile, output_vcf, hom_threshold=0.2, het_threshold=0.8):
    """
    Calculate allele depth (AD), genotype (GT), and variant allele frequency (VAF)
    for CNVs in a VCF file.

    Args:
        bamfile (str): Path to the input BAM file.
        vcffile (str): Path to the input VCF file (gzipped if ending with .gz).
        output_vcf (str): Path for the output VCF file with updated fields.
        hom_threshold (float): Threshold for calling homozygous reference. Default is 0.2.
        het_threshold (float): Threshold for calling homozygous alternate. Default is 0.8.
    """
    try:
        bam = pysam.AlignmentFile(bamfile, "rb")
    except Exception as e:
        logging.error(f"Error opening BAM file: {bamfile}. {e}")
        return

    try:
        vcf = pysam.VariantFile(vcffile, "r")
    except Exception as e:
        logging.error(f"Error opening VCF file: {vcffile}. {e}")
        return

    # Add VAF to the VCF header if not present
    if 'VAF' not in vcf.header.formats:
        vcf.header.formats.add("VAF", "1", "Float", "Variant Allele Frequency")
    
    output = pysam.VariantFile(output_vcf, "w", header=vcf.header)

    for record in vcf:
        chrom = record.chrom
        start = record.pos - 1  # Convert to 0-based for pysam
        end = record.stop  # VCF stop position is inclusive, pysam is 0-based exclusive
        
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
        
        total_count = ref_count + alt_count
        region_length = end - start
        avg_depth = total_count / region_length if region_length > 0 else 0
        vaf = alt_count / total_count if total_count > 0 else 0

        # Update AD, GT, and VAF fields for each sample
        for sample in record.samples:
            record.samples[sample]['AD'] = [ref_count, alt_count]

            # Determine genotype based on depth thresholds
            if avg_depth <= hom_threshold:
                record.samples[sample]['GT'] = (0, 0)  # Homozygous reference
            elif avg_depth >= het_threshold:
                record.samples[sample]['GT'] = (1, 1)  # Homozygous alternate
            else:
                record.samples[sample]['GT'] = (0, 1)  # Heterozygous
            
            # Add VAF field
            record.samples[sample]['VAF'] = vaf
        
        output.write(record)
    
    bam.close()
    vcf.close()
    output.close()
    
    # Compress the output VCF with bgzip
    try:
        subprocess.run(['bgzip', '-f', output_vcf], check=True)
        logging.info(f"Output written to {output_vcf} and compressed successfully.")
    except subprocess.CalledProcessError as e:
        logging.error(f"Error compressing the output VCF: {e}")


def main():
    """
    Main function to parse command-line arguments and call the calculation function.
    """
    parser = argparse.ArgumentParser(description='Calculate AD, GT, and VAF for CNVs in a VCF file.')
    parser.add_argument('-b', '--bamfile', required=True, help='Input BAM file')
    parser.add_argument('-v', '--vcffile', required=True, help='Input VCF file, gz is autodetected')
    parser.add_argument('-o', '--output_vcf', required=True, help='Output VCF gz file with updated AD, GT, and VAF fields')
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
