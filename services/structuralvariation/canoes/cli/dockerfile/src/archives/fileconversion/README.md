# FileConversion

It is a python script that allows the conversion between different file formats.
There are four possible conversions : tsv d'AnnotSV to vcf, tsv d'AnnotSV to bed, bed to vcf, vcf to bed

## Installation

This script works with python 3.
You can install python 3 like this :
`apt-get install python3`

Then, you have to install the pip command:
`apt-get install python3-pip`
`pip3 install --upgrade pip`

It is then necessary to install pandas
`pip3 install pandas`

## How to use the script ?
The script takes different arguments as input:

Required arguments
- `-c or --convert` : Conversion perform by the script (`tsv2vcf/tsv2bed/vcf2bed/bed2vcf`)
- `-i or --infile` : Input file
- `-o or --outfile` : Output file

Optional arguments (according to the type of conversion)
- `-p or --profil` : Profile of tsv file (`bed/AnnotSV`)
- `-col or --column` : Column numbers corresponding to the chromosome, start, end, type of variation and genotype columns (<u>in that order</u>)
- `-ref or --reference` : Fasta file of the reference genome
- `-v or --version` : Version of the AnnotSV file (`default: '3.0'`)
- `-d or --default` : Score and Strand by default (`default: False`)

Depending on the type of conversion (--convert), some optional arguments become mandatory. Below is the syntax for each conversion :

- Conversion of AnnotSV to vcf : `python3 file_conversion.py -c tsv2vcf -p AnnotSV -v 3.0 -i <path_to_input> -o <path_to_output> -r <path_to_reference>`
- Conversion of AnnotSV to bed : `python3 file_conversion.py -c tsv2bed -p AnnotSV -v 3.0 -i <path_to_input> -o <path_to_output> -d True`
- Conversion of bed to vcf :
  - With header : `python3 file_conversion.py -c tsv2vcf -p bed -i <path_to_input> -o <path_to_output> -r <path_to_reference>`
  - Without header : `python3 file_conversion.py -c tsv2vcf -p bed -i <path_to_input> -col 0,1,2,3,4 -o <path_to_output> -r <path_to_reference>`
- Conversion vcf to bed : `python3 file_conversion.py -c vcf2bed -i <path_to_input> -o <path_to_output> -d True`

We will detail each conversion

### AnnotSV to vcf

This conversion supports AnnotSV version 3.0 and 2.3. To take into account future versions, it will be enough to add a configuration file in .txt format in the script folder. The output file is a VCF4.3 file.
The output file is a VCF file which keeps all the information of the annotation. The genotype of the samples is also created.  

<u>**IMPORTANT</u>: The genotype in the output file is 0/0 when the variation is not present, and by convention 1/1 when it is present. This does not mean that the patient is homozygous, only that he/she is a carrier of the variation (without knowing whether he/she is homozygous or heterozygous)**


### AnnotSV to bed

This conversion supports AnnotSV version 3.0 and 2.3. To take into account future versions, it will be enough to add a configuration file in .txt format in the script folder. The output file is a bed file with 4 columns (chromosome, start, end, name) or 6 columns (chromosome, start, end, name, score, strand) id default argument are activate.

### bed to vcf

This conversion transforms a bed file. The output file is a VCF4.3 (chromosome, start, end, name).
The conversion requires a header or column information given when the script is called.

### vcf to bed

This conversion transforms a VCF4.3 file. The output file is a bed file with 4 columns (chromosome, start, end, name)  or 6 columns (chromosome, start, end, name, score, strand) id default argument are activate.
