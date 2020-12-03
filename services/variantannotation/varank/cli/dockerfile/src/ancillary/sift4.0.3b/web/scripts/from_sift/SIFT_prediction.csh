#!/bin/csh
#		SIFT.csh

# Argument 1: the protein sequence file in fasta format
# Argument 2: the pathname to the protein sequence database
# Argument 3: the substitution file (file containing amino acid substitutions to be predicted on. 

setenv BLASTFILTER /usr/local/filter
setenv NCBI /usr/local/bin
setenv BLIMPS_DIR /usr/local/projects/GOSII/pkumar/sift3.0/blimps-3.9
set bindir = "/usr/local/projects/GOSII/pkumar/sift3.0/bin/"
#       Location of SIFT
#setenv SIFT_DIR /usr/local/projects/GOSII/pkumar/sift3.0/


set query_loc = $1
set subst = $2
set outdir = $3
set pid = $4
set tmpdir = $outdir/$pid"_result"

#	SIFT's output files are written here
set alignment = $tmpdir/$pid.selected.alignedfasta
 

$bindir/info_on_seqs $alignment $subst $tmpdir/$pid.prediction
 

exit
