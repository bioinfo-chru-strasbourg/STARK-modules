#!/bin/csh
#		SIFT.csh
# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software
#

# Argument 1: the protein sequence file in fasta format
# Argument 2: the pathname to the protein sequence database
# Argument 3: the substitution file (file containing amino acid substitutions to be predicted on. 

### Set these for your installation
#	Location of blastpgp
setenv NCBI /usr/local/packages/blast-2.2.18/bin/


#	Location of BLIMPS
setenv BLIMPS_DIR /opt/www/sift/sift3.0/blimps/

#	Location of SIFT 
setenv SIFT_DIR /opt/www/sift/sift3.0

#	SIFT's output files are written here
set tmpdir = /opt/www/sift/sift3.0/tmp

### Shouldn't need to make any more changes, look for output in $tmpsift
set bindir = $SIFT_DIR/bin
set root_file = $1:r
set tail_of_root = $root_file:t
set tmpfasta = $tmpdir/$tail_of_root.alignedfasta
set tmpsift = $tmpdir/$tail_of_root.SIFTprediction

$bindir/seqs_chosen_via_median_info.csh $1 $2 $4
$bindir/info_on_seqs $tmpfasta $3 $tmpsift
echo "Output in $tmpsift"

exit
