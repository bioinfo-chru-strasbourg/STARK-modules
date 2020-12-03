#!/bin/csh
# SIFTING2.csh
# args are  pid opt seq out info poly ? idfilt db address
#           1   2   3   4   5    6    7 8      9  10

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

# 3/8/04 If $address, return predictions output & put pid in subject

#	Need group write on all files created
umask 002

# Only integer arithmetic is possible in csh.  the intolerance probability threshold
# called in seqs_to_matrixweb is set at 0.05, must change this manually in this script

setenv BLIMPS_DIR /usr/local/blimps
setenv NCBI /usr/local/blast

set tmp = "../tmp"
set commentscsh = $tmp/$1.commentscsh
# no results in commentscsh -- just for debugging
set srcdir = "../src/src_for_csh/established_and_tested_src_for_csh"
set pid = $1 #pid is used to generate a lot of the files 
set option = $2 #seq, related but unaligned sequences, aligned sequences
set output = $4
set info = $5 # sequence median info
set polymorphism_file = $6 # file with substitutions e.g. P2S, SIFT will
			   # print predictions for these. 
set gap_option = 0 # 1 for to account for gaps, 0 to ignore gaps.  permanently set to 0. not modelling gaps 
set exp_option = 1 # 1 for exponential weight, 0 for linear.  use 1.
set seq_identity_filter = $8

# Psiblast parameters
set iterations = 2

# make a .err file.  This will contain 
# errors and warnings that don't stop the programs 
# from running (output from stderr) 
# in contrast to .error files
cp /dev/null $tmp/$pid.err

if ($option == query_seq) then
	set queryseq = $3 
else if ($option == related_seqs) then
	set unalignedseqs = $3
else if ($option == alignedseqs) then
	set alignedseqs = $3
endif

if ($option == query_seq) then
	# $9 holds the database
#echo "queryseq=$queryseq"
#echo "db=$9"
#echo "info=$info"
	./seqs_chosen_via_median_info.csh $queryseq $9 $info 
	if ($status != 0) then
		exit (-1)
	endif
#	mv $tmp/$pid.selectedfasta $tmp/$pid.alignedfasta
	set alignedseqs = $tmp/$pid.alignedfasta
	set option = "alignedseqs"
endif
#echo "option=$option"

if ($option == related_seqs) then
# related sequences already known. getting the alignment from PSI-BLAST, quick
# and dirty
	set alignedseqs = $tmp/$pid.alignedfasta
	$srcdir/separate_query_from_rest_of_seqs $unalignedseqs $tmp/$pid.queryseq $tmp/$pid.database
	$NCBI/formatdb -i $tmp/$pid.database -o T -p T 
# extremely large evalues and multipass threshold because want to make sure get all the
# sequences the user submits
	$NCBI/blastpgp -d $tmp/$pid.database -i $tmp/$pid.queryseq -o $tmp/$pid.psiblastout -m 0 -j 4 -e 10 -h 1 -b 399
	echo QUERY > $tmp/$pid.listseq
	grep ">" $tmp/$pid.database | cut -d" " -f1 | cut -c2- >> $tmp/$pid.listseq 
	$srcdir/seqs_from_psiblast_res $tmp/$pid.psiblastout $tmp/$pid.listseq 4 $tmp/$pid.queryseq $alignedseqs $pid
#	($srcdir/psiblast_res_to_fasta_dbpairwise $tmp/$pid.psiblastout $alignedseqs $iterations $tmp/$pid.queryseq) >>& $tmp/$pid.err
	($srcdir/seqs_to_msf $alignedseqs $tmp/$pid.msf) >>& $tmp/$pid.err
	set option = "alignedseqs" # goes on to next step
endif

if ($option == alignedseqs) then
#echo in alignedseqs
echo comments printed in $commentscsh
	# get the query sequence in $pid.queryseq (needed for logo)
	# change directory to tmp because if pass in a msf or clustal
	# alignment, will make a file called "mablock.*" in the home directory
	cd $tmp

	($srcdir/process_alignment $alignedseqs $alignedseqs.gapsremoved)
#echo "after process_alignment"
	($srcdir/seqs_to_matrixweb $alignedseqs.gapsremoved $polymorphism_file $output 0 0.05 $gap_option $exp_option $seq_identity_filter > $commentscsh) >>& $tmp/$pid.err
	$srcdir/allowed_subst_html $tmp/$pid.aatable $alignedseqs.gapsremoved 0.05 $gap_option $exp_option $seq_identity_filter
	unalias rm
	rm -f $tmp/mablock.$pid.*
	rm -f $tmp/$pid.domain*
	rm -f $tmp/$pid.seq.query.unfiltered
	rm -f $tmp/$pid.clumped*
	rm -f $tmp/$pid.block.*
	rm -f $tmp/$pid.blks
	rm -f $tmp/$pid.startingblock
	# MAIL RESULTS 
	# removed mail option 12/13/00. instead imposed cpu time limit
	if ($#argv > 9) then # there's an address in $10
		set return_address = "sift@fhcrc.org"
		/usr/bin/mailx -s "SIFT $pid Alignment" -r $return_address $10 < $alignedseqs
		/usr/bin/mailx -s "SIFT $pid Conditional Probabilities" -r $return_address $10 < $tmp/$pid.siftresults.matrix
                if (-e "$tmp/$pid.siftresults.predictions") then
		   /usr/bin/mailx -s "SIFT $pid Predictions" -r $return_address $10 < $tmp/$pid.siftresults.predictions
                endif
		set table_files = ($tmp/$pid.aatable*)
		@ index = 0
		foreach table_file ($table_files) 
			@ start = $index * 100 + 1
			@ end = ($index + 1) * 100
		# do not change a thing when writing subject!  
		# this was tricky and took me whole day!!
		set subject = "SIFT $pid table positions $start to $end"
			/usr/bin/mailx -s "$subject" -r $return_address $10 < $table_file
			@ index++
		end #foreach
		set error_files = ($tmp/$pid.*.error)
		foreach error_file ($error_files) {
			/usr/bin/mailx -s "SIFT $pid Errors" -r $return_address $10 < $error_file
		end #foreach errorfile
	endif # end of sending mail

	if ($status != 0) then
		echo Error in execution.
	endif
endif

exit (0)
