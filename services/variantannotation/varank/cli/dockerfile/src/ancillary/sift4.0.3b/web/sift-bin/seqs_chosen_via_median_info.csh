#!/bin/csh
# seqs_chosen_via_median_info.csh
#12/14/07 nice blastpgp  JGH
#	Group write required

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
umask 002

#setenv BLASTFILTER /usr/local/filter
setenv BLASTFILTER /usr/local/bin
setenv NCBI /usr/local/packages/blast/bin/ 
setenv BLASTMAT /usr/local/packages/blast/data/
set srcdir = "/opt/www/sift/src_from_howard/src/src_for_csh/established_and_tested_src_for_csh"
set tmpdir = "/opt/www/sift/tmp/"
set settingsdir = "../settings"

limit coredumpsize 1k
limit datasize 1024m

if ($2 == sp_tr) then
	#set seq_database = "/d3/databases/swiss/sp_tr.uni"
	set seq_database = "/opt/www/sift/databases/swiss/uniprot_trembl"
else if ($2 == sp) then
	#set seq_database = "/d3/databases/swiss/swiss.uni"
	set seq_database = "/opt/www/sift/databases/swiss/uniprot_sprot"
else if ($2 == nr) then
	#set seq_database = "/d3/databases/ncbi/nr"
	set seq_database = "/opt/www/sift/databases/ncbi/nr"
endif

unalias cp
unalias rm

set median_threshold = $3;

# 	Cleanup query sequence.  NCBI stuff is intolerant of format
#	deviations.

set tail = $1:t
set pid = $tail:r
#echo tail is $tail
set query = $tmpdir/$tail.query
#echo query is $query

./fastaseqs $1 $query	>& /dev/null
##echo database is $seq_database

##if ($#argv == 0) then
##	echo iterative_impala.csh.  fasta sequence
##	exit (-1);
##endif


# (1) PSIBLAST

 set iterations = 2
# the query is called $$.fasta, get pid from that

#######   FILTER ############################
echo "query is $query"
set outquery = "$query.imp"
#echo "outquery is $outquery"
$BLASTFILTER/seg $query -x > $query.imp
mv $query $query.unfiltered
mv -f $query.imp $query # this has the filtered query 

# b option, # of alignments to show (else it is 250)
# m option, type of alignment 0-pairwise, 6-flat master slave w/out identites

# check that query is a protein sequence
$srcdir/check_protein $query.unfiltered

echo $seq_database Query $query Out $query.out

# make error file & remove only if PSI-BLAST is finished within time allotted
# not filtering query because Jorja says new statistics take care of that

echo "PSI-BLAST time limit exceeded." > $tmpdir/$pid.time.error
limit cputime 30m

nice $NCBI/blastpgp -d $seq_database -i $query -o $query.out -I T -m 0 -j $iterations -e .0001 -h 0.002 -b 399

if ($status != 0) then
        echo "exiting because status not equal to 0"
        exit (-1)
endif
unlimit cputime

# finished within 10 minutes, remove error file
rm -f $tmpdir/$pid.time.error

echo done psiblast

# get the psiblast alignment of all the sequences found.

$srcdir/psiblast_res_to_fasta_dbpairwise $query.out $query.globalX $iterations $query.unfiltered
if ($status != 0) then
	exit (-1)
endif

$srcdir/clump_output_alignedseq $query.globalX $tmpdir/$pid.clumped .9 0 >& /dev/null

echo "choose seqs time limit exceeded." > $tmpdir/$pid.time.error
limit cputime 30m

cd $tmpdir   # change because PSI-BLAST used in next command puts output 
	    # in home directory 

($srcdir/choose_seqs_via_psiblastseedmedian $query.unfiltered $tmpdir/$pid.clumped $query.selectedclumped 1.0 $pid $median_threshold) >>& $tmpdir/$pid.err

if ($status != 0) then
       echo "exiting because stauts not equal to 0"
       exit (-1)
endif
unlimit cputime
rm -f $tmpdir/$pid.time.error

$srcdir/consensus_to_seq $query.selectedclumped $tmpdir/$pid.clumped.consensuskey $tmpdir/$pid.selected >& /dev/null

 
#perl5 $srcdir/get_sequences.pl $tmpdir/$pid.selected $seq_database  

#mv $tmpdir/$pid.selected.seqs $tmpdir/$pid.unaligned

set alignedseqs = $tmpdir/$pid.alignedfasta
#        $srcdir/separate_query_from_rest_of_seqs $tmpdir/$pid.unaligned $tmpdir/$pid.queryseq $tmpdir/$pid.database
#        $NCBI/formatdb -i $tmpdir/$pid.unaligned -p T
# extremely large evalues and multipass threshold because want to make sure 
# get all the sequences the user submits 
#        $NCBI/blastpgp -d $tmpdir/$pid.unaligned -i $query.unfiltered -o $tmpdir/$pid.psiblastout -m 0 -j $iterations -e 0.0001 -h 1 -b 399

#	rm $tmpdir/$pid.unaligned.* # remove formatted database formats
       
 $srcdir/seqs_from_psiblast_res $query.out $tmpdir/$pid.selected $iterations $query.unfiltered $alignedseqs $pid >& /dev/null 
  $srcdir/seqs_to_msf_web $alignedseqs $tmpdir/$pid.msf

	rm $tmpdir/$pid.selected
	rm $tmpdir/$pid.unaligned
	rm $query.globalX
	rm $tmpdir/$pid.psiblastout	
	rm $tmpdir/$pid.TEMP*
	rm $query.out
	rm $query
	rm $query.unfiltered
exit (0);


