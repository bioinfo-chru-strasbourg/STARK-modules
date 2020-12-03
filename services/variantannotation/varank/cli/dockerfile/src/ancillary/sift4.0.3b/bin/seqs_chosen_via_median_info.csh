#!/bin/csh
# seqs_chosen_via_median_info.csh
# Arg1 = protein sequence in fasta format
# Arg2 = blastpgp database

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software

set tmpdir = $SIFT_DIR/tmp
set bindir = $SIFT_DIR/bin

unalias cp
unalias rm

set seq_database = $2
set median_threshold = 2.75;

# 	Cleanup query sequence.  NCBI stuff is intolerant of format
#	deviations.

set tail = $1:t
set pid = $tail:r
echo tail is $tail
set query = $tmpdir/$tail.query
echo query is $query

$bindir/fastaseqs $1 $query	>& /dev/null



# (1) PSIBLAST

 set iterations = 1 
# the query is called $$.fasta, get pid from that

#######   FILTER ############################
#$BLASTFILTER/seg $query -x > $query.imp
cp $query $query.unfiltered
#mv -f $query.imp $query # this has the filtered query 

# b option, # of alignments to show (else it is 250)
# m option, type of alignment 0-pairwise, 6-flat master slave w/out identites

# make error file & remove only if PSI-BLAST is finished within time allotted
# not filtering query because Jorja says new statistics take care of that

echo "PSI-BLAST time limit exceeded." > $tmpdir/$pid.time.error
limit cputime 30m

$NCBI/blastpgp -d $seq_database -i $query -o $query.out -m 0 -j $iterations -e .0001 -h 0.002 -b 399

if ($status != 0) then
        echo "exiting because stauts not equal to 0"
        exit (-1)
endif
unlimit cputime

# finished within 10 minutes, remove error file
rm -f $tmpdir/$pid.time.error

##echo done psiblast

# get the psiblast alignment of all the sequences found.

$bindir/psiblast_res_to_fasta_dbpairwise $query.out $query.globalX $iterations $query.unfiltered
if ($status != 0) then
	exit (-1)
endif

$bindir/clump_output_alignedseq $query.globalX $tmpdir/$pid.clumped .9 0 >& /dev/null

echo "choose seqs time limit exceeded." > $tmpdir/$pid.time.error
limit cputime 30m

$bindir/choose_seqs_via_psiblastseedmedian $query.unfiltered $tmpdir/$pid.clumped $query.selectedclumped 1.0 $pid $median_threshold >& /dev/null

if ($status != 0) then
       echo "exiting because stauts not equal to 0"
       exit (-1)
endif
unlimit cputime
rm -f $tmpdir/$pid.time.error

$bindir/consensus_to_seq $query.selectedclumped $tmpdir/$pid.clumped.consensuskey $tmpdir/$pid.selected >& /dev/null


$bindir/seqs_from_psiblast_res $query.out $tmpdir/$pid.selected $iterations $query.unfiltered $tmpdir/$pid.alignedfasta $pid

	rm $tmpdir/$pid.clumped*
	rm $query.selectedclumped* 
	rm $tmpdir/$pid.selected
	rm $tmpdir/$pid.unaligned
	rm $query.globalX
	rm $tmpdir/$pid.psiblastout	
	rm $pid.TEMP*
	rm $query
	rm $query.unfiltered
exit (0)


