#!/bin/csh
# seq_to_subfamily_block.csh
#12/14/07 nice blastpgp JGH

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.


#  (1)  Search database using psiblast for similar sequences to a query 
#       sequence ($1).
#         (a)  a fasta database of these sequences is made.
#
#  (2)  Gibbs & MOTIF is executed on the query sequence ($1) to find a conserved
# 	domain.
#  (3)  Sequences found in (1) with this domain are obtained.
#  (4)	A seed of sequences that are 90% identical is made.
#  (5)  Iterative psiblast to join the next closest sequence to the current
#	PSSM.
#  (6)  Reexecute Gibbs and MOTIF to find the conserved regions from these
#	closely related sequences

#12/13/00 added 10 min. time limits for PSI-BLAST and MOTIF
# (6) is not executed to save time & not necessary for SIFT


setenv BLASTFILTER /usr/local/filter
setenv NCBI /usr/local/blast

set srcdir = "../src/src_for_csh/established_and_tested_src_for_csh"
set tmpdir = "../tmp"
set settingsdir = "../settings"

if ($2 == sp_tr) then
	set seq_database = "/d3/databases/swiss/sp_tr.uni"
else if ($2 == sp) then
	set seq_database = "/d3/databases/swiss/swiss.uni"
else if ($2 == nr) then
	set seq_database = "/d3/databases/ncbi/nr"
endif

unalias cp
unalias rm

# 	Cleanup query sequence.  NCBI stuff is intolerant of format
#	deviations.
set tail = $1:t
set pid = $tail:r
echo tail is $tail
set query = $tmpdir/$tail.query
echo query is $query
#echo $query
./fastaseqs $1 $query	>& /dev/null
echo database is $seq_database

if ($#argv == 0) then
	echo iterative_impala.csh.  fasta sequence
	exit 0
endif


# (1) PSIBLAST

set iterations = 2
# the query is called $$.fasta, get pid from that
echo i am here

#######   FILTER ############################
$BLASTFILTER/seg $query -x > $query.imp
echo la de da
mv $query $query.unfiltered
echo whrer
mv -f $query.imp $query # this has the filtered query 
echo is the
# b option, # of alignments to show (else it is 250)
# m option, type of alignment 0-pairwise, 6-flat master slave w/out identites
echo $seq_database Query $query Out $query.out

# 12/12/00 commented out checkpoint file & pssm file from blast output

# make error file & remove only if PSI-BLAST is finished within time allotted
echo "PSI-BLAST time limit exceeded." > $tmpdir/$pid.time.error
limit cputime 30m

nice $NCBI/blastpgp -d $seq_database -i $query -o $query.out -m 0 -j $iterations -e .0001 -h 0.002 -b 399

if ($status != 0) then
        echo "exiting because stauts not equal to 0"
        exit (-1)
endif
unlimit
# finished within 10 minutes, remove error file
rm -f $tmpdir/$pid.time.error

echo done psiblast

# get the psiblast alignment of all the sequences found.

$srcdir/psiblast_res_to_fasta_dbpairwise $query.out $query.globalX $iterations $query.unfiltered
if ($status != 0) then
	exit (-1)
endif
echo done psiblast res_to_fasta
echo query is $query

#option 1 for clumping to minimize work gibbs has to do
@ clustered_seq = (`$srcdir/clump $query.globalX $tmpdir/$pid.clumped .9 1 | awk '{print $5}'`)
echo "clumpint $clustered_seq"
# Clumping too severe, too few sequences left so use option 0
if ($clustered_seq < 5) then 
	echo "WARNING!  After clumping, only $clustered_seq was left."
	echo "Reclumping..."
	@ clustered_seq = (`$srcdir/clump $query.globalX $tmpdir/$pid.clumped .9 0 | awk '{print $5}'`)
	echo "GIBBS and motif on $clustered_seq sequences"
endif
if ($status != 0) then
	exit (-1)
endif

echo done clumping in $tmpdir/$pid.clumped
# (2) MOTIF 
#  MUST CHANGE tmp tmpdir because protomat puts output in directory
cd $tmpdir
echo "Time to execute MOTIF exceeded maximum time limit." > $tmpdir/$pid.time.error
limit cputime 30m
./protomat_motif.csh $tmpdir/$pid.clumped
if (-e $tmpdir/$pid.mblks) then
        unlimit
	rm -f $tmpdir/$pid.time.error # finished motif, remove error file with
					   # time limit 
else
        echo "exceeded time limit"
        exit (-1)
endif

#mv $pid.* $tmpdir
#AUG.17, 2000 used both motifs & gibbs regions to select sequences
# Sept. 5, 2000 use just motif to select regions.  saves time & doesn't
# make that much of a difference 
$srcdir/gibbs_to_regions $tmpdir/$pid.mblks $tmpdir/$pid.regions
if ($status != 0) then
	exit (-1)
endif

echo status of gibbs_motifs_to_regions $status
# (3)  Now that the conserved region has been identified by MOTIF
#      Get the conserved regions from all of the sequences
$srcdir/psiblast_res_to_domain $query.out $tmpdir/$pid.domain $tmpdir/$pid.regions $iterations $query.unfiltered
if ($status != 0) then
        exit (-1)
endif
echo $status for psiblast_res_to_domain

# clump conserved regions to reduce redundancy.  also affects error correction
# factor
#option 0 for clumping so that iterative search against all consensus sequences
# gets those that are in their own cluster as well 
$srcdir/clump $tmpdir/$pid.domain $tmpdir/$pid.domainclumped .9 0
if ($status != 0) then
	exit (-1)
endif

$NCBI/formatdb -i $tmpdir/$pid.domainclumped -o T -p T

#echo done psiblast_res_tofasta

# from sequences that are > 90% identical with query sequence in the conserved
# regions, this is the seed. put in block format 
$srcdir/seqs_to_block $tmpdir/$pid.domainclumped $tmpdir/$pid.startingblock $tmpdir/$pid.queryfasta
if ($status != 0) then
	exit (-1)
endif

# CHOOSE SEED BY PERCENT IDENTITY
# iterative procedure that keeps adding sequences from $3 until information
# decreases.  $2 is the output 
$srcdir/choose_seqs_via_psiblastseedidentity $tmpdir/$pid.startingblock $tmpdir/$pid.block.clumped9 $tmpdir/$pid.domainclumped $tmpdir/$pid.queryfasta .9 $pid

echo done choose_seqs_via_psiblastseedidentity $status
# now that consensus sequences have been chosen, get the  original sequences
# that made up the consensus sequences
$srcdir/consensus_to_seq $tmpdir/$pid.block.clumped9.seq $tmpdir/$pid.domainclumped.consensuskey $tmpdir/$pid.block.9.seq

# get the alignment of these selected sequences from PSI-BLAST.
$srcdir/seqs_from_psiblast_res $query.out $tmpdir/$pid.block.9.seq $iterations $query.unfiltered $tmpdir/$pid.selectedfasta $pid  
if ($status != 0) then
        exit (-1)
endif

# turn into blocks
$srcdir/seqs_to_blockdomains $tmpdir/$pid.selectedfasta $tmpdir/$pid.subblock $tmpdir/$pid.regions 
if ($status != 0) then
	exit (-1)
endif

# get alignment
$srcdir/seqs_to_msf_web $tmpdir/$pid.selectedfasta $tmpdir/$pid.msf

# get sequences from psiblast file
$srcdir/seqs_unaligned_from_psiblast_res $query.out $tmpdir/$pid.block.clumped9.seq $iterations $query.unfiltered $query.unalignedsubfamily $pid 
 
#cleanup files from previous protomat.csh
rm -f $pid.outprotomat $pid.cf $pid.gblks $pid.mblks $pid.mcob $pid.mot $pid.motifj.pros $pid.motifj.pros.sn


#rm -f $pid.domain*
#rm -f $pid.block.*
#rm -f $pid.clumped.*
rm -f $pid.motifj.pros.sn
#rm -f $pid.regions
#rm -f $pid.seq.query.*
#rm -f $pid.startingblock
rm -f $pid.TEMP_PSIBLAST

exit (0)


