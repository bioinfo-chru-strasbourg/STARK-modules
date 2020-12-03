#!/bin/csh
# gibbs_then_iteratie_psiblast.csh

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

setenv BLASTFILTER /usr/local/filter
setenv NCBI /usr/local/blast

set srcdir = "../src/src_for_csh"
set tmpdir = "../tmp"
set settingsdir = "../settings"

if ($2 == sp_tr) then
	set seq_database = "/d3/databases/swiss/sp_tr.uni"
else if ($2 == sp) then
	set seq_database = "/d3/databases/swiss/swiss38.uni"
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

set iterations = 4
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

$NCBI/blastpgp -d $seq_database -i $query -o $query.out -m 0 -j $iterations -C $query.chk -Q $query.pssm -e .0001 -b 399
echo done psiblast

$srcdir/psiblast_res_to_fasta_dbpairwise $query.out $query.globalX $iterations $query.unfiltered
echo done psiblast res_to_fasta
echo query is $query

#option 1 for clumping to minimize work gibbs has to do
@ clustered_seq = (`$srcdir/clump $query.globalX $tmpdir/$pid.clumped .9 1 | awk '{print $5}'`)
echo "clumpint $clustered_seq"
# Clumping too sever, option 0
if ($clustered_seq < 5) then 
	echo "WARNING!  After clumping, only $clustered_seq was left."
	echo "Reclumping..."
	@ clustered_seq = (`$srcdir/clump $query.globalX $tmpdir/$pid.clumped .9 0 | awk '{print $5}'`)
	echo "GIBBS and motif on $clustered_seq sequences"
endif

echo done clumping in $tmpdir/$pid.clumped
# (2) GIBBS 
echo stomethin gwrong
#  MUST CHANGE tmp tmpdir because protomat puts output in directory
cd $tmpdir
set d = (`./protomat.csh $tmpdir/$pid.clumped`)
echo $d
echo status of protomat.csh $status
echo gibbs

echo $tmpdir/$pid.clumped

#mv $pid.* $tmpdir
#AUG.17, 2000 used both motifs & gibbs regions to select sequences
$srcdir/gibbs_motifs_to_regions $tmpdir/$pid.gblks $tmpdir/$pid.motifj.pros.sn $tmpdir/$pid.regions

echo status of gibbs_motifs_to_regions $status

$srcdir/psiblast_res_to_domain $query.out $tmpdir/$pid.domain $tmpdir/$pid.regions $iterations $query.unfiltered
echo $status for psiblast_res_to_domain

#option 0 for clumping so that iterative search against all consensus sequences
# gets those that are in their own cluster as well  NOT CLUMPING AT THIS TIME
#$srcdir/clump $tmpdir/$pid.domaingibbs1 $tmpdir/$pid.domaingibbsclumped 1 0

$NCBI/formatdb -i $tmpdir/$pid.domain -o T -p T

#echo done psiblast_res_tofasta

#
$srcdir/seqs_to_block $tmpdir/$pid.domain $tmpdir/$pid.startingblock $tmpdir/$pid.queryfasta
echo done seqs_to_block $status

# CHOOSE SEED BY PERCENT IDENTITY
$srcdir/choose_seqs_via_psiblastseedidentity $tmpdir/$pid.startingblock $tmpdir/$pid.block.9 $tmpdir/$pid.domain $tmpdir/$pid.queryfasta .9 $pid

echo done choose_seqs_via_psiblastseedidentity $status

#$srcdir/robust $tmpdir/$pid.startingblockmotif $tmpdir/$pid.blockmotif.9.seq $tmpdir/$pid.domainmotif .9 $pid $tmpdir/$pid.robustmotif

$srcdir/seqs_from_psiblast_res $query.out $tmpdir/$pid.block.9.seq $iterations $query.unfiltered $tmpdir/$pid.selectedfasta $pid  

#cleanup files from previous protomat.csh
rm -f $pid.outprotomat $pid.cf $pid.gblks $pid.mblks $pid.mcob $pid.mot $pid.motifj.pros $pid.motifj.pros.sn

./protomat.csh $tmpdir/$pid.selectedfasta
echo done with second protomat $status

exit (0)


