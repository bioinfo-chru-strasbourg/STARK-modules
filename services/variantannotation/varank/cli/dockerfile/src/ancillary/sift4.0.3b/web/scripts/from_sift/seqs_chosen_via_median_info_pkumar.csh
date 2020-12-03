#!/bin/csh
# seqs_chosen_via_median_info.csh
#12/14/07 nice blastpgp  JGH
#	Group write required
umask 002

setenv BLASTFILTER /usr/local/filter
setenv NCBI /usr/local/bin
setenv BLIMPS_DIR /usr/local/projects/GOSII/pkumar/sift3.0/blimps-3.9
set query_loc = $1
set aln_loc = $2
set subst_loc = $3;
set pid = $4 # unique identifier
set outdir = $5 #output directory 
set srcdir = "/usr/local/projects/GOSII/pkumar/sift3.0/bin/"
set tmpdir = $outdir/$pid"_result"
mkdir -p $tmpdir

limit coredumpsize 1k
limit datasize 1024m



unalias cp
unalias rm

set median_threshold = 2.75; # if the median information content ends up being btween 2.75 and 3.25, that's good. 

#set tail = $1:t
 # the protein sequence in fasta format. the one the user is interested in.  

$srcdir/clump_output_alignedseq $aln_loc $tmpdir/$pid.clumped .9 0 >& $tmpdir/$pid.clump_output_alignedseq.out 

cd $tmpdir   # change because PSI-BLAST used in next command puts output 
	    # in home directory 

$srcdir/choose_seqs_via_psiblastseedmedian $query_loc $tmpdir/$pid.clumped $pid.selectedclumped 1.0 $pid $median_threshold >& $tmpdir/$pid.choose_seqs_via_psiblastseedmedian.out
$srcdir/consensus_to_seq $pid.selectedclumped $tmpdir/$pid.clumped.consensuskey $tmpdir/$pid.selected
set alignedseqs = $tmpdir/$pid.alignedfasta
set fileSize=`stat -c %s $tmpdir/$pid.selected`  

if ($fileSize == 0) then
        $srcdir/choose_seqs_via_psiblastseedmedian_original $query_loc $tmpdir/$pid.clumped $pid.selectedclumped 1.0 $pid $median_threshold >& $tmpdir/$pid.choose_seqs_via_psiblastseedmedian.out
        $srcdir/consensus_to_seq $pid.selectedclumped $tmpdir/$pid.clumped.consensuskey $tmpdir/$pid.selected
endif
     
 
/usr/local/projects/SIFT/scripts/extract_selected.pl $pid $aln_loc $outdir 
/usr/local/projects/SIFT/scripts/from_sift/SIFT_prediction.csh $query_loc $subst_loc $outdir $pid 

exit (0);


