#!/bin/csh
#
#	Compile the SIFT programs
#	>>>Edit the cccb script for your installation of BLIMPS<<<
#
unalias mv

./cccb choose_seqs_via_psiblastseedmedian
./cccb clump_output_alignedseq
./cccb consensus_to_seq
./cccb info_on_seqs
./cccb psiblast_res_to_fasta_dbpairwise
./cccb seqs_from_psiblast_res

mv choose_seqs_via_psiblastseedmedian ../bin
mv clump_output_alignedseq ../bin
mv consensus_to_seq ../bin
mv info_on_seqs ../bin
mv psiblast_res_to_fasta_dbpairwise ../bin
mv seqs_from_psiblast_res ../bin

exit
