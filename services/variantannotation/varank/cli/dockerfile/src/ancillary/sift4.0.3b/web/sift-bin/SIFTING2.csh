#!/bin/csh
# SIFTING2.csh
# args are  pid opt seq out info poly ? idfilt db address
#           1   2   3   4   5    6    7 8      9  10
# 3/8/04 If $address, return predictions output & put pid in subject
#12/11/07 Added nice to blastpgp execution. JGH
#3/6/08 Remove .psiblastout and .database* at end JGH
# 10/2/08 filter out first line ONTENT-TYPE generated from website
# which was throwing off following programs when submitting alignment

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

#       Need group write on all files created
umask 002
#limit filesize 10m

# Only integer arithmetic is possible in csh.  the intolerance probability threshold
# called in seqs_to_matrixweb is set at 0.05, must change this manually in this script

setenv BLIMPS_DIR /opt/www/sift/blimps-3.8
setenv NCBI /usr/local/packages/blast/bin
#Added by pkumar Malloc_check set to 0 to avoid glibc double free or curruption error
setenv MALLOC_CHECK_ 0

set tmp = "/opt/www/sift/tmp/"
set bin = "/opt/www/sift/htdocs/sift-bin/"
set commentscsh = $tmp/$1.commentscsh

# no results in commentscsh -- just for debugging
set srcdir = "/opt/www/sift/src_from_howard/src/src_for_csh/established_and_tested_src_for_csh"
set scriptsdir = "/opt/www/sift/htdocs/scripts"
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
        $bin/seqs_chosen_via_median_info.csh $queryseq $9 $info
        if ($status != 0) then
                exit (-1)
        endif
#       mv $tmp/$pid.selectedfasta $tmp/$pid.alignedfasta
        set alignedseqs = $tmp/$pid.alignedfasta
        set option = "alignedseqs"
endif
#echo "option=$option"
set option_links_in_msf =  0

if ($option == NCBI_blink) then
#echo "in here $pid seq numb"
echo " seq $3 4 $4 5 $5 6 $6  7 $7 8 $8 optin 9  $9 10 $10<BR>"
        echo "perl $srcdir/get_BLINK_seq.pl $3 $tmp/$pid.unaligned $9 " >> $tmp/log;
        $srcdir/get_BLINK_seq.pl $3 $tmp/$pid.unaligned $9
#       (perl $srcdir/get_BLINK_seq.pl $3 $tmp/$pid.unaligned $9) >>& /dev/null
        # copy sequences to unaligned seqs
        set unalignedseqs = $tmp/$pid.unaligned
        set option = "related_seqs"
        set option_links_in_msf = 1
endif
if ($option == related_seqs) then
# related sequences already known. getting the alignment from PSI-BLAST, quick
# and dirty
        set alignedseqs = $tmp/$pid.alignedfasta
        $srcdir/separate_query_from_rest_of_seqs $unalignedseqs $tmp/$pid.queryseq $tmp/$pid.database
	#following is a new script replacing the above to avoid segfault
	#$scriptsdir/separate_query_from_database.pl $unalignedseqs $tmp/$pid.queryseq $tmp/$pid.database  
        $NCBI/formatdb -i $tmp/$pid.database -o T -p T
# extremely large evalues and multipass threshold because want to make sure get all the
# sequences the user submits
        nice $NCBI/blastpgp -d $tmp/$pid.database -i $tmp/$pid.queryseq -o $tmp/$pid.psiblastout -m 0 -j 4 -e 10 -h 1 -b 399
        echo QUERY > $tmp/$pid.listseq
        grep ">" $tmp/$pid.database | cut -d" " -f1 | cut -c2- >> $tmp/$pid.listseq
#       perl $srcdir/clean_up_converged_alignment.pl $tmp/$pid.psiblastout > $tmp/$pid.tmmp
#       mv $tmp/$pid.tmmp $tmp/$pid.psiblastout
#       grep -v "^CONVERGED" $tmp/$pid.psiblastout > $tmp/$pid.tmmp
#       mv $tmp/$pid.tmmp $tmp/$pid.psiblastout
        $srcdir/seqs_from_psiblast_res $tmp/$pid.psiblastout $tmp/$pid.listseq 4 $tmp/$pid.queryseq $alignedseqs $pid 
#       ($srcdir/psiblast_res_to_fasta_dbpairwise $tmp/$pid.psiblastout $alignedseqs $iterations $tmp/$pid.queryseq) >>& /dev/null #$tmp/$pid.err
        if ($option_links_in_msf == 0) then
                ($srcdir/seqs_to_msf $alignedseqs $tmp/$pid.msf) >>& /dev/null #$tmp/$pid.err
        else
        # to print links for BLink alignment
                ($srcdir/seqs_to_msf_web $alignedseqs $tmp/$pid.msf)
        endif

        set option = "alignedseqs" # goes on to next step
endif

if ($option == alignedseqs) then
#echo in alignedseqs
#echo comments printed in $commentscsh
        # get the query sequence in $pid.queryseq (needed for logo)
        # change directory to tmp because if pass in a msf or clustal
        # alignment, will make a file called "mablock.*" in the home directory
        cd $tmp
# 10-02-08
        (cat $alignedseqs | perl -pe 's/ONTENT\-TYPE\: APPLICATION\/OCTET-STREAM//' > $alignedseqs.2)
        ($srcdir/process_alignment $alignedseqs.2 $alignedseqs.gapsremoved)
#echo "after process_alignment"
	echo seqs_to_matrixweb $alignedseqs.gapsremoved $polymorphism_file $output 0 0.05 $gap_option $exp_option $seq_identity_filter 
        ($srcdir/seqs_to_matrixweb $alignedseqs.gapsremoved $polymorphism_file $output 0 0.05 $gap_option $exp_option $seq_identity_filter > $commentscsh) >>& /$tmp/$pid.err
        $srcdir/allowed_subst_html $tmp/$pid.aatable $alignedseqs.gapsremoved 0.05 $gap_option $exp_option $seq_identity_filter >& $tmp/remove 
        unalias rm
        rm -f $tmp/mablock.$pid.*
        rm -f $tmp/$pid.domain*
        rm -f $tmp/$pid.seq.query.unfiltered
        rm -f $tmp/$pid.clumped*
        rm -f $tmp/$pid.block.*
        rm -f $tmp/$pid.blks
        rm -f $tmp/$pid.startingblock
        rm -f $tmp/$pid.comments
#        rm -f $tmp/$pid.psiblastout
        rm -f $tmp/$pid.database*
        # MAIL RESULTS
        # removed mail option 12/13/00. instead imposed cpu time limit
	set table_files = ($tmp/$pid.aatable*)
        foreach table_file ($table_files)
        	cp $table_file $table_file.html
	end
	

        if ($#argv > 9) then # there's an address in $10
                set return_address = "sift@fhcrc.org"
                cp $alignedseqs $alignedseqs.txt
                cp $tmp/$pid.siftresults.matrix $tmp/$pid.siftresults.matrix.html
		cp $tmp/$pid.siftresults.predictions $tmp/$pid.siftresults.predictions.html
	        $bin/email_sift_aatable_results.pl $pid $10 

#               /usr/bin/mailx -s "SIFT $pid Alignment" -r $return_address $10 < $alignedseqs
#               /usr/bin/mailx -s "SIFT $pid Conditional Probabilities" -r $return_address $10 < $tmp/$pid.siftresults.matrix
# changed to muttx.  Subject and To: must be flushed right
# because piping to mutt
# and header must be in proper format

# mailing alignment ##
#echo "From: SIFT <sift@fhcrc.org>\
#Subject: SIFT $pid Alignment \
#To: $10 \
#"  | /opt/sfw/bin/mutt  -a $alignedseqs.txt -H /dev/stdin
#
## mailing probabilities
#echo "From: SIFT <sift@fhcrc.org> \
#Subject: SIFT $pid Conditional Probabilities \
#To: $10 \
#" | /opt/sfw/bin/mutt -a $tmp/$pid.siftresults.matrix.html -H /dev/stdin
#
#                if (-e "$tmp/$pid.siftresults.predictions") then
##                  /usr/bin/mailx -s "SIFT $pid Predictions" -r $return_address $10 < $tmp/$pid.siftresults.predictions
#                  cp $tmp/$pid.siftresults.predictions $tmp/$pid.siftresults.predictions.html
#echo "From: SIFT <sift@fhcrc.org> \
#Subject: SIFT $pid Predictions \
#To: $10\
#" |  /opt/sfw/bin/mutt -a $tmp/$pid.siftresults.predictions.html -H /dev/stdin
#                endif
#                set table_files = ($tmp/$pid.aatable*)
#                @ index = 0
#                foreach table_file ($table_files)
#                        @ start = $index * 100 + 1
#                        @ end = ($index + 1) * 100
#                # do not change a thing when writing subject!
#                # this was tricky and took me whole day!!
#                set subject = "SIFT $pid table positions $start to $end"
#                cp $table_file $table_file.html
##                       /usr/bin/mailx -s "$subject" -r $return_address $10 < $table_file.html
#echo "From: SIFT <sift@fhcrc.org> \
#Subject: $subject \
#To: $10\
#" | /opt/sfw/bin/mutt  -a $table_file.html -H /dev/stdin
#                        @ index++
#                end #foreach
#                set error_files = ($tmp/$pid.*.error)
#                foreach error_file ($error_files) {
#                        /usr/bin/mailx -s "SIFT $pid Errors" -r $return_address $10 < $error_file
#
#                end #foreach errorfile
#        endif # end of sending mail

#        if ($status != 0) then
#                echo Error in execution.
#        endif
#endif

exit (0)

