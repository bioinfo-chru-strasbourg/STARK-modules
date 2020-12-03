#!/bin/csh
#  logos.csh <pid> <aligned sequences> <polymorphism file>
# polymorphism file = '-' if there is no polymorphism file
#  will take the first sequence of aligned sequences as the referece sequence
#  file
#  The files colors, marks, wave and default.amino.frq must be in
#  the current directory; also makelogob (executable version of C program)

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

set pid = $1
set polymorphismfile = $3
unalias rm
unalias cp
set srcdir = "../src/src_for_csh/established_and_tested_src_for_csh"
set tmp = "../tmp"
echo entered logos.csh
set path = ($path /opt/sfw/bin)

touch marks
touch wave
$srcdir/separate_query_from_rest_of_seqs $2 $tmp/$pid.queryseqlogo $tmp/$pid.tmp
if ($status != 0) then
	exit -1
endif
cp ~btest/bin/colors $tmp/colors.$pid
$srcdir/seqs_to_blockMAXWIDTH $2 $tmp/$pid.blk
if ($status != 0) then
	exit -1
endif

#cp ~btest/bin/makelogob .  # use own makelogob
#cp ~btest/bin/colors $tmp/colors.$pid
cp ~btest/bin/default.amino.frq $tmp/
cp marks $tmp
cp wave $tmp
#   This creates one logo PS file for each block in the file $1
#       matrix_logob executes makelogob

#~btest/bin/matrix_logob $tmp/$pid.blk - $pid > /dev/null
# must be in working directory for symvec and makelogop files to be written
cd $tmp
$srcdir/matrix_logobPN $tmp/$pid.blk - $pid > /dev/null

rm colors.*
rm makelogop.*
rm symvec.*

set logodir = $tmp/logo.$pid
echo $logodir
$srcdir/separate_blocks $tmp/$pid.blk 
			 # block files are now separated with suffix 1, 2,...
set blockdir = $tmp/$pid.blk
echo $blockdir

foreach postscript_infile ($logodir.*)
        set extension = $postscript_infile:e
	set outfile = $pid.$extension.ps
        set blockfile = $blockdir$extension
        $srcdir/process_logo_ps $postscript_infile $outfile $blockfile $polymorphismfile $tmp/$pid.queryseqlogo
	echo $status of process_logo
	rm $postscript_infile
	echo $status of for loop
        echo $postscript_infile $outfile $blockfile
end
rm $tmp/$pid.blk*
rm $tmp/$pid.queryseqlogo

exit(0)

