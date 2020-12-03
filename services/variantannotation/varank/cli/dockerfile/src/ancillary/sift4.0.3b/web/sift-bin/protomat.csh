#!/bin/csh
#               protomat.csh <file of proteins in fasta format>
#       Final blocks in $name.mblks and .gblks
#	from Jorja 06-06-00

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
unalias rm
unalias mv
limit coredumpsize 1k
limit datasize 256m

set bindir = ~ptest/bin
set path = ($bindir $path)

set pros = $1
set temp = $pros:t; set name = $temp:r
touch $name.out
#       For dups:  motifj 4 -$pros 0 dups 17 (0 seqs => n/2 to start)
$bindir/motifj 4 -$pros >& /dev/null
#  Run motomat again so sequences aren't clumped/re-ordered
$bindir/motomat $name.mot 1 1 -10 >& /dev/null
if ( -e $name.blks ) then
#  Add weights to blocks
   $bindir/blweight $name.blks $name.mblks P M >& /dev/null
#       Produce a "multiple alignment"
   $bindir/blalign $name.mblks >> $name.out
#       Make cobbler sequence = $name.mcob
#   echo "TY     2"           > $name.cf
#   echo "BL     $name.mblks"   >> $name.cf
#   echo "DB     $pros"    >> $name.cf
#   echo "OU     $name.mcob"    >> $name.cf
#   echo "SU     $bindir/default.iij" >> $name.cf
#   $bindir/cobbler $name.cf >& /dev/null
#               Show the cobbler sequence now
#   echo "" >> $name.out
#   echo "" >> $name.out
#   echo "      **COBBLER sequence from MOTIF**" >> $name.out
#   cat $name.mcob | tr '\015' ' ' >> $name.out
else
   echo "ERROR: No blocks produced by MOTIF" >> $name.out
endif
#
#       Clean up
rm $name.mot $name.blks >& /dev/null

#
#===========================================================
#       Make blocks using Gibbs now  
#
echo " " >> $name.out
echo " " >> $name.out
echo "              **BLOCKS from GIBBS**" >> $name.out
echo " " >> $name.out
echo getting minlen from $name.motifj.pros

set minlen = (`grep " MINIMUM" $name.motifj.pros |  awk '{print $2}'`)
if ($minlen < 8) then
   echo "ERROR: Minimum sequence length is 8 for GIBBS" >> $name.out
   echo "ERROR: No blocks produced by GIBBS" >> $name.out
else
   #    Model heuristic; needs MINIMUM sequence length from $name.motifj.pros
   set model = (`$bindir/model.run $name.motifj.pros`)
   $bindir/gibbs $name.motifj.pros $model -f -s3 >& $name.gblks

   #    Convert gibbs output for motomat
#   grep ">" $name.lis > $name.temp
   cp $name.motifj.pros.sn $name.temp
   $bindir/blk2mot $name.motifj.pros $name.temp $name.mot
   $bindir/motomat $name.mot 1 1 -15 >& /dev/null
if (-e $name.cf) then
	unalias rm
	rm -f $name.cf
endif
   if ( -e $name.blks ) then
   #  Add weights to blocks
      $bindir/blweight $name.blks $name.gblks P M >& /dev/null
   #    Make multiple alignment
      $bindir/blalign $name.gblks >> $name.out
   #    Make cobbler sequence = $name.gcob
#      echo "TY  2"           > $name.cf
#      echo "BL  $name.gblks"   >> $name.cf
#      echo "DB  $pros"    >> $name.cf
#      echo "OU  $name.gcob"    >> $name.cf
#      echo "SU  $bindir/default.iij" >> $name.cf
#      $bindir/cobbler $name.cf >& /dev/null
   #            Show the cobbler sequence now
#      echo "" >> $name.out
#      echo "" >> $name.out
#      echo "      **COBBLER sequence from GIBBS**" >> $name.out
#      cat $name.gcob | tr '\015' ' '>> $name.out
   else
      #  $name.gblks contains the Gibbs output now 
      cat $name.gblks >> $name.out
      echo "ERROR: No blocks produced by GIBBS" >> $name.out
      rm $name.gblks
   endif
endif

#       Clean up
rm  motomat.err $name.mot $name.blks >& /dev/null
#rm $name.cf
#rm $name.motifj.pros.sn
rm  $name.temp $name.motifj.pros 
rm $name.out
#
exit 0


