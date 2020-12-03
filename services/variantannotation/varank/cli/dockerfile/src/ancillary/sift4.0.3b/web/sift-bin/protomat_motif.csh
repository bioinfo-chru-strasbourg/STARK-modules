#!/bin/csh
#               protomat.csh <file of proteins in fasta format>
#       Final blocks in $name.mblks and .gblks
#	from Jorja 06-06-00
#Dec.13, 2000 make sure I don't call status anywhere because I have
# a timer in seq_to_subfamily_block.csh and if I use something
# like echo, then echo call will be successful although protomat
# call may not be

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
unalias rm
unalias mv
limit coredumpsize 1k
limit datasize 256m

set protomat_bindir = ~ptest/bin
set path = ($protomat_bindir $path)

set pros = $1
set temp = $pros:t; set name = $temp:r
touch $name.outprotomat
#       For dups:  motifj 4 -$pros 0 dups 17 (0 seqs => n/2 to start)
#echo $pros
$protomat_bindir/motifj 4 -$pros >& /dev/null
#  Run motomat again so sequences aren't clumped/re-ordered
$protomat_bindir/motomat $name.mot 1 1 -10 >& /dev/null
#echo motifs in $name.mot
if ( -e $name.blks ) then
#  Add weights to blocks
   $protomat_bindir/blweight $name.blks $name.mblks P M >& /dev/null
#       Produce a "multiple alignment"
   $protomat_bindir/blalign $name.mblks >> $name.outprotomat
#       Make cobbler sequence = $name.mcob
   echo "TY     2"           > $name.cf
   echo "BL     $name.mblks"   >> $name.cf
   echo "DB     $pros"    >> $name.cf
   echo "OU     $name.mcob"    >> $name.cf
   echo "SU     $protomat_bindir/default.iij" >> $name.cf
#   $protomat_bindir/cobbler $name.cf >& /dev/null  # commented out 9/5/2000
#               Show the cobbler sequence now
#   echo "" >> $name.outprotomat
#   echo "" >> $name.outprotomat
#   echo "      **COBBLER sequence from MOTIF**" >> $name.outprotomat
#   cat $name.mcob | tr '\015' ' ' >> $name.outprotomat
else
   echo "ERROR: No blocks produced by MOTIF" >> $name.outprotomat
endif
exit 0  # sept. 5 2000, jsut motif

#
#       Clean up
# commented out aug.17 2000
#rm $name.mot $name.blks >& /dev/null
#
#===========================================================
#       Make blocks using Gibbs now  
#
echo in gibbs
echo " " >> $name.outprotomat
echo " " >> $name.outprotomat
echo "              **BLOCKS from GIBBS**" >> $name.outprotomat
echo " " >> $name.outprotomat

set minlen = (`grep " MINIMUM" $name.motifj.pros |  awk '{print $2}'`)
if ($minlen < 8) then
   echo "ERROR: Minimum sequence length is 8 for GIBBS" >> $name.outprotomat
   echo "ERROR: No blocks produced by GIBBS" >> $name.outprotomat
else
   #    Model heuristic; needs MINIMUM sequence length from $name.motifj.pros
   set model = (`$protomat_bindir/model.run $name.motifj.pros`)
   $protomat_bindir/gibbs $name.motifj.pros $model -f -s3 >& $name.gblks
echo status of gibbs $status

#exit # PN -- I don't have a lis file

   #    Convert gibbs output for motomat
#   grep ">" $name.lis > $name.temp
   cat $name.motifj.pros.sn >> $name.temp
   $protomat_bindir/blk2mot $name.motifj.pros $name.temp $name.mot
   $protomat_bindir/motomat $name.mot 1 1 -15 >& /dev/null
   if ( -e $name.blks ) then
   #  Add weights to blocks
      $protomat_bindir/blweight $name.blks $name.gblks P M >& /dev/null
   #    Make multiple alignment
      $protomat_bindir/blalign $name.gblks >> $name.out
   #    Make cobbler sequence = $name.gcob
      echo "TY  2"           > $name.cf
      echo "BL  $name.gblks"   >> $name.cf
      echo "DB  $pros"    >> $name.cf
      echo "OU  $name.gcob"    >> $name.cf
      echo "SU  $protomat_bindir/default.iij" >> $name.cf
      $bindir/cobbler $name.cf >& /dev/null
   #            Show the cobbler sequence now
      echo "" >> $name.out
      echo "" >> $name.out
      echo "      **COBBLER sequence from GIBBS**" >> $name.out
      cat $name.gcob | tr '\015' ' '>> $name.out
   else
      #  $name.gblks contains the Gibbs output now 
      cat $name.gblks >> $name.out
      echo "ERROR: No blocks produced by GIBBS" >> $name.out
      rm $name.gblks
   endif
endif

#       Clean up
# commented out aug. 17, 2000
#rm  motomat.err $name.mot $name.blks >& /dev/null
#rm $name.cf $name.temp $name.motifj.pros $name.motifj.pros.sn
#rm $name.out
#
exit


