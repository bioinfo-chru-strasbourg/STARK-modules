#!/bin/csh
#
#PHAT_submission.csh
# step (1) pro_to_prf turn sequence ($1, $2) into profile with matrix values
#		      ($4, $5) as specified by user
#		      Output $.prf automatically created by program
# step (2) swat	      swat on sequence ($1) and database ($3) given by user	

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

setenv BLIMPS_DIR /usr/local/blimps

#set tmp = ../tmp
set bin = .
#set data = ../data

# >& /dev/null sends both standard output and standard error to /dev/null
($bin/pro_to_prf $4 $5 $1 $2) >& /dev/null 
(/usr/local/bin/swat $1 $3 -M $1.prf > $1.out ) >& $1.err 

cat $1.out
rm $1*
rm $2
rm $3

exit(0)

