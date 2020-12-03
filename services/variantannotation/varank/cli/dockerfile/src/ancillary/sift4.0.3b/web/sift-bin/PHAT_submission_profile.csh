#!/bin/csh
#
#PHAT_submission.csh
# step (1) pro_to_prf turn sequence ($1, $2) into profile with matrix values
#		      ($3, $4) as specified by user
#		      Output $.prf automatically created by program

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

setenv BLIMPS_DIR /usr/local/blimps

set bin = .
set tmp = ../tmp
set data = ../data

# >& /dev/null sends both standard output and standard error to /dev/null
($bin/pro_to_prf $3 $4 $1 $2) >& /dev/null 

#cat $1.prf
#rm $1

exit(0)

