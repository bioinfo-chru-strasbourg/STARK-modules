#!/bin/csh
#	SIFT.start.csh and 10 other arguments.. all for SIFTING.csh
#	SIFTING.csh & mails results to get around
#	SIFT_seq_submission.pl waiting for SIFTING.csh to finish
#

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
set bin = .
if ($#argv < 10) then
	echo "USAGE SIFT.start.csh <SIFTING arguments> <mail addr>"
	exit(-1)
endif
echo "SIFT results will be mailed to $10"
$bin/SIFTING.csh $1 $2 $3 $4 $5 $6 $7 $8 $9 $10 >& /dev/null & 
# running in background 
exit(0)
