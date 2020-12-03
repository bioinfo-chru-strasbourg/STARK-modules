#!/bin/csh
#		clean.csh
#	Removes temporary sift files older than $hours old
#	Change hours to whatever you want
#	Makes a csh file which is then executed
#

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
set hours = 12
#	arg to clean.pl is number of hours ago
set bin = /home/blocks/bin
set dir = /home/sift/tmp
set out = $dir/$$.csh
echo $out
ls -lgt $dir | awk '{print $5 " " $6 " " $7 " " $8}' | $bin/clean.pl $hours > $out
cd $dir >& /dev/null
chmod a+x $out
./$out >& /dev/null

exit
