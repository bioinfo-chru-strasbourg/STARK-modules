#!/bin/csh
#		Clear out ~blocks/tmp/ & ~sift/tmp/ directories
#		Remove all files that aren't today's (run at 20:00)
#

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
unalias rm
unalias ls
#	Get today's date
#	Problem here:  date dd = "06", ls dd = " 6", but zdate is just "6"
#	So squeeze extra spaces out of ls output with "tr -s" below
set zdate = `date '+%h %d' | sed "s/ 0/  /g"`

#	May be too many files, so break it up
#---------------------------------------------------------------------------
set tmpdir = /home/sift/tmp
@ c = 9999
while ($c >= 0)
#	Get tmp file dates - awk is $8 for howard, $9 for morgan
   set files = (`ls -ldt $tmpdir/$c* | tr -s ' ' ' ' | grep -v "$zdate" | awk '{print $9}'`)
   foreach file ($files)
 	rm -f $file >& /dev/null
   end

   set files = (`ls -ldt $tmpdir/NP_${c}* | tr -s ' ' ' ' | grep -v "$zdate" | awk '{print $9}'`)
   foreach file ($files)
 	rm -f $file >& /dev/null
   end

   @ c -= 1
end

set files = (`ls -ldt $tmpdir/[A-Z]* | tr -s ' ' ' ' | grep -v "$zdate" | awk '{print $9}'`)
foreach file ($files)
	rm -f $file >& /dev/null
end
rm -f $tmpdir/formatdb.log

exit(0)

