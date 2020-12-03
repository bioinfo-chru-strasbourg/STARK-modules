#!/bin/csh
#
#	Remove files at least 12 hours old after noon
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
#set zdate = `date '+%h %d %H' | sed "s/ 0/  /g"`
set hour = `date '+%H'`
@ cut = $hour - 12
echo "hour=$hour cut=$cut"
if ($cut < 0) then
   exit
endif

set tmpdir = /home/sift/tmp
@ c = 9999
while ($c >= 1)
echo "c=$c"
#	Get today's tmp file dates - awk is $8 for howard, $9 for morgan
   set files = (`ls -ldt $tmpdir/$c* | tr -s ' ' ' ' | grep "$zdate" | awk '{print $9}'`)
   foreach file ($files)
	set four = `ls -lgt $file | awk '{print $7}' | cut -d":" -f1`
#echo "four=$four"
	if ($four < $cut) then
 		rm -f $file >& /dev/null
	endif
   end
   @ c -= 1
end
exit

set files = (`ls -ldt $tmpdir/[A-Z]* | tr -s ' ' ' ' | grep "$zdate" | awk '{print $9}'`)
foreach file ($files)
	set four = `ls -lgt $file | awk '{print $7}' | cut -d":" -f1`
	if ($four < $cut) then
		rm -f $file >& /dev/null
	endif
end

exit(0)

