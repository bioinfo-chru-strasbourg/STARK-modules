#!/bin/csh
#		Clear some files out of ~/tmp/ directory
#  Remove today's *.err files earlier than this hour
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

set hdate = `date '+%h %d %H'`
echo $hdate
set hour = `echo $hdate | awk '{print $3}'`
#echo $hour

#---------------------------------------------------------------------------
@ c = 299
while ($c >= 100)
#	set files = (`ls -ldt ~/tmp/$c* | tr -s ' ' ' ' | grep "$zdate" | awk '{print $8}'`)
#	set times = (`ls -ldt ~/tmp/$c* | tr -s ' ' ' ' | grep "$zdate" | awk '{print $7}'`)
 	set files = (`ls -ldt ~/tmp/$c*.err ~/tmp/$c*.seq* | tr -s ' ' ' ' | grep "$zdate" | awk '{print $8}'`)
 	set times = (`ls -ldt ~/tmp/$c*.err ~/tmp/$c*.seq* | tr -s ' ' ' ' | grep "$zdate" | awk '{print $7}'`)
	@ i = 0
	while ($i <= $#files)
#	echo "$i $files[$i] $times[$i]"
		set h = `echo $times[$i] | awk -F: '{print $1}'`
		if ($h < $hour) then
			rm -fr $files[$i] >& /dev/null
		endif
		@ i += 1
	end
   	@ c -= 1
end
#---------------------------------------------------------------------------
@ c = 99
while ($c >= 1)
	set files = (`ls -ldt ~/tmp/$c*.err ~/tmp/$c*.seq* | tr -s ' ' ' ' | grep "$zdate" | awk '{print $8}'`)
	set times = (`ls -ldt ~/tmp/$c*.err ~/tmp/$c*.seq* | tr -s ' ' ' ' | grep "$zdate" | awk '{print $7}'`)
	@ i = 0
	while ($i <= $#files)
#	echo "$i $files[$i] $times[$i]"
		set h = `echo $times[$i] | awk -F: '{print $1}'`
		if ($h < $hour) then
			rm -fr $files[$i] >& /dev/null
		endif
		@ i += 1
	end
   	@ c -= 1
end

#----------------------------------------------------------------------------
exit(0)

