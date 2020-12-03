#!/bin/csh
#		catfile <filename> <title> <PRE>
# centers title, has white background.

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

echo Content-type: text/html
echo ""
echo '<body bgcolor=white>'
echo '<H1><center>'
echo $2
echo '</center></H1>'
if ($3 == "PRE") then
	echo '<PRE>'
endif
if (-e $1) then
	cat $1
else 
	echo "Cannot open file $1"
endif

if ($3 == "PRE") then
	echo '</PRE>'
endif
</body>

#done
exit 0
