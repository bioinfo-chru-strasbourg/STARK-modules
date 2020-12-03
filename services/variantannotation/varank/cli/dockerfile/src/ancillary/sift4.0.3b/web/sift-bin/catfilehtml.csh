#!/bin/csh
#		catfile <filename> 
# file is in html format

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

echo Content-type: text/html
echo ""
if (-e $1) then
	cat $1
else 
	echo "Cannot open file $1"
endif


#done
exit 0
