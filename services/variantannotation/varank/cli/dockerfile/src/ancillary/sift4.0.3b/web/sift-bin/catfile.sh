#!/bin/sh
#		catfile <filename>

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

echo Content-type: text/html
echo

echo "<PRE>"
if [ $1 ]; then
	cat $1
fi
echo "</PRE>"
