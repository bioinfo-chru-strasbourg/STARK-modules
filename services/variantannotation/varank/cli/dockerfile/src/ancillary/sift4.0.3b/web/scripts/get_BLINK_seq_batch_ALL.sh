
# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
for pid in `cat $1`
do
	echo Searching best hits for $pid
	/home/pkumar/SIFT/scripts/get_BLINK_seq.pl $pid out
done

