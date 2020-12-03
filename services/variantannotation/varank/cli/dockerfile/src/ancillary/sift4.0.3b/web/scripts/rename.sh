
# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
for part in `cat $1`
do
	cd $part
	ls -l | grep homologs | awk '{print $NF}' > list
	mkdir renamed
	for file in `cat list`; do pid=`head -n 2 $file | grep ">" | awk -F\| '{print $4}'`;cp $file renamed/$pid.alignedfasta; done
	rm list;
	cd ../
done
