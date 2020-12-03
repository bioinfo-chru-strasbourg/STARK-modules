
# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
for batch in `cat $1`
do
	cd $batch/renamed
	mkdir query
	ls -l | grep alignedfasta | awk '{print $NF}' > aln.list
	for file in `cat aln.list`
	do	
		pid=`echo $file |sed s/\.alignedfasta//`
		/usr/local/projects/GOSII/pkumar/tmp/scripts/get_first_seq.pl $file > query/$pid.seq
	done
	cd ../../
done
