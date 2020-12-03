
# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
for batch in `cat $1`
do
	cd $batch/renamed
	for file in `ls *.alignedfasta`; do pid=`echo $file | sed s/.alignedfasta//`;path=`pwd`;echo -e $pid "\t" $path; done >> ../../../aln_locations_$2.txt
	cd ../../
done
