
# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
for batch in `cat $1`
do
        cd $batch/renamed/query
        for file in `ls *.seq`; do pid=`echo $file | sed s/.seq//`;path=`pwd`;echo -e $pid "\t" $path; done >> ../../../../query_location_$2.txt
        cd ../../../
done

