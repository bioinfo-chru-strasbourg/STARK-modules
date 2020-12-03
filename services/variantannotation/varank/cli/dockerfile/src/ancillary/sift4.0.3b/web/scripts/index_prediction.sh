
# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
for batch in `cat $1`
do
        cd $batch
	ls -l | grep result | awk '{print $NF }' > rlist
	for result in `cat rlist`
	do
		pid=`echo $result | awk -F_result '{print $1}'`
		file=`ls $result/*.prediction`
		path=`pwd`
		echo -e "$pid\t$path/$file"	
	done
        cd ../
done

