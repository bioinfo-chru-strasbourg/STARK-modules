
# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
for partition in `cat $1`
do
	/usr/local/projects/GOSII/pkumar/tmp/scripts/submitMuscle_batch_partition.sh $partition
done
