
# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
inDir_name=$1;
inDir=`pwd`/$inDir_name
outDir_name=`echo $1 | sed s/partition/aligned/`
outDir=`pwd`/$outDir_name
mkdir $outDir
for file in `ls $1`
do
	/usr/local/projects/GOSII/pkumar/tmp/scripts/submitMuscle.sh $inDir/$file $outDir
done
