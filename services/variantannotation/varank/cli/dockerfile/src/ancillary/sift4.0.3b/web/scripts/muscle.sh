
# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
muscle="/home/pkumar/SIFT/software/muscle3.6/muscle"
format_align="/usr/local/projects/GOSII/pkumar/tmp/scripts/format_align.pl"
inFile=$1
inFile_name=`echo $inFile | awk -F\/ '{print $NF}'`
outFile_name=$inFile_name.aln.out
outDir=$2
outFile=$outDir/$outFile_name

$muscle -in $inFile -out $outFile_name -maxiters 2 -stable -quiet
mv $outFile_name $outDir
$format_align $outFile > $outDir/$inFile_name.aln.formatted.out
rm $outFile
