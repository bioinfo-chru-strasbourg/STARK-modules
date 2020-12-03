
# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

qsub -P 08010 -j y -V -cwd  -l "arch=lx*64"  /usr/local/projects/GOSII/pkumar/tmp/scripts/muscle.sh $1 $2
