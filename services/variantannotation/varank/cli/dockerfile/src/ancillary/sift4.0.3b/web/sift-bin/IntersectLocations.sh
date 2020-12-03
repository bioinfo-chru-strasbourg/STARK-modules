#!/bin/sh

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
bin="/opt/www/sift/htdocs/sift-bin"
java -Xmx500m -jar $bin/IntersectFeatures.jar $1 $2 $3 $4 $5
