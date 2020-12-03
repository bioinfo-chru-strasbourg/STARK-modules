#!/bin/sh
bin="/usr/local/projects/SIFT/sift4.0/bin/"
java -Xmx500m -jar $bin/IntersectFeatures.jar $1 $2 $3 $4 $5
