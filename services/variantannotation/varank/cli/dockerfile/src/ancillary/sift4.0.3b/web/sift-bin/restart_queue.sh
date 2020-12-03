#!/bin/sh

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

#	Stop the current queue and wait

su - sift -c "touch /home/sift/bin/quitSIFT" > /dev/null 2>&1
wait 5 
rm -f /home/sift/bin/quitSIFT
su - sift -c "/home/sift/bin/start_queue.sh" > /dev/null 2>&1

exit
