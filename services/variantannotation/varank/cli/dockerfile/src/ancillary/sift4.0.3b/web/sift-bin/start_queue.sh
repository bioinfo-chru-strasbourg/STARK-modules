#! /bin/sh

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

# start the queues

SERVER="root"
bin="/opt/www/sift/htdocs/sift-bin/"
tmp="/opt/www/sift/tmp"
queue="/opt/www/sift/queue"
# -----------------------------------------------------------------
#  Check if the user is $SERVER. If not then exit.
# -----------------------------------------------------------------
#if [ $USER != $SERVER ]; then
#  echo ""
#  echo Sorry, $USER you have to be $SERVER to run this script.
#  exit
#fi
#
cd $bin

#	Make sure the queue file exists & has group write permission
#	Don't start a queue if one is already running
qs="SIFT SIFT2"
for q in $qs; do
#	nq="`ps -fu sift | grep -c ${q}_queue`"
        nq="`ps fu | grep ${q}_queue | grep -v grep| grep quit|wc -l`"
	#echo "ps fu | grep ${q}_queue | grep -v grep| wc -l"
	#echo nq=$nq
	if [ $nq -gt 0 ]; then
		echo ${q}_queue is already running
	else
		rm -f quit$q > /dev/null 2>&1
		touch $queue/${q}_queue
		chmod 777 $queue/$q\_queue
 		$bin/run_queue.pl $queue/${q}_queue 2 quit $q > /dev/null 2>&1 &
		echo $queue/${q}_queue started
	fi
done

exit
