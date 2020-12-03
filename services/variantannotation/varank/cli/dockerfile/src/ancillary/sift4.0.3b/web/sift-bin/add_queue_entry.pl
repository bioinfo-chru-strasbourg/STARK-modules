#!/usr/local/bin/perl
#	add_queue.pl:	Add an entry to a queue file = $1
#	Other arguments are written to the queue file
#	EG "add_queue.pl LAMA_queue 'LAMA block_file database'"
#	Executes lockfile

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
$bin = "/opt/www/sift/htdocs/sift-bin/";
$tmp = "/opt/www/sift/tmp";
$lockfile = "/usr/bin/lockfile";
$queue = "/opt/www/sift/queue/";
# usage
if (@ARGV < 2) {
    print "usage: $0 queuefile command ...\n";
    exit(1);
}
&check_queue_running;
$queuefile = shift(@ARGV);
$lockqueuefile = "$queuefile.lock";

# lockfile will wait until it can lock the file
`$lockfile $lockqueuefile`;
# append the address and command to the queue file
open(FILE, ">>$queuefile") || die ("Cannot open queue file");
print FILE "@ARGV\n";
#system("more $queuefile > $tmp/remove");
#system ("echo @ARGV > $tmp/remove");
#system ("more $queuefile > $tmp/remove1");
close(FILE);
chmod(0777, $queuefile);

# remove the lock file
unlink($lockqueuefile);
exit(0);

sub check_queue_running {
	$num_queues = 0;
	$num_queues = `ps -ef | grep SIFT | grep queue|grep quit| grep -v grep| wc -l `;
	chomp $num_queues;
	if ($num_queues != 2){
		system ("$bin/start_queue.sh");	
	}
}
