#!/usr/bin/perl
#	add_queue.pl:	Add an entry to a queue file = $1
#	Other arguments are written to the queue file
#	EG "add_queue.pl LAMA_queue 'LAMA block_file database'"
#	Executes lockfile

# 3/13/08 Stop queuing after maxq entries

$maxq = 1000;

# usage
if (@ARGV < 2) {
    print "usage: $0 queuefile command ...\n";
    exit(1);
}
$queuefile = shift(@ARGV);
$lockqueuefile = "$queuefile.lock";

$nq = 0;
open (IN, "<$queuefile");
while (<IN>)
{  $nq++; }
close(IN);
print STDERR "nq=$nq\n";
if ($nq > $maxq) { exit(-1); }

# lockfile will wait until it can lock the file
`./lockfile $lockqueuefile`;

# append the address and command to the queue file
open(FILE, ">>$queuefile");
print FILE "@ARGV\n";
close(FILE);

chmod(0664, $queuefile);

# remove the lock file
unlink($lockqueuefile);
exit(0);
