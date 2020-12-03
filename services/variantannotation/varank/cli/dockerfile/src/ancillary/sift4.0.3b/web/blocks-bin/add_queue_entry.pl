#!/usr/bin/perl
#	add_queue.pl:	Add an entry to a queue file = $1
#	Other arguments are written to the queue file
#	EG "add_queue.pl LAMA_queue 'LAMA block_file database'"
#	Executes lockfile

# usage
if (@ARGV < 2) {
    print "usage: $0 queuefile command ...\n";
    exit(1);
}

$queuefile = shift(@ARGV);

$lockqueuefile = "$queuefile.lock";

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
