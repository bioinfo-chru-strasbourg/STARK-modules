#!/usr/bin/perl
#	run_queue.pl:  Process requests from a queue			
#		$1 = name of queue file
#	Executes lockfile

# IMPORTANT VARIABLE
$stopfile = "quit";
$sleeptime = 60; # in seconds



# usage
if (@ARGV < 3)
{
    print "usage: $0 queuefile sleeptime stopfile\n";
    exit(1);
}

# get the args
$queuefile = $ARGV[0];
$sleeptime = $ARGV[1];
$stopfile = $ARGV[2];
print "$queuefile $sleeptime $stopfile\n";
#exit;


$tmpqueuefile = "$queuefile.tmp";
$lockqueuefile = "$queuefile.lock";


# test to see if the queue file exists
if (!(-e $queuefile)) {
    `touch $queuefile`;
}

# test to see if the queue file is readable
if (!(-r $queuefile)) {
    print "Queuefile $queuefile is not readable.  Exiting.\n";
    exit(1);
}


# the file exists.  Start processing the queue
# stop if the stopfile exists
while (!(-e $stopfile)) {
    $command = &readlineorsleep();
    system($command);
}


# the stopfile exists
exit(0);




sub readlineorsleep {
    local($line);

    open(FILE, $queuefile);
    
    # while there is nothing in the file, loop
    while(eof(FILE)) {
	# check to see if the stopfile says to stop
	if (-e $stopfile) {
	    exit(1);
	}
	# close file
	close(FILE);
	# sleep
	sleep $sleeptime;
	# open file
	open(FILE, $queuefile);
    }

    # make the temp file and the lockfile
    unlink($tmpqueuefile);
    # lockfile will wait until it can lock the file
    `./lockfile $lockqueuefile`;
    open(TMPFILE, ">$tmpqueuefile");

    $line = <FILE>;
    
    # copy the rest of the queuefile into the temp queue
    while (<FILE>) {
	print TMPFILE;
    }

    close(FILE);
    close(TMPFILE);
    unlink($queuefile);
    rename($tmpqueuefile, $queuefile);
    unlink($tmpqueuefile);
    unlink($lockqueuefile);

    chmod(0664, $queuefile);

    # return the line without the return
    chop($line);
    $line;
}
