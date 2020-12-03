#!/usr/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
	$program_call = "/usr/bin/ls";
	$tmpdir = "../tmp";
	$pid = $$;
        system("./add_queue_entry.pl SIFT_queue $program_call $tmpdir/$pid.* > $tmpdir/$pid.test 2>/dev/null");

exit(0)
