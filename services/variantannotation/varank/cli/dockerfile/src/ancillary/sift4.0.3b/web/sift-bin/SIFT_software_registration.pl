#!/usr/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
$| = 1;

use Fcntl qw/:flock/;
use CGI qw/:standard/;

$filename = "../www/registered_users";
open (LOGFILE, ">>$filename");

select LOGFILE;
$1=|; # flush out buffer
flock (LOGFILE, &LOCK_EX);
print "*************************\n";
$date = `date`;
print "$date";
print "Name:\t$first_name\t$last_name\t$email\n";
print "Institution:\t$institution\n";
print "Address:\t$add1\n";
print "Address:\t$add2\n";
print "Address:\t$city\t$state\t$country\n";
if ($institution_type eq "A") {
	print "Academic\n";
} elsif ($institution_type eq "C") {
	print "Company\n";
}

flock (LOGFILE, &LOCK_UN);
select STDOUT;
close (LOGFILE);
print redirect ("sift.tar.Z");

exit;
