#!/usr/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
$summaryFile = @ARGV[0];
open( SUMMARY, $summaryFile ) || die("not Possible");
while (<SUMMARY>) {
        if ( $_ =~ m/(rs\d+)/ ) {
                $rsid = $1;
        }
        if ( $_ =~ m/VALIDATED=(.+)/ ) {
                $val = $1;
        	print "$rsid\t$val\n";
	}
}

