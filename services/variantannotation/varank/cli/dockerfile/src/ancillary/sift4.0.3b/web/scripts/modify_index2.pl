#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

$file = @ARGV[0];

open (FILE,$file) or die ("Can't open");

while (<FILE>){
	chomp;
	@elts = split "\t",$_;
	$pid = @elts[0];
	$loc1 = @elts[1];
	$pid =~ s/\s+$//;
	$pid =~ s/^\s+//;
	$loc1 =~ s/\s+$//;
        $loc1 =~ s/^\s+//;
	$loc1 =~ /(\/usr\/local\/projects\/GOSII\/pkumar\/tmp\/result\/\/part.+?\/.+?result)\/.+?_result/;
	$loc2 = $1;
	print "$pid\t$loc2\n";
}
