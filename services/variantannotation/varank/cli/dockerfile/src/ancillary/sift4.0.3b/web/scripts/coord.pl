#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

$file = @ARGV[0];
$coord = @ARGV[1];

open (FILE, $file);
$len = 0;
while (<FILE>){
	chomp;
	if ($_ =~ />/){
		next;
	}
	$len = $len + length($_);
	if ($coord <= $len){
	
		$local_coord = $coord - ($len - length($_));
		print "$coord\t$len\t$local_coord\n";
		$nu = substr $_, $local_coord-1,1;
		print $nu,"\n";
		last;
	}
}

