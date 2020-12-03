#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

$file = @ARGV[0];

open (PRED,$file) || die ("cannot open");

while (<PRED>){
	chomp;
	@elts = split "\t", $_;
	$pid = @elts[0];
	$loc = @elts[1];
	$fasta = "$loc\/$pid\.selected\.alignedfasta";
	$filesize = -s $fasta;
	if ($filesize < 20){
		print "$pid\t$loc\n";	
	}
}
