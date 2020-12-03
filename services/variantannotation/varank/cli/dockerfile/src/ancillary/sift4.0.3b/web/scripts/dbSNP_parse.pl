#!/usr/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
$flatFile    = @ARGV[0];
$summaryFile = @ARGV[1];
open( SUMMARY, $summaryFile ) || die("not Possible");
while (<SUMMARY>) {
	if ( $_ =~ m/(rs\d+)/ ) {
		$rsid = $1;
	}
	if ( $_ =~ m/ACC=.*(NP_\d+)/ || $_ =~ m/ACC=.*(XP_\d+)/ ) {
		$proteinid = $1;
		$index{$rsid} = $proteinid;
	}
}


$rsid          = "";
$proteinid     = "";
$aa_position   = "";
$aa1           = "";
$aa2           = "";
open( FLAT, $flatFile ) || die("can't open");
while (<FLAT>) {
	chomp;
	if ( $_ =~ m/(^rs.\d+)/ ) {
		print "$rsid\t$proteinid\t$aa1\t$aa2\t$pos\n";
		$rsid          = $1;
		$proteinid     = @index{$rsid};
		$aa_position   = "";
		$aa1           = "";
		$aa2           = "";		
	}

	if ( $_ =~ m/fxn-class=coding-nonsynonymous/ ) {
		if ($aa2 eq ""){
			if ( $_ =~ m/aa_position=(\d+)/ ) {
                        	$pos = $1;
                	}
                	if ( $_ =~ m/residue=(.+?) \|/ ) {
                        	$aa2 = $1;
                	}
		}
		else{
			next;
		}
	
	}
	if ( $_ =~ m/fxn-class=reference/ ) {
		if ($aa1 eq ""){
			if ( $_ =~ m/residue=(.+?) \|/ ) {
                        	$aa1 = $1;
   	             	}
		}
		else{
			next;
		}
	}
}
print "$rsid\t$proteinid\t$aa1\t$aa2\t$pos\n";
