#! /usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

$file = @ARGV[0];
open( OUTFILE1, ">errored.list" );
open( OUTFILE2, ">noblink.list" );
open( FASTA, $file ) || die("can't open");

$flag = 0;
while (<FASTA>) {
	$errored = "";
	$noblink = "";
	chomp;
	if ($_ =~ /Searching/){
		$flag = 1;
		$search_line = $_;
		$search_line_flag = 0;
		next;
	}
	elsif ($_ =~ /Unable to/i || $_ =~ /uninitialized value/i){
		$flag = 0;
		if ($_ =~ /Unable/i){
			$search_line =~ /Searching best hits for (.+)/;
			$errored = "$1";
			print OUTFILE1 "$errored\n";
		}	
		else{
			$search_line =~ /Searching best hits for (.+)/;
                        $noblink = "$1";
                        print OUTFILE2 "$noblink\n";

		}	
	}
	if ($flag == 1 ){
		if ($search_line_flag == 0){
			print "$search_line\n";
	                print "$_\n";
        	        $search_line_flag = 1;
		}
		else{
			print "$_\n";
		}
	}
}
