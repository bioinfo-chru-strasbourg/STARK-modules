#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
$file = @ARGV[0];
$queryseq = @ARGV[1];
$database = @ARGV[2];

open (QUERYSEQ,">$queryseq");
open (DATABASE,">$database");

open (FILE,$file);
$count = 0;
while (<FILE>){
	if ($flag == 1){
		if ($prev ne ""){
			print DATABASE "$prev";
		}
		print DATABASE "$_";	
		$prev = "";
	}
	else{
		if ($_ =~ /\>/){
			$count ++;
			if ($count == 2){
				$prev = $_;
				$flag = 1;
				next;
			}
			else{
				$prev = "";
			}
		}	
		print QUERYSEQ "$_";
	}
	
}	
