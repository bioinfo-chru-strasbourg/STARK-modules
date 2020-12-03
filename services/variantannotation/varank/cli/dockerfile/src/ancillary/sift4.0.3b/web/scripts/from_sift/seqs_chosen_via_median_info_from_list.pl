#!/usr/local/bin/perl

$file = @ARGV[0];

open (FILE,$file) || die ("Cannot Open");

while (<FILE>){
	system("/usr/local/projects/SIFT/scripts/from_sift/seqs_chosen_via_median_info_pkumar.csh $_");
}
unlink ($file);
