#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

$file = @ARGV[0];

open (FILE,$file);

while (<FILE>){
	if ($_ =~ /homologs/ && $flag==1){
		print $_;
        }

	if ($_ =~ /homologs/){
		$flag = 1;	
		$homolog = $_;
		
	}
	if ($_ =~ /\>gi/ && $flag==1){
		$gi = $_;
		$flag=0;

	}
}
