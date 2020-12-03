#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

$pred = @ARGV[0];
open (PRE,$pred) || die ("Can't open");

while (<PRE>){
        chomp;
        @p_elts = split "\t",$_;
	splice @p_elts, 6,1;
	$line = join("\t",@p_elts);
	print "$line\n";
 
}

