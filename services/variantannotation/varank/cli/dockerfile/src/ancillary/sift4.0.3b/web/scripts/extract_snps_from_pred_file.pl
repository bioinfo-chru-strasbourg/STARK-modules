#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

$pred = @ARGV[0];
$extract = @ARGV[1];
open (PRED,$pred) || die ("Can't open");

while (<PRED>){
        chomp;
        @p_elts = split "\t",$_;
        $rsid = @p_elts[0];
        $p_index{$rsid} = $_;
}

open (EXT,$extract) || die ("Can't open");

while (<EXT>){
	chomp;
	$rsid = $_;
	$detail_rec = $p_index{$rsid};
	print "$detail_rec\n";
}
