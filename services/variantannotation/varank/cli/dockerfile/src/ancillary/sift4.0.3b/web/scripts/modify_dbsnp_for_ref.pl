#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
$dbsnp_revised = @ARGV[0];




open (DBSNP,$dbsnp_revised) || die ("cannot open");

while (<DBSNP>){
        chomp;
        @elts = split "\t",$_;
        $rsid = @elts[0];
        $pid = @elts[1];
        $aa1 = @elts[2];
        $aa2 = @elts[3];
        $pos = @elts[4];
	print "$rsid\t$pid\t$aa1\t$aa1\t$pos\n";
}

