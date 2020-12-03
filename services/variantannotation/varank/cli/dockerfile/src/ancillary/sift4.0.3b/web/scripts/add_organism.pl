#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

$pred = @ARGV[0];
$org = @ARGV[1];

open (ORG,$org) || die ("Can't open");
while (<ORG>){
        chomp;
        @o_elts = split /\t/,$_;
        $rsid = @o_elts[0];
        $organism = @o_elts[1];
	$rsid =~ s/^\s+//;
	$rsid =~ s/\s+$//;
	$organism =~ s/^\s+//;
        $organism =~ s/\s+$//;
	#print "$rsid$organism\n";
        $o_index{$rsid} = $organism;
}

open (PRE,$pred) || die ("Can't open");

while (<PRE>){
        chomp;
        @p_elts = split "\t",$_;
        $rsid = @p_elts[0];
        $organism = $o_index{$rsid};
        print "$organism\t$_\n";
}



