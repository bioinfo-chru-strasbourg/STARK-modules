#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
$dbsnp = @ARGV[0];
$ver = @ARGV[1];


open (VER,$ver) || die ("cannot open");

while (<VER>){
	chomp;
	@elts = split "\t", $_;
	$pid = @elts[0];
	$rsid = @elts[1];	
	$flag1 = @elts[4];
	$flag2 = @elts[5];
	$index1{$rsid} = $flag1;
	$index2{$rsid} = $flag2;
}


open (DBSNP,$dbsnp) || die ("cannot open");

while (<DBSNP>){
        chomp;
        @elts = split "\t",$_;
        $rsid = @elts[0];
        $pid = @elts[1];
        $aa1 = @elts[2];
        $aa2 = @elts[3];
        $pos = @elts[4];
	$flag1 = $index1{$rsid};
	$flag2 = $index2{$rsid};
	if ($flag1 =~/TRUE/ && $aa1 !~ /\*/ && $aa1 !~ /\-/ && $aa1 !~ /X/ && $aa1 ne "" ){
		if ($aa2 !~ /\*/ && $aa2 !~ /\-/ && $aa2 !~ /X/ && $aa2 ne ""){
			print "$rsid\t$pid\t$aa1\t$aa2\t$pos\n";
		}
		else{
			print "$rsid\t$pid\t$aa1\t$aa1\t$pos\n";
		}
	}
	else{
		if ($flag2 =~ /TRUE/ && $aa2 !~ /\*/ && $aa2 !~ /\-/ && $aa2 !~ /X/ && $aa2 ne "" ){
			if ($aa1 !~ /\*/ && $aa1 !~ /\-/ && $aa1 !~ /X/ && $aa1 ne ""){
				print "$rsid\t$pid\t$aa2\t$aa1\t$pos\n";
			}
			else{
				print "$rsid\t$pid\t$aa2\t$aa2\t$pos\n";
			}
		}
	}
}

