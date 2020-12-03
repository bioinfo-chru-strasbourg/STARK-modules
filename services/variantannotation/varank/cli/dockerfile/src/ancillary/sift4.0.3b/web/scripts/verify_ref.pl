#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

$query_locs = @ARGV[0];
$dbsnp = @ARGV[1];

open (QUERY,$query_locs) || die ("Can't open");

while (<QUERY>){
        chomp;
        @q_elts = split "\t",$_;
        $pid = @q_elts[0];
        $pid =~ s/\..+?//g;
        $q_loc = @q_elts[1];
        $q_index{$pid} = $q_loc;
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
	$q_loc = $q_index{$pid};
	$seq = "";
	open (QUERY,$q_loc) || next; 
	while (<QUERY>){
		chomp;
		#print "$_\n";
		if ($_ eq ""){next;}
		if ($_ !~ /\>/){
			$seq=$seq.$_;
		}
	}	
	$seq =~ s/^\s+//;
	$seq =~ s/\s+$//;
	$aa_actual = substr($seq,$pos-1,1);
	if ($aa1 eq $aa_actual){
		$flag1 = "TRUE"
	}
	else {
		$flag1 = "FALSE";
	}
	
	if ($aa2 eq $aa_actual){
                $flag2 = "TRUE"
        }
        else {
                $flag2 = "FALSE";
        }

	print "$pid\t$rsid\t$aa1\t$aa_actual\t$flag1\t$flag2\n";
}
