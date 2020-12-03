#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
$best = @ARGV[0];
$best_ref = @ARGV[1];
$all = @ARGV[2];
$all_ref = @ARGV[3];

open (BEST,$best) || die ("cannot open");
undef @elts;
while (<BEST>){
	chomp;
        @elts = split "\t",$_;
        $org = @elts[0];
        $rsid = @elts[1];
        $pid = @elts[2];
        $aa1 = @elts[3];
        $aa2 = @elts[4];
        $pos = @elts[5];
        $val = @elts[6];
        $sub = @elts[7];
        $pred = @elts[8];
        $score = @elts[9];
        $median = @elts[10];
        $num_pos = @elts[11];
        $num_aln = @elts[12];
        $best_rec = "$rsid\t$pid\t$sub\t$aa2\tBEST_HIT\t$pred\t$score\t$median\t$num_pos\t$num_aln";
	$best_index{$rsid} = $best_rec;
	push @rsid_arr,$rsid;
}



open (BEST_REF,$best_ref) || die ("cannot open");
undef @elts;
while (<BEST_REF>){
	chomp;
	@elts = split "\t",$_;
	$org = @elts[0];
	$rsid = @elts[1];
	$pid = @elts[2];
	$aa1 = @elts[3];
	$aa2 = @elts[4];
	$pos = @elts[5];
	$val = @elts[6];
	$sub = @elts[7];
	$pred = @elts[8];
	$score = @elts[9];
	$median = @elts[10];
	$num_pos = @elts[11];
	$num_aln = @elts[12];

	$best_ref_rec = "$aa1\tBEST_HIT\t$pred\t$score\t$median\t$num_pos\t$num_aln";
	$best_ref_index{$rsid} = $best_ref_rec;

}

open (ALL,$all) || die ("cannot open");
undef @elts;
while (<ALL>){
	chomp;
        @elts = split "\t",$_;
        $org = @elts[0];
        $rsid = @elts[1];
        $pid = @elts[2];
        $aa1 = @elts[3];
        $aa2 = @elts[4];
        $pos = @elts[5];
        $val = @elts[6];
        $sub = @elts[7];
        $pred = @elts[8];
        $score = @elts[9];
        $median = @elts[10];
        $num_pos = @elts[11];
        $num_aln = @elts[12];

        $all_rec = "$aa2\tALL_HITS\t$pred\t$score\t$median\t$num_pos\t$num_aln";
        $all_index{$rsid} = $all_rec;
}

open (ALL_REF,$all_ref) || die ("cannot open");
undef @elts;
while (<ALL_REF>){
	chomp;
        @elts = split "\t",$_;
        $org = @elts[0];
        $rsid = @elts[1];
        $pid = @elts[2];
        $aa1 = @elts[3];
        $aa2 = @elts[4];
        $pos = @elts[5];
        $val = @elts[6];
        $sub = @elts[7];
        $pred = @elts[8];
        $score = @elts[9];
        $median = @elts[10];
        $num_pos = @elts[11];
        $num_aln = @elts[12];

        $all_ref_rec = "$aa1\tALL_HITS\t$pred\t$score\t$median\t$num_pos\t$num_aln";
        $all_ref_index{$rsid} = $all_ref_rec;
}

$best_rec = "";
$best_ref_rec = "";
$all_rec = "";
$all_ref_rec = "";
foreach $rsid (@rsid_arr){
	$best_rec = $best_index{$rsid};
	$best_ref_rec = $best_ref_index{$rsid};
	$all_rec = $all_index{$rsid};
	$all_ref_rec = $all_ref_index{$rsid};
	print "$best_rec\t$best_ref_rec\t$all_rec\t$all_ref_rec\n";
}

