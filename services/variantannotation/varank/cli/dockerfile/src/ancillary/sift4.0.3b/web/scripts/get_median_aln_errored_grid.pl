#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

$query_locs = @ARGV[0];
$aln_locs = @ARGV[1];
$sub_locs = @ARGV[2];
$outdir = @ARGV[3];
$err_locs = @ARGV[4];
open (QUERY,$query_locs) || die ("Can't open");

while (<QUERY>){
	chomp;
	@q_elts = split "\t",$_;
	$pid = @q_elts[0];
	$pid =~ s/\..+//g;
	$q_loc = @q_elts[1];
	$q_index{$pid} = $q_loc;
	push @pid_arr,$pid;
}

open (ALN,$aln_locs) || die ("Can't open");

while (<ALN>){
        chomp;
        @a_elts = split "\t",$_;
        $pid = @a_elts[0];
	$pid =~ s/\..+//g;
        $a_loc = @a_elts[1];
        $a_index{$pid} = $a_loc;
}

open (SUB,$sub_locs) || die ("Can't open");

while (<SUB>){
        chomp;
        @s_elts = split "\t",$_;
        $pid = @s_elts[0];
        $pid =~ s/\..+//g;
        $s_loc = @s_elts[1];
        $s_index{$pid} = $s_loc;
}

open (ERR,$err_locs) || die ("Can't open");
while (<ERR>){
        chomp;
        @e_elts = split "\t",$_;
        $pid = @e_elts[0];
        $pid =~ s/\..+?//g;
        $e_loc = @e_elts[1];
	$e_loc =~ /(\/usr\/local\/projects\/GOSII\/pkumar\/tmp\/result\/\/partition.+?)\/.+?_result/;
	$e_loc_new = $1;
        $e_index{$pid} = $e_loc;
	$q_loc = $q_index{$pid};
        $a_loc = $a_index{$pid};
        $s_loc = $s_index{$pid};
	print "$pid\n";
	system("qsub -P 08010 -j y -V -cwd  -o $e_loc_new -l \"arch=lx*64\" /usr/local/projects/GOSII/pkumar/tmp/scripts/from_sift/seqs_chosen_via_median_info_pkumar.csh $q_loc $a_loc s_loc $pid $e_loc_new");
}



