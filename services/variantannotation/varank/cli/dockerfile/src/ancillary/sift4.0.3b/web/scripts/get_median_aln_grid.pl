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
$num_split = @ARGV[4];
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

$num_files = -1;
$num_dirs  = -1;

@pid_arr = sort @pid_arr;
foreach $pid (@pid_arr){
	$num_files++;
        if ( $num_files % $num_split == 0 ) {
                $num_dirs = $num_dirs + 1;
                $outdir_partition = "$outdir/partition_$num_dirs";
		$result_dir = "$outdir_partition\/$pid\_result";
		$q_loc = $q_index{$pid};
        	$a_loc = $a_index{$pid};
		 if (exists $s_index{$pid}){
                        $s_loc = $s_index{$pid};
			system("qsub -P 08010 -j y -V -cwd -l \"arch=lx*64\" /usr/local/projects/SIFT/scripts/from_sift/seqs_chosen_via_median_info_pkumar.csh $q_loc $a_loc $s_loc $pid $outdir_partition");	
                }
                else{
                        next;
                }

        }
	else{
		$outdir_partition = "$outdir/partition_$num_dirs";
		$result_dir = "$outdir_partition\/$pid\_result";
                $q_loc = $q_index{$pid};
                $a_loc = $a_index{$pid};
		if (exists $s_index{$pid}){
                        $s_loc = $s_index{$pid};
                        system("qsub -P 08010 -j y -V -cwd -l \"arch=lx*64\" /usr/local/projects/SIFT/scripts/from_sift/seqs_chosen_via_median_info_pkumar.csh $q_loc $a_loc $s_loc $pid $outdir_partition");
                }
                else{
                        next;
                }

	}
}

