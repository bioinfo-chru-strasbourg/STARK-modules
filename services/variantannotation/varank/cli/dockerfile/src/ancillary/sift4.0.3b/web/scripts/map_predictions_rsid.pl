#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

$sub_locs = @ARGV[0];
$pre_locs = @ARGV[1];
$dbsnp = @ARGV[2];
$val = @ARGV[3];

open (PRE,$pre_locs) || die ("Can't open");

while (<PRE>){
	chomp;
	@p_elts = split "\t",$_;
	$pid = @p_elts[0];
	$pid =~ s/\..+?//g;
	$p_loc = "@p_elts[1]";
	$p_index{$pid} = $p_loc;
	push @pid_arr,$pid;
}

open (SUB,$sub_locs) || die ("Can't open");
while (<SUB>){
        chomp;
        @s_elts = split "\t",$_;
        $pid = @s_elts[0];
        $pid =~ s/\..+?//g;
        $s_loc = @s_elts[1];
	$s_index{$pid} = $s_loc;
}

open (DBSNP,$dbsnp) || die ("Can't open");
while (<DBSNP>){
        chomp;
        @d_elts = split "\t",$_;
        $rs = @d_elts[0];
        $d_index{$rs} = $_;
	#print "$_\n";
}

open (VAL,$val) || die ("Can't open");
while (<VAL>){
        chomp;
        @v_elts = split "\t",$_;
        $rs = @v_elts[0];
	$val_status = @v_elts[1];
        $v_index{$rs} = $val_status;
        #print "$_\n";
}




foreach $pid (@pid_arr){
	undef @subst_arr;
	undef @pred_arr;
	undef @subst_arr_sorted;
	undef @F;
	undef @sort_col;
	undef @elts2;
        chomp;

	$s_loc = $s_index{$pid};
	$p_loc = $p_index{$pid};
	#print "$p_loc\n";
	#print "$pid\t$s_loc\n";
	open (subst,$s_loc) || die ("cannot open $s_loc");
	while (<subst>){
		chomp;
		@elts2 = split "\t", $_;
		$code = @elts2[0];
		$rsid = @elts2[1];
		if ( $code =~ /(\d+)/ ) {
			$pos = $1;
		}
		push @subst_arr, "$rsid\t$code\t$pos";
		#print "$rsid\t$code\t$pos\n";
		$column = 2;
	}
	
	foreach (@subst_arr) {
        	@F = split /\t/, $_;
                push @sort_col, $F[$column];
        }
        @subst_arr_sorted = @subst_arr[ sort { $sort_col[$a] <=> $sort_col[$b] } 0 .. $#sort_col ];

	
	open (pred,$p_loc) || die ("cannot open");
	while (<pred>){
                chomp;
		if ($_ !~ /WARNING/i){
			push @pred_arr,$_;
		} 
        }
	for ($x = 0; $x < scalar @subst_arr_sorted; $x++){
		@elts = split "\t",@subst_arr_sorted[$x];
		$rsid_temp = @elts[0];
		$rsid_temp =~ /\#(rs.+)/;
		$rsid = $1;
		#print "$rsid\n";
		$dbsnp_rec = $d_index{$rsid};
		$val_status = $v_index{$rsid};
		print "$dbsnp_rec\t$val_status\t@pred_arr[$x]\n";
	}

}


