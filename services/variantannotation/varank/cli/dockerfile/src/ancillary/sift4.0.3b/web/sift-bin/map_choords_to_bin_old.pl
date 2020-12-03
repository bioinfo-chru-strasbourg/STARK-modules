#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
$coding_info_dir = "/usr/local/projects/SIFT/siftdev/coding_info";
my $tmp             = "/opt/www/sift/tmp";
$bins_list = @ARGV[0];
$coords_list = @ARGV[1];
$pid = @ARGV[2];
if (scalar @ARGV != 3){
	print "Usage: perl map_choords_to_bin.pl bin_file coords_file pid\n";
	exit;
}
open (BINS, $bins_list) || die ("Cannot open bins list");
$first_bin = 1;
$prev_chr = 1;
while(<BINS>) {
	chomp;
	@elts = split /\t/, $_;
	$bin = @elts[0];
	$chr = @elts[1];
	$beg = @elts[2];
	$end = @elts[3];
	$index_bin{$bin} = "$chr\t$beg\t$end";
        if ($chr ne $prev_chr){
		$index_chr{$prev_chr} = "$first_bin\t$last_bin";
		$prev_chr = $chr;
		$first_bin = $bin;
		$last_bin = $bin;
	}
	else{
		$last_bin = $bin;	
	}
}
$index_chr{$prev_chr} = "$first_bin\t$last_bin"; 

open(COORDS,$coords_list) ||  die("Cannot open coords list");
open (COORDS_MAP, ">$tmp/$pid\_old_new_coords_map.txt");
$count = 0;
while (<COORDS>){
	$count++;
	chomp;
	$_ =~ s/^\s+//g;
	$_ =~ s/\s+$//g;
	@elts = split /\,/, $_;
	if (scalar @elts != 5){
		next;
	}
	$chr = @elts[0];
	$beg = @elts[1];
	$end = @elts[2];
	$orn = @elts[3];
	$allele = @elts[4];
	if ($end -$beg != 1){
		next;
	}
	$first_bin = (split /\t/,$index_chr{$chr})[0];
	$last_bin = (split /\t/,$index_chr{$chr})[1];
	for ($i = $first_bin; $i <= $last_bin; $i++){
		@elts = split /\t/, $index_bin{$i};
		$bin_start =@elts[1];
		$bin_end = @elts[2];
		if ($beg >= $bin_start && $beg <= $bin_end){
			$new_beg = $beg-$bin_start+1;
			$new_end = $end-$bin_start+1;
			$outfile = "$pid\_snps_chr$chr\_bin$i\_$bin_start\_$bin_end.txt";
			if (-e "$tmp/$outfile"){
				open (OUTFILE, ">>$tmp/$outfile");	
				print OUTFILE "$count\t$new_beg\t$new_end\t$orn\t$allele\n";
				print COORDS_MAP "$chr,$beg,$end,$orn,$allele\t$chr,$new_beg,$new_end,$orn,$allele\n";
				
			}
			else{
				open(OUTFILE, ">$tmp/$outfile") || die ("Cannot open new bin file");
				print OUTFILE "$count\t$new_beg\t$new_end\t$orn\t$allele\n";
				print COORDS_MAP "$chr,$beg,$end,$orn,$allele\t$chr,$new_beg,$new_end,$orn,$allele\n";
				$coding_info_file = "$chr.$bin_start-$bin_end.coding_info";
				$fasta_file = "$chr.$bin_start-$bin_end.fasta";
				push(@exec_arr,"$chr\t$i\t$outfile\t$coding_info_file\t$fasta_file");
			}
			last;
		}  
		else{
			next;
		}
	}
}

foreach (keys %index_bin){
	#print "$_\t$index_bin{$_}\n";
}
open (MAPFILE, ">$tmp/$pid\_snp_chr_map_file.txt");
foreach (@exec_arr){
	print MAPFILE "$_\n";	
}
