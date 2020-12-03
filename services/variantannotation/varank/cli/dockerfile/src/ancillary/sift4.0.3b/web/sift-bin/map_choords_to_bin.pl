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
open (OUTPAGE,">>$tmp/debug.txt") || die ("cannot open outpage");
print OUTPAGE "BINS LIST: $bins_list";
close(OUTPAGE);

open (BINS, $bins_list) || die ("Cannot open bins list: $bins_list");
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
$count = 0;
while (<COORDS>){
	$count++;
	chomp;
	$_ =~ s/^\s+//g;
	$_ =~ s/\s+$//g;
	@elts = split /\,/, $_;
	if (scalar @elts < 5 || scalar @elts > 6){
		next;
	}
	$chr = @elts[0];
	$beg = @elts[1];
	$end = @elts[2];
	$orn = @elts[3];
	$allele = @elts[4];
	$comment = @elts[5];
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
			$outfile = "$pid\_snps_chr$chr\_bin$i\_$bin_start\_$bin_end.txt";
			if (-e "$tmp/$outfile"){
				open (OUTFILE, ">>$tmp/$outfile");	
				print OUTFILE "$count\t$beg\t$end\t$orn\t$allele\t$comment\n";
				
			}
			else{
				open(OUTFILE, ">$tmp/$outfile") || die ("Cannot open new bin file");
				print OUTFILE "$count\t$beg\t$end\t$orn\t$allele\t$comment\n";
				$database= "Human_CHR$chr.sqlite";
				$table = "chr$chr\_$bin_start\_$bin_end";
				push(@exec_arr,"$chr\t$i\t$outfile\t$database\t$table");
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
