#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
my $tmp             = "/opt/www/sift/tmp";
$coords_list = @ARGV[0];
$pid = @ARGV[1];
if (scalar @ARGV != 2){
	print "Usage: perl map_choords_to_bin.pl coords_file pid\n";
	exit;
}

open(COORDS,$coords_list) ||  die("Cannot open coords list");
open (MAP,">>$tmp/$pid.snp_chr_map_file.txt") || die ("Cannot open $pid.snp_chr_map_file.txt");
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
	$bin_file = "$tmp/$pid.snps_chr$chr.txt";
	if (!exists ($index_bin_seen{$chr})){
		$index_bin_seen{$chr} = 1;
		print MAP "$chr\t$bin_file\n";
	}
	open (BIN_FILE,">>$bin_file") || die("Cannot open $bin_file");
	print BIN_FILE "$count\t$beg\t$end\t$orn\t$allele\n";
}
close(BIN_FILE);
close (MAP);
