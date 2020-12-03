#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
my $bin             = "/opt/www/sift/htdocs/sift-bin";
my $tmp             = "/opt/www/sift/tmp";

$batch_file = @ARGV[0];

open (BATCH_FILE ,"$batch_file");
while (<BATCH_FILE>){
	chomp;
	@elts = split /\t/, $_;
	$master_pid = $elts[0];
	$pid = $elts[1];
	$chrfile = $elts[2];
	$org = $elts[3];
	$seq_identity_filter = $elts[4];
	$COORD_SYSTEM = $elts[5];
	$address = $elts[8];
	$last_partition = $elts[6];
	$output_options = $elts[7];
	
	open (OUTPAGE,">>$tmp/debug.txt") || die ("cannot open outpage");
	print OUTPAGE "$bin/SIFT_chr_coords_submit.pl $master_pid $pid $chrfile $org $seq_identity_filter $last_partition $COORD_SYSTEM $output_options $address\n";
	close(OUTPAGE);
	system ("$bin/SIFT_chr_coords_submit.pl $master_pid $pid $chrfile $org $seq_identity_filter $last_partition $COORD_SYSTEM $output_options $address");	
}
close(BATCH_FILE);
