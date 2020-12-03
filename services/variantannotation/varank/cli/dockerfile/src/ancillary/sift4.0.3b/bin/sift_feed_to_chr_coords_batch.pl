#!/usr/local/bin/perl
my $SIFT_HOME = $ENV{'SIFT_HOME'};
my $bin             = "$SIFT_HOME/bin";

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software

$batch_file = $ARGV[0];
$tmp = $ARGV[1];
$Variation_db_dir = $ARGV[2];
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

	system ("$bin/snv_db_engine.pl $tmp $Variation_db_dir $master_pid $pid $chrfile $org $seq_identity_filter $last_partition $COORD_SYSTEM $output_options $address");	
}
close(BATCH_FILE);
