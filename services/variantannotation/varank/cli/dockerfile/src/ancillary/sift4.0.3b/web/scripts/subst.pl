#! /usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

$file = @ARGV[0];
$num_split = @ARGV[1];
open( SUB, $file ) || die("cannot open!");

undef %index;
while (<SUB>) {
	undef @entries;
	$entry = "";
	chomp;
	@elts = split "\t", $_;
	$rsid = @elts[0];
	$pid  = @elts[1];
	$aa1  = @elts[2];
	$aa2  = @elts[3];
	$pos  = @elts[4];
	$entry = "$aa1$pos$aa2\t\#$rsid";
	$index{$pid} = [] unless exists $index{$pid};
	push @{ $index{$pid} }, $entry;
}

$num_files = -1;
$num_dirs  = -1;
open( OUTFILE2, ">subst_locations.txt" );
foreach $pid ( sort keys %index ) {	
	$num_files++;
	if ( $num_files % $num_split == 0 ) {
		$num_dirs = $num_dirs + 1;
		mkdir "partition_$num_dirs";
		open( OUTFILE, ">partition_$num_dirs/$pid.subst" );
		@entries = @{ $index{$pid} };
		foreach $entry (@entries) {
			print OUTFILE "$entry\n";
		}
		print OUTFILE2 "$pid\t\/usr\/local\/projects\/SIFT\/datasets\/subst_files\/partition_$num_dirs\/$pid.subst\n";
		close OUTFILE;
	}
	else {
		open( OUTFILE, ">partition_$num_dirs/$pid.subst" );
		@entries = @{ $index{$pid} };
		foreach $entry (@entries) {
			print OUTFILE "$entry\n";
		}
		print OUTFILE2 "$pid\t\/usr\/local\/projects\/SIFT\/datasets\/subst_files\/partition_$num_dirs\/$pid.subst\n";
		close OUTFILE;
	}
}
close OUTFILE2;


