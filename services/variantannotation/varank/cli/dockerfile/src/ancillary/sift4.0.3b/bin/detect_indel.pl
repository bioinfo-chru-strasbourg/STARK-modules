#!/usr/local/bin/perl -w
use strict;

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software

my $original_fasta = $ARGV[0];
my $modified_fasta = $ARGV[1];
my $outfile        = $ARGV[2];
my $enst           = "";
my $first          = "";
my $second         = "";
if ( $original_fasta =~ /.+?\d+_(ENST\d+)\.fa/ ) {
	$enst = $1;
}

#if ($original_fasta =~ //)
#my $enst =
my $seq1;
my $seq2;
my $length1;
my $length2;
open( ORIGINAL, "$original_fasta" ) || die("Cannot open original fasta");

while (<ORIGINAL>) {
	chomp;
	if ( $_ =~ /^\>/ ) { next }
	$seq1 .= $_;
	$length1 = length($seq1);
}

open( MODIFIED, "$modified_fasta" ) || die("Cannot open modified fasta");
while (<MODIFIED>) {
	chomp;
	if ( $_ =~ /^\>/ ) { next }
	$seq2 .= $_;
	$length2 = length($seq2);
}

#print "$seq1\n\n$seq2\n";
open( OUTFILE, "> $outfile" )
  || die("Cannot open outfile for original/modified seqs");
my $outstr = &detect_indel( $seq1, $seq2 );
print "$outstr\n";

sub detect_indel {
	my @outarr;
	my $coords_change_original;
	my $coords_change_modified;
	my $aa_change_original;
	my $aa_change_modified;
	my $seq1 = $_[0];
	my $seq2 = $_[1];
	our $offset = 0;
	our $order_reverse;
	my $original_seq;
	my $modified_seq;
	if ( length $seq2 < length $seq1 ) {
		$original_seq  = $seq1;
		$modified_seq  = $seq2;
		$order_reverse = 1;
	}
	else {
		$original_seq  = $seq2;
		$modified_seq  = $seq1;
		$order_reverse = 0;

	}
	my $index = 0;
	my $mismatch_start;
	my $mismatch_stop;
	if ($original_seq ne $modified_seq){
		while (
			substr( $original_seq, $index, 1 ) eq substr( $modified_seq, $index, 1 ) 
	  	)
		{
			$index++;
			next;
		}
		$mismatch_start = $index;
	}
	else{$mismatch_start = 0;}

	$index = 1;
	if ($original_seq ne $modified_seq){
		while (
			substr( $original_seq, -$index, 1 ) eq
			substr( $modified_seq, -$index, 1 ) )
		{
			$index++;
			next;
		}
		$mismatch_stop = length($modified_seq) - $index + 1;
	}
	else{$mismatch_stop = 0;}
	if ( $mismatch_stop < $mismatch_start ) {
		$offset = $mismatch_start - $mismatch_stop;
		$mismatch_stop += $offset;

	}
	my $mismatch            = "MISMATCH = $mismatch_start-$mismatch_stop";
	my $left_flank_modified = uc substr( $modified_seq, 0, $mismatch_start );
	my $indel_modified      =
	  lc substr( $modified_seq, $mismatch_start,
		( $mismatch_stop - $mismatch_start ) );
	my $right_flank_modified = uc substr( $modified_seq, $mismatch_stop );

	my $modified = "$left_flank_modified$indel_modified$right_flank_modified";
	$modified =~ s/.{60}(?=.)/$&\n/g;
	$modified =~ s/[a-z]+/<font color=red>$&<\/font>/g;
	if ( $order_reverse == 0 ) {
		$first = ">$enst; $mismatch\n$modified\n\n";
		$coords_change_original = "$mismatch_start-$mismatch_stop";
		$aa_change_original = join ('',substr($left_flank_modified,-5,5),$indel_modified,substr($right_flank_modified,0,5));
		if (length $right_flank_modified == 0){
			$aa_change_original = substr($aa_change_original,0,20)."*";
		}
		else{
			$aa_change_original = substr($aa_change_original,0,20);
		}
	}
	else {
		$second = ">$enst; $mismatch\n$modified\n";
		$coords_change_modified = "$mismatch_start-$mismatch_stop";
		$aa_change_modified = join ('',substr($left_flank_modified,-5,5),$indel_modified,substr($right_flank_modified,0,5));
		if (length $right_flank_modified == 0){
                        $aa_change_modified = substr($aa_change_modified,0,20)."*";
                }
                else{
                        $aa_change_modified = substr($aa_change_modified,0,20);
                }


	}

############################################################################################################

	my $tmp = $original_seq;
	$original_seq   = $modified_seq;
	$modified_seq   = $tmp;
	$index          = 0;
	$mismatch_start = -1;
	$mismatch_stop  = -1;
	if ($original_seq ne $modified_seq){
		while (
			substr( $original_seq, $index, 1 ) eq substr( $modified_seq, $index, 1 )
	 	)
		{
			$index++;
			next;
		}
		$mismatch_start = $index;
	}
	else{$mismatch_start = 0;}
	$index = 1;
	if ($original_seq ne $modified_seq){
		while (
			substr( $original_seq, -$index, 1 ) eq
			substr( $modified_seq, -$index, 1 ) )
		{
			$index++;
			next;
		}
		$mismatch_stop = length($modified_seq) - $index + 1;
	}
	else{$mismatch_stop = 0;}
	$mismatch_stop += $offset;
	$mismatch            = "MISMATCH = $mismatch_start-$mismatch_stop";
	$left_flank_modified = uc substr( $modified_seq, 0, $mismatch_start );
	$indel_modified      =
	  lc substr( $modified_seq, $mismatch_start,
		( $mismatch_stop - $mismatch_start ) );
	$right_flank_modified = uc substr( $modified_seq, $mismatch_stop );

	$modified = "$left_flank_modified$indel_modified$right_flank_modified";
	$modified =~ s/.{60}(?=.)/$&\n/g;
	$modified =~ s/[a-z]+/<font color=red>$&<\/font>/g;
	if ( $first eq "" ) {
		$first = ">$enst; $mismatch\n$modified\n\n";
		$coords_change_original = "$mismatch_start-$mismatch_stop";
		$aa_change_original = join ('',substr($left_flank_modified,-5,5),$indel_modified,substr($right_flank_modified,0,5));
		if (length $right_flank_modified == 0){
                        $aa_change_original = substr($aa_change_original,0,20)."*";
                }
                else{
                        $aa_change_original = substr($aa_change_original,0,20);
                }


	}
	else {
		$second = ">$enst; $mismatch\n$modified\n";
		$coords_change_modified = "$mismatch_start-$mismatch_stop";
		$aa_change_modified = join ('',substr($left_flank_modified,-5,5),$indel_modified,substr($right_flank_modified,0,5));
		if (length $right_flank_modified == 0){
                        $aa_change_modified = substr($aa_change_modified,0,20)."*";
                }
                else{
                        $aa_change_modified = substr($aa_change_modified,0,20);
                }

	}

	print OUTFILE "$first$second";
	$outstr = "$coords_change_original\t$coords_change_modified\t$aa_change_original\t$aa_change_modified\t$length1\t$length2";
	return $outstr;

}
close(OUTFILE);

