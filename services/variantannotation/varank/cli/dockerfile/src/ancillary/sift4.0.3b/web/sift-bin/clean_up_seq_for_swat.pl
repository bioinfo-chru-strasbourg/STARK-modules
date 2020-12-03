#!/opt/local/bin/perl5

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use strict;

#clean_up_seq_for_swat.pl
#removes leading whitespaces from input file
# arguments : Input file, OUtput file
# 
my ($infile, $outfile) = @ARGV;

print "infile: $infile\n";

open (IN, $infile) || die "can't open $infile";
open (OUT, ">$outfile") || die "can't open $outfile";

my $line;
while ($line = <IN>) {
	chomp ($line);
	if ($line ne "") {
		last;
	}
}

print OUT "$line\n";
while ($line = <IN>) {
	print OUT "$line";
} 
