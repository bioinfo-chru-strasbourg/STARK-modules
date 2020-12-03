#!/usr/bin/perl
use strict;

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software

my $infile  = $ARGV[0];
my $query_outfile = $ARGV[1];
my $rest_outfile = $ARGV[2];

#my $outfile = $infile . ".firstseqremoved";

open (DB, $infile) || die "can't open $infile";
my $seekpos = tell DB;
my @titles; my @old_sequences;

my ($title, $sequence) = read_a_sequence (\$seekpos);
open (OUT2, ">$query_outfile") || die "can't open $query_outfile";
my $tmp = format_seq_for_printing ($title, $sequence);
print OUT2 $tmp;
close (OUT2);

while (!eof (DB)) {
        ($title, $sequence) = read_a_sequence (\$seekpos);
        push (@titles, $title);
        push (@old_sequences, $sequence);
}
close (DB);

open (OUT, ">$rest_outfile") || die "can't open $rest_outfile";
for (my $seq_index= 0; $seq_index < @old_sequences ; $seq_index++) {
        my $tmp_string = format_seq_for_printing ($titles[$seq_index],
                                $old_sequences[$seq_index]);
        print OUT $tmp_string;
}
close (OUT);


sub format_seq_for_printing{
        my ($id, $seq) = @_;
        my $newseq = ">$id\n";
        $_ = $seq;
        my $seq_length = tr/[A-Za-z\-]//;
#        print "seqlength $seq_length\n";
        my $loc = 0;
        while ($loc < $seq_length) {
                $newseq .=  substr ($seq, $loc, 60) . "\n";
                $loc += 60;
        }
        return ($newseq);
} # end format_seq_for_printing


sub
read_a_sequence
{

        my ($seekpos) = @_;

        seek DB, ${$seekpos}, 0;
        if (!eof (DB)) {
                $_ = <DB>;
                $$seekpos = tell DB;
        }
        chomp;
        my $firstline = substr ($_, 1, 10000);
        my $seq = "";
        do {
                if (!eof (DB)) {
                        $_ = <DB>;
                        chomp;
                        unless ($_ =~ /^>/) {
                                $seq .= $_;
                                $$seekpos = tell DB;
                        }
                }
        } while (!eof (DB) && !($_ =~ /^>/)) ;
        $seq =~ s/\s+//g;
        return ($firstline, $seq);
}


