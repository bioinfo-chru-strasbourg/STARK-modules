#!/usr/local/bin/perl -w

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use strict;
use LWP::Simple;

# retrieves sequence from NCBI database

my $gi_number = $ARGV[0];
my $outfile = $ARGV[1];

my $best_hits_or_all_hits_option = "";
if (@ARGV >= 3) {
        $best_hits_or_all_hits_option = $ARGV[2];
}

my $URL = "http://www.ncbi.nlm.nih.gov/sutils/blink.cgi?pid=" . $gi_number;
if ($best_hits_or_all_hits_option eq "BEST") {
        $URL .= "&cut=100&org=1";
}
my $NCBI_content = get ($URL);

if (!defined ($NCBI_content) || $NCBI_content eq "") {
        print "ERROR!  Unable to contact NCBI right now.  Please try later.\n";
        exit (-1);
}
if ($NCBI_content =~ /No hits found/) {
        # ERROR
        print "ERROR!  No BLAST hits were found for gi:$gi_number by NCBI BLink .  Please check that this sequence is a protein.<BR>You can also submit your pr
otein sequence directly to SIFT and have SIFT choose the protein sequences  <A HREF=\"http://blocks.fhcrc.org/sift/SIFT_seq_submit2.html\">here.</A>\n";

        exit (-1);
}
if ($NCBI_content =~ /ERROR/) {
        print "ERROR trying to retrieve sequences from NCBI!  Please follow this <A HREF=$URL>link</A HREF> to NCBI to see the error. \n";
        exit (-1);
}

my @lines = split ("\n", $NCBI_content);
if (@lines == 0) {
        print "ERROR! Unable to retrieve sequences from NCBI. Please try later.\n";
        exit (-1);

}

my $lineindex = 0;
while ($lineindex < @lines &&
        !($lines[$lineindex]  =~ /TABLE/)) {
        $lineindex++;
}
$lineindex++;

do {
        $lineindex++;
} until ($lines[$lineindex] =~ /img src/ || $lineindex >= @lines);

my %seen;
my $gi_hits = $gi_number . ",";
while ($lineindex < @lines ) {
        my $loc = index ($lines[$lineindex], "val=");
        while ($loc >= 0) {
                my $gi = substr ($lines[$lineindex], $loc + 4);
                $loc = index ($gi, " ");
                $gi = substr ($gi, 0, $loc);
                if (!exists ($seen{$gi})) {
                        $gi_hits .= $gi . "," ;
                        $seen{$gi} = 1;
                }
                $lines[$lineindex] = substr ($lines[$lineindex], $loc+1, 100000);
                $loc = index ($lines[$lineindex], "val=");
        }
        if ($loc < 0) {
                $loc = index ($lines[$lineindex], "value=");
                if ($loc >= 0 && $lines[$lineindex] =~ /list_uids/) {
                        my $list = substr ($lines[$lineindex], $loc + 6);
                        $list =~ s/>//;
                        my @gis = split (/\,/, $list);
                        for (my $i = 0; $i < @gis; $i++) {
                                $gi_hits .= $gis[$i] . ",";
                                $seen{$gis[$i]} = 1;
                        }
                }
        }
        $lineindex++;
}
$gi_hits =~ s/\,$//;
print "$gi_hits\n";

my $seq_URL = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&email=sift\@fhcrc.org&id=$gi_hits&rettype=fasta&retmode=text";
my $seq = get ($seq_URL);
print $seq;
open (OUT, ">$outfile") || die "can't open $outfile";
print OUT $seq ;
close (OUT);

exit (0);

