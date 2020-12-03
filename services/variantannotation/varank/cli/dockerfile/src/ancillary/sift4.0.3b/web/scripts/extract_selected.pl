#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

$pid = @ARGV[0]; #protein ID
$a_loc = @ARGV[1]; #alignment file path
$outdir = @ARGV[2];
$tmpdir = "$outdir/$pid\_result"; #from seqs_chosen_via_median_info_pkumar.csh script
$selected = "$tmpdir/$pid.selected"; #contains sequence GIs to extract
$selected_aligned = "$tmpdir/$pid.selected.alignedfasta";
#print "$selected\n$selected_aligned\n";
open(OUTFILE,">$selected_aligned");
open (ALN_FILE,$a_loc) || die ("cannot open!");
while (<ALN_FILE>){
                if ($_ =~ />/ || /\|/ ){
                        chomp;
                        $header.=$_;
                        if($seq){
                                push(@seqs,$seq);
                        }
                        $seq = "";
                }
                else{
                        chomp;
                        $seq.= $_;
                        if ($header){
                                push(@headers,$header);
			}
                        $header="";
                }
}
push(@seqs,$seq);

for ($x=0; $x < scalar @seqs; $x++){
	$header = @headers[$x];
	$seq = @seqs[$x];
        @elements = split /\|/,$header;
	$gi = @elements[1];
	$gi_seq_index{$gi}=$seq;
	$gi_header_index{$gi}=$header;	
}

open (SELECTED,$selected) || die ("cannot open");
while (<SELECTED>){
	chomp;
	@elts = split /\|/,$_;	
	$gi = @elts[1];
	#print "$_\t$gi\n";
	push @gilist,$gi;
}

foreach $gi (@gilist){
	$header = $gi_header_index{$gi};
	$seq = $gi_seq_index{$gi};
	$seq =~ s/(.{1,60})/$1\n/g;
	print OUTFILE "$header\n$seq\n";
}



