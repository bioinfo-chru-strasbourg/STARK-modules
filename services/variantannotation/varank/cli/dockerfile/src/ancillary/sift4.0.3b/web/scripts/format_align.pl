#! /usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

$file = @ARGV[0];

open (ALN,$file) || die ("cannot open!");

#put all seqs and header in array
while (<ALN>){
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
$head_query = @headers[0];
$seq_query = @seqs[0];

#################################

# remove/remember leading and trailing gaps from query

$temp = $seq_query;
$reversed = reverse ($temp);

@chars = split ("",$temp);
$beg = 0; $end = scalar @chars - 1;
$len = scalar @chars;
for ($i = 0; $i < scalar(@chars); $i++){
        if (@chars[$i] eq '-'){
                $beg = $beg + 1;
        }
        else{
                last;
        }
}

@chars = split ("",$reversed);
for ($i = 0; $i < scalar(@chars); $i++){
        if (@chars[$i] eq '-'){
                $end = $end - 1;
        }
        else{
                last;
        }
}

$seq_query_copy= substr($seq_query,$beg,$end-$beg+1);

#################################

#index indels from query seq

@chars = split("",$seq_query_copy);
for ($pos = 0; $pos < length($seq_query_copy); $pos++){
        if (@chars[$pos] eq  '-'){
                $indel_table{$pos}=1;
        }
        else {
                $indel_table{$pos}=0;
        }
}

################################

#process all seqs

for ($i = 0; $i < scalar (@seqs); $i++){
        $header = @headers[$i];
        $seq = @seqs[$i];
        $seq = substr($seq,$beg,$end-$beg+1);
        @chars = split ("",$seq);
        for ($pos = 0; $pos < scalar @chars; $pos++){
                if ($indel_table{$pos}  == 0){
                        $new_seq.=@chars[$pos];
                }                
        }
        
        ($a,$b)= $new_seq =~ m/^(\-*)\S.*?(\-*)$/;
        	
        ($ax = $a) =~ s/\-/X/g;
		($bx = $b) =~ s/\-/X/g;
		
		$new_seq =~ s/^$a/$ax/;
		$new_seq =~ s/$b$/$bx/;
        
        #$new_seq =~ s/\-+$/X/;
        $new_seq =~ s/(.{1,60})/$1\n/g;
        print "$header\n$new_seq\n";
        
        $new_seq = "";
}

