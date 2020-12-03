# COPYRIGHT 2008. This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software

use POSIX;
use Class::Struct;

struct Region => {
        chr=> '$',
        beg => '$',
        end => '$', 
	orientation => '$',
	id => '$',
	cum_length => '$',
	count => '$', 
	sequence => '$'
	
};

sub
merge_regions
{
        my (@old_regions) = @_;
        if (@old_regions == 0) { return;}
        my $chr;
                my @a = @old_regions;
                my @newlist;
                my $added_last = 0;
                for (my $i = 1; $i < @a; $i++) {
                        my $current_r = $a[$i-1];
                        my $next_r = $a[$i];
                        my $first_merge = 1;
                        while ($next_r->beg <= ($current_r->end  + 1) &&
                                                $i < @a)
                        {
                                if ($next_r->end > $current_r->end) {
                                        $current_r->end ($next_r->end);
                                }
                                $i++;
                                if ($i < @a) {
                                        $next_r = $a[$i];
                                }
                           }
                        push (@newlist, $current_r);
                        if ($i == @a) { $added_last = 1;}
                }
                if (!$added_last) {
                        push (@newlist, $a[$#a]);
                }
                my @new_regions = @newlist;
        return (@new_regions);
}

sub
trim_3prime_of_Met {
	my ($seq) = @_;

	my $loc = index ($seq, "M");
	if ($loc == -1) {
#	       print "trying towrite an error\n";
	#	my $error_file = $seq_file . ".error";
        #        open (ERROR, ">$error_file");
		print ERROR "trimmed away all Met's and still didn't match\n";
		print ERROR "$seq\n";
		close (ERROR);
		die "trim Met unsuccesfule";
		exit (-1);
	}
	my $trimmed_seq = substr ($seq, $loc + 1, 1000);
print "trimeed seq $trimmed_seq\n";
	return $trimmed_seq;
}

sub
matches {
	my ($substring, $string) = @_;
	my $loc = index ($string, $substring);
	if ($loc == -1) { #substring noti in string return false
		return 0; 
	} else {
		return 1;
	}
}

sub
dna_sequence_contains_stop
{
	my ($seq) = @_;
	my $aaseq = translate_dna_to_aa ($seq);
	if ($aaseq =~ /\*/) {
		return 1;
	} else {
		return 0;
	}

}

sub
trim_last_codon_if_stop
{
        my ($seq1_ref, $seq2_ref) = @_;
        my $seq1 = $$seq1_ref;
        my $seq2 = $$seq2_ref;

        $_ = $seq1;
        my $length = tr/A-Za-z//;
        my $seq1lastcodon = substr ($seq1, $length - 3, 3);
        my $seq2lastcodon = substr ($seq2, $length - 3, 3);
        $seq1lastcodon = uc ($seq1lastcodon);
        $seq2lastcodon = uc ($seq2lastcodon);
        print "last codon is $seq1lastcodon $seq2lastcodon\n";
        if ($seq1lastcodon eq "TGA" || $seq2lastcodon eq "TGA" ||
            $seq1lastcodon eq "TAG" || $seq2lastcodon eq "TAG" ||
            $seq1lastcodon eq "TAA" || $seq2lastcodon eq "TAA") {
                $$seq1_ref = substr ($seq1, 0, $length - 3);
                $$seq2_ref = substr ($seq2, 0, $length - 3);
        }

}



sub
translate_dna_to_aa  {

	my ($dna_sequence) = @_;

	$dna_sequence = uc($dna_sequence);
my %dna_to_aa = ( "TTT", "F" ,
        "TTC",  "F",
	"TTA", "L",
	"TTG", "L",

	"TCT", "S", 
	"TCC", "S",
	"TCA", "S",
	"TCG", "S", 
	
	"TAT", "Y", 
	"TAC", "Y",
	"TAA", "*",
	"TAG", "*", 

	"TGT", "C", 
	"TGC", "C",
	"TGA", "*",
	"TGG", "W",

	"CTT", "L",
	"CTC", "L", 
	"CTA", "L",
	"CTG", "L",

	"CCT", "P",
	"CCC", "P",
	"CCA", "P",
	"CCG", "P",

	"CAT", "H", 
	"CAC", "H",
	"CAA", "Q",
	"CAG", "Q",

	"CGT", "R",
	"CGC", "R",
	"CGA", "R",
	"CGG", "R",

	"ATT", "I",
	"ATC", "I",
	"ATA", "I",
	"ATG", "M",

	"ACT", "T",
	"ACC", "T",
	"ACA", "T",
	"ACG", "T",

	"AAT", "N",
	"AAC", "N",
	"AAA", "K",
	"AAG", "K",

	"AGT", "S",
	"AGC", "S",
	"AGA", "R",
	"AGG", "R",

	"GTT", "V",
	"GTC", "V",
	"GTA", "V",
	"GTG", "V", 
		
	"GCT", "A",
	"GCC", "A",
	"GCA", "A",
	"GCG", "A",

	"GAT", "D",
	"GAC", "D",
	"GAA", "E",
	"GAG", "E",

	"GGT", "G",
	"GGC", "G",
	"GGA", "G",
	"GGG", "G"
 );
 
$_ = $dna_sequence;
my $dna_length = tr/A-Z//;
my $dna_translated = 0;
my $prot_seq = "";
while ($dna_translated < $dna_length ) {
	my $codon = substr ($dna_sequence, $dna_translated, 3);
#	$_ = $codon; my $codon_length = tr/A-Z//;
#	print "codon length $codon_length\n";
#	if ( ($codon_length % 3) == 0) {
		my $aa;
		if (exists ($dna_to_aa{$codon})) {
			$aa = $dna_to_aa{$codon};
		} else {
			$aa = "X";
		}
		#print "$codon: $aa\n";
		$prot_seq .= $aa;
#	}
	$dna_translated += 3;

}
#print "\n\n";
$_ = $prot_seq;
my $prot_length = tr/[A-Z*\-]//;
if (substr ($prot_seq, $prot_length -1, 1) eq "*") {
	$prot_seq = substr ($prot_seq, 0, $prot_length -1);
}
return $prot_seq;

}

sub
loc_of_seq_diff {
	my ($allele1, $allele2) = @_;
	$_ = $allele1;
	my $seq_length = tr/A-Z//;
	my $i = 0;
	while ($i < $seq_length) {
		if (substr ($allele1, $i, 1) ne substr ($allele2, $i, 1)) {
			return $i;
		}
		$i++;
	}
}


sub get_id
{
       my ($filename) = @_;

        open (INB, $filename) || die "can't open $filename";

        my $seq ="";

        my $line = <INB>;
	chomp ($line);
	close (INB);
	return (substr ($line, 1, 1000));
}

sub
read_normal_sequence_with_name
{

	my ($filename) = @_;

        open (INA, $filename) || die "can't open $filename";

        my $seq ="";

        my $firstline = <INA>;
	chomp ($firstline); 
	$firstline = substr ($firstline, 1, 10000);
        my $line;
	while ($line = <INA>) {
                chomp $line;
                $line =~ s/\s+//g; # remove spaces
                $seq .= $line;
        }
        close (INA);
    return ($firstline, $seq);
}


sub read_normal_sequence
{
        my ($filename) = @_;

        open (INA, $filename) || die "can't open $filename";

        my $seq ="";

        my $line = <INA>;
        while ($line = <INA>) {
                chomp $line;
                s/\s+//g; # remove spaces
		$seq .= $line;
        }
	close (INA);
    return ($seq);
}


sub reading_frame_of_subst 
# where is the amino acid change??
{

	my ($codon1, $codon2) = @_;
	print "$codon1 and $codon2\n";
	if (substr ($codon1, 0, 1) ne substr ($codon2, 0, 1)) {
	        $reading_frame = 0; # mutation is at 1rst position of codon
	} elsif (substr ($codon1, 1, 1) ne substr ($codon2, 1, 1)) {
       		 $reading_frame = -1;
	} elsif (substr ($codon1, 2, 1) ne substr ($codon2, 2, 1)) {
       		 $reading_frame = -2;
	}
	return $reading_frame;
} 

sub
retrieve_seq_from_chr_file
{
	my ($gz_chr_file, $start_pos, $end_pos) = @_;
	my $dna_seq;
	my $totalbases = $end_pos - $start_pos + 1;

	my $seqtmp  = `zcat $gz_chr_file | head -2  | tail -1`;
	$_ = $seqtmp;
	my $linewidth = tr/[a-zA-Z]//;

	if ($start_pos < 0) {$start_pos = 0;}

	my $linestart = floor ($start_pos / $linewidth) + 1;  # add 1 to take
                                # into account  the header > line;
	if ($start_pos % $linewidth == 0) {
       	 $linestart = $linestart - 1;
	}

	my $numlines =  ceil ($totalbases / $linewidth) + 1;
	$linestart += $numlines;

	my @seq = `zcat $gz_chr_file | head -$linestart | tail -$numlines`;
	my $basestart = $start_pos % $linewidth;
	if ($basestart == 0) { $basestart = $linewidth;}

	my $basecount = 0; my $line = 0;
	while ($basestart <= $linewidth && $basecount < $totalbases) {
       	 $dna_seq .=  substr ($seq[$line], $basestart - 1, 1);
       	 $basestart++;
       	 $basecount++;
	}
	$line++;
	while ( $basecount < $totalbases){
       	 if ($basecount + $linewidth <= $totalbases) {
                chomp ($seq[$line]);
		$dna_seq .= $seq[$line];
                $basecount += $linewidth;
                $line++;
        } else {
                my $index = 0;
                while ($basecount < $totalbases) {
                        $dna_seq .=  substr($seq[$line], $index, 1);
                        $index++;
                        $basecount++;
                }
        }
	}

return $dna_seq;
}

sub print_sequence {

my ($seq_file, $id, $prot_seq) = @_;

        open (OUT, ">$seq_file") || die "can't open $seq_file";
        print OUT ">$id\n";

        $_ = $prot_seq;
        my $seq_length = tr/[A-Za-z]//;
        print STDERR "seqlength $seq_length\n";
        my $loc = 0;
        while ($loc < $seq_length) {
                my $line =  substr ($prot_seq, $loc, 60);
                print OUT "$line\n";
                $loc += 60;
        }
        close (OUT);
} # end of print_sequence
	
sub format_seq_for_printing{
	my ($id, $seq) = @_;
	$newseq = ">$id\n";
        $_ = $seq;
        my $seq_length = tr/[A-Za-z0-9\-]//;
    #    print "seqlength $seq_length\n";
        my $loc = 0;
        while ($loc < $seq_length) {
                $newseq .=  substr ($seq, $loc, 60) . "\n";
                $loc += 60;
        }
	return ($newseq);
} # end format_seq_for_printing

sub reverse {
	my ($dna_seq) = @_;
        $_ = $dna_seq;
        my $new_dna_seq;
	my $length = tr/[A-Za-z\-]//;
        for (my $i = $length - 1; $i >= 0 ; $i--) {
                my $old = substr ($dna_seq, $i, 1);
                $new_dna_seq .= $old;
        }
        return $new_dna_seq;

}

	
sub complement {
	my ($dna_seq ) = @_;

	my $new_dna_seq = "";
	my %comp = ( "a" => "t",
		     "A" => "T",
		     "t" => "a", 
		     "T" => "A",
		     "G" => "C",
		     "g" => "c",
		     "C" => "G",
		     "c" => "g",
		     "-" => "-",
		     "N" => "N",
		     "/" => "/",
	       	     "n" => "n");

	$_ = $dna_seq;
	my $length = tr/[A-Za-z\-\/]//;
	for (my $i = 0; $i < $length; $i++) {
		my $old = substr ($dna_seq, $i, 1);
		my $new;
		if (exists ($comp{$old})) {
			$new = $comp{$old};
		} else {
			print "\n$old not changed\n"; 
			$new = $old; };
		$new_dna_seq .= $new;
	}
	return $new_dna_seq;
}	

sub get_complete_sequence
{
my ($seqDBFileName, $seqName) = @_;

my $sequence;
open(SEQDB, "$seqDBFileName") ||
  die "Error: Could not open file $seqDBFileName\n";
my $line = "empty";
my @seq;

# Extract sequences
  while (defined($line)) {
    if ($line =~ m/^>(\S+)(\|\S+)?\s/) {
      if ($1 =~ m/$seqName/i) {
        push(@seq, $line);
        while (defined($line = <SEQDB>) && $line !~ m/^>/) {
          push(@seq, $line);
        }

        for (my $i = 0; $i < @seq; $i++) {
                $sequence .= $seq[$i];
        }

        last;
      }
    }

    unless (defined($line = <SEQDB>)) {
      print "Error: Sequence $seqName not found\n";
      last;
    }

}

close(SEQDB);
return $sequence;

} # end get_complete_sequence

sub 
get_sequence_with_name
{
# copied from Nick's get_sequence.pl, given query name and database
# return sequence

my ($seqDBFileName, $seqName) = @_;

my $sequence;
open(SEQDB, "$seqDBFileName") ||
  die "Error: Could not open file $seqDBFileName\n";
my $line = "empty";
my @seq;

# Extract sequences
  while (defined($line)) {
    if ($line =~ m/^>(\S+)(\|\S+)?\s/) {
      if ($1 =~ m/$seqName/i) {
        push(@seq, $line);
        while (defined($line = <SEQDB>) && $line !~ m/^>/) {
          push(@seq, $line);
        }

       $sequence = ">$seqName\n";
        for (my $i = 1; $i < @seq; $i++) {
                chomp ($seq[$i]);
                $sequence .= $seq[$i];
        }
#	$seq[$i] . "\n";	
        last;
      }
    }

    unless (defined($line = <SEQDB>)) {
      print "Error: Sequence $seqName not found\n";
      last;
    }

}

close(SEQDB);
return $sequence;
}

sub
get_sequence
{
# copied from Nick's get_sequence.pl, given query name and database
# return sequence

my ($seqDBFileName, $seqName) = @_;

my $sequence;
open(SEQDB, "$seqDBFileName") ||
  die "Error: Could not open file $seqDBFileName\n";
my $line = "empty";
my @seq;

# Extract sequences
  while (defined($line)) {
    if ($line =~ m/^>(\S+)(\|\S+)?\s/) {
      if ($1 =~ m/$seqName/i) {
        push(@seq, $line);
        while (defined($line = <SEQDB>) && $line !~ m/^>/) {
          push(@seq, $line);
        }

#	$sequence = ">$seqName\n";
	for (my $i = 1; $i < @seq; $i++) {
		chomp ($seq[$i]);
		$sequence .= $seq[$i];
	}

        last;
      }
    }

    unless (defined($line = <SEQDB>)) {
      print "Error: Sequence $seqName not found\n";
      last;
    }

}

close(SEQDB);
return $sequence;
}

sub
percent_id_without_Ns
{
        my ($alignment_seq1, $alignment_seq2, $beg, $end) = @_;

        my $identical = 0; my $total_length = 0;
        for (my $i = $beg; $i < $end; $i++) {
                my $letter1 = substr ($alignment_seq1, $i, 1);
                my $letter2 = substr ($alignment_seq2, $i, 1);
                if ($letter1 eq $letter2 && $letter1 ne "N" && $letter2 ne "N")
                {
                        $identical++;
                }
		if ($letter1 ne "N" && $letter2 ne "N") {
			$total_length++;
		}
	}
	return ($identical, $total_length);
}

sub
percent_id_with_gaps
{
        my ($alignment_seq1, $alignment_seq2, $beg, $end) = @_;

        my $identical = 0; my $total_length = 0;
        for (my $i = $beg; $i < $end; $i++) {
                my $letter1 = substr ($alignment_seq1, $i, 1);
                my $letter2 = substr ($alignment_seq2, $i, 1);
                if ($letter1 eq $letter2)
                {
                        $identical++;
                }
                $total_length++;
        }

        return ($identical, $total_length);

} # end percent_id_of_region
sub
percent_id_of_region_ungapped
{
        my ($alignment_seq1, $alignment_seq2, $beg, $end) = @_;

        my $identical = 0; my $total_length = 0;
        for (my $i = $beg; $i < $end; $i++) {
                my $letter1 = substr ($alignment_seq1, $i, 1);
                my $letter2 = substr ($alignment_seq2, $i, 1);
                if ($letter1 eq $letter2 && $letter1 ne "N" && $letter1 ne "-"
			&& $letter2 ne "-" )
                {
                        $identical++;
                }
                unless ($letter1 eq "-" || $letter2 eq "-") {
                        $total_length++;
                }
	}			
	return ($identical, $total_length);
}

sub
percent_id_of_region
{
	my ($alignment_seq1, $alignment_seq2, $beg, $end) = @_;

	my $identical = 0; my $total_length = 0;	
	for (my $i = $beg; $i < $end; $i++) {
		my $letter1 = substr ($alignment_seq1, $i, 1);
		my $letter2 = substr ($alignment_seq2, $i, 1);
		if ($letter1 eq $letter2 && $letter1 ne "N" && $letter1 ne "-")
		{
			$identical++;
		}
		unless ($letter1 eq $letter2 && $letter1 eq "-") {
			$total_length++;
		}	
	}

	return ($identical, $total_length);	

} # end percent_id_of_region

sub
region_length
{
        my ($reg) = @_;
        my $end = $reg->end;
        my $beg = $reg->beg;
        print "region length is \n";
        print "here" .  ($end - $beg + 1) . "\n";
        return ($end - $beg + 1);
}

sub
within
{
        my ($loc, $start, $end) = @_;
        if ($loc >= $start && $loc <= $end) {
                return 1;
        }
        return 0;
}

sub
read_a_sequence2
{
       my ($seekpos) = @_;

        seek DB2, ${$seekpos}, 0;
        if (!eof (DB2)) {
                $_ = <DB2>;
                $$seekpos = tell DB2;
        }
        chomp;
        my $firstline = substr ($_, 1, 10000);
        my $seq = "";
        do {
                if (!eof (DB2)) {
                        $_ = <DB2>;
                        chomp;
                        unless ($_ =~ /^>/) {
                                $seq .= $_;
                                $$seekpos = tell DB2;
                        }
                }
        } while (!eof (DB2) && !($_ =~ /^>/)) ;
#       if ($seq eq "") {
                # read next line so can proceed
#                $$seekpos = tell DB2;
#       }
        $seq =~ s/\s+//g;
        return ($firstline, $seq);
}

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
#	if ($seq eq "") {
		# read next line so can proceed 
#		 $$seekpos = tell DB;
#	}	
	$seq =~ s/\s+//g;
        return ($firstline, $seq);
}

sub
read_a_sequence_keep_space
{
        my ($seekpos) = @_;

        seek DB, ${$seekpos}, 0;
        if (!eof (DB)) {
                $_ = <DB>;
                $$seekpos = tell DB;
        }
        my $firstline = substr ($_, 1, 10000);
        my $seq = "";
        do {
                if (!eof (DB)) {
                        $_ = <DB>;
                    #    chomp;
                        unless ($_ =~ /^>/) {
                                $seq .= $_;
                                $$seekpos = tell DB;
                        }
                }
        } while (!eof (DB) && !($_ =~ /^>/)) ;
#       if ($seq eq "") {
                # read next line so can proceed
#                $$seekpos = tell DB;
#       }
        $seq =~ s/\n/ /g;
        return ($firstline, $seq)


}

sub 
intersect
{
	my ($beg1, $end1, $beg2, $end2) = @_;
	if (overlapping ($beg1, $end1, $beg2, $end2)) {
		if ($beg1 < $beg2) { $start = $beg2; } else { $start = $beg1;}
		if ($end1 < $end2) { $end = $end1; } else { $end = $end2;}
		return ($start, $end);
	}
	return (-1,-1);	

}

sub
overlapping
{
        my ($beg1, $end1, $beg2, $end2) = @_;

        if (within ($beg1, $beg2, $end2)) {
                return 1;
        } elsif (within ($end1, $beg2, $end2)) {
                return 1;
        } elsif (within ($beg2, $beg1, $end1)) {
                return 1;
        } elsif (within ($end2, $beg1, $end1)) {
                return 1;
        } else {
                return 0;
        }
}

# returns hash for a file, 2nd field is the key and the 3rd field
# is the value 4th field, is the delimiter
sub
make_hash
{
	my ($file, $keyindex, $valueindex, $delimiter) = @_;
	my %hash;	
	open (HASH, $file) || die "can't open $file";
	my $line; 
	while ($line = <HASH>) {
		chomp ($line);;
		my @fields ;
		if ($delimiter eq "") {
			@fields = split (/\s+/, $line);
		} else {
			@fields = split (/$delimiter/, $line);
		}
#print "fields key $fields[$keyindex] value $fields[$valueindex]\n";
		if ($valueindex eq "") {
			 $hash{$fields[$keyindex]} = $line;
		} else {
			$hash{$fields[$keyindex]} = $fields[$valueindex];
		}
	}
	close (HASH);
	return (%hash);	
}

sub
by_number {
	if ($a < $b) {
		return -1;
	} elsif ($a == $b) {
		return 0;
	} elsif ($a > $b) {
		return 1;
	}
}

sub
find_sequence 
{
        my ($title, $db) = @_;
        open (SEQ, $db) || die "can't open $db";
        my $done = 0;

        while (!$done) {
                if (!eof (SEQ)) {
                        $_ =<SEQ>;
                } else {
                        close (SEQ);
                        return "";
                }
                my @fields= split;
                if ($fields[0] eq ">" . $title ) {
                        my $seq = "";
                        while (!eof (SEQ)) {
                                $_ = <SEQ>;
                                chomp;
                                if (/^>/) {
                                        close (SEQ);
                                        return ($seq);
                                } else {
                                        $seq .= $_;
                                }
                        }
                        if (eof (SEQ)) {
                                close (SEQ);
                                return "";
                        }
                }
        }
} # end find_sequence 

1;
