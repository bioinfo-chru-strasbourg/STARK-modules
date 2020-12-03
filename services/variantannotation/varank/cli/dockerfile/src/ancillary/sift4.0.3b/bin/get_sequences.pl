#! /opt/local/bin/perl5 -w

#####################################################################
#
# Filename: get_sequences.pl - get specific sequences from a database
#
# Input:
#	File containing names of sequences
#	Name of a sequence
# Action:
#	Retrieves the sequences from the sequence database
# Output:
#	File called $seqName.seqs
#
# By: Bob Chan, bobchan@howard.fhcrc.org
# Last modified: 12/13/00
#
# COPYRIGHT 2000. This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software
#
#
#####################################################################

## Constants --------------------------------------------------------

my %seq_names_hash;

## Main

# Handle command-line arguments
unless (scalar(@ARGV) == 2) {
  die "Usage: $0 seq_names_file seq_db_file\n";
}

my $seq_names_filename = shift;
my $seq_db_filename = shift;

open(SEQS, "$seq_names_filename") ||
  die "Error: Could not open file $seq_names_filename for reading!\n";

while (<SEQS>) {
  chomp;
  my $line = $_;
  if ($line =~ m/(\S+)/) {
    $seq_names_hash{$1} = 1;
  }
}

close(SEQS);

# Read sequence database file
my $line = "empty";
my @seqs;

open(SEQDB, "$seq_db_filename") ||
  die "Error: Could not open file $seq_db_filename\n";

# Extract the sequences
while (defined($line)) {
  if ($line =~ m/^>\s*([\w_]+)(\|\S+)?\s/) {
    if (exists($seq_names_hash{$1})) {
      push(@seqs, $line);
      while (defined($line = <SEQDB>) && $line !~ m/^>/) {
        push(@seqs, $line);
      }
    } else {
      $line = <SEQDB>;
    }
  } else {
    $line = <SEQDB>;
  }
}

close(SEQDB);

  my $seqOutFile = $seq_names_filename . ".seqs";
  open(OUTFILE, ">$seqOutFile") ||
    die "Error: Could not open file $seqOutFile for writing\n";
  print OUTFILE @seqs;
  close(OUTFILE);

exit(0);
