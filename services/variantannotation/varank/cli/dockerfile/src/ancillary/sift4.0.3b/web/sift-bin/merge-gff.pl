#!/usr/local/bin/perl
#
# File supplied - merge-gff.pl -  a program to merge together overlapping GFF. This
# code assumes that the GFF in question is sorted by UID, then begin position
#
#

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use Getopt::Long;
#use lib $ENV{COMP_PERL_LIB};
#use utils;

sub is_overlap{
  my $g1 = shift;
  my $g2 = shift;
  my $result = 0;

  if($g1->[0] eq $g2->[0] && $g1->[2] eq $g2->[2]){
    if($g1->[3] >= $g2->[3] && $g1->[3] <= $g2->[4]){
      $result=1;
    }
  }
  return $result
}

my $gff1_file,$gff2_file,$ft1,$ft2;
my $g1_last;
undef $gff1_file,$gff2_file,$ft1,$ft2;
my @g1;
my $ol_count;

GetOptions("G1=s" => \$gff1_file,
	   "G2=s" => \$gff2_file,
	   "F1=s" => \$ft1,
	   "F2=s" => \$ft2
	   );

$usage="
       $0 <options>

         options:

             -G1 GFF file to merge

             -G2 GFF file that will contain merged output

             -F1 Feature type in file to be merged 

             -F2 Feature type of output. If not present the value in F1 will be assumed 

Please ensure that the file is sorted in order in the 1st and 3rd columns.
";

if(!defined($gff1_file) || !defined($gff2_file) || !defined($ft1)){
  print "$usage\n";
  exit(1);
}

if(!defined($ft2)){
  $ft2 = $ft1;
}


open(GFF1, $gff1_file);
open(GFF2, ">$gff2_file");

$g1_last = <GFF1>;
chomp $g1_last;
push @g1,$g1_last;

$ol_count=1;

while(<GFF1>){
  chomp $_;
  $g1_curr = $_;
  # get current record
  @gc= split(/\t/,$g1_curr);
  # get last record
  @gl= split(/\t/,pop(@g1));
  # now compare
  if(is_overlap(\@gc,\@gl) == 1){
    # then overlaps, update the gff end position and
    # store the update value in array.
    if($gc[4] > $gl[4]){
      $gl[4] = $gc[4];
    }
    push @g1,(join "\t",@gl);
    $ol_count++;
  }
  else{
    # doesnt overlap therefore output the last record and
    # store the current value for next comparison
    printf GFF2 join "\t",@gl,$ol_count,"\n";
    push @g1,(join "\t",@gc);
    $ol_count=1;
  }
}
# ensure to pop any record still in @gc
#my $remain = pop(@gc);
if($#gc > 0){
  printf GFF2 join "\t",@gc,$ol_count,"\n";
}

close(GFF1);
close(GFF2);
