#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.


#ARG1 snp file

$snpfile = @ARGV[0];
open (SNPFILE,"$snpfile") || die("Cannot open input SNP file. Please check argument 2");
open (GFFFILE,">$snpfile.journal") || die ("Cannot write to location. Please check permissions");
while (<SNPFILE>){
	chomp;
	@elts = split /\,/, $_;
	$chr = @elts[0];
	$coord1 = @elts[1];
	$coord2 = @elts[2];
	if ($coord1 == $coord2){
		$coord2 = $coord1 + 1;
	}
	$orn = @elts[3];
	if ($orn eq "1"){
		$orn = "+";
	}
	else{$orn = "-"}
	$alleles = @elts[4];
	print GFFFILE "$chr\tSource\tType\t$coord1\t$coord2\t\.\t$orn\t\.\n";

}
close (SNPFILE);
close (GFFFILE);
system ("sort -k1,1 -k4,4n -k5,5n $snpfile.journal > $snpfile.gff");
unlink("$snpfile.journal");
