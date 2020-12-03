#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
$denormal = @ARGV[0];
if ($denormal =~ /(\d+?)\_snps_classified_chr(.+)_.+\.denormal/){
	$pid = $1;
	$chr = $2;
}
else {
	print "Illegal file name";
	exit;
}

open (DENORMAL, "$denormal") || die ("cannot open denormal file");

while (<DENORMAL>){
	chomp;
	$coord1 ="" ;
	$coord2="" ;
	$orn="" ;
	$ntpos1="" ;
	$ntpos2="" ;
	$aapos1="" ;
	$aapos2="" ;
	$snp="" ;
	$codon1="" ;
	$codon2="" ;
	$aa1="" ;
	$aa2="" ;
	$nt1 = "";
	$nt2 = "";
	$conservation="" ;
	$CDS = 0;
	$seq = "";
	
	
	@elts = split /\t/, $_;
	$id = $elts[0];
	$coords = $elts[1];
	$seq  = $elts[2];
	$rsid = $elts[3];		
	$ensg = $elts[5];
	$enst = $elts[6];
	$region = $elts[7];
	$positions = $elts[8];
	$subst = $elts[9];
	$aa_pos = $elts[10];
	
	#if ($_ =~ /INTRON/){
	#	next;
	#}
	
	if ($seq =~ /.+? (.+?) .+?/){
		$indel = $1;
		
	}	
	
	if ($_ =~ /CDS/){
		$CDS = 1;
	}
	
	if ($coords =~ /(\d+?)\-(\d+?)\:(-?1)\:.*?/){
		$coord1 = $1;
		$coord2 = $2;
		$coord1 = $bin_start+$coord1 -1;
		$coord2 = $bin_start+$coord2 -1;
		$orn = $3;
	}
	if ($positions =~ /\[(.+?)-(.+?) (.+?)-(.+?)]/){
		$ntpos1 = $1;
		$ntpos2 = $2;
		$aapos1 = $3;
		$aapos2 = $4;
	}
	#if ($subst =~ /AA_DELETION/
	if ($subst =~ /subst_synonymous\[(.+?)\:(.+?) (.+?)\:(.+?) (.+?)\]/i){
		$snp = "Synonymous";
		$codon1 = $1;
		$aa1 = $2;
		if ($aa1 eq " "){
			$aa1 ="*"
		}
		$codon2 = $3;
		$aa2 = $4;
		$conservation = $5;
		for ($i = 0; $i < 3; $i ++){
			if (substr($codon1, $i,1) ne substr($codon2, $i,1)){				
				$nt2 = uc(substr($codon2, $i,1));
			}
		}
		if ($orn == -1){
			if ($nt2 =~ /A/i){
				$nt2 = "T";
			}
			elsif ($nt2 =~ /T/i){
				$nt2 = "A";
			}
			elsif ($nt2 =~ /G/i){
				$nt2 = "C";
			}
			elsif ($nt2 =~ /C/i){
				$nt2 = "G";
			}
		}
		
	}
	if ($subst =~ /subst_NONSYNONYMOUS\[(.+?)\:(.+?) (.+?)\:(.+?) (.+?)\]/i){
		$snp = "Nonsynonymous";
		$codon1 = $1;
		$aa1 = $2;
		if ($aa1 eq " "){
			$aa1 ="*"
		}
		$codon2 = $3;
		$aa2 = $4;
		$conservation = $5;		
		for ($i = 0; $i < 3; $i ++){
			if (substr($codon1, $i,1) ne substr($codon2, $i,1)){				
				$nt2 = uc(substr($codon2, $i,1));
			}
		}
		if ($orn == -1){
			if ($nt2 =~ /A/i){
				$nt2 = "T";
			}
			elsif ($nt2 =~ /T/i){
				$nt2 = "A";
			}
			elsif ($nt2 =~ /G/i){
				$nt2 = "C";
			}
			elsif ($nt2 =~ /C/i){
				$nt2 = "G";
			}
		}
	}
	
	
	print "$chr\t$id\t$coord1\t$coord2\t$orn\t$rsid\t$ensg\t$enst\t$subst\t$region\t$snp\t$ntpos1\t$ntpos2\t$aapos1\t$aapos2\t$CDS\n";
	#print "$coords\n";
	
}
