#!/usr/local/bin/perl
$SIFT_HOME = $ENV{'SIFT_HOME'};
$bin = "$SIFT_HOME/bin";
 
###############################################################################
#                                                                             #
#       Copyright (c) 2009 J. Craig Venter Institute.                         #
#       All rights reserved.                                                  #
#                                                                             #
###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.    #
#                                                                             #
###############################################################################



require "$bin/DNA_PROT_SUBROUTINES.pl";

######################################################################

use strict;

use Getopt::Std;

use vars qw($opt_s $opt_c $opt_n $opt_o);

getopts("s:c:n:o:");
my $usage = "usage: 
$0 

	-s <Query SNPs Filename>
	-c <coding info file>
	-n <reference nucleotide sequence FASTA>
	-o <output file root>

	All values should be in local 0 space based coordinates.
	
";
my $tmpdir = "/opt/www/sift/tmp/";
my $blossum_hash_ref=loadBlossumMatrix();

if(!(
	defined($opt_s) && 
	defined($opt_c) && 
	defined($opt_o) && 
	defined($opt_n))){
	print STDERR $usage;
	die;
}

###############################################################################

my $UPSTREAM_PROMOTER_REGION=1000;
my $DOWNSTREAM_GENOMIC_REGION=1000;
my $UPSTREAM_SS_REGION=10;
my $DOWNSTREAM_SS_REGION=6;
my $FLANKING_LENGTH=100;

###############################################################################
my $output_file_root=$opt_o;
my $sift_pid = (split /_/,$output_file_root)[0];
my $normalize_out="$output_file_root\.normal";
my $denormalize_out="$output_file_root\.denormal";

###############################################################################

open(QUERY_SNP_FH, "<$opt_s") || die "Could not open Query SNPs File: $opt_s\n";
open(CODING_INFO_FH, "<$opt_c") || die "Could not open Coding Info File: $opt_c\n";
open(REF_NUC_FH, "<$opt_n") || die "Could not open Reference Nucleotide Sequence FASTA File: $opt_n\n";

###############################################################################

open(NORM_OUT_FH, ">$normalize_out") || die "Could not open $normalize_out\n";
open(DENORM_OUT_FH, ">$denormalize_out") || die "Could not open $denormalize_out\n";

###############################################################################

my @cds_arr;
my @genes_arr;
my @snps_arr;
my %gene_begin_hash;
my %gene_end_hash;
my %gene_ori_hash;
my %gene_transcript_hash;
my %transcript_exon_hash;
my %transcript_exon_id_hash;
my %transcript_ori_hash;
my %transcript_cds_hash;
my %protein_domain_hash; # keyed by transcript_id
my %gene_name_hash;
my %repeat_hash;

###############################################################################
# Read in coding information

print STDERR "Reading coding information...\n";
while(<CODING_INFO_FH>){

	my @data_line=split "\t";

	if($data_line[0] eq "EXON"){
		my $gene_id=$data_line[1];
		my $transcript_id=$data_line[2];
		my $exon_id=$data_line[3];
		my $begin=$data_line[4];
		my $end=$data_line[5];
		my $orientation=$data_line[6];

		#######################################################################

		# Track gene's begin
		if(!defined($gene_begin_hash{$gene_id})){
			$gene_begin_hash{$gene_id}=$begin;
		}else{
			if($gene_begin_hash{$gene_id}>$begin){
				$gene_begin_hash{$gene_id}=$begin;
			}
		}

		# Track gene's end
		if(!defined($gene_end_hash{$gene_id})){
			$gene_end_hash{$gene_id}=$end;
		}else{
			if($gene_end_hash{$gene_id}<$end){
				$gene_end_hash{$gene_id}=$end;
			}
		}

		#######################################################################

		${$gene_transcript_hash{$gene_id}}{$transcript_id}=1;

		# These three arrays are parallel.
		push @{$transcript_exon_hash{$transcript_id}}, "$begin#$end";
		push @{$transcript_exon_id_hash{$transcript_id}}, $exon_id;
		push @{$transcript_ori_hash{$transcript_id}}, $orientation;

		#######################################################################

	}elsif($data_line[0] eq "CDS"){
		my $gene_id=$data_line[1];
		my $transcript_id=$data_line[2];
		my $begin=$data_line[3];
		my $end=$data_line[4];

		#######################################################################

		$transcript_cds_hash{$transcript_id}="$begin#$end";
		
	}elsif($data_line[0] eq "SNP"){
		my %snp_struct;
		$snp_struct{dbsnp_id}=$data_line[1];
		$snp_struct{begin}=$data_line[2];
		$snp_struct{end}=$data_line[3];

		if($data_line[4]==-1){
			$snp_struct{allele}=revcomp($data_line[5]);
		}else{
			$snp_struct{allele}=$data_line[5];
		}

		push @snps_arr, \%snp_struct;

	}elsif($data_line[0] eq "PROTEIN"){
		my %prot_struct;
### Pauline -- designate which protein domains to store
		unless ($data_line[3] =~ /Low complexity/ || 
			$data_line[3] =~ /Low_complexit/) {
			$prot_struct{transcript_id}=$data_line[1];
			$prot_struct{protein_id}=$data_line[2];
			$prot_struct{description}=$data_line[3];
			$prot_struct{begin}=$data_line[4];
			$prot_struct{end}=$data_line[5];

			push @{$protein_domain_hash{$data_line[1]}}, "$data_line[4]\t$data_line[5]\t$data_line[3]";
		}
	}elsif($data_line[0] eq "HUGO_ID"){
		$gene_name_hash{$data_line[1]}=$data_line[2];
	}elsif($data_line[0] eq "REPEAT"){
		$repeat_hash{join "#", ($data_line[1], $data_line[2], $data_line[3])}=1;
	}

}


# Compute the begin/end of the transcripts
my %transcript_extent_hash;
foreach my $transcript_id(keys %transcript_exon_hash){
	my ($min_begin, $max_end)=split /#/, ${$transcript_exon_hash{$transcript_id}}[0];
	foreach my $exon(@{$transcript_exon_hash{$transcript_id}}){
		my ($begin, $end)=split /#/, $exon;
		if($begin<$min_begin){
			$min_begin=$begin;
		}
		if($end>$max_end){
			$max_end=$end;
		}	
	}
	$transcript_extent_hash{$transcript_id}="$min_begin#$max_end";
}

###############################################################################
# Output status

print STDERR "Genes Read:\n";
foreach my $gene_id(keys %gene_end_hash){
	print STDERR "\t$gene_id: $gene_begin_hash{$gene_id}-$gene_end_hash{$gene_id}\n";
}

print STDERR "Transcripts Read:\n";
foreach my $transcript_id(keys %transcript_extent_hash){
	my $transcript_extent=$transcript_extent_hash{$transcript_id};
	$transcript_extent=~s/#/-/;
	print STDERR "\t$transcript_id: $transcript_extent\n";
}

###############################################################################

print STDERR "Reading FASTA File...\n";
my $defline=0;
my $sequence="";
while(<REF_NUC_FH>){
	chomp;
	if(/^>(.+)/){
		print STDERR "\tDefline: $1\n";
		$defline++;
		if($defline>1){
			die "This FASTA file should only have one sequence record in it.\n";
		}
	}else{
		$sequence.=uc($_);	
	}
}

print STDERR "\tLength of reference sequence read: " . length($sequence) . "\n";
###############################################################################

while(<QUERY_SNP_FH>){

	chomp;

	# Parse query snp line
	my ($snp_id, $begin, $end, $orientation, $qry_allele)=split /\t/;

	if(!(defined($begin) && defined($end) && defined($orientation) && defined($qry_allele))){
		die "Your query SNP file is incomplete.\n";
	}

	# Extract reference allele on the reference
	my $reference_allele=substr($sequence, $begin, $end-$begin);
	my $upstream_flank=substr($sequence, $begin-$FLANKING_LENGTH, $FLANKING_LENGTH);
	my $downstream_flank=substr($sequence, $end, $FLANKING_LENGTH);

	$reference_allele=($reference_allele eq "")?"-":$reference_allele;

	# Check to see if snp overlaps with known dbSNP
	my $dbSNP_matches_arr_ref=CheckSNPs($begin, $end, $reference_allele);
	my $dbSNP_matches_str=join ";", @{$dbSNP_matches_arr_ref};
	$dbSNP_matches_str=($dbSNP_matches_str eq "")?"novel":$dbSNP_matches_str;

	# If orientation is reversed, convert it to forward strand of reference
	my @nonref_forw_qry_alleles;
	if($orientation==-1){
		my @qry_alleles_arr=split /\//, $qry_allele;
		for(my $i=0; $i<=$#qry_alleles_arr; $i++){
			$qry_alleles_arr[$i]=revcomp($qry_alleles_arr[$i]);
		}
		$qry_allele=join  "\/", @qry_alleles_arr;
	}

	# See if SNP is in repeat region
	my $repeat_list_ref=CheckRepeats($begin, $end);
	my $repeat_list_str=join ";", @{$repeat_list_ref};

	# This is a per SNP analysis
	my $snp_loc_info="$snp_id\t$begin-$end:$orientation:$qry_allele\t$upstream_flank $reference_allele $downstream_flank\t$dbSNP_matches_str\t$repeat_list_str";

	# Get the list of genes we have read in
	my @gene_array=sort keys %gene_transcript_hash;

	# Output the snp location
	print NORM_OUT_FH "$snp_loc_info\n";

	my $associated_with_gene=0;
	foreach my $gene_id(@gene_array){
		my $gene_printed=0;

		foreach my $transcript_id(sort keys %{$gene_transcript_hash{$gene_id}}){
			my ($hit_type, $hit_info)=CheckTranscript($begin, $end, $transcript_id);
#		print "$snp_id is in $transcript_id $hit_type $hit_info\n";
			my $transcript_info="$transcript_id\t$hit_info";
			my $gene_info="$gene_name_hash{$gene_id} $gene_id";

			if($hit_type ne "IRRELEVANT"){
print STDERR "entered in here $snp_id is in $transcript_id $hit_type $hit_info\n";
				$associated_with_gene=1;
				if(!($gene_printed)){
					print NORM_OUT_FH "\t$gene_info\n";
					$gene_printed=1;
				}

				if($hit_type eq "CDS"){

					my ($snp_cds_begin, $snp_cds_end, $snp_ori, $ref_cds_seq)=
						getCDSSequence($begin, $end, $transcript_id);
					my $mut_arr_ref=
						GetDonorMutation($snp_cds_begin, $snp_cds_end, $snp_ori, $qry_allele, $ref_cds_seq, $transcript_id, $snp_id);

			
			 my ($distance_to_upstream_intron,
                                $adjacent_upstream_intron_size,
                         $distance_to_downstream_intron,
                        $adjacent_downstream_intron_size) =
                                                CheckFlankingIntrons ($begin, $end, $transcript_id);
                        print NORM_OUT_FH "$snp_id adjacent introns $distance_to_upstream_intron  $adjacent_upstream_intron_size  $distance_to_downstream_intron  $adjacent_downstream_intron_size\n";

	
					my $prot_begin=int($snp_cds_begin/3);
                                        # tbs - 2007/01/10 - made fix by subtracting 1
					my $prot_end=int(($snp_cds_end-1)/3)+1;

					#my $domain_hits_arr_ref=CheckProtein($prot_begin, $prot_end, $transcript_id);
					#my $domain_hits_str=join ";", @{$domain_hits_arr_ref};
				
					#my $cds_info="$gene_id $transcript_id $hit_info\t$snp_cds_begin-$snp_cds_end $prot_begin-$prot_end $domain_hits_str";
					my $cds_info="[$snp_cds_begin-$snp_cds_end $prot_begin-$prot_end]";
					print NORM_OUT_FH "\t\t$transcript_info $cds_info\n";

					if($#{$mut_arr_ref}==-1){
						print DENORM_OUT_FH "$snp_loc_info\t$gene_info\t$transcript_info\t$cds_info\n";
					}else{
						foreach my $mut(@{$mut_arr_ref}){
							print NORM_OUT_FH "\t\t\t$mut\n";
							print DENORM_OUT_FH "$snp_loc_info\t$gene_info\t$transcript_info\t$cds_info\t$mut\n";
						}
					}

				}else{
					print NORM_OUT_FH "\t\t$transcript_info\n";
					print DENORM_OUT_FH "$snp_loc_info\t$gene_info\t$transcript_info\n";
				}
			}
			
		}
	}

	# If there were no genes or SNP was not asscoiated with any transcripts
	if(!$associated_with_gene){
		print DENORM_OUT_FH "$snp_loc_info\n";
	}
}

###############################################################################

sub GetDonorMutation{
	my ($begin, $end, $snp_ori, $allele_str, $refer_cds_seq, $transcript_id, $id) = @_;

	my $target_length=$end-$begin;
	
	# Pauline, if it only has one allele, assume the other is a blank
	unless ($allele_str =~ /\//) {
		$allele_str .= "/ ";
	}
	my @allele_arr=split ( '/', $allele_str);
	my $ref_allele=substr($refer_cds_seq, $begin, $end-$begin);
	my @mut_array;
#	 print NORM_OUT_FH "GetDonorMutation $id allele ref: $ref_allele mutated_allele: $allele_str end\n";

	foreach my $allele(@allele_arr){
# print NORM_OUT_FH "Foreach $id allele ref: $ref_allele mutated_allele: $allele end $snp_ori\n";

		if($snp_ori==-1){
			$allele=revcomp($allele);
		}

		 $allele=~ s/\-//g; # moved out here so that insertions not compared 
		$allele =~ s/\s+//g;
	
		if($ref_allele ne $allele){
			# Don't report if the allele is the same as the reference
			print NORM_OUT_FH "$id allele ref: $ref_allele mutated_allele: $allele end\n";
			# Compute donor allele length
			$allele=~s/-//g;
			my $allele_length=length($allele);

			# Compute if there will be a frame shift mutation
			my $frameshift_mut=((abs($allele_length-$target_length)%3)!=0)?1:0;

			# Compute which codons will be affected
			my ($target_codon_begin, $target_codon_end);
			my ($affected_aa_begin, $affected_aa_end);

			$target_codon_begin=$begin-($begin%3);
			$target_codon_end=($end-1)-(($end-1)%3)+3;

			$affected_aa_begin=$target_codon_begin/3;

			# Apply the donor allele to the reference seqeuence
			my $donor_cds_seq=$refer_cds_seq;

			my $allele_length = $allele =~ tr/A-Za-z//;
#			print NORM_OUT_FH "Reference\n$refer_cds_seq\n";
			substr($donor_cds_seq, $begin, $end-$begin)=lc($allele);
			$donor_cds_seq =~ s/\s+//g;
#		print NORM_OUT_FH "Donor\n$donor_cds_seq\n";	
			my ($donor_protein)=split / /, translate($donor_cds_seq);
			# check if there is a downstream Met
### Pauline - should be more intelligent about this
			if ($begin <= 50) {
				$donor_protein = translate_with_downstream_Met
					($donor_cds_seq, $end);	
			} 
			my $donor_prot_len=length($donor_protein);

			my ($refer_protein)=split / /, translate($refer_cds_seq);
			$_ = $refer_protein;
			my $refer_protein_length =tr/[A-Za-z\-\ ]//;
# Pauline's edits
			print_sequence ("$sift_pid\_$transcript_id"  . ".fa", 
				$transcript_id, $refer_protein) ; 
			print_sequence ("$sift_pid\_$transcript_id" . "_" . $id . ".fa",
				$transcript_id . "_" . $id , $donor_protein);
#			print NORM_OUT_FH "pairwise between original and mutated\n";
#Prateek - disabling clustalw for now
			#my ($percent_id, $percent_aligned, $aligned_donor_protein) = 
				#run_clustal ($refer_protein, $donor_protein);


#			my $blast_result = blast_pairwise ($refer_protein, $donor_protein);
#			my ($e_value, $percent_id, $percent_aligned) = 
#				get_blast_scores ($blast_result, $refer_protein_length);
	#		print "$refer_protein\n";
#			my $aligned_donor_protein = 
#				get_blast_alignment ($blast_result, $refer_protein_length);
#Prateek - back to original format
			#print NORM_OUT_FH $id . "\t" . "Identity: $percent_id Aligned: $percent_aligned\%" . "\n";
			#CheckProteinDomains ($transcript_id, $refer_protein, $aligned_donor_protein, $id);
			my ($distance_to_upstream_intron,  
				$adjacent_upstream_intron_size,
	                 $distance_to_downstream_intron,  
			$adjacent_downstream_intron_size) = 
#						CheckFlankingIntrons ($begin, $end, $transcript_id);
#			print NORM_OUT_FH "adjacent introns $distance_to_upstream_intron  $adjacent_upstream_intron_size  $distance_to_downstream_intron  $adjacent_downstream_intron_size\n";
	
			my $refer_prot_len=length($refer_protein);

			my $mutation_type;

			if($frameshift_mut){
				# affected region will be from substitution down since
				# the translation frame will be different.  In this case,
				# we need to compute the rest of the protein
				$affected_aa_end=length($refer_cds_seq)/3;

				$mutation_type="FRAMESHIFT[$refer_prot_len $donor_prot_len]";

			}elsif($allele_length == $target_length){
				# affected region will only be substitution region so
				# we only need to determine if it is substitution or deletion.
				$affected_aa_end=$target_codon_end/3;

				my $donor_codon=substr($donor_cds_seq, $target_codon_begin, $target_codon_end-$target_codon_begin);
				my $refer_codon=substr($refer_cds_seq, $target_codon_begin, $target_codon_end-$target_codon_begin);
				my $donor_aa=translate($donor_codon);
				my $refer_aa=translate($refer_codon);
				$donor_aa=($donor_aa eq " ")?"*":$donor_aa;

				my $conserved_modification;
				if(isAminoAcidChangeConservative($blossum_hash_ref, $refer_aa, $donor_aa)>=0){
					$conserved_modification="c";	
				}else{
					$conserved_modification="U";	
				}

				my $allele_info="$refer_codon:$refer_aa $donor_codon:$donor_aa $conserved_modification";

				# See if substitution codes for different amino acids
				if($donor_aa eq $refer_aa){
					$mutation_type.="subst_synonymous[$allele_info]";
				}else{
					$mutation_type.="subst_NONSYNONYMOUS[$allele_info]";
				}	

			}else{
				# only affected region will be affected and will result in either
				# an extended or shortened protein, but non targeted regions	
				# will be unaffected.
				$affected_aa_end=$target_codon_end/3;

				if($allele_length<$target_length){
					$mutation_type="AA_DELETION[$refer_prot_len $donor_prot_len]";
				}else{
					$mutation_type="AA_INSERTION[$refer_prot_len $donor_prot_len]";	
				}

			}

			my $domain_hits_arr_ref=CheckProtein($affected_aa_begin, $affected_aa_end, $transcript_id, $transcript_id);
			my $dom_hits_str=join ";", @{$domain_hits_arr_ref};

			$mutation_type.="\t$affected_aa_begin-$affected_aa_end $dom_hits_str";

			push @mut_array, $mutation_type;
		}
	}

	return (\@mut_array);
}

###############################################################################

sub getCDSSequence{
	my ($begin, $end, $transcript_id)=@_;
	my ($snp_cds_begin, $snp_cds_end, $snp_ori, $cds_seq)=(undef, undef, undef, "");

	my $prev_exons=0;
	my $stop_exon_begin;
	my $stop_exon_end;

	my $is_insert=($begin==$end)?1:0;
	my $ins_pos=$begin;

	# Get the CDS begin/end in genomic coordinates for the transcript
	# cds_begin < cds_end, even if transcript is on opposite strand
	my ($cds_begin, $cds_end)=split /#/, $transcript_cds_hash{$transcript_id};

	# For each exon in the transcript
	for(my $i=0; $i<=$#{$transcript_exon_hash{$transcript_id}}; $i++){

		my $exon=${$transcript_exon_hash{$transcript_id}}[$i];
		my $exon_ori=${$transcript_ori_hash{$transcript_id}}[$i];

		my ($ex_begin, $ex_end)=split /#/, $exon;

		if(overlaps($ex_begin, $ex_end, $cds_begin, $cds_end)){

			# Limit only returning CDS
			if(isIn($cds_begin, $ex_begin, $ex_end)){
				$ex_begin=$cds_begin;
			}
			if(isIn($cds_end, $ex_begin, $ex_end)){
				$ex_end=$cds_end;
			}

			# Compute exon length
			my $cds_length=$ex_end-$ex_begin;

			# Reverse complement if the exon is on opposite strand
			if($exon_ori==1){
				$cds_seq.=substr($sequence, $ex_begin, $cds_length);
				if(overlaps($begin, $end, $ex_begin, $ex_end) || ($is_insert && isIn($ins_pos, $ex_begin, $ex_end))){
					$snp_cds_begin=$begin-$ex_begin+$prev_exons;
					$snp_cds_end=$end-$ex_begin+$prev_exons;
					$snp_ori=1;
				}
			}elsif($exon_ori==-1){
				$cds_seq.=revcomp(substr($sequence, $ex_begin, $cds_length));
				if(overlaps($begin, $end, $ex_begin, $ex_end) || ($is_insert && isIn($ins_pos, $ex_begin, $ex_end))){
					$snp_cds_end=$ex_end-$begin+$prev_exons;
					$snp_cds_begin=$ex_end-$end+$prev_exons;
					$snp_ori=-1;
				}
			}else{
				die "$exon_ori not 1|-1\n";
			}

			$prev_exons+=$cds_length;
		}
	}

	return ($snp_cds_begin, $snp_cds_end, $snp_ori, $cds_seq);
}

###############################################################################

sub CheckGenes{
	my $begin=shift;
	my $end=shift;
	my @hit_genes;

	my $is_insert=($begin==$end)?1:0;
	my $ins_pos=$begin;
		
	foreach my $gene(keys %gene_begin_hash){
		if(overlaps($begin, $end, $gene_begin_hash{$gene}, $gene_end_hash{$gene}) || 
			($is_insert && isIn($ins_pos, $gene_begin_hash{$gene}, $gene_end_hash{$gene}))){
			push @hit_genes, $gene;
		}
	}

	return(\@hit_genes);

}
#########################################################
# Pauline added 12-17-08
# for those indels occuring in exons,see how close to nearest intron, and
# the size of the itron

sub CheckFlankingIntrons {
	my $begin=shift;
        my $end=shift;
        my $transcript_id=shift;

        my $is_insert=($begin==$end)?1:0;
        my $ins_pos=$begin;
	
	 my $adjacent_downstream_intron_size = -999;
         my $distance_to_downstream_intron = -999;
	 my $adjacent_upstream_intron_size = -999;
         my $distance_to_upstream_intron = -999;

	my $trans_ori=$transcript_ori_hash{$transcript_id}[0];
#	print NORM_OUT_FH "orientation $trans_ori\n";

	# order exons by chromosome order, not orientation order
	my @ordered_exon_array;
	my @tmp_exon_array = @{$transcript_exon_hash{$transcript_id}};
	if ($trans_ori == -1) {
		for (my $i = $#tmp_exon_array; $i >= 0; $i--) {
			push (@ordered_exon_array, $tmp_exon_array[$i]);
		}
	} else {
		@ordered_exon_array = @tmp_exon_array;
	}
   # Check which exon it hits
       for(my $i=0; $i<=$#ordered_exon_array; $i++){
           my ($ex_begin, $ex_end)=split /#/, $ordered_exon_array[$i];
#		print NORM_OUT_FH "CheckFlankingIntrons exon $i $begin $end $ex_begin $ex_end\n";
          if(overlaps($begin, $end, $ex_begin, $ex_end) ||
                        ($is_insert && isIn($ins_pos, $ex_begin, $ex_end))){
                        # Check exons
#print NORM_OUT_FH "CheckFlankingIntrons It's in an exon intron\n";
# It's in the Exon PCN 
			if ($i >= 1) {
				my ($up_exon_begin, $up_ex_end) = split /#/, 
					$ordered_exon_array[$i-1];
				$distance_to_upstream_intron = $begin - $ex_begin + 1;
				 $adjacent_upstream_intron_size = $ex_begin - $up_ex_end ;
#			print NORM_OUT_FH "up stats $up_exon_begin $up_ex_end\n";
			}

			if ($i < $#ordered_exon_array ) {
				my ($down_ex_begin, $down_ex_end) = split /#/, $ordered_exon_array[$i+1];
				$distance_to_downstream_intron = $ex_end - $end + 1;
				$adjacent_downstream_intron_size = $down_ex_begin - $ex_end ; 
			} # end for downstream exon
	last;	
	} # end if overlaps
	}	# end check all exons
	return ($distance_to_upstream_intron,  $adjacent_upstream_intron_size,
		 $distance_to_downstream_intron,  $adjacent_downstream_intron_size);
} # end subroutine 
###############################################################################

sub CheckTranscript{
	my $begin=shift;
	my $end=shift;
	my $transcript_id=shift;
	my $hit_type;
	my $hit_info;

	# Get transcript extents
	my ($trans_begin, $trans_end)=split /#/, $transcript_extent_hash{$transcript_id};
	# Get transcript orientation
	my $trans_ori=$transcript_ori_hash{$transcript_id}[0];

	my $is_insert=($begin==$end)?1:0;
	my $ins_pos=$begin;

	if(overlaps($begin, $end, $trans_begin, $trans_end) || 
		($is_insert && isIn($ins_pos, $trans_begin, $trans_end))){
	# Check if it hits a transcript

		my $prev_ex_end=undef;
		my $prev_ex_begin=undef;

		# Check which exon it hits
		for(my $i=0; $i<=$#{$transcript_exon_hash{$transcript_id}}; $i++){
			my ($ex_begin, $ex_end)=split /#/, ${$transcript_exon_hash{$transcript_id}}[$i];

			if(overlaps($begin, $end, $ex_begin, $ex_end) ||
				($is_insert && isIn($ins_pos, $ex_begin, $ex_end))){
			# Check exons

				$hit_info="EXON." . ($i+1) . " $transcript_exon_id_hash{$transcript_id}[$i]";
				$hit_type=CheckCDS($begin, $end, $transcript_id);
				$hit_info.=" $hit_type";

				last;
			}elsif(defined($prev_ex_end) && collapsable($begin, $end, $prev_ex_end, $ex_begin)){
			# Check introns

				$hit_type="INTRON";
				$hit_info="INTRON." . ($i);
				
				#Compute distance to 5'ss and 3'ss
				my ($intron_begin, $intron_end);
				if($prev_ex_end<$ex_begin){
					# Forward
					($intron_begin, $intron_end)=($prev_ex_end, $ex_begin);
				}else{
					# Reverse
					($intron_begin, $intron_end)=($ex_end, $prev_ex_begin);
				}

				my ($five_ss, $three_ss);
				if($trans_ori==1){
					$five_ss=$begin-$intron_begin;
					$three_ss=$intron_end-$end;
				}else{
					$five_ss=$intron_end-$end;
					$three_ss=$begin-$intron_begin;
				}

				$hit_info.=" $transcript_exon_id_hash{$transcript_id}[$i-1] +$five_ss $transcript_exon_id_hash{$transcript_id}[$i] -$three_ss";

				if($five_ss<=$UPSTREAM_SS_REGION){
				# Check 5' splice site
					$hit_type="5'SS";
					$hit_info.=" $hit_type";
				}

				if($three_ss<=$DOWNSTREAM_SS_REGION){
				# Check 3' splice site
					$hit_type="3'SS";
					$hit_info.=" $hit_type";
				}
	
				last;
			}

			$prev_ex_end=$ex_end;
			$prev_ex_begin=$ex_begin;
		}

	}else{

	# Check if it is in promoter or downstream
		if($trans_ori==1){
			if(overlaps($begin, $end, $trans_begin-$UPSTREAM_PROMOTER_REGION, $trans_begin)){
				$hit_type="PROMOTER";
				$hit_info=sprintf("PROMOTER -%i", $trans_begin-$end);
			}elsif(overlaps($begin, $end, $trans_end, $trans_end+$DOWNSTREAM_GENOMIC_REGION)){
				$hit_type="DOWNSTREAM";
				$hit_info=sprintf("DOWNSTREAM +%i", $begin-$trans_end);
			}else{
				$hit_type="IRRELEVANT";
			}
		}elsif($trans_ori==-1){
			if(overlaps($begin, $end, $trans_begin-$DOWNSTREAM_GENOMIC_REGION, $trans_begin)){
				$hit_type="DOWNSTREAM";
				$hit_info=sprintf("DOWNSTREAM +%i", $trans_begin-$end);
			}elsif(overlaps($begin, $end, $trans_end, $trans_end+$UPSTREAM_PROMOTER_REGION)){
				$hit_type="PROMOTER";
				$hit_info=sprintf("PROMOTER -%i", $begin-$trans_end);
			}else{
				$hit_type="IRRELEVANT";
			}
		}else{
			die "transcript orientation for $transcript_id is not (1|-1)\n";
		}
	}
	

	return($hit_type, $hit_info);
}

###############################################################################

sub CheckCDS{
	my $begin=shift;
	my $end=shift;
	my $transcript_id=shift;

	my $is_insert=($begin==$end)?1:0;
	my $ins_pos=$begin;

	my ($cds_begin, $cds_end)=split /#/, $transcript_cds_hash{$transcript_id};

	if(overlaps($begin, $end, $cds_begin, $cds_end) || ($is_insert && isIn($ins_pos, $cds_begin, $cds_end))){
		return("CDS");
	}else{
		my ($tr_begin, $tr_end)=split /#/, $transcript_extent_hash{$transcript_id};
		if(overlaps($begin, $end, $tr_begin, $tr_end) || ($is_insert && isIn($ins_pos, $tr_begin, $tr_end))){
			if($transcript_ori_hash{$transcript_id}[0]==1){
				if($end<$cds_begin){
					return("5'UTR");
				}elsif($begin>$cds_end){
					return("3'UTR");
				}
			}elsif($transcript_ori_hash{$transcript_id}[0]==-1){
				if($end<$cds_begin){
					return("3'UTR");
				}elsif($begin>$cds_end){
					return("5'UTR");
				}
			}
		}else{
			return("Not in transcript.");
		}
	}

}
########################################################
# Pauline #12-05-08
sub blast_pairwise {
	my ($query_seq, $subject_seq) = @_;
	my $qtmpfile = "tmp007";
	my $stmpfile = "tmp008";
	system ("echo \">query\n$query_seq\n\" > $qtmpfile");	
	system ("echo \">subject\n$subject_seq\n\" > $stmpfile");
	my $result = `blast -p blastp -i $qtmpfile -j $stmpfile`; # | grep Identities | head -n 1`;

#	print NORM_OUT_FH "$result\n";
	return ($result);
}

# returns e-values, percent identity, percent aligned
sub 
get_blast_scores
{
	my ($blast_result, $query_protein_length) = @_;
	my @blast_lines = split (/\n/, $blast_result);
        my $i = 0; my $done = 0;
	while (!$done && $i < @blast_lines) {
	if ($blast_lines[$i] =~ /Score/ &&
                                $blast_lines[$i] =~ /Expect/) {
                        $done = 1;
                } else {
                        $i++;
                }
	}
	$blast_lines[$i] =~ s/^\s+//;
	my @tmp_a = split (/\s+/, $blast_lines[$i]);
	my $e_value = $tmp_a[7];
	$i++;
	$blast_lines[$i] =~ s/^\s+//;
	@tmp_a = split (/\s+/, $blast_lines[$i]);
	my $pid = $tmp_a[3];
	$pid =~ s/\(//; $pid =~ s/\)//;
	my ($id, $length_aligned) = split (/\//, $tmp_a[2]);
#print NORM_OUT_FH "tmp  $blast_lines[$i]\n";
#print NORM_OUT_FH  " id $pid length aligned $length_aligned\n";
#print NORM_OUT_FH  "e-value $e_value pid $pid \n";
	return ($e_value, $pid, 
		sprintf ("%.1f", $length_aligned/$query_protein_length * 100)  );
}

#################
# Pauline #12-08-08
sub
get_blast_alignment {

	my ($blast_result, $query_protein_length) = @_;

	my @blast_lines = split (/\n/, $blast_result);
	my $i = 0; my $done = 0;
	print NORM_OUT_FH $blast_result;
	while (!$done && $i < @blast_lines) {
#print $blast_lines[$i];
		if ($blast_lines[$i] =~ /Score/ && 
				$blast_lines[$i] =~ /Expect/) {
			$done = 1;
		} else { 
			$i++;
		}
	} 
	if (!$done) { return "";}
	# on score line
	do {
		$i++;
	}until ($blast_lines[$i] =~ /Query/) ;
#print "current line $blast_lines[$i]\n";
	my $alignedseq = "-" x $query_protein_length;	
#	print "befor  $alignedseq \n";
	my ($t, $start, $seq, $end, $s, $sstart, $sseq, $send);

	while ($blast_lines[$i] =~ /Query/) {
#		print "in blast line $blast_lines[$i]\n";
		($t, $start, $seq, $end) = split (/\s+/, $blast_lines[$i]);
		$i++; $i++;
		($s, $sstart, $sseq, $send) = split (/\s+/, $blast_lines[$i]);
		while ($seq =~ /\-/) {
# there are extra amino acids in the subject sequence. remove these from the
# subjct sequence
print NORM_FH_OUT "entering extra gap $seq $sseq\n"; 
			my $gap_loc = index ($seq, "-");
			substr ($seq, $gap_loc,1) = "";
			substr ($sseq, $gap_loc, 1) = "";
print NORM_FH_OUT "leaving extra gap $seq $sseq\n"
		}
		$i++;
		substr ($alignedseq, $start - 1, $end - $start + 1) = $sseq;
		$i++;
#		print "current blast line $blast_lines[$i]\n";
	}
#	print "\nstart $alignedseq end\n" ;
	return ($alignedseq);	
}
sub CheckProteinDomains {
# mutated protein must be same length as reference protein
# and aligned to it. checks which domains are missing
	my ($transcript_id, $ref_protein, $mutated_aligned_protein, $id) = @_;
#	print NORM_OUT_FH "entered CheckProteinDomains\n";	
	$_ = $ref_protein;	
	my $ref_protein_length = tr/[A-Za-z\-\ ]//;
	$_ = $mutated_aligned_protein;
	my $mutated_aligned_length =  tr/[A-Za-z\-\ ]//;
	if ($mutated_aligned_length != $ref_protein_length ) {
		print "ERROR in CheckProteinDomains! Ref $ref_protein_length mutated $mutated_aligned_length\n";
		exit (-1);
	}
	my $protein_info_arr_ref=$protein_domain_hash{$transcript_id};
	my @problems;
	if (!exists ($protein_domain_hash{$transcript_id} )) {
		return;
	}

# change to read blast alignment
# see that all the domains are covered
	my $domains_OK = 0;
	my $total_domains = 0; # scalar (@{$protein_info_arr_ref}); 

	 foreach my $protein_info(@{$protein_info_arr_ref}){
                my ($domain_begin, $domain_end, $domain_info)=split /\t/, $protein_info;
		unless ($domain_info =~ /Low complexity/) {
# don't count low complexity domains 
		$total_domains++;
		chomp ($domain_info);
		my $domain_length =  $domain_end -$domain_begin +1;
		my $ref_subseq = substr ($ref_protein, $domain_begin, 
					$domain_length);
# domains start with 0
		my $mut_subseq = substr ($mutated_aligned_protein, 
					$domain_begin,
					$domain_length);

		if ($ref_subseq eq $mut_subseq) {
#			print NORM_OUT_FH "$transcript_id pfam domain $domain_info OK  $domain_begin $domain_end\n";
#			print NORM_OUT_FH $ref_subseq . "\n$mut_subseq\n";
			$domains_OK++;
			# domain is present, everything is fine
		} else {
#			push (@problems, $protein_info);
#			print  NORM_OUT_FH  "enterd potential problems\n$ref_subseq\n$mut_subseq\n";
			my $match = 0; my $missing = 0; my $mismatch = 0;
			for (my $j = 0; $j < $domain_length; 
				$j++) {
				if (substr ($ref_subseq, $j, 1) eq 
				substr ($mut_subseq, $j, 1) ) {
					$match++;
				} elsif ( substr ($mut_subseq, $j, 1) eq "-") {
					$missing++;
				} elsif (substr ($ref_subseq, $j, 1) ne  substr ($mut_subseq, $j,1) ) {
					$mismatch++;	
				}
			} # end for gone through domain length 
			#print NORM_OUT_FH "matches $match missing " . 
			#$missing/$domain_length * 100 . 
			#" ( $missing / $domain_length )  mismatches " .
			#$mismatch / $domain_length * 100 . 
			#" ( $mismatch / $domain_length )\n";
			if ($match < $domain_length) {
				# something's mismatched or missing 
				my $percent_missing = sprintf ("%1.f", 
					$missing /  $domain_length * 100);
				my $percent_mismatch =  sprintf ("%1.f", 
					$mismatch / $domain_length * 100);
			print NORM_OUT_FH "$id $transcript_id Affected pfam domain ";
			chomp ($domain_info);
			chomp ($protein_info);
			my $protein_info_t =  $protein_info; $protein_info_t =~ s/\n/\t/g;
			print NORM_OUT_FH $protein_info_t .   			 " Missing: $percent_missing\% ($missing\/$domain_length\) Mismatch $percent_mismatch\% ($mismatch\/$domain_length\) \n";
			}		 
		} # end else , some part of domain is missing
	} # end exclude lo complexity domains
	}# end for each domain 	
	print NORM_OUT_FH "$id $transcript_id $domains_OK  / $total_domains domains OK " . ($domains_OK  / $total_domains) . "\n"; 
} # end CheckProteinDomains 
###############################################################################

sub CheckProtein{
	my $snp_begin=shift;
	my $snp_end=shift;
	my $transcript_id=shift;
	my @domain_hits;

	my $is_insert=($snp_begin==$snp_end)?1:0;
	my $ins_pos=$snp_begin;

	my $protein_info_arr_ref=$protein_domain_hash{$transcript_id};

	foreach my $protein_info(@{$protein_info_arr_ref}){
		my ($domain_begin, $domain_end, $domain_info)=split /\t/, $protein_info;
		if(overlaps($snp_begin, $snp_end, $domain_begin, $domain_end) || 
			($is_insert && isIn($ins_pos, $domain_begin, $domain_end))){
			push @domain_hits, $domain_info;
		}
	}

	return \@domain_hits;
}

###############################################################################

sub CheckSNPs{
	my $begin=shift;
	my $end=shift;
	my $reference_allele=shift;
	my @matches;

	foreach my $snp_struct (@snps_arr){
		if(overlaps($begin, $end, ${$snp_struct}{begin}, ${$snp_struct}{end}) ||
			($begin==${$snp_struct}{begin} && $end==${$snp_struct}{end})){

			my @dbSNP_allele_arr=split /\//, ${$snp_struct}{allele};
			my @nonref_allele_arr;
			foreach my $dbSNP_allele(@dbSNP_allele_arr){
				if($reference_allele ne $dbSNP_allele){
					push @nonref_allele_arr, $dbSNP_allele;
				}
			}
			my $nonref_allele_arr=join "/", @nonref_allele_arr;

			push @matches, "${$snp_struct}{dbsnp_id}:$nonref_allele_arr";
		}
	}
	return(\@matches);
}

###############################################################################

sub CheckRepeats{
	my $begin=shift;
	my $end=shift;
	my @hit_repeats;

	my $is_insert=($begin==$end)?1:0;
	my $ins_pos=$begin;
		
	foreach my $repeat(keys %repeat_hash){
		my ($r_id, $r_begin, $r_end)=split /#/, $repeat;
		if(overlaps($begin, $end, $r_begin, $r_end) || 
			($is_insert && isIn($ins_pos, $r_begin, $r_end))){
			my $repeat_size=$r_end-$r_begin;
			push @hit_repeats, "$r_id:$repeat_size";
		}
	}

	return(\@hit_repeats);

}

###############################################################################

#	sub distance_to_exon{
#		my $begin=shift;
#		my $end=shift;
#		my $gene_id=shift;
#		my @results;

#		foreach my $transcript_id(keys %{$gene_transcript_hash{$gene_id}}){
#			my $exons_region_arr_ref=$transcript_exon_hash{$transcript_id};
#			my $exons_id_arr_ref=$transcript_exon_id_hash{$transcript_id};
#			my $transcript_orientation=${$transcript_ori_hash{$transcript_id}}[0];

#			my $exon;
#			my ($exon_upstream_dist, $exon_downstream_dist);
#			my ($exon_upstream_exon_id, $exon_downstream_exon_id);

#			for(my $i=0; $i<=$#{$exons_region_arr_ref}; $i++){
#				my ($ex_begin, $ex_end)=split /-/, ${$exons_region_arr_ref}[$i];
#				my $downstream_dist=$ex_begin-$end;
#				my $upstream_dist=$begin-$ex_end;

#				if($downstream_dist>=0 && ($downstream_dist<$exon_downstream_dist || !defined($exon_downstream_dist))){
#					$exon_downstream_dist=$downstream_dist;
#					$exon_downstream_exon_id=${$exons_id_arr_ref}[$i];
#				}
#				if($upstream_dist>=0 && ($upstream_dist<$exon_upstream_dist || !defined($exon_upstream_dist))){
#					$exon_upstream_dist=$upstream_dist;
#					$exon_upstream_exon_id=${$exons_id_arr_ref}[$i];
#				}
#			}
#			push @results, "$transcript_id,$exon_upstream_dist,$exon_upstream_exon_id,$exon_downstream_dist,$exon_downstream_exon_id";
#		}

#		return(\@results);
#		
#	}

###############################################################################

sub isIn{
	my ($pos, $begin, $end)=@_;

	if($begin>$end){
		($begin, $end)=($end, $begin);
	}
	if($pos>$begin && $pos<$end){
		return(1);
	}else{
		return(0);
	}
}

sub collapsable{
	my ($abegin, $aend, $bbegin, $bend)=@_;

	if($abegin>$aend){
		($abegin, $aend)=($aend, $abegin);
	}
	if($bbegin>$bend){
		($bbegin, $bend)=($bend, $bbegin);
	}

	if($abegin>=$bbegin && $abegin<=$bend){
		return(1);
	}	
	if($aend>=$bbegin && $aend<=$bend){
		return(1);
	}	

	if($bbegin>=$abegin && $bbegin<=$aend){
		return(1);
	}	
	if($bend>=$abegin && $bend<=$aend){
		return(1);
	}	

}

sub overlaps{
        my $abegin = shift;
        my $aend= shift;
        my $bbegin = shift;
        my $bend= shift;

		#if($abegin==$aend){
		#	print STDERR "$abegin==$aend, should you be using 'collapsable'?\n";
		#}
		#if($bbegin==$bend){
		#	print STDERR "$abegin==$aend, should you be using 'collapsable'?\n";
		#}

        # Make sure begins are less than ends
        if($abegin>$aend){
                ($abegin, $aend)=($aend, $abegin);
        }
        if($bbegin>$bend){
                ($bbegin, $bend)=($bend, $bbegin);
        }

        my $alength=$aend-$abegin;
        my $blength=$bend-$bbegin;

        # Compute possible overlaps
        my $region1=$aend-$bbegin;
        my $region2=$bend-$abegin;

        # The smaller of the two regions is the overlap
        my $overlap=($region1<$region2)?$region1:$region2;

        # If the overlap is negative, there is no overlap
        $overlap=($overlap<0)?0:$overlap;
        $overlap=($overlap>$alength)?$alength:$overlap;
        $overlap=($overlap>$blength)?$blength:$overlap;

        return($overlap);
}

###############################################################################

sub revcomp{

=head2 str revcomp(str sequence);

This function will return the sequence reverse complemented.

=cut
    my $sequence = shift;

    $sequence =~ tr/atugcyrswkmbdhvnxATUGCYRSWKMBDHVNX/taacgryswmkvhdbnxTAACGRYSWMKVHDBNX/;
    my $revcomp = "";
    for (my $j=length($sequence)-1;$j>-1;$j--) {
        $revcomp .= substr($sequence,$j,1);
    }
    return $revcomp;

}

###############################################################################

sub translate{

	my %codon_table=(
		"CCC"=>"P", "CCT"=>"P", "CCA"=>"P", "CCG"=>"P",
		"TCC"=>"S", "TCT"=>"S", "TCA"=>"S", "TCG"=>"S",
		"ACC"=>"T", "ACT"=>"T", "ACA"=>"T", "ACG"=>"T",
		"GCC"=>"A", "GCT"=>"A", "GCA"=>"A", "GCG"=>"A",

		"CTC"=>"L", "CTT"=>"L", "CTA"=>"L", "CTG"=>"L",
		"TTC"=>"F", "TTT"=>"F", "TTA"=>"L", "TTG"=>"L",
		"ATC"=>"I", "ATT"=>"I", "ATA"=>"I", "ATG"=>"M",
		"GTC"=>"V", "GTT"=>"V", "GTA"=>"V", "GTG"=>"V",

		"CAC"=>"H", "CAT"=>"H", "CAA"=>"Q", "CAG"=>"Q",
		"TAC"=>"Y", "TAT"=>"Y", "TAA"=>" ", "TAG"=>" ",
		"AAC"=>"N", "AAT"=>"N", "AAA"=>"K", "AAG"=>"K",
		"GAC"=>"D", "GAT"=>"D", "GAA"=>"E", "GAG"=>"E",

		"CGC"=>"R", "CGT"=>"R", "CGA"=>"R", "CGG"=>"R",
		"TGC"=>"C", "TGT"=>"C", "TGA"=>" ", "TGG"=>"W",
		"AGC"=>"S", "AGT"=>"S", "AGA"=>"R", "AGG"=>"R",
		"GGC"=>"G", "GGT"=>"G", "GGA"=>"G", "GGG"=>"G"
	);

	my $cds=uc(shift);
	$cds=~s/U/T/g;

	my $protein="";
	my $length=length($cds);
	if($length%3){
		print STDERR "CDS is not a multiple of 3.\n";
	}
	
	for(my $i=0; $i<$length; $i+=3){
		my $codon=substr($cds, $i, 3);
		my $aa=$codon_table{$codon};
		if(!defined($aa)){
			print STDERR "Undefined codon?  '$codon'\n";
		}
		$protein.=$aa;
	}
	return($protein);
}

###############################################################################

sub isAminoAcidChangeConservative{
	my $blossum_hash_ref=shift;
	my $aa1_str=shift;
	my $aa2_str=shift;

	if(!defined($blossum_hash_ref)){
		die "Error on isAminoAcidChangeConservative. BLOSSUM Hash is undefined.\n";
	}

	if(length($aa1_str) ne length($aa2_str)){
		return(-100);
	}

	my $aa_str_len=length($aa1_str);

	for(my $i=0; $i<$aa_str_len; $i++){
		my $aa1=substr($aa1_str,$i,1);
		my $aa2=substr($aa2_str,$i,1);

		if(($aa1 eq " " && $aa2 ne " ") || ($aa2 eq " " && $aa1 ne " ")){
			return(-1);
		}else{
			if(${${$blossum_hash_ref}{$aa1}}{$aa2}<0){
				return(-1);
			}
		}
	}
	return(1);
}

sub loadBlossumMatrix{

# This is the blossum30 matrix;
	my $blossum_matrix=

"   A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
A  4 -1  0  0 -3  1  0  0 -2  0 -1  0  1 -2 -1  1  1 -5 -4  1  0  0  0 -7 
R -1  8 -2 -1 -2  3 -1 -2 -1 -3 -2  1  0 -1 -1 -1 -3  0  0 -1 -2  0 -1 -7 
N  0 -2  8  1 -1 -1 -1  0 -1  0 -2  0  0 -1 -3  0  1 -7 -4 -2  4 -1  0 -7 
D  0 -1  1  9 -3 -1  1 -1 -2 -4 -1  0 -3 -5 -1  0 -1 -4 -1 -2  5  0 -1 -7 
C -3 -2 -1 -3 17 -2  1 -4 -5 -2  0 -3 -2 -3 -3 -2 -2 -2 -6 -2 -2  0 -2 -7 
Q  1  3 -1 -1 -2  8  2 -2  0 -2 -2  0 -1 -3  0 -1  0 -1 -1 -3 -1  4  0 -7 
E  0 -1 -1  1  1  2  6 -2  0 -3 -1  2 -1 -4  1  0 -2 -1 -2 -3  0  5 -1 -7 
G  0 -2  0 -1 -4 -2 -2  8 -3 -1 -2 -1 -2 -3 -1  0 -2  1 -3 -3  0 -2 -1 -7 
H -2 -1 -1 -2 -5  0  0 -3 14 -2 -1 -2  2 -3  1 -1 -2 -5  0 -3 -2  0 -1 -7 
I  0 -3  0 -4 -2 -2 -3 -1 -2  6  2 -2  1  0 -3 -1  0 -3 -1  4 -2 -3  0 -7 
L -1 -2 -2 -1  0 -2 -1 -2 -1  2  4 -2  2  2 -3 -2  0 -2  3  1 -1 -1  0 -7 
K  0  1  0  0 -3  0  2 -1 -2 -2 -2  4  2 -1  1  0 -1 -2 -1 -2  0  1  0 -7 
M  1  0  0 -3 -2 -1 -1 -2  2  1  2  2  6 -2 -4 -2  0 -3 -1  0 -2 -1  0 -7 
F -2 -1 -1 -5 -3 -3 -4 -3 -3  0  2 -1 -2 10 -4 -1 -2  1  3  1 -3 -4 -1 -7 
P -1 -1 -3 -1 -3  0  1 -1  1 -3 -3  1 -4 -4 11 -1  0 -3 -2 -4 -2  0 -1 -7 
S  1 -1  0  0 -2 -1  0  0 -1 -1 -2  0 -2 -1 -1  4  2 -3 -2 -1  0 -1  0 -7 
T  1 -3  1 -1 -2  0 -2 -2 -2  0  0 -1  0 -2  0  2  5 -5 -1  1  0 -1  0 -7 
W -5  0 -7 -4 -2 -1 -1  1 -5 -3 -2 -2 -3  1 -3 -3 -5 20  5 -3 -5 -1 -2 -7 
Y -4  0 -4 -1 -6 -1 -2 -3  0 -1  3 -1 -1  3 -2 -2 -1  5  9  1 -3 -2 -1 -7 
V  1 -1 -2 -2 -2 -3 -3 -3 -3  4  1 -2  0  1 -4 -1  1 -3  1  5 -2 -3  0 -7 
B  0 -2  4  5 -2 -1  0  0 -2 -2 -1  0 -2 -3 -2  0  0 -5 -3 -2  5  0 -1 -7 
Z  0  0 -1  0  0  4  5 -2  0 -3 -1  1 -1 -4  0 -1 -1 -1 -2 -3  0  4  0 -7 
X  0 -1  0 -1 -2  0 -1 -1 -1  0  0  0  0 -1 -1  0  0 -2 -1  0 -1  0 -1 -7 
* -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7 -7  1";

	my %blossum_hash;
	my @blossum_rows=split /\n/, $blossum_matrix;
	my @col_aa_arr=split /\s+/, $blossum_rows[0];
	for(my $i=1; $i<=$#blossum_rows; $i++){
		my @row_aa_arr=split /\s+/, $blossum_rows[$i];
		my $row_aa=$row_aa_arr[0];
		for(my $j=1; $j<=$#col_aa_arr; $j++){
			my $col_aa=$col_aa_arr[$j];
			if(length($col_aa)==1 && ($col_aa=~/[A-Z\*]/)){
				if($row_aa_arr[$j]=~/\d+/){
					${$blossum_hash{$col_aa}}{$row_aa}=$row_aa_arr[$j];
				}else{
					print STDERR "Error Reading BLOSSUM Matrix.  Score expected.\n";
					print STDERR "Received '$row_aa_arr[$j]'.\n";
					die;
				}
			}else{
				print STDERR "Error Reading BLOSSUM Matrix. Single Letter Amino Acid Expected.\n";
				print STDERR "Received '$col_aa'.\n";
				die;
			}
		}
	}
	return(\%blossum_hash);
}


###############################################################################

sub
translate_with_downstream_Met
{
    my ($donor_cds_seq, $mut_location) = @_;

    my $new_cds_seq = substr ($donor_cds_seq, $mut_location, 30000000);
    $new_cds_seq = uc ($new_cds_seq);
	my $loc = index ($new_cds_seq, "ATG");
   	$new_cds_seq = substr ($new_cds_seq, $loc, 30000000);
	my $new_aa = translate ($new_cds_seq);
	return ($new_aa);
} 

##############################################################
# Pauline 12-18-08
sub 
run_clustal
{
	my ($refer_protein, $donor_protein) = @_;
	open (CLUSTAL_OUT, ">tmp.clustal") || die "can't open tmp.clustal\n";
	print CLUSTAL_OUT ">ReferenceProtein\n$refer_protein\n>DONOR\n$donor_protein\n";
	close (CLUSTAL_OUT);
	system ("clustalw -INFILE=tmp.clustal -OUTPUT=PIR");
	open (DB, "tmp.pir") || die "can't open clustal file";
	my $seekpos = tell DB;
	my ($d, $ref_sequence) =  read_a_sequence (\$seekpos);
	my ($d1, $donor_sequence) =  read_a_sequence (\$seekpos);
	$ref_sequence =~ s/\*//; $ref_sequence =~ s/\s+//;
	$donor_sequence =~ s/\*//; $donor_sequence =~ s/\s+//;

	my $identical = 0; my $rseq =""; my $dseq = "";
	my $aligned_length = $ref_sequence =~ tr/[A-Za-z\-]//;
	for (my $i = 0; $i < $aligned_length; $i++) {
		my $rbase = substr ($ref_sequence, $i, 1);
		my $dbase = substr ($donor_sequence, $i,1);

		if ($rbase ne "-" ) {
			$rseq = $rseq .  $rbase;
			$dseq = $dseq . $dbase;

			if ($rbase eq $dbase) {
				$identical++;
			}
		}
	} # gone through all of alignement
	my $num_gaps = $dseq =~ tr/\-//;
	my $ref_length = $refer_protein =~ tr/[A-Za-z]//;
 
	my $percent_aligned = sprintf ("%.2f", ( 1 - ($num_gaps/$ref_length)) * 100);	
	my $percent_identical = sprintf ("%.2f", $identical / $ref_length * 100); 
	return ( $percent_identical,  $percent_aligned ,  $dseq);

}
