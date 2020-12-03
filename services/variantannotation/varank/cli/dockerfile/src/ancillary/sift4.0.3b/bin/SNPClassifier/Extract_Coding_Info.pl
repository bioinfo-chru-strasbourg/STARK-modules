#!/usr/local/bin/perl

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

######################################################################
#
#	Revision: $Revision: 4813 $
#	Date:     $Date: 2009-10-21 09:46:23 -0400 (Wed, 21 Oct 2009) $
#	Author:   $Author: kli $
#
######################################################################

my $SPLASH="Copyright (c) 2005 J. Craig Venter Institute.\nAll rights reserved.";
my $VERS_INFO=qw$Revision: 4813 $;

print STDERR "\n\"$0\"\n$SPLASH\nRevision: $VERS_INFO\n\n";

######################################################################

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::TranslationAdaptor;
use Bio::EnsEMBL::Variation::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Variation::Variation;

use Getopt::Std;

use vars qw($opt_c $opt_b $opt_e $opt_O $opt_B $opt_A $opt_x $opt_f $opt_n);

getopts("c:b:e:O:B:A:xf:n");
my $usage = "usage: 
$0 
	[required]
	-c chromosome
	-b begin (absolute ensemble 1-residue based coordinates)
	-e end (absolute ensemble 1-residue based coordinates)

	-O \"<organism>\", for example \"homo sapiens\"
	-B <ncbi build version>, for example 21
	-A <ncbi assembly version>, for example 34d

	[optional]
	-x automatically expand extraction range to include
		all the exons overlapping the range specified.
	-f preferred output filename
	-n skip SNP extraction (because lite database is not available)

	This program will:

	Read in gene, intron, exon, CDS, and SNP information that was
	found in the range of begin and end at the specified chromosome.
	
";

if(!(
	defined($opt_c) && 
	defined($opt_b) && 
	defined($opt_e) && 
	defined($opt_O) && 
	defined($opt_B) && 
	defined($opt_A))){
	print STDERR $usage;
	die;
}

my $organism=lc($opt_O);
my $build_ver=lc($opt_B);
my $assembly_ver=lc($opt_A);
my $skip_snps=$opt_n;

###############################################################################

my $host;
my $user;
my $core_dbname;
my $snp_dbname;
my $password;

$host='tcagdev0.jtc.jcvsf.org';
$user='webuser';
$core_dbname='homo_sapiens_core_19_34a';
$password='webpass';

# SNP & SNP Database
$host='ensembldb.ensembl.org';
$user='anonymous';

$organism=~s/ /_/;

$core_dbname="$organism\_core_$build_ver\_$assembly_ver";

if(!$skip_snps){
	$snp_dbname="$organism\_variation_$build_ver\_$assembly_ver";
}

###############################################################################

print STDERR "Using $core_dbname for core data...\n";

if($skip_snps){
	print STDERR "Skipping dbSNP extraction.\n";
}else{
	print STDERR "Using $snp_dbname for SNP data...\n";
}

print STDERR "Are you sure these are the most recent???\n";

###############################################################################

print STDERR "Connection to core database...\n";
my $core_db=new Bio::EnsEMBL::DBSQL::DBAdaptor(
	-host => $host, 
	-user => $user, 
	-dbname => $core_dbname,
	-port => 5306
);
print STDERR "Done.\n";

my $snp_db;
if(!$skip_snps){
	print STDERR "Connection to snp database...\n";
	$snp_db=new Bio::EnsEMBL::Variation::DBSQL::DBAdaptor(
		-host => $host, 
		-user => $user, 
		-dbname => $snp_dbname,
		-port =>5306
	);
	print STDERR "Done.\n";
}

###############################################################################

if(!$skip_snps){
	$core_db->add_db_adaptor('variation', $snp_db);
	#$snp_db->add_db_adaptor('core', $core_db);
	$snp_db->dnadb($core_db);
}

###############################################################################

my ($chrom, $begin, $end)=($opt_c, $opt_b, $opt_e);

###############################################################################

my $slice_adaptor=$core_db->get_SliceAdaptor();

if($opt_x){

	my $survey_slice=$slice_adaptor->fetch_by_region('chromosome', $chrom, $begin, $end);

	print STDERR "Getting all transcripts in the area of $chrom.$begin-$end\n";

	# Keep track of outer extents of exons regions
	my $min_begin;
	my $max_end;

	my $survey_genes_ref=$survey_slice->get_all_Genes();
	foreach my $gene(@{$survey_genes_ref}){

		my $transcripts_ref=$gene->get_all_Transcripts();

		foreach my $transcript(@{$transcripts_ref}){
			my $exons_ref=$transcript->get_all_Exons();
			my $display_id=$transcript->display_id();

			foreach my $exon(@{$exons_ref}){
				if($min_begin>$exon->start()){
					$min_begin=$exon->start();
				}
				if($max_end<$exon->end()){
					$max_end=$exon->end();
				}
			}

		}

	}

	my $orig_begin=$begin;
	my $orig_end=$end;
	my $orig_length=$end-$begin+1;

	my $begin_extension=($min_begin<0)?(-$min_begin+1):0;
	my $end_extension=($max_end-$orig_length);
	$end_extension=($end_extension<0)?0:$end_extension;

	$begin-=$begin_extension;
	$end+=$end_extension;

	printf STDERR "Extending begin by $begin_extension bps\n";
	printf STDERR "Extending end by $end_extension bps\n";

	print STDERR "Resetting extraction range to $chrom.$begin-$end\n";

}

###############################################################################


my $output_filename;

if(defined($opt_f)){
	$output_filename=$opt_f;
}else{
	$output_filename="$chrom\.$begin-$end";
}


open(OUTPUT_FH, ">$output_filename\.coding_info") || 
	die "Could not open $output_filename\.coding_info\n";
print STDOUT "$output_filename\.coding_info\n";

###############################################################################

print OUTPUT_FH join "\t", ("EXTRACTION", $chrom, $begin-1, $end, "\n");

print STDERR "Getting coding information on $chrom.$begin-$end\n";

my $slice=$slice_adaptor->fetch_by_region('chromosome', $chrom, $begin, $end);
my $sequence=$slice->seq();
my @transcripts_arr;

print STDERR "Getting Gene information...\n";
my $genes_ref=$slice->get_all_Genes();
foreach my $gene(@{$genes_ref}){
	
	my $gene_name=$gene->display_id();

	# Get Transcript information for this gene
	my $transcripts_ref=$gene->get_all_Transcripts();

	my $human_readable_name=$gene->external_name();
	my $display_id=$gene->display_id();
        my $full_gene_description=$gene->description();
	if($full_gene_description =~ /(.*)\s\[/){
	        $full_gene_description = $1;
	}
	my @gene_description_arr = split(/\(/,$full_gene_description);
	my $gene_description = $gene_description_arr[0];
	print OUTPUT_FH join "\t", ("HUGO_ID", $display_id, $human_readable_name, "\n");
	print OUTPUT_FH join "\t", ("PROT_DESC", $display_id, $gene_description, "\n");
	foreach my $transcript(@{$transcripts_ref}){
		push @transcripts_arr, $transcript;

		my $cds_start=$transcript->coding_region_start();
		my $cds_end=$transcript->coding_region_end();
		my $transcript_name=$transcript->display_id();
		my $exons_ref=$transcript->get_all_Exons();

		foreach my $exon(@{$exons_ref}){
			my $exon_name=$exon->display_id();
			my $exon_start=$exon->start();
			my $exon_end=$exon->end();
			my $exon_strand=$exon->strand();
			print OUTPUT_FH join "\t", 
				("EXON", $gene_name,$transcript_name,$exon_name,$exon_start-1,$exon_end,$exon_strand, "\n");
		}

		print OUTPUT_FH join "\t",
			("CDS", $gene_name,$transcript_name,$cds_start-1,$cds_end,"\n");
	}

}

###############################################################################

if(!$skip_snps){
	print STDERR "Getting SNP information...\n";
	my $snp_adaptor=$snp_db->get_VariationFeatureAdaptor();
	my $variation_feature_array_ref=$snp_adaptor->fetch_all_by_Slice($slice);

	foreach my $snp(@{$variation_feature_array_ref}){
		my ($dbsnp_id, $start, $stop, $strand, $score, $source, $alleles)=(
			$snp->display_id(),
			$snp->start(),
			$snp->end(),
			$snp->strand(),
			"", #$snp->score(),
			$snp->source(),
			$snp->allele_string()
		);

		my ($tcag_start, $tcag_stop);
		if($start <= $stop){
			# Substitution or deletion
			($tcag_start, $tcag_stop)=($start-1, $stop);
		}elsif(($start-$stop)==1){
			# Insertion
			($tcag_start, $tcag_stop)=($start, $start);
		}else{
			die "Could not understand coordinate system for SNP.\n";
		}

		print OUTPUT_FH join "\t",
			("SNP", $dbsnp_id, $tcag_start, $tcag_stop, $strand, $alleles, "\n");
	}
}

###############################################################################
# Extract Protein Information
print STDERR "Getting Protein information...\n";
my $translation_adaptor = $core_db->get_TranslationAdaptor();

foreach my $transcript(@transcripts_arr){

	my $transcript_id = $transcript->display_id();
	my $translation = $translation_adaptor->fetch_by_Transcript($transcript);

	if(!defined($translation)){
		next;
	}
	my $protein_id = $translation->display_id();

	my %protein_features=(
		"Prints"=>"prints", 
		"pfscan"=>"?", 
		"scanprosite"=>"prosite", 
		"Signalp"=>"Signal_peptide", 
		"Seg"=>"Low_complexity", 
		"ncoils"=>"Coiled_coil", 
		"tmhmm"=>"Transmembrane", 
		"Pfam"=>"Pfam"
	);

	my @results;

	foreach my $prot_feat(keys %protein_features){

		my $prot_feat_ref=$translation->get_all_ProteinFeatures($prot_feat);

		foreach my $result(@{$prot_feat_ref}){
			my ($prot_feat_id, $start, $end, $hstart, $hstend, $desc, $interpro_ac)=
				($result->{hname},
				$result->{start}-1,
				$result->{end},
				$result->{hstart},
				$result->{hend},
				$result->idesc(),
				$result->interpro_ac()
			);

			my $full_desc=$protein_features{$prot_feat};
			$full_desc.=($desc ne "")?":$desc":"";

			my $outstring=join "\t", ("PROTEIN", $transcript_id, $protein_id, $full_desc, $start, $end);
			print OUTPUT_FH "$outstring\n";
		}
	}
}

###############################################################################
# Fetch Repeat regions
print STDERR "Getting Repeat information...\n";

my $repeat_features_ref=$slice->get_all_RepeatFeatures();
foreach my $repeat_feature(@{$repeat_features_ref}){

	my ($start, $stop, $strand, $display_id, $score, $hstart, $hend)=(
		$repeat_feature->start(),
		$repeat_feature->end(),
		$repeat_feature->strand(),
		$repeat_feature->display_id(),
		$repeat_feature->score(),
		$repeat_feature->hstart(),
		$repeat_feature->hend()
	);

	my $outstring=join "\t", ("REPEAT", $display_id, $start-1, $stop, "\n");
	print OUTPUT_FH "$outstring";
}


OutputFASTAFile("$output_filename\.fasta", $chrom, $begin, $end);

###############################################################################

print STDERR "done.\n";

###############################################################################

sub OutputFASTAFile{

# This function outputs the file specified by filename at the coordinates
# given by chrom, begin and end.

    my $filename=shift;
    my $chrom=shift;
    my $begin=shift;
    my $end=shift;

    print STDERR "Outputting: $filename\n";
    my $slice_adaptor=$core_db->get_SliceAdaptor();
    my $slice=$slice_adaptor->fetch_by_region('chromosome', $chrom, $begin, $end);
    my $sequence=$slice->seq();

    open (FASTA_FH, ">$filename") || die "Could not open $filename for writing sequence.\n";

    # Output defline
    print FASTA_FH ">$chrom\.$begin-$end\n";

    # Output sequence.
    my $length=length($sequence);
    my $width=50;
    my $pos=0;
    do{
        my $out_width=($width>$length)?$length:$width;
        print FASTA_FH substr($sequence, $pos, $width) . "\n";
        $pos+=$width;
        $length-=$width;
    }while($length>0);

    close(FASTA_FH);
}


