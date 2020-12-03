#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software
#
# October 11, 2009 created as a script option to SIFT indels
#-------------------------------------------------------------------------
use DBI;
use Tie::IxHash;
use File::Copy;
use Getopt::Std;
$ENV{'SIFT_HOME'} = '/usr/local/projects/SIFT/sift4.0.1';
$SIFT_HOME = $ENV{'SIFT_HOME'};
use vars qw($opt_i $opt_c $opt_d $opt_o);
getopts("i:c:d:o:");
my $usage = "usage: 
$0 
        -i <Query indels filename with complete path>
        -c <coding info directory path>
	-d <Variation db directory path>
        -o <Optional: output file with complete path - default=$SIFT_HOME/tmp>

        All values should be in local 0 space based coordinates.
        
";
$| = 1;
if(!(
        defined($opt_i) && 
	defined($opt_d) &&
        defined($opt_c))){
        print STDERR $usage;
        die;
}
#       Set file permissions to rw-rw----
system("umask 006");
my $bin             = "$SIFT_HOME/bin";
my $pid             = $$;
my $tmp = defined($opt_o) ? "$opt_o/$pid" : "$SIFT_HOME/tmp/$pid";
require "$bin/SIFT_subroutines.pm";
mkdir $tmp;
chmod 0777, $tmp;
my $coding_info_dir = "$opt_c";
my $Variation_db_dir = "$opt_d";
$snp_classifier     = "$bin/Classify_SNPs_Indels.pl";
$reformat_chrfile   = "$bin/reformat_chrfile.pl";
$detect_indel 	    = "$bin/detect_indel.pl";
$model_transcript   = "$bin/model_transcript.pl";
$small_indels_error_intersect_file = "$tmp/$pid.small_indels_error_intersect.gff";
$small_indels_error_gff = "$SIFT_HOME/db/small_indel_errors_ref.gff";
$repeat_db_file = "$SIFT_HOME/db/repeat_db.gff";
$indelfile_to_gff = "$bin/indelfile_to_gff.pl";
$intersect_locations = "$bin/IntersectLocations.sh";
$detect_repeat = "$bin/detect_repeat.pl";

#get user output options
$all_transcripts = 1;
my $output_options = "";

$chr_file = $opt_i;
## Check for validity of user inputs
my $all_chr_file = $tmp . "/$pid.chrfile";

copy($chr_file, $all_chr_file) or die "File cannot be copied.";

$bin_file = "$coding_info_dir/bins.list";
chmod 0777, $all_chr_file;
#check input validity and reformat chrfile
system("perl $reformat_chrfile $all_chr_file $tmp");

#display wait message
print
"Your job id is $pid and is currently running.\n";


#convert chrfile to gff for intersecting with small indel errors file.
system("perl $indelfile_to_gff $all_chr_file");
$all_chr_gff_file = "$all_chr_file.gff";

#Now intersect gff chrfile with small indels ref file
system("$intersect_locations $all_chr_gff_file gff $small_indels_error_gff gff simple > $small_indels_error_intersect_file");

#create an index %user_small_indels_error_index with key as "chr#startindex"  stop index redundant also problem with Intersect script not detecting intersect when start = stop i.e. insertion.
my %user_small_indels_error_index;
if (-e $small_indels_error_intersect_file){
	open (SMALL_INDELS,"$small_indels_error_intersect_file"); 
	while (<SMALL_INDELS>){
		next if ($_ =~ /^\d+|^Total/i);
		chomp;
		my @elts = split /\t/, $_;
		my $chr = $elts[1];
		my $start = $elts[2];
		my $key = "$chr\t$start";
		$user_small_indels_error_index{"$key"} = 1;	
	}
}

#Creat binned SNP files for SNP Classifier
system("$bin/map_coords_to_bin_indels.pl $bin_file $all_chr_file $pid $tmp");
#create index of old coords to new coords using pid_old_new_coords_map.txt
open( COORDS_MAP, "$tmp/$pid\_old_new_coords_map.txt" )
  || die("Cannot open old-new coords map file");
while (<COORDS_MAP>) {
	chomp;
	@elts       = split /\t/, $_;
	$old_coords = @elts[0];
	$new_coords = @elts[1];
	@elts2 = split /,/, $new_coords;
	$chr = $elts2[0];
	$new_start = @elts2[1];
	$new_stop = @elts2[2];
	$key = "$chr\t$new_start\t$new_stop";
	#print "$key *** $old_coords<br>";
	$index_old_new_coords{$key} = $old_coords;
}
#Run snp classifier on files in snp_chr_map_file
open( SNP_CHR_MAP_FILE, "$tmp/$pid\_snp_chr_map_file.txt" )
  || die("Cannot open SNP_CHR_MAP_FILE");
while (<SNP_CHR_MAP_FILE>) {
	chomp;
	@cols             = split /\t/, $_;
	$chromosome       = @cols[0];
	$bucket           = @cols[1];
	$snp_file         = "$tmp/@cols[2]";
	$coding_info_file = "$coding_info_dir/@cols[3]";
	$chr_fasta_file   = "$coding_info_dir/@cols[4]";
	#print "$snp_classifier -s $snp_file -c $coding_info_file -n $chr_fasta_file -o $tmp/$pid\_snps_classified_chr$chromosome\_bin$bucket<br>";
	system(
"$snp_classifier -s $snp_file -c $coding_info_file -n $chr_fasta_file -o $tmp/$pid\_snps_classified_chr$chromosome\_bin$bucket"
	);

}
close(SNP_CHR_MAP_FILE);
#Parse denormal files to extract coding info for input coords
open( SNP_CHR_MAP_FILE, "$tmp/$pid\_snp_chr_map_file.txt" )
  || die("Cannot open SNP_CHR_MAP_FILE");
my %coord_result_hash;
my $en_trancript_id_string;
while (<SNP_CHR_MAP_FILE>) {
	chomp;
	@cols             = split /\t/, $_;
	$chromosome       = @cols[0];
	$bucket           = @cols[1];
	$snp_file         = "$tmp/@cols[2]";
	$coding_info_file = "$coding_info_dir/chr$chromosome/@cols[3]";
	$chr_fasta_file   = "$coding_info_dir/chr$chromosome/@cols[4]";
	$denormal_file    =
	  "$tmp/$pid\_snps_classified_chr$chromosome\_bin$bucket.denormal";
	$normal_file    =
          "$tmp/$pid\_snps_classified_chr$chromosome\_bin$bucket.normal";

	open( DENORMAL, $denormal_file ) || die "(denormal file does not exist)";
	while (<DENORMAL>) {
		@elts        = split /\t/, $_;
		$coord       = @elts[1];
		@coord_elts  = split /\:/, $coord;
		$allele      = $coord_elts[scalar @coord_elts - 1];
		$coord       = join( ',', @coord_elts );
		@coord_elts2 = split( '-', $coord, 2 );
		$coord       = join( ',', @coord_elts2 );
		@elts2 = split /,/, $coord;
		$new_start = @elts2[0];
		$new_stop = @elts2[1];
		$key = "$chromosome\t$new_start\t$new_stop";
		$old_coord= $index_old_new_coords{$key};
		$coord       = $old_coord;
		$nu_change = $elts[2];
                if ($nu_change =~ /([ATGC]+) ([ATGC]+) ([ATGC]+)/){
                        $left_flank = uc substr ($1,-5,5);
                        $deletion = lc $2;
                        $right_flank = uc substr($3,0,5);
                        $nu_change_modified = "$left_flank-$deletion-$right_flank";
                }
                elsif($nu_change =~ /([ATGC]+) - ([ATGC]+)/){
                        $left_flank = lc substr ($1,-5,5);
                        $insertion = uc $allele;
                        $right_flank = lc substr($2,0,5);
                        $nu_change_modified = "$left_flank-$insertion-$right_flank";
                }
		$nu_change = $nu_change_modified;

		if ( $_ =~ /AA_DELETION | AA_INSERTION | FRAMESHIFT/x ) {
			$subst = $&;
			$region = "CDS";
                        if ( $_ =~ /^(\d+).+?ENSG(\d+).+?ENST(\d+).+/) {
				$hit_id = $1;
                                $en_transcript_id = "ENST$3";
				$ensg = "ENSG$2";
                                $result = "$hit_id\t$ensg\t$en_transcript_id\t$subst\t$region";
				#append seq (original and modified) to result
				$original_seq = "$tmp/$pid\_$en_transcript_id.fa";
				$modified_seq = "$tmp/$pid\_$en_transcript_id\_$hit_id.fa";
				if (-e $modified_seq && -e $original_seq){
					$result.="\t$original_seq\t$modified_seq\t$nu_change";
				}
				push @{ $coord_result_hash{$coord} }, $result;		#hash of arrays: multiple ENSTs per coord
                        }
                }
		elsif ( $_ =~ /INTRON|PROMOTOR|DOWNSTREAM|UPSTREAM/ix ) {
                        $subst = "";
			$region = "$&";
                        if ( $_ =~ /^(\d+).+?ENSG(\d+).+?ENST(\d+).+/) {
                                $hit_id = $1;
                                $en_transcript_id = "ENST$3";
                                $ensg = "ENSG$2";
                                $result = "$hit_id\t$ensg\t$en_transcript_id\t$subst\t$region\t\t\t$nu_change";
                                #append seq (original and modified) to result
                                push @{ $coord_result_hash{$coord} }, $result;          #hash of arrays: multiple ENSTs per coord
                        }
                }
		
		elsif ( $_ =~ /CDS/i ) {
                        $subst = "";
                        $region = "CDS";
                        if ( $_ =~ /^(\d+).+?ENSG(\d+).+?ENST(\d+).+/) {
                                $hit_id = $1;
                                $en_transcript_id = "ENST$3";
                                $ensg = "ENSG$2";
                                $result = "$hit_id\t$ensg\t$en_transcript_id\t$subst\t$region\t\t\t$nu_change";
                                #append seq (original and modified) to result
                                push @{ $coord_result_hash{$coord} }, $result;          #hash of arrays: multiple ENSTs per coord
                        }
                }
			
	}
	close (DENORMAL);
	unlink ("$denormal_file");
	unlink ("$normal_file");
}

#create original - modified protein files
my $repeat_detect_input_file = "$tmp/$pid\_repeat_detect_input_file.txt";
open (REPEAT_INPUT,">$repeat_detect_input_file");

foreach $coord (keys %coord_result_hash){
	my @elts;
	my $repeat_input_row = "";
	my @indel_result;
	my @indel_location_features;
	my $chr = (split /,/, $coord)[0];
	my $indel_start = (split /,/, $coord)[1];
	my $indel_stop = (split /,/, $coord)[2];
	@result_set = @{ $coord_result_hash{$coord} };
	foreach $result (@result_set){
		chomp $result;
		@elts = split /\t/,$result;
		my $hit_id = $elts[0];
		my $ensg = $elts[1];
		my $enst = $elts[2];
		my $subst = $elts[3];
		my $region = $elts[4];
		my $original_seq = $elts[5];
		my $modified_seq = $elts[6];
		@indel_location_features = split /\t/, `$model_transcript $enst $indel_start $indel_stop $Variation_db_dir`;
		$indel_start_CDS = $indel_location_features[0];
		$indel_stop_CDS = $indel_location_features[1];
		$indel_loc_pct =  $indel_location_features[2];
		$indel_span =  $indel_location_features[3];
		$num_5UTR = $indel_location_features[4];
		$num_exon = $indel_location_features[5];
		$num_intron = $indel_location_features[6];
		$num_3UTR = $indel_location_features[7];
		$tx_length = $indel_location_features[8];
		$CDS_length = $indel_location_features[9];
		$protein_length = $indel_location_features[10];
		$last_exon_length = $indel_location_features[11];
		$visual_tx = $indel_location_features[12];
		$NMD_protein_region = int(($last_exon_length + 50)/3);
		$result.="\t$indel_start_CDS\t$indel_stop_CDS\t$indel_loc_pct\t$indel_span\t$visual_tx";
		#print "$last_exon_length $NMD_protein_region $original_length $modified_length<br>";
		my $outfile = "";
		if (-e $original_seq && -e $modified_seq){
			$outfile = "$tmp/$pid\_$enst\_$hit_id\_comparison.fasta";
			@indel_result = split /\t/, `$detect_indel $original_seq $modified_seq $outfile`;
			$aa_coords_original =  $indel_result[0];
			$aa_coords_modified =  $indel_result[1];
			$aa_change_original =  $indel_result[2];
			$aa_change_modified =  $indel_result[3];
			$original_length = $indel_result[4];
			$modified_length = $indel_result[5];
			if ($num_exon < 2) {
				$causes_NMD = "N\/A";
			}
			elsif ($original_length - $modified_length > $NMD_protein_region){
				$causes_NMD = "YES";
			}
			else{
				$causes_NMD = "NO";
			}
			chomp $aa_coords_original,$aa_coords_modified,$aa_change_original,$aa_change_modified;	
			#print "@indel_result<br>";
			$result.="\t$aa_coords_original\t$aa_change_original->$aa_change_modified\t$causes_NMD";
			$repeat_input_row = "$chr\t$enst\t$indel_start\t$indel_stop\t".join ("\t", split(/\-/,$aa_coords_original) )."\t$original_seq\t$modified_seq";
			print REPEAT_INPUT "$repeat_input_row\n";	
		}
	}	
	$coord_result_hash{$coord} = [ @result_set ];
	
}
close (REPEAT_INPUT);

#Run the detect_repeat script on the repeat_input file created in previous loop
my $repeat_detect_output_file = "$tmp/$pid\_repeat_detect_output_file.txt";
system("$detect_repeat $repeat_detect_input_file $repeat_db_file > $repeat_detect_output_file");

my %repeat_index;
if (-e $repeat_detect_output_file){
	open (REPEAT_OUTPUT,"$repeat_detect_output_file");
	while (<REPEAT_OUTPUT>){
		chomp;
		my @elts = split /\t/, $_;
		my $chr = $elts[0];
		my $indel_start  = $elts[2];
		if ($elts[10] ne "" && $elts[11] ne ""){	
			my $repeat_change = $elts[10]."-->".$elts[11];
			$repeat_index{"$chr\t$indel_start"} = $repeat_change;
		}
	}	
}


## SIFT prediction operations

$exp_option = 1;

#$info = $names{info};
$comments            = "$tmp/$pid.comments";
$out                 = $tmp . "/$$.siftresults";


#build combined classification (no predictions yet for indels)  file from  %coord_result_hash, hash of arrays. Also get gene infor from Human_supp db
$db_supp = DBI->connect( "dbi:SQLite:dbname=$Variation_db_dir/Human_Supp.sqlite","", "", { RaiseError => 1, AutoCommit => 1 } );
$db_supp->do('PRAGMA synchronous=1');
$sth_db_supp_geneinfo = $db_supp->prepare("select * from GENE_INFO where ENST = ?");


$combined_siftresults_file = "$pid.combined_siftresults";
open (CLASSIFICATION_FILE,">$tmp/$pid.classification_file");
foreach $coord (keys %coord_result_hash) {
	my @coord_elts = split /,/, $coord;
	
	#warning if user coordinates intersect with small indel errors
	my $chr = $coord_elts[0];
	my $coord_start = $coord_elts[1];
	my $key_for_small_indels = "$chr\t$coord_start";
	if ($user_small_indels_error_index{$key_for_small_indels} == 1){
		$small_indel_warning = "Warning: NCBI reference miscall!";
	}
	else{
		$small_indel_warning = "";
	}
	
	my $repeat_change = $repeat_index{$key_for_small_indels};  
	#Populate CLASSIFICATION_FILE with results obtained (multiple results per coord for multiple ENSTs hit)
	@result_set = @{ $coord_result_hash{$coord} };
	foreach $result (@result_set){
		chomp $result;
		@elts             = split /\t/, $result;
		$comparison_file = "";
		$hit_id = @elts[0];
		$ensg = @elts[1];
		$enst = @elts[2];
		$subst = @elts[3];
		$region = @elts[4];
		$original_seq = @elts[5];
		$modified_seq = @elts[6];
		$nu_change = $elts[7];
		$indel_start_CDS = $elts[8];
		$indel_stop_CDS = $elts[9];
		$indel_loc_pct = $elts[10];
		$indel_span = $elts[11];
		$visual_tx = $elts[12];
		$coords_change = $elts[13];
		$aa_change = $elts[14];
		$causes_NMD = $elts[15];
		$sth_db_supp_geneinfo->execute($enst);
        	@rows1 = $sth_db_supp_geneinfo->fetchrow_array();
			$gene_name = $rows1[3];
	        	$gene_desc = $rows1[4];
        		$ensfm = $rows1[5];
	        	$fam_desc = $rows1[6];
        		$gene_status = $rows1[7];
	        	$fam_size = $rows1[8];
	        	$kaks_mouse = $rows1[9];
		        $kaks_macaque = $rows1[10];
        		$mim_status = $rows1[11];

		unless($original_seq eq ""){
			$comparison_file = "$tmp\/$pid\_$enst\_$hit_id\_comparison.fasta";
		}
		print CLASSIFICATION_FILE "$coord\t$hit_id\t$ensg\t$enst\t$subst\t$coords_change\t$nu_change\t$aa_change\t$region\t$comparison_file\t$indel_span\t$indel_loc_pct\t$visual_tx\t$causes_NMD\t$small_indel_warning\t$repeat_change\t$gene_name\t$gene_desc\t$ensfm\t$fam_desc\t$gene_status\t$fam_size\t$kaks_mouse\t$kaks_macaque\t$mim_status\n";
		unlink ("$original_seq");
		unlink ("$modified_seq");
	}
}
close(CLASSIFICATION_FILE);

#create condensed Classification file (with only one transcript) in case user selects that option.
open( CLASSIFICATION_FILE, "$tmp/$pid.classification_file" )
  || die("Cannot open classification file");
open (CONDENSED_CLASSIFICATION_FILE, ">$tmp/$pid.condensed_classification_file");
my %seen_index;
while (<CLASSIFICATION_FILE>){
	chomp;
	my @elts = split /\t/, $_;
	my $coords = $elts[0];
	my $num_elts = num_elts_in_row($_);
	if ($seen_index{$coords} ne ""){
		my $prev_row = $seen_index{$coords};
		my $prev_num_elts = num_elts_in_row($prev_row);
		if ($num_elts > $prev_num_elts){
			$seen_index{$coords} = $_;
		}
	}
	else{
		$seen_index{$coords} = $_;
	}
}

foreach my $key (keys %seen_index){
	print CONDENSED_CLASSIFICATION_FILE $seen_index{$key},"\n";
}
close (CONDENSED_CLASSIFICATION_FILE);
close (CLASSIFICATION_FILE);

#read from combined file and print table
if ($all_transcripts eq "1"){
	open( CLASSIFICATION_FILE, "$tmp/$pid.condensed_classification_file" )|| die("Cannot open classification file");
}
else{
	open( CLASSIFICATION_FILE, "$tmp/$pid.classification_file" )|| die("Cannot open classification file");
}

$heading_tsv = get_output_heading_tsv($output_options);
$heading_html = get_output_heading_html($output_options);
open( OUTFILETABLE, ">$tmp/$pid\_predictions.html" );
open( OUTFILETSV,   ">$tmp/$pid\_predictions.tsv" );
print OUTFILETSV "$heading_tsv\n";
print OUTFILETABLE '<table border="1" cellspacing="0" cellpadding="4">';
print OUTFILETABLE '<thead>';
print OUTFILETABLE '<tr><th>';
$heading_html =~ s?\t?</th><th>?g;
print OUTFILETABLE "$heading_html</th></tr>\n";
my $warningflag = 0;

while (<CLASSIFICATION_FILE>) {
	chomp;
	@elts = split /\t/, $_;
	$coord =  @elts [0];
	$hit_id = @elts[1];
	$ensg = @elts[2];
	$enst = @elts[3];
	$subst = @elts[4];
	$coords_change = @elts[5];
	$nu_change = $elts[6];
	$aa_change = @elts[7];
	$region = @elts[8];
	$comparison_file = @elts[9];
	$indel_span = $elts[10];
	$indel_loc_pct = $elts[11];
	$visual_tx = $elts[12];
	$causes_NMD = $elts[13];
	$small_indel_warning = $elts[14]; 
	$repeat_change = $elts[15];
	$table_row_html = get_output_row_html($output_options,$_);
	$table_row_tsv = get_output_row_tsv($output_options,$_);
        print OUTFILETSV "$table_row_tsv\n";
	print OUTFILETABLE "<tr>\n";
	my @fields = split( '\t', $table_row_html );
		for $cell (@fields) {
			#if small indel error detected, sho warning in coords column
			if ($cell =~ /\d+,\d+/){
				$cell.="<br><a href=\"\/www\/chr_coords_example_indels.html#miscall\" target=\"_blank\"><font color=red>$small_indel_warning<\/font><\/a>";
				print OUTFILETABLE "<td>$cell</td>";
			}
			#link to original and modfied sequences file
			elsif ($cell =~ /comparison/i){
				print OUTFILETABLE "<td><a href=\"/sift-bin/catfile.csh?$cell+Alignment+PRE\" target=_blank> $enst change</td>";
			}
			#color indels red
			elsif ($cell =~ /([a-z]|[A-Z])+\*?->([a-z]|[A-Z])+\*?/){
				$cell =~ s/[a-z]+/<font color=red>$&<\/font>/g;
				print OUTFILETABLE "<td>$cell</td>";
			}
			#hover text for visual transcript
			elsif($cell =~ /\%/ || $cell =~ /\(\d+-\d+\)-\(\d+-\d+\)/ || $cell =~ /CDS|INTRON/){
				$hover_cell = "<a title=\'$visual_tx\' class=body_con>$cell</a>";
				print OUTFILETABLE "<td>$hover_cell</td>";
			}
			#color indel(nucleotide) red
			elsif ($cell =~ /[atgc]+-[ATGC]+-[atgc]+/){
                                $cell =~ s/[ATGC]+/<font color=red>$&<\/font>/g;
                                print OUTFILETABLE "<td>$cell</td>";
                        }
			elsif ($cell =~ /[ATGC]+-[atgc]+-[ATGC]+/){
                                $cell =~ s/[atgc]+/<font color=red>$&<\/font>/g;
                                print OUTFILETABLE "<td>$cell</td>";
                        }


			else{
				print OUTFILETABLE "<td>$cell</td>";
			}
		}
		print OUTFILETABLE "</tr>\n";
		$coord     = "";
                $ensg       = "";
                $enst  = "";
                $subst = "";
                $region     = "";
                $hit_id      = "";
                $comparison_file      = "";
		$visual_tx = "";
		$indel_loc_pct = "";
		$indel_span = "";
}
print OUTFILETABLE "</tr>\n</tbody>\n</table>\n<BR>";
close(FILE);
close(OUTFILETABLE);
close(OUTFILETSV);

sub get_output_heading_tsv{
        $oo = @_[0] ;
        @elts = split /,/, $oo;
	$heading_tsv = "Coordinates\tGene ID\tTranscript ID\tSubstitution Type\tRegion\tAmino acid position change\tIndel Location\tNucleotide change\tAmino acid change\tCauses NMD\tRepeat detected\tTranscript visulization";
        @options = ("Gene Name","Gene Desc","Protein Family ID","Protein Family Desc","Transcript Status","Protein Family Size","Ka/Ks (Mouse)","Ka/Ks (Macaque)","OMIM Disease");

        return $heading_tsv;
}

sub get_output_heading_html{
        $oo = @_[0] ;
        @elts = split /,/, $oo;
	$heading_html =
"Coordinates\tGene ID\tTranscript ID\tSubstitution Type\tRegion\t<a href=\"\/www\/chr_coords_example_indels.html#Amino_acid_position_change\" target=\"_blank\">Amino acid position change</a>\t<a href=\"\/www\/chr_coords_example_indels.html#Indel_location\" target=\"_blank\">Indel location</a>\t<a href=\"\/www\/chr_coords_example_indels.html#Nucleotide_change\" target=\"_blank\">Nucleotide change</a>\t<a href=\"\/www\/chr_coords_example_indels.html#Amino_acid_change\" target=\"_blank\">Amino acid change</a>\t<a href=\"\/www\/chr_coords_example_indels.html#Protein_sequence_change\" target=\"_blank\">Protein Sequence Change</a>\t<a href=\"\/www\/chr_coords_example_indels.html#NMD\" target=\"_blank\">Causes Nonsense Mediated Decay (NMD)</a>\t<a href=\"\/www\/chr_coords_example_indels.html#repeat\" target=\"_blank\">Repeat detected</a>";
	@options = ("Gene Name","Gene Desc","Protein Family ID","Protein Family Desc","Transcript Status","Protein Family Size","Ka/Ks (Mouse)","Ka/Ks (Macaque)","OMIM Disease");

        return $heading_html;
}


sub get_output_row_tsv {
        $oo = @_[0];
        $classification_line = @_[1];
	chomp $classification_line;
        @elts = split /\t/, $classification_line;
        $coord =  $elts [0];
        $hit_id = $elts[1];
        $ensg = $elts[2];
        $enst = $elts[3];
        $subst = $elts[4];
        $coords_change = $elts[5];
        $nu_change = $elts[6];
        $aa_change = $elts[7];
        $region = $elts[8];
        $comparison_file = $elts[9];
        $indel_span = $elts[10];
        $indel_loc_pct = $elts[11];
        $visual_tx = $elts[12];
        $causes_NMD = $elts[13];
        $small_indel_warning = $elts[14];
	$repeat_change = $elts[15];
	$table_row_tsv = "$coord $small_indel_warning\t$ensg\t$enst\t$subst\t$region\t$coords_change\t$indel_loc_pct\t$nu_change\t$aa_change\t$causes_NMD\t$repeat_change\t$visual_tx";
        return $table_row_tsv;

}

sub get_output_row_html {
        $oo = @_[0];
        $classification_line = @_[1];
        chomp $classification_line;
        @elts = split /\t/, $classification_line;
        $coord =  $elts [0];
        $hit_id = $elts[1];
        $ensg = $elts[2];
        $enst = $elts[3];
        $subst = $elts[4];
        $coords_change = $elts[5];
        $nu_change = $elts[6];
        $aa_change = $elts[7];
        $region = $elts[8];
        $comparison_file = $elts[9];
        $indel_span = $elts[10];
        $indel_loc_pct = $elts[11];
        $visual_tx = $elts[12];
        $causes_NMD = $elts[13];
        $small_indel_warning = $elts[14];
	$repeat_change = $elts[15];
        $table_row_html ="$coord\t$ensg\t$enst\t$subst\t$region\t$coords_change\t$indel_loc_pct\t$nu_change\t$aa_change\t$comparison_file\t$causes_NMD\t$repeat_change";
        return $table_row_html;

}


sub num_elts_in_row{
	my $row = @_[0];
	my $num_elts;
	chomp $row;
	my @elts = split /\t/, $row;
	foreach my $elt(@elts){
		if ($elt ne ""){
			$num_elts++;
		}
	}	
	return $num_elts;
}


#-------------------------------------------------------------------------
exit(0);

#-------------------------------------------------------------------------
#
# parameter: a string that is the html QUERY_STRING environment
#variable
# returns: an associative array of name/value pairs.  The name is the
#key.
sub parse_query {
	local ($query_string) = @_;
	local ( %ans, @q, $pair );

	#print $query_string;
	# break up into individual name/value lines
	@q = split( /&/, $query_string );

	foreach $pair (@q) {

		# break the name/value pairs up
		# use split rather than regular expressions because the value may
		# have
		#  newlines in it
		split( /=/, $pair, 2 );

		# change '+' to ' '
		$_[1] =~ s/\+/ /g;

		# change the escaped characters (has to be after the split on '&'
		# and '=')
		$_[1] =~ s/%(..)/pack("c",&hextodec("$1"))/eg;

		$ans{ $_[0] } = $_[1];
	}

	return %ans;
}

#-------------------------------------------------------------------------
# parameter: a hex representation of a number (doesn't need to be a
# string)
# returns: the decimal representation of the number
sub hextodec {
	unpack( "N", pack( "H8", substr( "0" x 8 . shift, -8 ) ) );
}

#-------------------------------------------------------------------------
# $names = &general_parse($ENV{CONTENT_TYPE}, $QUERY_STRING);
# parameters:   CONTENT_TYPE
#               QUERY_STRING
# returns: an associative array of name/value pairs.  The name is the
# key.

# WARNING:  Some of this routine is program-dependent!!!

# CONTENT_TYPE: application/x-www-form-urlencoded
# QUERY_STRING: key1=val1&key2=val2

# CONTENT_TYPE: multipart/form-data; boundary=<boundary>
# QUERY_STRING: <boundary>
#               Content-Disposition: form-data; name="key1"
#               <blank line>
#               val1
#               <boundary>
#               Content-Disposition: form-data; name="key2"
#               <blank line>
#               val2
#               <boundary>

sub general_parse {
	local ( $content_type, $query_string ) = @_;
	local ( %ans, @q, $pair, $loc, $boundary, $temp, $temp1 );

	if ( $content_type eq "application/x-www-form-urlencoded" ) {

		# break up into individual name/value lines
		@q = split( /&/, $query_string );

		foreach $pair (@q) {

			# break the name/value pairs up
			# use split rather than regular expressions because the value
			# may have
			#  newlines in it
			split( /=/, $pair, 2 );

			# change '+' to ' '
			$_[1] =~ s/\+/ /g;

			# change the escaped characters (must be after the split on '&'
			# and '=')
			$_[1] =~ s/%(..)/pack("c",&hextodec("$1"))/eg;

			$ans{ $_[0] } = $_[1];
		}    #end of foreach $pair

	}    #end of if ($content_type)
	else {
		$loc = index( $content_type, "boundary=" );
		if ( $loc > 0 ) {
			$temp = substr( $content_type, $loc + 9 );

		 #               Why is this necessary? (boundary= doesn't match actual)
			$boundary = "--" . $temp;

			# break up into individual name/value lines
			@q = split( /$boundary/, $query_string );

			foreach $pair (@q) {

				# break the name/value pairs up
				$loc = index( $pair, "name=" );
				$temp = substr( $pair, $loc + 5 );

				#         $loc = index($temp, "\n\n");
				$loc = index( $temp, "\n" );
				$temp1 = substr( $temp, $loc + 2 );

				#   Get rid of stuff after the name; including semicolon if any
				$loc_semi = index( $temp, ";" );
				$loc_eol  = index( $temp, "\n" );
				$loc      = $loc_eol;
				if ( $loc_semi > 0 && $loc_semi < $loc ) {
					$loc = $loc_semi;
				}
				if ( $loc > 0 ) { $temp = substr( $temp, 0, $loc ); }

				#               Get rid of quotes around the name
				$temp =~ s/\"//g;

				#               Still has a trailing whitespace character ...
				$temp =~ s/\s//g;

		  #               Substitute spaces with nothing
		  #               Need to strip leading/ending whitespace off of $temp1,
		  #               but be careful not to strip off internal CRs
		  #               MAC file lines end in just \r, no \n, so makelis won't
		  # find all
		  #               of the sequences; DOS file lines end in \r\n, UNIX in
		  #\n.
		  #               Change \r\n to \n and then \r to \n
#######PROGRAM -SPECIFIC!!!!!!!######################
		 #In my case, I want to keep the newlines in fields which have "file" or
		 # 'seq"
		 # and remove newlines everywhere else.
		 #if ( $temp =~ /file/ || $temp =~ /seq/ || $temp =~ /subst/ ) {
				$temp1 =~ s/\r\n/\n/g;
				$temp1 =~ s/\r/\n/g;

				#}

			 # for other variables that are NOT files or file-like, remove extra
			 #whitespace
			 #else { $temp1 =~ s/\s//g; }
				if ( $temp ne "" ) { $ans{$temp} = $temp1; }
			}    # end of foreach
		}    #end of if loc > 0
		else {
			print "Cannot parse\n";
			print "content_type=$content_type\n";
			print "query_string=$query_string\n";
		}
	}
	return %ans;

	#print "</PRE>";
}    # end of general_parse

# returns hash for a file, 2nd field is the key and the 3rd field
# is the value 4th field, is the delimiter
sub make_hash {
	my ($file) = @_;
	my %hash;
	open( HASH, $file ) || die "can't open $file";
	my $line;
	while ( $line = <HASH> ) {
		chomp($line);
		if ( exists( $hash{$line} ) ) {
			$hash{$line}++;
		}
		else {
			$hash{$line} = 1;
		}
	}
	close(HASH);
	return (%hash);
}

sub update_IP_logfile {
	my ( $queuefile, $IP_address ) = @_;

	$lockqueuefile = "$queuefile.lock";

	# lockfile will wait until it can lock the file
	`./lockfile $lockqueuefile`;

	# append the address and command to the queue file
	open( FILE, ">>$queuefile" );
	print FILE "$IP_address\n";
	close(FILE);

	chmod( 0664, $queuefile );

	# remove the lock file
	unlink($lockqueuefile);

}


