#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use Class::Struct;
#
# Nov 28, 2000 added option to email results
# 7/4/01 gap option is turned off permanently
#11/11/03 Added SIFT_queue stuff  JGH
# 3/8/04 Don't wait for search to finish, send them format.pl URL  JGH
#-------------------------------------------------------------------------
use DBI;
use Tie::IxHash;
#my $db =
#  DBI->connect( "dbi:SQLite:dbname=/opt/www/sift/db/1.5/SIFT_DB_1.5.sqlite",
#	"", "", { RaiseError => 1, AutoCommit => 1 } );
#$sth_aln_enst = $db->prepare("select address_alignment from enst_address where ens_transcript_id = ?");


$| = 1;
require 'SIFT_subroutines.pm';

#       Set file permissions to rw-rw----
system("umask 006");
my $bin             = "/opt/www/sift/htdocs/sift-bin";
my $tmp             = "/opt/www/sift/tmp";
my $datadir         = "/opt/www/sift/data";
$program_call       = $bin . "/SIFTING2_for_subst_only.csh";
#$return_address     = "sift\@fhcrc.org";
$snp_classifier     = "$bin/Classify_SNPs.pl";
$sequence_extractor = '/usr/local/projects/SIFT/snpclassifier/temp/bio.pl';
my $master_pid = $ARGV[0];
my $pid = $ARGV[1];
my $all_chr_file = $ARGV[2];
my $organism = $ARGV[3];
my $Human_db_dir    = "/opt/www/sift/db/$organism";
my $output_options = $ARGV[7];
my $address = $ARGV[8];
my $seq_identity_filter = $ARGV[4];

check_ip_counts();

my $last_partition = 0;
if ($ARGV[5] !~ /NOT_LAST/i){
	$last_partition = 1;
}
my $COORD_SYSTEM = $ARGV[6];

struct VariantInfo => {
        rsid => '$',
        orn_db => '$',
        enst => '$',
        ensp => '$',
        region => '$',
        snp_type => '$',
        codon1 => '$',
        codon2 => '$',
        nt1 => '$',
        nt2 => '$',
        alleles_db => '$',
        AA1 => '$',
        AA2 => '$',
        AAPOS => '$',
        score => '$',
        median => '$',
        seqs_rep => '$',
        freq_av => '$',
        freq_ceu => '$',

# gene info
        ensg => '$',
        enst => '$',
        ensp => '$',
        gene_name => '$',
        gene_desc => '$',
        ensfm => '$',
        fam_desc => '$',
        gene_status => '$',
        fam_size => '$',
        kaks_mouse => '$',
        kaks_macaque => '$',
        mim_status => '$',
};


update_status($master_pid,$pid,"Running");
$bin_file = "$Human_db_dir/bins.list";

#Creat binned SNP files for querying human chr database- the script also creates snp_chr_map_file to map bin file to chr database and table
system("$bin/map_choords_to_bin.pl $bin_file $all_chr_file $pid");

open (SNP_CHR_MAP_FILE, "$tmp/$pid\_snp_chr_map_file.txt");
while (<SNP_CHR_MAP_FILE>){
        chomp;
        @elts = split /\t/, $_;
        $chr = @elts [0];
        $bin = @elts[1];
        $snp_chr_bin_file = @elts[2];
        $db_chr_str = @elts[3];
        $table_chr = @elts[4];
        $db_chr =DBI->connect( "dbi:SQLite:dbname=$Human_db_dir/Human_CHR$chr.sqlite","", "", { RaiseError => 1, AutoCommit => 1 } );
	$db_chr->do('PRAGMA synchronous=1');
	$db_chr->do('PRAGMA cache_size=4000');
	#$sth_db_chr = $db_chr->prepare("select * from $table_chr where chr = \'chr$chr\' AND COORD1 = ? AND COORD2 = ? ");
        open (BIN_FILE, "$tmp/$snp_chr_bin_file") || die ("Cannot open snp_chr_bin_file");
        $coord2_str = "";
	$coord1_str = "";
	while (<BIN_FILE>){
		$count ++;
                chomp;
                @elts2 = split '\t',$_;
                $coord1 = @elts2[1];
                $coord2 = @elts2[2];
                $orn_usr = @elts2[3];
                $alleles = @elts2[4];
		$usr_comment = @elts2[5];
                $alleles =~ /(.+?)\/(.+?)/;
                $nt2a = $1;
                $nt2b = $2;
	# PCN if given in - strand, then reverse complement
                if ($orn_usr == -1) {
                        $nt2a =~ tr/ACGTacgt/TGCAtgca/;
                        $nt2b =~ tr/ACGTacgt/TGCAtgca/;
                        $orn_usr = 1;
                }

		$coord2_str.="$coord2,";
		$coord1_str.="$coord1,";
		### comment this block if using in operator see below.
		#$sth_db_chr->execute($coord1,$coord2);
		#while (@rows = $sth_db_chr->fetchrow_array()){
        	#	push @query_result, join("\t",@rows);
		
	        #}
		#system("echo $snp_chr_bin_file $count >> $tmp/$pid.query_status.txt");		
		###
		
	}
	close (BIN_FILE);
	chop $coord1_str;
	chop $coord2_str;
	### uncomment following block if using in operator - see above
	$sth_db_chr = $db_chr->prepare("select * from $table_chr where chr = \'chr$chr\' AND COORD2 in ($coord2_str) AND COORD1 in ($coord1_str)");
	#system ("echo \"select * from $table_chr where chr = \'chr$chr\' AND COORD2 in \($coord2_str\) AND COORD1 in \($coord1_str\)\" > $tmp/remove");
	$sth_db_chr->execute();
        while (@rows = $sth_db_chr->fetchrow_array()){
		push @query_result, join("\t",@rows);
		
	}
	system("echo $snp_chr_bin_file $count $time>> $tmp/$pid.query_status.txt"); 
	###
}
close (SNP_CHR_MAP_FILE);
#index all results obtained into index_query_result. hash of arrays.
$db_supp = DBI->connect( "dbi:SQLite:dbname=$Human_db_dir/Human_Supp.sqlite","", "", { RaiseError => 1, AutoCommit => 1 } );
$db_supp->do('PRAGMA synchronous=1');
$sth_db_supp_geneinfo = $db_supp->prepare("select * from GENE_INFO where ENST = ?");
$sth_db_supp_allelefreq = $db_supp->prepare("select * from ALLELE_FREQ where RSID = ?"); 


## PCN - Jan 25, 2010 add transcript info
$db_tx = DBI->connect( "dbi:SQLite:dbname=$Human_db_dir/TRANSCRIPT_DB.sqlite",
         "",  "", { RaiseError => 1, AutoCommit => 1} );
$db_tx->do('PRAGMA synchronous=1');
$sth_db_strand_tx = $db_tx->prepare ("select orn from ENST_REGION where ENST = ?");

%tx_enst_hash;

foreach $row (@query_result){
	chomp $row;
	@elts = split /\t/, $row;	
	$chr = @elts[0];
	$coord1 = @elts[1];
	$coord2 = @elts[2];
	$orn = @elts[3];
	$nt1 = @elts[10];
	$nt2 = @elts[11];
	$enst = @elts[6];
	$rsid = @elts[4];
	if($rsid =~ /(rs\d+)\:.+?/){
		$rsid_formatted = $1;
	}
	else{
		$rsid_formatted = "";
	}
	$ensp = @elts[7];
	
	if ($ensp !~ /\d/){
		$ensp_exists = 0;
	}
	else{
		$ensp_exists = 1;
	}
	$CDS = @elts[20];
	$AA1_VALID = @elts[21];

## PCN 01-25-10 add strand orientation 
	$sth_db_strand_tx->execute ($enst);
	@rows_strand = $sth_db_strand_tx->fetchrow_array();
	$orn_tx = @rows_strand[0];
	$tx_enst_hash{$enst} = $orn_tx;

	$sth_db_supp_geneinfo->execute($enst);
        @rows1 = $sth_db_supp_geneinfo->fetchrow_array();
	$ensg = @rows1[2];
        $gene_name = @rows1[3];
	$gene_desc = @rows1[4];
	$ensfm = @rows1[5];
	$fam_desc = @rows1[6];
	$gene_status = @rows1[7];
	$fam_size = @rows1[8];
	$kaks_mouse = @rows1[9];
	$kaks_macaque = @rows1[10];
	$mim_status = @rows1[11];


	
	if ($rsid_formatted =~ /rs/){
		$sth_db_supp_allelefreq->execute($rsid_formatted);
	        @rows2 = $sth_db_supp_allelefreq->fetchrow();
        	$freq_av1 = @rows2[1];
		if ($freq_av1 =~ /(0\.\d+)/){
			$fav1 = $freq_av1;
		}
        	$freq_av2 = @rows2[2];
		if ($freq_av2 =~ /(0\.\d+)/){
                        $fav2 = $freq_av2;
                }

        	$freq_ceu1 = @rows2[3];
		if ($freq_ceu1 =~ /(0\.\d+)/){
                        $fceu1 = $freq_ceu1;
                }

        	$freq_ceu2 = @rows2[4];
		if ($freq_ceu2 =~ /(0\.\d+)/){
                        $fceu2 = $freq_ceu2;
                }
		$freq_av = "$fav1:$fav2";
		$freq_ceu = "$fceu1:$fceu2";

	}
	else{
		$freq_av1 = "";
                $freq_av2 = "";
                $freq_ceu1 = "";
                $freq_ceu2 = "";
		$freq_av="";
		$freq_ceu="";

	}
	$row2 = "$row\t$ensg\t$gene_name\t$gene_desc\t$ensfm\t$fam_desc\t$gene_status\t$fam_size\t$kaks_mouse\t$kaks_macaque\t$mim_status\t$freq_av\t$freq_ceu";
	#system("echo $row2 >> $tmp/remove1");
	$key = "$chr\t$coord1\t$coord2\t$CDS\t$ensp_exists";
	push @ {$index_query_result{$key}}, $row2;	
}
$rsid = "";$enst = "";$ensp="";$coord1="";$coord2 = "";$chr = "";$nt1 = "";$nt2="";

#parse snp_chr_map_file to fetch data from previously created $index_query_result.
open (SNP_CHR_MAP_FILE, "$tmp/$pid\_snp_chr_map_file.txt");
open (BAD_ENST_FILE, ">$tmp/$pid.bad_enst_file.txt");
while (<SNP_CHR_MAP_FILE>){
	chomp;
	@elts = split /\t/, $_;
	$chr = @elts [0];	
	$bin = @elts[1];
	$snp_chr_bin_file = @elts[2];
	$db_chr_str = @elts[3];
	$table_chr = @elts[4];
	open (BIN_FILE, "$tmp/$snp_chr_bin_file") || die ("Cannot open snp_chr_bin_file");
	while (<BIN_FILE>){
		chomp;
		@elts2 = split '\t',$_;
		$coord1 = @elts2[1];
		$coord2 = @elts2[2];
		$orn_usr = @elts2[3];
## PCN if orn_usr is -1, reverse complement alleles, because in database
# it's referenced to +1 orientation
		$alleles = @elts2[4];
		$alleles =~ /(.+?)\/(.+?)/;
		$usr_comment = @elts2[5];
		$nt2a = $1;
		$nt2b = $2;
		if ($orn_usr == -1) {
                        $nt2a =~ tr/ACGTacgt/TGCAtgca/;
                        $nt2b =~ tr/ACGTacgt/TGCAtgca/;
                        $orn_usr = 1;
                }


		$key_ensp_CDS = "chr$chr\t$coord1\t$coord2\t1\t1"; 			#CDS is true, ensp is true,
		@rows_ensp_CDS = @{ $index_query_result{$key_ensp_CDS} };
		$key_no_ensp_CDS = "chr$chr\t$coord1\t$coord2\t1\t0";                   #CDS is true, ensp is false,
                @rows_no_ensp_CDS = @{ $index_query_result{$key_no_ensp_CDS} };
#Prateek had below commented out, only looking at rows_ensp_CDS
# Pauline allowed rows_no_ensp_CDS to be seen
		@rows_CDS = (@rows_ensp_CDS,@rows_no_ensp_CDS);
		#replaced above with below - not using no-ensp rows any more
		#will also not output anything for not found in db
		#@rows_CDS = @rows_ensp_CDS;
# Pauline commented out above, print out rows with non-ensp id's
		$key_ensp_no_CDS = "chr$chr\t$coord1\t$coord2\t0\t1";			#CDS is false , ensp is true
                @rows_ensp_no_CDS = @{ $index_query_result{$key_ensp_no_CDS} };
		$key_no_ensp_no_CDS = "chr$chr\t$coord1\t$coord2\t0\t0";    		#CDS is false , ensp is false
                @rows_no_ensp_no_CDS = @{ $index_query_result{$key_no_ensp_no_CDS} };
		@rows_no_CDS = (@rows_ensp_no_CDS,@rows_no_ensp_no_CDS); 		# final combined - this and $rows_CDS will now be used.
		
			
		if (scalar @rows_CDS == 0){						# no CDS found for this coords
			if (@rows_no_CDS == 0){ 					#coord not in db
		#		next;
				#added above for not outputting anything not in db
				my $m = make_new_variant ();
                                $m->region ("NON-GENIC");

                                push @coord_arr_no_CDS, "$chr\t$coord1\t$coord2\t$orn_usr\t$alleles\t$usr_comment";
                                $coord_result_no_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_usr\t$alleles\t$usr_comment"} = $m;
			}
			else{			# rows outside of CDS found (promotor, Downstream or UTR)
				my $v = make_variant_with_elts (@rows_no_CDS[0]);
# PCN 012610 changed to print orientation of user, which should be +1 strand
# because alleles will be printed out with respect to +1 strand, so it makes more sense
                                push @coord_arr_no_CDS, "$chr\t$coord1\t$coord2\t$orn_usr\t$alleles\t$usr_comment";
                                $coord_result_no_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_usr\t$alleles\t$usr_comment"} = $v;

			}
		}
		else{			#CDS found - iterate through rows- choose one with right nt2 - if nt2 is ref then modify subst
			$selected_row = -1;
			$ref_allele_row = -1;
			$AA1_VALID_found = 0;
			for ($i = 0; $i < scalar @rows_CDS; $i++){
				$row = @rows_CDS[$i];
				chomp $row;
				#print "$row<BR>";
				#system("echo $row >> $tmp/remove");
				@elts_CDS = split /\t/, $row;
				$enst = @elts_CDS[6];
                                $ensg = @elts_CDS[26];
				$orn_db = @elts_CDS[3];
				$nt1 = @elts_CDS[10];
				$nt2 = @elts_CDS[11];
				$AA1_VALID = @elts_CDS[21];
				if ($AA1_VALID == 1){
					$AA1_VALID_found = 1;
					$ref_allele_row = $i;  #save last row so if seclected row==-1, this row is used assuming  ref allele analysis reqd.
#					if ($nt2 eq $nt2b || $nt2 eq $nt2a ){
# Pauline commented out Prateek's code 01-22-10
# homozygous variants not processed, heterozygous only
# don't care which orientation, reference doesn't have to be first
# add check for orientation (orn)
#print "nucleotides $nt2 $nt1  and $nt2a $nt2b\n";
#print "transcript orientation is $tx_enst_hash{$enst} $enst $ensg orn db $orn_db\n";
					if ( (($nt2 eq $nt2b && $nt1 eq $nt2a)
							 ||
					     ($nt2 eq $nt2a && $nt1 eq $nt2b)) 
						&& $orn_db eq 
						 $tx_enst_hash{$enst}) { 
#print "orn_tx $orn_tx orn_db $orn_db\n";
                                               	$selected_row = $i;
                                               	last;
                                	}

				}
				else{				#this CDS has bad enst. go to next 
					 $ref_allele_row = $i;
					#save this coord in  BAD_ENST_FILE for later
					next;
				}
			}
			if($selected_row != -1){			#valid CDS
                                $row = @rows_CDS[$selected_row];        
                                chomp $row;
				my $v = make_variant_with_elts ($row);
				my $alleles_db = $v->alleles_db;

                                push @coord_arr_CDS, "$chr\t$coord1\t$coord2\t$orn_usr\t$alleles_db\t$usr_comment";
                                $coord_result_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_usr\t$alleles_db\t$usr_comment"} = $v;
                        }
	
			elsif($selected_row == -1 && $AA1_VALID_found ==1){		#this means user input requires pred on ref allele
				if ($nt2a eq $nt2b) {
#this means user input requires pred on ref allele
#print  "chr $chr $coord1 $coord2 in ab98ua\n";
                                $row = @rows_CDS[$ref_allele_row];
                                chomp $row;
                                my $v = make_variant_with_elts ($row);
# reset everything to synonmyous
                                $v->alleles_db ($v->nt1 . "\/" . $v->nt1);
                                $v->snp_type ("Synonymous");
                                $v->codon2 ($v->codon1);
                                $v->AA2 ($v->AA1);
                                $v->score ("N/A");
                                $v->median ("N/A");
                                $v->seqs_rep ("N/A");
                                my $s = $v->nt1 . "\/" . $v->nt1;
                                push @coord_arr_CDS, "$chr\t$coord1\t$coord2\t$orn_usr\t$s\t$usr_comment";
                                $coord_result_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_usr\t$s\t$usr_comment"} = $v;
                                } else {
# not a reference allele, don't know how to process
					 $row = @rows_CDS[$ref_allele_row];
                                        chomp $row;
                                 my $v = make_variant_with_elts ($row);
# reset everything to synonmyous
                                $v->snp_type ("Unknown");
                                $v->codon1 ("");
                                $v->codon2 ("");
                                $v->AA1 ("NA");
                                $v->AA2 ("");
                                $v->score ("");
                                $v->median ("N/A");
                                $v->seqs_rep ("N/A");
                                 $v->alleles_db ($nt1 . "\/" . $nt2);
                                my $alleles_db = $v->alleles_db;
                                push @coord_arr_CDS, "$chr\t$coord1\t$coord2\t$orn_usr\t$alleles_db\t$usr_comment";
                                $coord_result_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_usr\t$alleles_db\t$usr_comment"} = $v;


                                        print BAD_ENST_FILE "$chr\t$coord1\t$coord2\t$orn_usr\t$nt1\/$nt2\t$usr_comment\n";
                                next;

                                }

			}
			elsif($AA1_VALID_found == 0){		#All ensts found bad. save this coord in a different file.
				 $row = @rows_CDS[$ref_allele_row];
                                        chomp $row;
                                 my $v = make_variant_with_elts ($row);
# reset everything to synonmyous
				my $error_message = "<A HREF=\"http://sift.jcvi.org/www/SIFT_help.html#SIFT_OUTPUT_GENOME_Error\" TARGET=_BLANK>Gene Annotation Error</A>";
                                $v->snp_type ("Unknown");
                                $v->codon1 ("");
                                $v->codon2 ("");
                                $v->AA1 ("NA");
                                $v->AA2 ("");
				$v->ensp ($error_message);
                                $v->score ("");
                                $v->median ("N/A");
                                $v->seqs_rep ("N/A");
                                 $v->alleles_db ($nt1 . "\/" . $nt2);
                                my $alleles_db = $v->alleles_db;
                                push @coord_arr_CDS, "$chr\t$coord1\t$coord2\t$orn_usr\t$alleles_db\t$usr_comment";
                                $coord_result_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_usr\t$alleles_db\t$usr_comment"} = $v;

				print BAD_ENST_FILE "$chr\t$coord1\t$coord2\t$orn_usr\t$nt1\/$nt2\t$usr_comment\n";
				next;
			}
		}
	}
	close (BIN_FILE)
}
close (SNP_CHR_MAP_FILE);
close (BAD_ENST_FILE);


#create enst_pred_list for coords in CDS .
open( ENSTFILE, ">>$tmp/$pid.enstfile" );
foreach $coord (@coord_arr_CDS) {
	my $variant_result    = $coord_result_CDS_hash{$coord};
        my $line = return_variant_line ($variant_result);
#print "line in here ENST is $line for coord $coord\n";
        my $alleles = (split /\t/, $coord)[4];
#print "coord after is $coord alleles is $alleles\n";
        my $key = $variant_result->enst . "\t" . $variant_result->AAPOS;

# Pauline commented out, so another prediction printed out even ifvariant
# is at same location

#	if (exists ($hash_enst_aapos_seen{$key})){
#		next;	
#	}	
#	else{
		$hash_enst_aapos_seen{$key} = "$alleles";
		$enst_line          = "$coord\t$line";
		push( @enst_pred_list, $enst_line );
		print ENSTFILE "$enst_line\n";	
#	}
}
close(ENSTFILE);

#create enst_subst_list for coords not in CDS - these will not be predicted by SIFT.
open( ENSTFILE, ">>$tmp/$pid.enstfile" );
foreach $coord (@coord_arr_no_CDS) {
        $variant_result    = $coord_result_no_CDS_hash{$coord};
	my $line = return_variant_line ($variant_result);

        $enst_line          = "$coord\t$line";
        push( @enst_pred_list, $enst_line );
        print ENSTFILE "$enst_line\n";
}
close(ENSTFILE);

#read from combined file and print table
open( ENST_PRED_FILE, "$tmp/$pid.enstfile" )
  || die("Cannot open predictions file");

$heading_tsv = get_output_heading($output_options);
$heading_html = $heading_tsv;
open( OUTFILETABLE, ">$tmp/$pid\_predictions.html" );
open( OUTFILETSV,   ">$tmp/$pid\_predictions.tsv" );
print OUTFILETSV "$heading_tsv\n";
print OUTFILETABLE '<table border="1" cellspacing="0" cellpadding="4">';
print OUTFILETABLE '<thead>';
print OUTFILETABLE '<tr><th>';
$heading_html =~ s?\t?</th><th>?g;
print OUTFILETABLE "$heading_html</th></tr>\n";
my $warningflag = 0;

while (<ENST_PRED_FILE>) {
	chomp;
	$table_row = get_output_row($output_options, $_);
	print OUTFILETABLE "<tr>\n";
	if ($table_row =~ /Synonymous/ || $table_row =~ /\w\d+\*/){
		$is_synonymous = 1;
	}
	else{
		$is_synonymous = 0;
	}
	my @fields = split( '\t', $table_row );
	$num_cells = scalar @fields;
	$count = 0;
	for $cell (@fields) {
		chomp $cell;
		$count ++;
		$cell_tsv = $cell;
		$cell_html = $cell;
		if ( $cell =~ /DAMAGING/ || $cell =~ /Warning/ ) {
			$cell_html =  "<font color=red>$cell</font>";
			$cell_tsv = $cell;
		}
		elsif ( $cell =~ /NOT PREDICTED/i || $cell =~ /not found/i ) {
			$cell_html = "<font color=blue>$cell</font>";
			$cell_tsv = $cell;
		}
		elsif ($cell =~ /(rs\d+):(.+?)/i){
			$dbsnp_id = $1;$mutation = $2;
			$cell_html =  "<A HREF=\"http:\/\/www.ncbi.nlm.nih.gov\/sites\/entrez?db=snp&cmd=search&term=\+$dbsnp_id\" target=\"_blank\">$dbsnp_id:$mutation</A>";
			$cell_tsv = $cell;
		}
		if($count >= 9 && $count <=12 && $is_synonymous == 1){
			$cell_html =  "N\/A";
			$cell_tsv = "N\/A";
			
		}
		print OUTFILETABLE "<td>$cell_html</td>";
		if ($count == $num_cells){
			print OUTFILETSV "$cell_tsv";
		}
		else{
			print OUTFILETSV "$cell_tsv\t";
		}

	}
	print OUTFILETABLE "</tr>\n";
	print OUTFILETSV "\n";

}

# print errors so people know what to look up
#open (BAD_ENST_FILE, "$tmp/$pid.bad_enst_file.txt");
#my $line;
my $error_message_invalid = "Error! Ensembl gene annotation did not match genome sequence. Please try ncbi37, or annotate by hand.";
#while ($line = <BAD_ENST_FILE>) {
#        chomp ($line);
#        my @fields = split (/\t/, $line);
#        if ($COORD_SYSTEM eq "SPACE"){
#                # join all coordinates, do nothing
#        } else {
                # splice out beg coordinate
#                splice (@fields,1,1);
#        }
#        my $variant = join (",",@fields);
#        $variant =~ s/\,$//;
#        print OUTFILETSV $variant . "\t" . $error_message_invalid . "\n";
#        print OUTFILETABLE "<tr><td>$variant</td><td colspan=\"11\">$error_message_invalid</td></tr>\n";
#}
#close (BAD_ENST_FILE);


print OUTFILETABLE "</tr>\n</tbody>\n</table>\n<BR>";
if ( $warningflag == 1 ) {

	print OUTFILETABLE
"<font color=red>* Low confidence means that the protein alignment does not have enough sequence diversity. Because the position artifically appears to be conserved, an amino acid may incorrectly predicted to be damaging.</font><BR><BR>";
}
#close(FILE);
close(OUTFILETABLE);
close(OUTFILETSV);

update_status($master_pid,$pid,"Complete");
system("rm -f $tmp/$pid.*.siftresults.predictions");

#if final batch - combine all prediction files and update status and stats
open(TSV_COMBINED,">$tmp/$master_pid\_predictions.tsv")|| die ("Cannot open prediction tsv files for combining");
open(HTML_COMBINED,">$tmp/$master_pid\_predictions.html")||die("Cannot open prediction tsv files for combining");
print TSV_COMBINED "$heading_tsv\n";
print HTML_COMBINED '<table border="1" cellspacing="0" cellpadding="4">';
print HTML_COMBINED '<thead>';
print HTML_COMBINED '<tr><th>';

print HTML_COMBINED "$heading_html</th></tr>\n";


if ($last_partition == 1){
	for ($i = $master_pid + 1; $i <= $pid; $i ++){
		open(TSV,"$tmp/$i\_predictions.tsv")|| die ("Cannot open prediction tsv files for combining");
		open(HTML,"$tmp/$i\_predictions.html")|| die ("Cannot open prediction tsv files for combining");
		while (<TSV>){
			#combined TSV file
			chomp;
			if ($_ =~ /Coordinates/i){
				next;
			}
			else{
				print TSV_COMBINED "$_\n";
			}
			$total_num++;		
			if ($_ =~ /CDS/){
				$CDS_num ++;
			}
			if ($_ =~ /TOLERATED/){
                                $tolerated_num ++;
                        }
			if ($_ =~ /DAMAGING/){
                                $damaging_num ++;
                        }
			if ($_ =~ /Nonsynonymous/){
                                $nonsynonymous_num ++;
                        }
			if ($_ =~ /Synonymous/){
                                $synonymous_num ++;
                        }

			if ($_ =~ /novel/){
                                $novel_num ++;
                        }

			#combined html file
			@fields = split( '\t', $_ );
			if ($_ =~ /Synonymous/ || $_ =~ /\w\d+\*/){
				$is_synonymous = 1;
			}
			else{
				$is_synonymous = 0;
			}
			print HTML_COMBINED "<tr>\n";
			$count = 0;
                	for $cell (@fields) {
				$count ++;
                        	if ( $cell =~ /DAMAGING/ || $cell =~ /Warning/ ) {
                                	$cell =  "<font color=red>$cell</font>";
                        	}
                        	elsif ( $cell =~ /NOT PREDICTED/i || $cell =~ /not found/i ) {
                                	$cell  = "<font color=blue>$cell</font>";
                        	}
				elsif ($cell =~ /(rs\d+):(.+?)/i){
                                	$dbsnp_id = $1;$mutation = $2;
                                	$cell = "<A HREF=\"http:\/\/www.ncbi.nlm.nih.gov\/sites\/entrez?db=snp&cmd=search&term=\+$dbsnp_id\" target=\"_blank\">$dbsnp_id:$mutation</A>";
                        	}
				if($count >= 9 && $count <=12 && $is_synonymous == 1){
		                        $cell =  "N\/A";
                		}
#				if ($cell =~ /Error/) {
#					print HTML_COMBINED "<td colspan=\"11\">$cell</td>";
#				} else {
                               	print HTML_COMBINED "<td>$cell</td>";
#				}			
                	}
                	print HTML_COMBINED "</tr>\n";
		}
                print HTML_COMBINED "</tr>\n";
	}
	print HTML_COMBINED "</tr>\n</tbody>\n</table>\n<BR>";
	update_status($master_pid,$master_pid,"Complete");
	update_stats($master_pid,$total_num,$CDS_num,$tolerated_num,$damaging_num,$nonsynonymous_num,$synonymous_num,$novel_num);
}
close (TSV_COMBINED);
close (HTML_COMBINED);
#print the table - we are not doing this now. instead the user can download or view partitioned results via status table
#open( OUTFILETABLE, "$tmp/$combined_siftresults_file.predictions.table.html" )
#  || die("Cannot open predictions table");
#while (<OUTFILETABLE>) {
#	print;
#}
#print
#"Click <A HREF=\"\/tmp\/$combined_siftresults_file.predictions.tsv\">here</A> to download the following table in tab separated format. You can open it using excel with delimiter set as TAB";
#print "<BR>";

#print "<BR>";
#print
#"<i>If you received a warning that the sequences were not diverse enough, you can have SIFT choose more diverse sequences <A HREF=\"/www/SIFT_seq_submit2.html\">here.</A></i><BR><BR>";

#email the results
if ( $address ne "" && $last_partition == 1) {
	open( MESSAGE, ">$tmp/$master_pid.email_message.txt" );
	print MESSAGE
"Dear User\n\nThank you for using SIFT.\n\nPlease find the results of your recent query attached with this message.\nRemember this job id \"$master_pid\" for any future correspondance.\nDo not hesitate to contact us if you have any questions about SIFT.\n\nThanks\nSIFT Team\nJ Craig Venter Institute (West Coast Campus)\n10355 Science Center Drive\nSan Diego, CA 92121\nUSA";
	close(MESSAGE);
	system(
"mutt -F /opt/www/sift/htdocs/.muttrc -a $tmp/$master_pid\_predictions.tsv -s \"SIFT Results for Job ID $master_pid\" $address <$tmp/$master_pid.email_message.txt"
	);
}
sub update_stats{
	$tmp = "/opt/www/sift/tmp";
	$master_pid = @_[0];
	$total_num = @_[1];
	$CDS_num = @_[2];
	$tolerated_num = @_[3];
	$damaging_num = @_[4];
	$nonsynonymous_num = @_[5];
	$synonymous_num = @_[6];
	$novel_num = @_[7];
	$predicted_num = $tolerated_num + $damaging_num;
	$CDS_pct = int($CDS_num * 100 / $total_num);
	$predicted_pct = int($predicted_num * 100 / $CDS_num) ;
	$tolerated_pct = int($tolerated_num *100 / $predicted_num);
	$damaging_pct = 100 - $tolerated_pct;
	$nonsynonymous_pct = int($nonsynonymous_num*100/$CDS_num);
	$synonymous_pct = 100 - $nonsynonymous_pct;
	$novel_pct = int($novel_num * 100 / $total_num);
	open (OUTPAGE,"$tmp/$master_pid.outpage.html") || die ("Cannot open outpage for updating status");
	open (OUTPAGE_MODIFIED,">$tmp/$master_pid.outpage.modified.html") || die ("Cannot open swap outpage file");
	while (<OUTPAGE>){
		chomp;
		if ($_ =~ /Number of input \(non-intronic\) variants/){
                        print OUTPAGE_MODIFIED "Number of input (non-intronic) variants: $total_num <BR>\n";
                }
		elsif ($_ =~ /Coding variants/ && $_ !~ /predicted/){
			print OUTPAGE_MODIFIED "Coding variants: $CDS_pct% ($CDS_num out of $total_num) <BR>\n";
		}
		elsif($_ =~ /Coding variants predicted/){
			print OUTPAGE_MODIFIED "Coding variants predicted: $predicted_pct% ($predicted_num out of $CDS_num) <BR>\n";
		}
		elsif($_ =~ /Tolerated/){
                        print OUTPAGE_MODIFIED "Tolerated: $tolerated_pct% ($tolerated_num out of $predicted_num) <BR>\n";
                }
		elsif($_ =~ /Damaging/){
                        print OUTPAGE_MODIFIED "Damaging: $damaging_pct% ($damaging_num out of $predicted_num) <BR>\n";
                }
		elsif($_ =~ /Nonsynonymous/){
                        print OUTPAGE_MODIFIED "Nonsynonymous: $nonsynonymous_pct% ($nonsynonymous_num out of $CDS_num) <BR>\n";
                }
		elsif($_ =~ /Synonymous/){
                        print OUTPAGE_MODIFIED "Synonymous: $synonymous_pct% ($synonymous_num out of $CDS_num) <BR>\n";
                }
		elsif($_ =~ /Novel/){
                        print OUTPAGE_MODIFIED "Novel: $novel_pct% ($novel_num out of $total_num) <BR>\n";
                }
		else{
			print OUTPAGE_MODIFIED "$_\n";
		}

	}
	close (OUTPAGE_MODIFIED);
        close (OUTPAGE);
        system ("chmod 755 $tmp/$master_pid.outpage.modified.html");
        rename "$tmp/$master_pid.outpage.modified.html", "$tmp/$master_pid.outpage.html"

}

sub update_status{
	$tmp = "/opt/www/sift/tmp";
	$master_pid = @_[0];
	$slave_pid = @_[1];
	$status = @_[2];
	#print "$master_pid $slave_pid $status<BR>";
	open (OUTPAGE,"$tmp/$master_pid.outpage.html") || die ("Cannot open outpage for updating status");
	open (OUTPAGE_MODIFIED,">$tmp/$master_pid.outpage.modified.html") || die ("Cannot open swap outpage file");
	while (<OUTPAGE>){
        	chomp;
        	if ($_ =~ /.+?(Partitioned set .+?)\<\/td\>\<td\>(.+?)\<.+?\>$slave_pid\<.+?/){
                	$partition = $1;
                	$job_size = $2;

	                if ($status =~ /running/i){
	                        print OUTPAGE_MODIFIED "<tr><td>$partition</td><td>$job_size</td><td>$slave_pid</td><td><font color=blue>Running..</font></td><td>Not available</td></td><td>Not available</td></tr>\n";
	                }
	                elsif($status =~ /complete/i){
	                        print OUTPAGE_MODIFIED "<tr><td>$partition</td><td>$job_size</td><td>$slave_pid</td><td><font color=green>Complete</font></td><td><a href=/sift-bin/catfile.csh?$tmp/$slave_pid\_predictions.html target=\"_blank\">$slave_pid table</a></td><td><A HREF=/tmp/$slave_pid\_predictions.tsv>$slave_pid results</A></td></tr>\n";
	                }
	        }
		elsif($_ =~ /.+?Complete set\<\/td\>\<td\>(.+?)\<.+?\>$slave_pid\<.+?/){
			$partition = "Complete set";
                        $job_size = $1;

                        if ($status =~ /running/i){
                                print OUTPAGE_MODIFIED "<tr BGCOLOR=\"\#DCDCDC\"><td>$partition</td><td>$job_size</td><td>$slave_pid</td><td><font color=blue>Running..</font></td><td>Not available</td></td><td>Not available</td></tr></table>\n";
                        }
                        elsif($status =~ /complete/i){
                                print OUTPAGE_MODIFIED "<tr BGCOLOR=\"\#DCDCDC\"><td>$partition</td><td>$job_size</td><td>$slave_pid</td><td><font color=green>Complete</font></td><td><a href=/sift-bin/catfile.csh?$tmp/$slave_pid\_predictions.html target=\"_blank\">$slave_pid table</a></td><td><A HREF=/tmp/$slave_pid\_predictions.tsv>$slave_pid results</A></td></tr></table>\n";
                        }
                }

	        else{
	                print OUTPAGE_MODIFIED "$_\n";
	        }
	}
	close (OUTPAGE_MODIFIED);
	close (OUTPAGE);
	system ("chmod 755 $tmp/$master_pid.outpage.modified.html");
	rename "$tmp/$master_pid.outpage.modified.html", "$tmp/$master_pid.outpage.html"
}
#-------------------------------------------------------------------------
exit(0);


sub get_output_heading{
        $oo = @_[0] ;
        @elts = split /,/, $oo;
        $heading =
"Coordinates\tCodons\tTranscript ID\tProtein ID\tSubstitution\tRegion\tdbSNP ID\tSNP Type\tPrediction\tScore\tMedian Info\t\# Seqs at position";
        @options = ("Gene ID","Gene Name","Gene Desc","Protein Family ID","Protein Family Desc","Transcript Status","Protein Family Size","Ka/Ks (Mouse)","Ka/Ks (Macaque)","OMIM Disease","Average Allele Freqs","CEU Allele Freqs");


        for ($i = 0 ; $i < scalar @elts; $i++){
                if (@elts[$i] eq "1"){
                        $heading.= "\t@options[$i]";
                }
        }
        $heading.="\tUser Comment";
        return $heading;
}


sub get_output_row {
        $oo = @_[0];
        $pred_line = @_[1];
        @elts = split /\t/, $pred_line;
        $chrom= @elts [0];
        $coord_begin = @elts[1];
        $coord_end = @elts[2];
        $orn = @elts[3];
        $alleles = @elts[4];
        $usr_comment = @elts[5];
        $proteinid = @elts[6];
        $subst = @elts[7];
        $ensp = @elts[8];
        $rsid = @elts[9];
        $region = @elts[10];
        if ($region =~ /CDS/){
                $region = "EXON CDS";
        }
        elsif($region =~ /3\'UTR/){
                $region = "3\' UTR";
        }
        elsif($region =~ /5\'UTR/){
                $region = "5\' UTR";
        }
        elsif($region =~ /3\'UTR/){
                $region = "3\' UTR";
        }


        $codons = @elts[11];
        $snp_type = @elts[12];
        $score = @elts[13];
        $median = @elts[14];
        $seqs_rep = @elts[15];
        if ($score <= 0.05){
                if ($median > 3.25){
                        $pred = "DAMAGING *Warning! Low confidence.";
                }
                else{
                        $pred = "DAMAGING";
                }
        }
        else{
                $pred = "TOLERATED";
        }
        if ($score eq ""){
                $pred = "Not scored";
                $score = "NA";
                $median = "NA";
                $seqs_rep = "NA";
        }
        if ($COORD_SYSTEM eq "SPACE"){
                $coords = "$coord_begin-$coord_end";
        }
        else{
                $coords = "$coord_end";
        }
        $table_row =
"$chrom,$coords,$orn,$alleles\t$codons\t$proteinid\t$ensp\t$subst\t$region\t$rsid\t$snp_type\t$pred\t$score\t$median\t$seqs_rep";
        @elts2 = split /,/, $oo;
	for ($i = 0 ; $i < scalar @elts2; $i++){
                if (@elts2[$i] eq "1"){
                	$table_row.= "\t@elts[16+$i]";
                }
        }
        $table_row.="\t$usr_comment\n";

        return $table_row;

}





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

sub round {
    my($number) = shift;
    return int($number + .5);
}

sub
make_new_variant
{
        my $v = VariantInfo->new();
        $v->rsid ("NA");
        $v->orn_db ("");
        $v->enst ("");
        $v->ensp ("");
        $v->region ("");
        $v->snp_type ("NA");
        $v->codon1 ("");
        $v->codon2 ("");
        $v->nt1 ("");
        $v->nt2 ("");
        $v->alleles_db => '$',
        $v->AA1 ("NA");
        $v->AA2 ("NA");
        $v->AAPOS ("NA");
        $v->score ("");
        $v->median  ("");
        $v->seqs_rep  ("");
        $v->freq_av  ("");
        $v->freq_ceu  ("");
        $v->ensg  ("");
        $v->enst  ("");
        $v->ensp  ("");
        $v->gene_name  ("");
        $v->gene_desc ("");
        $v->ensfm  ("");
        $v->fam_desc  ("");
        $v->gene_status  ("");
        $v->fam_size  ("");
        $v->kaks_mouse  ("");
        $v->kaks_macaque  ("");
        $v->mim_status  ("");

        return $v;
}

sub
make_variant_with_elts
{
        my ($line) = @_;

        my @elts = split (/\t/, $line);
        my $v = VariantInfo->new();

        $v->rsid (@elts[4]);
        $v->enst (@elts[6]);
        $v->ensp (@elts[7]);
        $v->region (@elts[8]);
        $v->snp_type (@elts[9]);
        $v->nt1 (@elts[10]);
        $v->nt2 (@elts[11]);
        $v->alleles_db (@elts[10] . "\/" . @elts[11]);
#print "alleles db " . $v->alleles_db . "\n";
        $v->codon1 ( @elts[14]);
        $v->codon2 ( @elts[15]);
        $v->orn_db ( @elts[3]);
        $v->AA1 ( @elts[16]);
        $v->AA2 ( @elts[17]);
        $v->AAPOS ( @elts[19]);
        $v->score ( @elts[23]);
        $v->median ( @elts[24]);
        $v->seqs_rep ( @elts[25]);
        $v->ensg (@elts[26]);
        $v->gene_name (@elts[27]);
       $v->gene_desc (@elts[28]);
        $v->ensfm (@elts[29]);
        $v->fam_desc (@elts[30]);
        $v->gene_status (@elts[31]);
        $v->fam_size (@elts[32]);
        $v->kaks_mouse (@elts[33]);
        $v->kaks_macaque (@elts[34]);
        $v->mim_status (@elts[35]);
        $v->freq_av (@elts[36]);
        $v->freq_ceu (@elts[37]);
        return $v;
}

sub
return_variant_line
{
        my ($v) = @_;

        my $subst ="";
        if ($v->AA1 eq "NA") {
                $subst = "NA";
        } else {
                $subst = $v->AA1 . $v->AAPOS . $v->AA2;
        }
        my $codonchange = $v->codon1 . "-" . $v->codon2;

        my @fields ;
         push (@fields, $v->enst, $subst, $v->ensp, $v->rsid, $v->region,
                $codonchange, $v->snp_type, $v->score, $v->median, $v->seqs_rep, $v->ensg,
                $v->gene_name, $v->gene_desc, $v->ensfm, $v->fam_desc, $v->gene_status,
                $v->fam_size, $v->kaks_mouse, $v->kaks_macaque, $v->mim_status,
                $v->freq_av, $v->freq_ceu);
        my $line = join ("\t", @fields);
        return ($line);


};

sub check_ip_counts
{
	
## Check that this IP address hasn't been used too much
my $IP_address = $ENV{REMOTE_ADDR};

#       print "<HR>" . $IP_address . "<BR></HR> ";
my $remote_host = $ENV{REMOTE_HOST};

	my $ip_counts =
`cat  /home/blocks/apache/logs/access_log  | grep POST | grep $IP_address | wc -l `;
chomp($ip_counts);
if ( $ip_counts == "" ) {
        $ip_counts =
`cat /home/blocks/apache/logs/access_log  | grep POST | grep $remote_host | wc -l`;
        chomp($ip_counts);
}
#       print $ip_counts. "<BR>";
my $upper_limit = 50;
if ( $address ne "" ) {
        $upper_limit = 1000;
}
if ( $ip_counts > $upper_limit ) {
        print "<H1>Your computer has exceeded its daily limit.</H1><BR>";
        print
"Please download <A HREF=\"/\">SIFT software</A HREF> directly to your computer or <A HREF=\"/sift-bin/contact.pl\">contact</A HREF> us so that we can help you.  Thank you for using SIFT. <BR>";
        exit;
}

}

