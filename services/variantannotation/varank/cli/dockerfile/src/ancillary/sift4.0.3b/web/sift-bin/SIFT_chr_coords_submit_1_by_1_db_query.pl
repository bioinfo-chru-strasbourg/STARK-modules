#!/usr/local/bin/perl
#
# Nov 28, 2000 added option to email results
# 7/4/01 gap option is turned off permanently
#11/11/03 Added SIFT_queue stuff  JGH
# 3/8/04 Don't wait for search to finish, send them format.pl URL  JGH
#-------------------------------------------------------------------------

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use DBI;
use Tie::IxHash;
my $db =
  DBI->connect( "dbi:SQLite:dbname=/opt/www/sift/db/1.5/SIFT_DB_1.5.sqlite",
	"", "", { RaiseError => 1, AutoCommit => 1 } );
$sth_aln_enst = $db->prepare("select address_alignment from enst_address where ens_transcript_id = ?");


$| = 1;
require 'SIFT_subroutines.pm';

#       Set file permissions to rw-rw----
system("umask 006");
my $bin             = "/opt/www/sift/htdocs/sift-bin";
my $tmp             = "/opt/www/sift/tmp";
my $pid             = $$;
my $datadir         = "/opt/www/sift/data";
my $Human_db_dir    = "/opt/www/sift/db/Human_db";	
my $coding_info_dir = "/opt/www/sift/coding_info";
$program_call       = $bin . "/SIFTING2_for_subst_only.csh";
$return_address     = "sift\@fhcrc.org";
$snp_classifier     = "$bin/Classify_SNPs.pl";
$sequence_extractor = '/usr/local/projects/SIFT/snpclassifier/temp/bio.pl';

# output the beginning text to be used on all pages
print "Content-type: text/html\n\n";
print "<body bgcolor=white>\n";

if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
	print "This script should be referenced with a METHOD of POST\n";
	exit;
}

read( STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"} );
%names = &general_parse( $ENV{"CONTENT_TYPE"}, $QUERY_STRING );
my $address;
chomp( $names{address} );
$names{address} =~ s/\s+//;
if ( $names{address} ne "" ) {
	$address = $names{address};
}

#foreach (keys %names) {
#        print "The key $_ contains $names{$_}\n";
#}

## Check for validity of user inputs
my $all_chr_file = $tmp . "/$pid.chrfile";
if ( $names{CHR} eq "" && $names{CHR_file} eq "" ) {
	print
"<H1>Error</H1> Please enter some chromosome coordinates with substitutions.<P>\n";
	exit;
}

my $organism = $names{organism};
$organism =~ s/^\s+//;
$organism =~ s/\s+$//;
$bin_file = "$coding_info_dir/$organism/bins.list";
if ( $organism =~ /Select Organism/i ) {
	print
"<H1>Error</H1> Please select organism after pressing back button on your browser.<P>\n";
	exit;
}

#Read input list of chromosome coordinates and add to all_chr_string
open( CHR_FILE, ">$all_chr_file" );
my $all_chr_string;
if( $names{CHR_file} !~ /\d/){
	$names{CHR_file} = "";
}
if( $names{CHR} !~ /\d/){
        $names{CHR} = "";
}

if ($names{CHR_file} ne "" && $names{CHR} ne ""){
	print "Please choose only one of the following input methods after clicking the back button on your browser<BR>";
	print "1. Paste the input in the relevant textbox<BR>";
	print "2. Upload the text file containing input data"; 
	exit;	
}
if ($names{CHR_file} eq "" && $names{CHR} eq ""){
        print "Please choose one of the following input methods after clicking the back button on your browser<BR>";
        print "1. Paste the input in the relevant textbox<BR>";
        print "2. Upload the text file containing input data";
        exit;
}


if ( $names{CHR_file} ne "" ) {
	$names{CHR_file} =~ s/\r/\n/g;
	$names{CHR_file} =~ tr/A-Z/a-z/;
	if ( $names{CHR_file} =~ /\d/ ) {
		print CHR_FILE uc( $names{CHR_file} );
	}
	$all_chr_string = uc( $names{CHR_file} ), "\t";
	$input_method = "FILE";
}
else{
	$names{CHR} =~ tr/A-Z/a-z/;
	$names{CHR} =~ s/^\n//;
	$names{CHR} =~ s/\r/\n/g;
	print CHR_FILE uc( $names{CHR} );
	$all_chr_string .= uc( $names{CHR} ), "\t";
	$input_method = "TEXT";
}
close(CHR_FILE);

#Creat binned SNP files for SNP Classifier - the script also creates snp_chr_map_file to map bin file to chr database and table
system("$bin/map_choords_to_bin.pl $bin_file $all_chr_file $pid");

#parse snp_chr_map_file to fetch data from human databases.
open (SNP_CHR_MAP_FILE, "$tmp/$pid\_snp_chr_map_file.txt");
while (<SNP_CHR_MAP_FILE>){
	chomp;
	@elts = split /\t/, $_;
	$chr = @elts [0];	
	$bin = @elts[1];
	$snp_chr_bin_file = @elts[2];
	$db_chr_str = @elts[3];
	$table_chr = @elts[4];
	$db_chr =DBI->connect( "dbi:SQLite:dbname=/opt/www/sift/db/Human_db/Human_CHR$chr.sqlite","", "", { RaiseError => 1, AutoCommit => 1 } );
	open (BIN_FILE, "$tmp/$snp_chr_bin_file") || die ("Cannot open snp_chr_bin_file");
	while (<BIN_FILE>){
		chomp;
		@elts2 = split '\t',$_;
		$coord1 = @elts2[1];
		$coord2 = @elts2[2];
		$orn_usr = @elts2[3];
		$alleles = @elts2[4];
		$alleles =~ /(.+?)\/(.+?)/;
		$nt2a = $1;
		$nt2b = $2;
		print "select * from $table_chr where chr = \'chr$chr\' AND COORD1 = $coord1 and COORD2 = $coord2 AND nt2 = \'$nt2b\' AND CDS = 1 AND ENSP <> \'\'";
		$sth_db_chr_CDS = $db_chr->prepare("select * from $table_chr where chr = \'chr$chr\' AND COORD1 = $coord1 and COORD2 = $coord2 AND nt2 = \'$nt2b\' AND CDS = 1 AND ENSP <> \'\'");
		$sth_db_chr_no_CDS = $db_chr->prepare("select * from $table_chr where chr = \'chr$chr\' AND COORD1 = $coord1 and COORD2 = $coord2 and (nt1 = \'$nt2a\' or nt1 = \'$nt2b\') AND CDS = 0 ");
        	$sth_db_chr_CDS->execute();
		@rows_CDS = $sth_db_chr_CDS->fetchrow_array() ;
		
		print "3";
		if (scalar @rows_CDS == 0){	#if CDS not found, then maybe UTR or downstream.. save as supp info for user
			print "4";
			$sth_db_chr_no_CDS->execute();
			@rows_no_CDS = $sth_db_chr_no_CDS->fetchrow_array() ;
			if (scalar @rows_no_CDS == 0){		#coordinate not valid of intron
				print "5";
				$rsid  = "NA";
                                $enst = "NA";
                                $ensp = "NA";
                                $region = "INTRON";
				$snp_type = "NA";
                                $codon1 = "NA";
                                $codon2 = "NA";
                                $AA1 = "NA";
                                $AA2 = "NA";
                                $AAPOS2 = "NA";
				push @coord_arr_no_CDS, "$chr\t$coord1\t$coord2\t$orn_usr\t$alleles";
				$result = "$rsid\t$enst\t$ensp\t$region\t$snp_type\t$codon1\t$codon2\t$AA1\t$AA2\t$AAPOS2";
                        	$coord_result_no_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_usr\t$alleles"} = $result;
				system("echo $chr\t$coord1\t$coord2\t$orn_usr\t$alleles >> $tmp/$pid.query_status.txt");	

			}
			else{					#if region is UTR or downstream etc.
				print "6";
				$rsid  = @rows_no_CDS[4];
                        	$enst = @rows_no_CDS[6];
                        	$ensp = @rows_no_CDS[7];
                        	$region = @rows_no_CDS[8];
				$snp_type = @rows_no_CDS[9];
                        	$codon1 = @rows_no_CDS[14];
                        	$codon2 = @rows_no_CDS[15];
				$orn_db = @rows_no_CDS[3];
                        	$AA1 = @rows_no_CDS[16];
                        	$AA2 = @rows_no_CDS[17];
                        	$AAPOS2 = @rows_no_CDS[19];
				push @coord_arr_no_CDS, "$chr\t$coord1\t$coord2\t$orn_db\t$alleles";
				$result = "$rsid\t$enst\t$ensp\t$region\t$snp_type\t$codon1\t$codon2\t$AA1\t$AA2\t$AAPOS2";
                                $coord_result_no_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_db\t$alleles"} = $result;
				system("echo $chr\t$coord1\t$coord2\t$orn_db\t$alleles >> $tmp/$pid.query_status.txt");
			}
			
		}
		else{				#CDS found.. SIFT will run on these.
			print "7";
			$rsid  = @rows_CDS[4];
			$orn_db = @rows_CDS[3];
			$enst = @rows_CDS[6];
			$ensp = @rows_CDS[7];
			$region = @rows_CDS[8];
			$snp_type = @rows_CDS[9];
			$codon1 = @rows_CDS[14];
			$codon2 = @rows_CDS[15];
			$alleles_db = "@rows_CDS[10]\/@rows_CDS[11]";
			$AA1 = @rows_CDS[16];
			$AA2 = @rows_CDS[17];
			$AAPOS2 = @rows_CDS[19];
			push @coord_arr_CDS, "$chr\t$coord1\t$coord2\t$orn_db\t$alleles_db";			
			$result = "$rsid\t$enst\t$ensp\t$region\t$snp_type\t$codon1\t$codon2\t$AA1\t$AA2\t$AAPOS2";
			$coord_result_CDS_hash{"$chr\t$coord1\t$coord2\t$orn_db\t$alleles_db"} = $result;
			system("echo $chr\t$coord1\t$coord2\t$orn_db\t$alleles_db >> $tmp/$pid.query_status.txt");
		}


	}

}

#create enst_subst_list for coords in CDS - these will be predicted by SIFT.
foreach $coord (@coord_arr_CDS) {
	$result    = $coord_result_CDS_hash{$coord};
	@elts             = split /\t/, $result;
	$rsid = @elts[0];
	$enst = @elts[1];
	$ensp = @elts[2];
	$region = @elts[3];
	$snp_type = @elts[4];
	$codon1 = @elts[5];
	$codon2 = @elts[6];
	$aa1              = @elts[7];
	$aa2              = @elts[8];
	$pos              = @elts[9];
	
	$enst_line          = "$coord\t$enst\t$aa1$pos$aa2\t$ensp\t$rsid\t$region\t$codon1-$codon2\t$snp_type";
	push( @enst_subst_list, $enst_line );
}

#create enst_subst_list for coords not in CDS - these will not be predicted by SIFT.
foreach $coord (@coord_arr_no_CDS) {
        $result    = $coord_result_no_CDS_hash{$coord};
        @elts             = split /\t/, $result;
        $rsid = @elts[0];
        $enst = @elts[1];
        $ensp = @elts[2];
        $region = @elts[3];
        $snp_type = @elts[4];
        $codon1 = @elts[5];
        $codon2 = @elts[6];
        $aa1              = @elts[7];
        $aa2              = @elts[8];
        $pos              = @elts[9];
	if ($aa1 eq "NA"){
		$subst = "NA";
	}
	else{
		$subst ="$aa1$pos$aa2";
	}
        $enst_line          = "$coord\t$enst\t$subst\t$ensp\t$rsid\t$region\t$codon1-$codon2\t$snp_type";
        push( @enst_subst_list, $enst_line );
}

#create subst files for each enst and a hash table with enst as key and subst file name as value
foreach $enst_line (@enst_subst_list) {
	@elts = split /\t/, $enst_line;
	$enst = @elts[5];
	$enst =~ s/^\s+//g;
	$enst =~ s/\s+$//g;
	if ( $enst !~ /\d+/ ) {
		next;
	}
	$region = @elts[9];
	if ($region =~ /INTRON/i){
		next;	
	}

	open( SUBST_FILE, ">>$tmp/$pid.$enst.substfile" );
	$subst = @elts[6];
	$subst =~ s/^\s+//g;
	$subst =~ s/\s+$//g;

	print SUBST_FILE $subst, "\n";

	close(SUBST_FILE);
	$enst_subst_file_hash{$enst} = "$pid.$enst.substfile";
}

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

## SIFT prediction operations

$exp_option = 1;

#$info = $names{info};
$comments            = "$tmp/$pid.comments";
$out                 = $tmp . "/$$.siftresults";
$seq_identity_filter = $names{seq_identity_filter};
$seq_identity_filter = trim($seq_identity_filter);
chomp($seq_identity_filter);

######  Calling the program #########

print "<A NAME=top><H1><center>S<font color=red>I</font>F<font
color=blue>T</font> results</center></H1></A><BR>\n";
$sequences_to_select = "BEST";
if ( $names{sequences_to_select} =~ /ALL/ ) {
	$sequences_to_select = "ALL";
}
if ( $address ne "" ) {
	print "Results will also be mailed to $address.<BR><BR>\n";
}

select(STDOUT);
$|       = 1;
$counter = -1;
foreach $enst_line (@enst_subst_list) {
	$counter++;
	@elts = split /\t/, $enst_line;
	$enst   = @elts[5];
	$subst_file = $enst_subst_file_hash{$enst};
	chomp;
	$region = @elts[9];
        if ($region =~ /INTRON/i){
                next;
        }

	$address_aln = "";
	$alignedfasta = "$tmp\/$pid\.$enst.alignedfasta";
	if ($index_enst_seen{$enst} != 1){
		$sth_aln_enst->execute($enst);
		@rows_enst_address = $sth_aln_enst->fetchrow_array() ;
		if ( scalar @rows_enst_address == 0 ) {
        		open (DUMMY, ">$tmp/$pid\.$enst.alignedfasta");
			close(DUMMY);
		}
		else{
			$address_aln_enst = "$datadir\/" . @rows_enst_address[0];
	                system("cp $address_aln_enst $alignedfasta");
		}
		$index_enst_seen{$enst} = 1;
	}
	else{
		next;
	}



# if aligned fasta not found - then add a tag to the gi_line and replace in gi_subst_file

	$comments = "$tmp/$pid.$enst.comments";
	$out      = $tmp . "/$$.$enst.siftresults";

### program call ###
	#print "$program_call $pid alignedseqs $alignedfasta $out 0 $subst_file 0 $seq_identity_filter $sequences_to_select $address<BR>";
	system(
"$program_call $pid alignedseqs $alignedfasta $out 0 $subst_file 0 $seq_identity_filter $sequences_to_select $address > $comments 2>&1"
	);
	print "</PRE>";

	# SIFT successful -- OUTPUT RESULTS
	$commentscsh = $tmp . "/$pid.commentscsh";

	$psiblastout = $alignedfasta;
	$seqno       = `grep ">" $psiblastout | wc -l`;
	$seqno--;    # subtract the QUERY sequence

	system("cat $tmp/$pid.*.error >> $tmp/$pid.err");

	#errors or warnings that don't disrupt program
	open( OUT, "<$tmp/$pid.err" );
	while ( $_ = <OUT> ) {

	   # don't want to print Jorja's messages, only my own.  so my error
	   # messages all have ERROR in the line, but don't print blastpgp error
	   # for reading checkpoing file
		if ( /ERROR/ || /WARNING/ || /\*\*\*/ ) {
			unless (/blastpgp/) {

				#print;
			}
		}
	}
	close(OUT);
}

# get results from individual prediction files into an index.
open( ENSTFILE, ">$tmp/$pid.enstfile" );
foreach $enst_line (@enst_subst_list) {
	#print $enst_line,"<BR>";
	print ENSTFILE "$enst_line\n";
	@elts             = split /\t/, $enst_line;
	$chromosome       = @elts[0];
	$coord1           = @elts[1];
	$coord2           = @elts[2];
	$orn              = @elts[3];
	$alleles	  = @elts[4];
	$enst = @elts[5];
	$subst            = @elts[6];
	$ensp = @elts[7];
	$rsid = @elts[8];
	$region =@elts[9] ;
	$codons = @elts[10];
	$snp_type =@elts[11];
	$enst =~ s/^\s+//g;
	$enst =~ s/\s+$//g;
        if ($region =~ /INTRON/i){
		next;
        }


	if ( $enst =~ /not found/i || $enst_seen{$enst} == 1 ) {  # not found present
		next;
	}

	else {                                                   #"not found" absent
		$enst_seen{$enst} = 1;
		$prediction_file = "$tmp/$pid.$enst.siftresults.predictions";
		open( PREDICTION_FILE, "$prediction_file" );
		while (<PREDICTION_FILE>) {
			chomp;
			if ( $_ =~
/Substitution at pos (\d+) from (.+?) to (.+?) is predicted to (.+?) with a score of (\d\.\d+)/i
			  )
			{
				$pos   = $1;
				$aa1   = $2;
				$aa2   = $3;
				$pred  = $4;
				$score = $5;
				$subst = "$aa1$pos$aa2";
				if ( $pred =~ /TOLERATED/ ) {
					$pred = "TOLERATED";
				}
				elsif ( $pred =~ /not predicted/i || $pred =~ /not scored/i ) {
					$pred = "NOT PREDICTED";
				}
				else {
					$pred = "DAMAGING";
				}
			}
			if ( $_ =~ /Median Sequence conservation: (.+)/i ) {
				$median = $1;
				if ( $median > 3.25 && $pred =~ /DAMAGING/i ) {
					$warningflag = 1;
					$pred        = $pred . " *Warning! Low confidence.";
				}
			}
			if ( $_ =~ /Sequences represented at this position:(.+)/i ) {
				$numseqs = $1;

				$output_row = "$pred\t$score\t$median\t$numseqs";
				$index_enst_pred{"$enst\t$subst"} = "$output_row";
				$chrom     = "";
				$orn       = "";
				$alleles  = "";
				$proteinid = "";
				$subst     = "";
				$rsid      = "";
				$pred      = "";
				$score     = "";
				$median    = "";
				$numseqs  = "";
			}

		}
	}
}
close(ENSTFILE);

print "<BR><BR>";

#build combined predictions file from individual predictions files.
$combined_siftresults_file = "$pid.combined_siftresults";

open (ENST_PRED_FILE,">$tmp/$pid.enst_pred_file");
foreach $enst_line (@enst_subst_list) {
	@elts             = split /\t/, $enst_line;
	$chromosome       = @elts[0];
	$coord1           = @elts[1];
	$coord2           = @elts[2];
	$orn              = @elts[3];
	$alleles = @elts[4];
	$enst = @elts[5];
	$subst            = @elts[6];
	$enst =~ s/^\s+//g;
	$enst =~ s/\s+$//g;
	if (exists($index_enst_pred{"$enst\t$subst"})){
		$output_row = $index_enst_pred{"$enst\t$subst"};
	}
	else{
		$output_row = "Not Scored\tNA\tNA\tNA";
	}
	print ENST_PRED_FILE "$enst_line\t$output_row\n";
		
}
close(ENST_PRED_FILE);

#read from combined file and print table
open( ENST_PRED_FILE, "$tmp/$pid.enst_pred_file" )
  || die("Cannot open predictions file");

$heading =
"Coordinates\tAlleles\tCodons\tProtein ID\tRegion\tSubstitution\tdbSNP ID\tSNP Type\tPrediction\tScore\tMedian Info\t\# Seqs at position";
open( OUTFILETABLE, ">$tmp/$combined_siftresults_file.predictions.table.html" );
open( OUTFILETSV,   ">$tmp/$combined_siftresults_file.predictions.tsv" );
print OUTFILETSV "$heading\n";
print OUTFILETABLE '<table border="1" cellspacing="0" cellpadding="4">';
print OUTFILETABLE '<thead>';
print OUTFILETABLE '<tr><th>';
$heading =~ s?\t?</th><th>?g;
print OUTFILETABLE "$heading</th></tr>\n";
my $warningflag = 0;

while (<ENST_PRED_FILE>) {
	chomp;
	@elts = split /\t/, $_;
	$chrom= @elts [0];
	$coord_begin = @elts[1];
	$coord_end = @elts[2];
	$orn = @elts[3];
	$alleles = @elts[4];
	$proteinid = @elts[5];
	$subst = @elts[6];
	$ensp = @elts[7];
	$rsid = @elts[8];
	$region = @elts[9];
	$codons = @elts[10];
	$snp_type = @elts[11];	
	$pred = @elts[12];
	$score = @elts[13];
	$median = @elts[14];
	$numseqs = @elts[15];	
	$table_row =
"$chrom,$coord_begin-$coord_end,$orn\t$alleles\t$codons\t$proteinid,$ensp\t$region\t$subst\t$rsid\t$snp_type\t$pred\t$score\t$median\t$numseqs\n";
		$chrom     = "";
		$orn       = "";
		$alleles = "";
		$codons = "";
		$ensp = "";	
		$proteinid = "";
		$region = "";
		$snp_type = "";
		$subst     = "";
		$rsid      = "";
		$pred      = "";
		$score     = "";
		$median    = "";
		$numseqs  = "";

		print OUTFILETSV "$table_row";
		print OUTFILETABLE "<tr>\n";
		my @fields = split( '\t', $table_row );
		for $cell (@fields) {
			if ( $cell =~ /DAMAGING/ || $cell =~ /Warning/ ) {
				print OUTFILETABLE "<td><font color=red>$cell</font></td>";
			}
			elsif ( $cell =~ /NOT PREDICTED/i || $cell =~ /not found/i ) {
				print OUTFILETABLE "<td><font color=blue>$cell</font></td>";
			}
			else {
				print OUTFILETABLE "<td>$cell</td>";
			}

		}
		print OUTFILETABLE "</tr>\n";
	

}
print OUTFILETABLE "</tr>\n</tbody>\n</table>\n<BR>";
if ( $warningflag == 1 ) {

	print OUTFILETABLE
"<font color=red>* Low confidence means that the protein alignment does not have enough sequence diversity. Because the position artifically appears to be conserved, an amino acid may incorrectly predicted to be damaging.</font><BR><BR>";
}
close(FILE);
close(OUTFILETABLE);
close(OUTFILETSV);

#print the table
open( OUTFILETABLE, "$tmp/$combined_siftresults_file.predictions.table.html" )
  || die("Cannot open predictions table");
while (<OUTFILETABLE>) {
	print;
}
print
"Click <A HREF=\"\/tmp\/$combined_siftresults_file.predictions.tsv\">here</A> to download the following table in tab separated format. You can open it using excel with delimiter set as TAB";
print "<BR>";

print "<BR>";
print
"<i>If you received a warning that the sequences were not diverse enough, you can have SIFT choose more diverse sequences <A HREF=\"/www/SIFT_seq_submit2.html\">here.</A></i><BR><BR>";

#email the results
if ( $address ne "" ) {
	open( MESSAGE, ">$tmp/$pid.email_message.txt" );
	print MESSAGE
"Dear User\n\nThank you for using SIFT.\n\nPlease find the results of your recent query attached with this message.\nRemember this job id \"$pid\" for any future correspondance.\nDo not hesitate to contact us if you have any questions about SIFT.\n\nThanks\nSIFT Team\nJ Craig Venter Institute (West Coast Campus)\n10355 Science Center Drive\nSan Diego, CA 92121\nUSA";
	close(MESSAGE);
	system(
"mutt -F /opt/www/sift/htdocs/.muttrc -a $tmp/$combined_siftresults_file.predictions.tsv -s \"SIFT Results for Job ID $pid\" $address <$tmp/$pid.email_message.txt"
	);
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


