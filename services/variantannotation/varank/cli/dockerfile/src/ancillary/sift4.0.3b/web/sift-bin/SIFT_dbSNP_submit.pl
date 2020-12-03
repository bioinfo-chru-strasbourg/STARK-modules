#!/usr/local/bin/perl
#
# February 12, 2005

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use DBI;
use Tie::IxHash;

my $db =
  DBI->connect(
	"dbi:SQLite:dbname=/opt/www/sift/db/1.4/SIFT_DB_1.4.sqlite",
	"", "", { RaiseError => 1, AutoCommit => 1 } );
# db_supp to retrieve allele freqs is hardcoded to ncbi36
my $db_supp = DBI->connect( "dbi:SQLite:dbname=/opt/www/sift/db/Human_db/Human_Supp.sqlite","", "", { RaiseError => 1, AutoCommit => 1 } );
$db_supp->do('PRAGMA synchronous=1');

$sth_best =
  $db->prepare(
	"select * from prediction where rsid = ? AND option = \'BEST_HITS\'");
$sth_all =
  $db->prepare(
	"select * from prediction where rsid = ? AND option = \'ALL_HITS\'");

$sth_aln_pid_all =
  $db->prepare(
        "select address_alignment_all from pid_address where pid = ?");

$sth_aln_pid_best =
  $db->prepare(
        "select address_alignment_best from pid_address where pid = ?");

$sth_db_supp_allelefreq = 
	$db_supp->prepare("select * from ALLELE_FREQ where RSID = ?");

#use strict;
$| = 1;
require 'SIFT_subroutines.pm';

#       Set file permissions to rw-rw----
system("umask 006");
my $bin            = "/opt/www/sift/htdocs/sift-bin";
my $tmp            = "/opt/www/sift/tmp";
my $pid            = $$;
my $return_address = "sift\@fhcrc.org";
my $datadir = "/opt/www/sift/data";
# be sure to take out "No PREDICTION" in $dbSNP_table file
my $dbSNP_table = "/opt/www/sift/htdocs/www/collated_predictions.tsv";

# output the beginning text to be used on all pages
print "Content-type: text/html\n\n";
print "<body bgcolor=white>\n";
if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
	print "This script should be referenced with a METHOD of POST\n";
	exit;
}

read( STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"} );
my %names = &general_parse( $ENV{"CONTENT_TYPE"}, $QUERY_STRING );
my $address;
chomp( $names{address} );
$names{address} =~ s/\s+//;
if ( $names{address} ne "" ) {
	$address = $names{address};
}
### substitution file #######
my $all_snp_file = $tmp . "/$pid.snpfile";
if ( $names{SNP} eq "" && $names{SNP_file} eq "" ) {
	print "<H1>Error</H1> Please enter some rs id's.<P>\n";
	exit;
}
open( SNP_FILE, ">$all_snp_file" );
my $all_snp_string;
if ( $names{SNP_file} ne "" ) {
	$names{SNP_file} =~ s/\r/\n/g;
	$names{SNP_file} =~ tr/A-Z/a-z/;
	print SNP_FILE $names{SNP_file};
	$all_snp_string = $names{SNP_file};
}
if ( $names{SNP} ne "" ) {
	$names{SNP} =~ tr/A-Z/a-z/;
	$names{SNP} =~ s/\r/\n/g;
	print SNP_FILE $names{SNP};
	$all_snp_string .= $names{SNP};
}
close(SNP_FILE);
#### output options
my $print_all_freq = 0; my $print_ceu_freq = 1;
if ($names{allfreq}) { $print_all_freq = 1;}
if ($names{ceufreq}) { $print_ceu_freq = 1;}

my @fields = split( /\n/, $all_snp_string );

#print "all snp string $all_snp_string<BR>\n";

my %snp_hash;
tie (%snp_hash, Tie::IxHash);

for ( my $i = 0 ; $i < @fields ; $i++ ) {
	$fields[$i] =~ s/\s+//;
	if ( $fields[$i] =~ /^rs/ ) {

		#	print "hashing $fields[$i]<BR>\n";
		$snp_hash{ $fields[$i] } = 1;
	}
}


## SIFT prediction operations

my $comments = "$tmp/$pid.comments";
my $out      = $tmp . "/$pid.html";
my $out_tsv = "$tmp/$pid.tsv";
open( HTML_OUT, ">$out" ) || die "Can't open $out";

######  Calling the program #########

print HTML_OUT "<A NAME=top><H1><center>S<font color=red>I</font>F<font 
color=blue>T</font> results</center></H1></A><BR>\n";
print HTML_OUT "Processing...\n";
print HTML_OUT "If your browser times out before results are shown, html results can be seen at http://sift.jcvi.org/tmp/$pid\.html and tsv results at http://sift.jcvi.org/tmp/$pid\.tsv Both files are stored for 1 hour before being deleted.<BR><BR>";

if ( $address ne "" ) {
	print "Results will also be mailed to $address.<BR><BR>\n";

	#system ("echo $address >> /home/sift/email_log");
}
select(STDOUT);
$| = 1;

# SIFT successful -- OUTPUT RESULTS
$commentscsh = $tmp . "/$pid.commentscsh";
$psiblastout = $tmp . "/$pid" . ".alignedfasta";
print_dbSNP_predictions( $dbSNP_table, \%snp_hash );
#close(HTML_OUT);
open( HTML_OUT, "<$out" );
while (<HTML_OUT>) {
	print;
}

close(HTML_OUT);
print "<BR>Click <A HREF=\"\/tmp\/$pid.tsv\">here</A> to download the table in tab separated format. You can open it using excel with delimiter set as TAB<br>";


#email the results
if ($address ne ""){
        open (MESSAGE, ">$tmp/$pid.email_message.txt");
        print MESSAGE "Dear User\n\nThank you for using SIFT.\n\nPlease find the results of your recent query attached with this message.\nRemember this job id \"$pid\" for any future correspondance.\nDo not hesitate to contact us if you have any questions about SIFT.\n\nThanks\nSIFT Team\nJ Craig Venter Institute (West Coast Campus)\n10355 Science Center Drive\nSan Diego, CA 92121\nUSA";
        close (MESSAGE);
        system("mutt -F /opt/www/sift/htdocs/.muttrc -a $tmp/$pid.tsv -s \"SIFT Results for Job ID $pid\" $address <$tmp/$pid.email_message.txt");
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
				if ( $loc_semi > 0 && $loc_semi < $loc ) { $loc = $loc_semi; }
				if ( $loc > 0 ) { $temp = substr( $temp, 0, $loc ); }

				#               Get rid of quotes around the name
				$temp =~ s/\"//g;

				#               Still has a trailing whitespace character ...
				$temp =~ s/\s//g;

		  #               Substitute spaces with nothing
		  #		Need to strip leading/ending whitespace off of $temp1,
		  #               but be careful not to strip off internal CRs
		  #               MAC file lines end in just \r, no \n, so makelis won't
		  # find all
		  #               of the sequences; DOS file lines end in \r\n, UNIX in
		  #\n.
		  #               Change \r\n to \n and then \r to \n
#######PROGRAM -SPECIFIC!!!!!!!######################
				$temp1 =~ s/\r\n/\n/g;
				$temp1 =~ s/\r/\n/g;
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

sub print_dbSNP_predictions {
	open( HTML_OUT, "<$out" );
	open (TSV_OUT,">$out_tsv");
        while (<HTML_OUT>) {
                print;
        }
        close (HTML_OUT);
	open( HTML_OUT, ">$out" );
	print HTML_OUT "Done.<br><br>\n";
	my $dbSNP_file = shift;
	my $snp_hash_ref = shift;

	open( DBSNP, $dbSNP_file ) || die "can't open $dbSNP_file";
	my $no_prediction_line = "<TD COLSPAN=4 ROWSPAN=2>No prediction.</TD>";
	$tsv_heading = "SNP\tAmino acid change\tProtein ID\tAmino acid\tPrediction (Orthologs)\tScore (Orthologs)\tMedian info (Orthologs)\tPrediction (Homologs)\tScore (Homologs)\tMedian info (Homologs)";
	print TSV_OUT "$tsv_heading\n"; 
	print HTML_OUT "<TABLE BORDER>\n";
	print HTML_OUT
"<TH ROWSPAN=2>SNP</TH><TH ROWSPAN=2>Amino acid change</TH><TH ROWSPAN=2>Protein ID</TH><TH ROWSPAN=2>Amino Acid</TH>\n";
	print HTML_OUT
	  "<TH COLSPAN=4>Using orthologues in the protein alignment</TH>\n";
	print HTML_OUT
	  "<TH COLSPAN=4>Using homologues in the protein alignment</TH></TR>\n";
	print HTML_OUT
"<TR><TH>Prediction</TH><TH>Score</TH><TH>Median Sequence IC</TH><TH>Alignment</TH>\n";
	print HTML_OUT
"<TH>Prediction</TH><TH>Score</TH><TH>Median info</TH><TH>Alignment</TH></TR>\n";
	my $dbsnp_line;
	


	foreach $snpid ( keys %snp_hash ) {
		$sth_best->execute($snpid);
		$rows_best = $sth_best->rows;
		my (
			$snpid,          $alignmenttype, $protein,
			$aasub,          $aa1_besthit,   $prediction1_besthit,
			$score1_besthit, $info1_besthit, $d7,
			$d8,             $aa2_besthit,   $prediction2_besthit,
			$score2_besthit, $info2_besthit, $d2,
			$d3
		);
		$snpid=$alignmenttype=$protein=$aasub=$aa1_besthit=$prediction1_besthit=$score1_besthit=$info1_besthit=$d7="";
                $d8= $aa2_besthit=$prediction2_besthit="";
                $score2_besthit= $info2_besthit= $d2=$d3 = "";

		$sth_best->bind_columns(
			undef,                 \$snpid,
			\$alignmenttype,       \$protein,
			\$aasub,               \$aa1_besthit,
			\$prediction1_besthit, \$score1_besthit,
			\$info1_besthit,       \$d7,
			\$d8,                  \$aa2_besthit,
			\$prediction2_besthit, \$score2_besthit,
			\$info2_besthit,       \$d2,
			\$d3
		);

		for ( $i = 0 ; $i < 1 ; $i++ ) {
			(
				$snpid,          $alignmenttype, $protein,
				$aasub,          $aa1_besthit,   $prediction1_besthit,
				$score1_besthit, $info1_besthit, $d7,
				$d8,             $aa2_besthit,   $prediction2_besthit,
				$score2_besthit, $info2_besthit, $d2,
				$d3
			  )
			  = $sth_best->fetchrow_array();
		}

		my $alignmenttype2_besthit = "BEST_HITS";

		$sth_all->execute($snpid);
		my (
			$snpid,         $alignmenttypei3_allhit, $protein,
			$aasub,         $aa3_allhit,             $prediction3_allhit,
			$score3_allhit, $info3_allhit,           $d3,
			$d4,            $aa4_allhit,             $prediction4_allhit,
			$score4_allhit, $info4_allhit,           $d5,
			$d6
		);
		$snpid=         $alignmenttypei3_allhit= $protein="";
                $aasub=         $aa3_allhit=             $prediction3_allhit="";
                $score3_allhit= $info3_allhit=           $d3="";
                $d4=            $aa4_allhit=             $prediction4_allhit="";
                $score4_allhit= $info4_allhit=           $d5=  $d6="";

		$sth_all->bind_columns(
			undef,                \$snpid,
			\$alignmenttype,      \$protein,
			\$aasub,              \$aa3_allhit,
			\$prediction3_allhit, \$score3_allhit,
			\$info3_allhit,       \$d3,
			\$d4,                 \$aa4_allhit,
			\$prediction4_allhit, \$score4_allhit,
			\$info4_allhit,       \$d5,
			\$d6
		);
		for ( $i = 0 ; $i < 1 ; $i++ ) {
			(
				$snpid,         $alignmenttypei3_allhit,
				$protein,       $aasub,
				$aa3_allhit,    $prediction3_allhit,
				$score3_allhit, $info3_allhit,
				$d3,            $d4,
				$aa4_allhit,    $prediction4_allhit,
				$score4_allhit, $info4_allhit,
				$d5,            $d6
			  )
			  = $sth_all->fetchrow_array();

		}
		if($sth_all->rows == 0 && $sth_best->rows == 0){
			next;
		}
		my $alignmenttype4_allhit = "ALL_HITS";
		
		if ($print_all_freq || $print_ceu_freq ) {
	
		$sth_db_supp_allelefreq->execute($snpid);
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
                $freq_av = "\t" . "$fav1:$fav2";
                $freq_ceu = "\t" . "$fceu1:$fceu2";
#print HTML_OUT "retrive freq $snpid freq $freq_av $freq_ceu<BR>\n";
# frequencies are being printed out, but not sure if Prateek is calculating
# freqs correctly. omit for now 
		if ($freq_av == 0) {
			$freq_av = "";
		}
		if ($freq_ceu == 0) {
			$freq_ceu = "";
		}
		} else { # no allele frequencies were requested 
			$freq_av = "";
			$freq_ceu = "";
		}
	
		$dbsnp_line = "$snpid\t$protein\t$aasub\t$aa1_besthit\t$alignmenttype\t$prediction1_besthit\t$score1_besthit\t$info1_besthit\t$d7\t$d8\t$aa2_besthit\t$alignmenttype2_besthit\t$prediction2_besthit\t$score2_besthit\t$info2_besthit\t$d2\t$d3\t$aa3_allhit\t$alignmentype3_allhit\t$prediction3_allhit\t$score3_allhit\t$info3_allhit\t$d3\t$d4\t$aa4_allhit\t$alignmenttype4_allhit\t$prediction4_allhit\t$score4_allhit\t$info4_allhit\t$d5\t$d6";
		$dbsnp_line = modify_with_html_tags ($dbsnp_line);
		my ($snpid, $protein, $aasub, $aa1_besthit, $alignmenttype,
                $prediction1_besthit, $score1_besthit, $info1_besthit, $d7, $d8,
                $aa2_besthit, $alignmenttype2_besthit,
                $prediction2_besthit, $score2_besthit, $info2_besthit, $d2, $d3,
                $aa3_allhit, $alignmentype3_allhit,
                $prediction3_allhit, $score3_allhit,
                $info3_allhit, $d3, $d4,
                $aa4_allhit, $alignmenttype4_allhit,
                $prediction4_allhit, $score4_allhit,
                $info4_allhit, $d5, $d6) = split (/\t/, $dbsnp_line);

		my ( $best_hit_prot_alignment_link, $all_hit_prot_alignment_link ) =
		  get_msf_links($protein);

		$snp_hash{$snpid} = "SEEN";
		$_ = $aasub;
		/([A-Z\*])([0-9]+)([A-Z\*])/;
		my $refaa = $1;
		my $pos   = $2;
		my $newaa = $3;

		#			print HTML_OUT "refaa $refaa new aa $newaa<BR>\n";
		print HTML_OUT "<TR><TD ROWSPAN=2>$snpid</TD>\n";
		print HTML_OUT "<TD ROWSPAN=2>$aasub</TD>\n";
		print HTML_OUT "<TD ROWSPAN=2>$protein</TD>\n";
		print TSV_OUT "$snpid\t$aasub\t$protein\t";	

		#			print HTML_OUT "<BR>score1 $score1_besthit\n";

		my $refaa_prediction_line = print_multiple_cells( ($refaa) );
		my $newaa_prediction_line = print_multiple_cells( ($newaa) );
		
		# print reference prediction (1rst row)
		my $prediction1_line =
		  print_multiple_cells( ( $prediction1_besthit, $score1_besthit ) );
		my $prediction2_line =
		  print_multiple_cells( ( $prediction2_besthit, $score2_besthit ) );
		if ( $aa1_besthit eq $refaa ) {
			$refaa_prediction_line .= $prediction1_line;
			$newaa_prediction_line .= $prediction2_line;
		}
		elsif ( $aa2_besthit eq $refaa ) {
			$refaa_prediction_line .= $prediction2_line;
			$newaa_prediction_line .= $prediction1_line;
		}

		if ( $info1_besthit > 0.0 ) {
			$refaa_prediction_line .=
			    "<TD ROWSPAN=2>$info1_besthit"
			  . "</TD><TD ROWSPAN=2>"
			  . $best_hit_prot_alignment_link . "</TD>";
		}
		else {
			$refaa_prediction_line .= $no_prediction_line;
		}
		my $prediction3_line =
		  print_multiple_cells( ( $prediction3_allhit, $score3_allhit ) );
		my $prediction4_line =
		  print_multiple_cells( ( $prediction4_allhit, $score4_allhit ) );
		if ( $aa3_allhit eq $refaa ) {
			$refaa_prediction_line .= $prediction3_line;
			$newaa_prediction_line .= $prediction4_line;
		}
		elsif ( $aa4_allhit eq $refaa ) {
			$refaa_prediction_line .= $prediction4_line;
			$newaa_prediction_line .= $prediction3_line;
		}

		#print HTML_OUT "$snpid  $info3_allhit<BR>";
		if ( $info3_allhit > 0 ) {
			$refaa_prediction_line .=
			    "<TD ROWSPAN=2>$info3_allhit"
			  . "</TD><TD ROWSPAN=2>"
			  . $all_hit_prot_alignment_link . "</TD>";
		}
		else {
			$refaa_prediction_line .= $no_prediction_line;

			#				$newaa_prediction_line .= $no_prediction_line;;
		}
		$tsv_row1 = "$refaa\t$prediction1_besthit\t$score1_besthit\t$info1_besthit\t$prediction3_allhit\t$score3_allhit\t$info3_allhit";
		$tsv_row2 = "\t\t\t$newaa\t$prediction2_besthit\t$score2_besthit\t$info2_besthit\t$prediction4_allhit\t$score4_allhit\t$info4_allhit";
		$tsv_row1 =~ s/<.+?>/ /g;
		$tsv_row2 =~ s/<.+?>/ /g;
		print TSV_OUT "$tsv_row1\n$tsv_row2\n";
		print HTML_OUT $refaa_prediction_line;
		print HTML_OUT "</TR><TR>";
		print HTML_OUT $newaa_prediction_line;
		print HTML_OUT "</TR>";

	}
	# print out SNPs not seen
	while ( my ( $id, $seen ) = each %snp_hash ) {
		if ( $seen ne "SEEN" ) {
			print HTML_OUT
			  "<TR><TD>$id</TD><TD colspan=11>Not found</TD></TR>\n";
			print TSV_OUT "$id Not found\n";
		}
	}

	print HTML_OUT "</TABLE><BR>\n";
	print HTML_OUT
"<font color=red>* Low confidence means that the protein alignment does not have enough sequence diversity. Because the position artifically appears to be conserved, an amino acid may incorrectly predicted to be damaging.  </font><BR>";
	print TSV_OUT "Low confidence means that the protein alignment does not have enough sequence diversity. Because the position artifically appears to be conserved, an amino acid may incorrectly predicted to be damaging.\n";
	close (TSV_OUT);

}

sub print_multiple_cells {

	my (@fields) = @_;
	my $line1;
	for ( my $i = 0 ; $i < @fields ; $i++ ) {
		$line1 .= "<TD>";
		$line1 .= "$fields[$i]";
		$line1 .= "</TD>\n";
	}
	return ($line1);
}

sub get_msf_links {
	my ($protein) = @_;
	$sth_aln_pid_all->execute($protein);
        my $all_hits_alignment_file = "$datadir\/".$sth_aln_pid_all->fetchrow_array();
	$sth_aln_pid_best->execute($protein);
        my $best_hits_alignment_file = "$datadir\/".$sth_aln_pid_best->fetchrow_array();
	#print "$all_hits_alignment_file<BR>$best_hits_alignment_file";

	#my $all_hits_alignment_file =
#`ls /home/sift/BLink_alignments/Predictions_ALL_HIT_Proteins/*/$protein*.msf`;

	#my $best_hits_alignment_file =
#`ls /home/sift/BLink_alignments/Predictions_BEST_HIT_Proteins/*/$protein*.msf`;
	chomp($all_hits_alignment_file);
	chomp($best_hits_alignment_file);

#	print HTML_OUT "find allhits alignment $all_hits_alignment_file best $best_hits_alignment_file<BR>\n";
	my $newfile1 = "$tmp/" . $protein . ".allhits.msf";
	my $newfile2 = "$tmp/" . $protein . ".besthits.msf";
	system("cp $all_hits_alignment_file  $newfile1");
	system("cp $best_hits_alignment_file $newfile2");
	system("chmod 660 $newfile2");
	system("chmod 660 $newfile1");

	#print HTML_OUT "newfile1 is $newfile1 newfile2 $newfile2<BR>\n";
	my $allhits_link =
	    "<A HREF=\"/sift-bin/catfile.csh?"
	  . $newfile1
	  . "+Alignment+PRE\" TARGET=_blank>Alignment of all hits</A>";

	my $besthits_link =
	    "<A HREF=\"/sift-bin/catfile.csh?"
	  . $newfile2
	  . "+Alignment+PRE\" TARGET=_blank>Alignment of best hits</A>\n";
	return ( $besthits_link, $allhits_link );
}

sub modify_with_html_tags {
	my ($dbsnp_line) = @_;

	#	print HTML_OUT "dbsnp line $dbsnp_line";
	my @fields = split( /\t/, $dbsnp_line );
	for ( my $i = 0 ; $i < @fields ; $i++ ) {
		if ( $fields[$i] =~ /DELETERIOUS/ ) {
			$fields[$i] = "<font color=red>DAMAGING</font>";
		}
	}
	my $index;
	foreach $index ( ( 7, 21 ) ) {
		if (
			$fields[$index] > 3.25
			&& (   $fields[ $index - 2 ] =~ /DAMAGING/
				|| $fields[ $index + 5 ] =~ /DAMAGING/ )
		  )
		{
			$fields[$index] .=
			  "<BR><i><font color=red>Warning! Low confidence.*</font></i>";
		}
	}
	my $newline;
	for ( my $i = 0 ; $i < @fields ; $i++ ) {
		$newline .= $fields[$i] . "\t";
	}

	#	print HTML_OUT $newline;
	return ($newline);

}


