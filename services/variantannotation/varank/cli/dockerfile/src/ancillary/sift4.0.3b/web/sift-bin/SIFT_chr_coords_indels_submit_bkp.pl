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
$sth_aln_pid_all =
  $db->prepare("select address_alignment_all from pid_address where pid = ?");

$sth_aln_pid_best =
  $db->prepare("select address_alignment_best from pid_address where pid = ?");

$sth_aln_gi_all =
  $db->prepare(
"select address_alignment_all from pid_address where pid = (select pid from pid_gi where gi = ?)"
  );

$sth_aln_gi_best =
  $db->prepare(
"select address_alignment_best from pid_address where pid = (select pid from pid_gi where gi = ?)"
  );

$sth_pid_subst = $db->prepare("select substitution from dbsnp where pid = ?");
$sth_gi_subst  =
  $db->prepare(
"SELECT substitution FROM dbsnp where pid  = (select pid from pid_gi where gi = ?)"
  );

$| = 1;
require 'SIFT_subroutines.pm';

#       Set file permissions to rw-rw----
system("umask 006");
my $bin             = "/opt/www/sift/htdocs/sift-bin";
my $tmp             = "/opt/www/sift/tmp";
my $pid             = $$;
my $datadir         = "/opt/www/sift/data";
my $coding_info_dir = "/opt/www/sift/coding_info";
$program_call       = $bin . "/SIFTING2_for_subst_only.csh";
$return_address     = "sift\@fhcrc.org";
$snp_classifier     = "$bin/Classify_SNPs_Indels.pl";
$reformat_chrfile   = "$bin/reformat_chrfile.pl";
$detect_indel 	    = "$bin/detect_indel.pl";
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
$organism = "Homo_sapien_partitioned";
$bin_file = "$coding_info_dir/$organism/bins.list";
if ( $organism =~ /Select Organism/i ) {
	print
"<H1>Error</H1> Please select organism after pressing back button on your browser.<P>\n";
	exit;
}

##set paths to coding info and chromosome fasta files for selected organism and chromosome.
$coding_info_file = "$coding_info_dir/$organism/$chromosome.*.coding_info";
$chr_fasta_file   = "$coding_info_dir/$organism/$chromosome.*.fasta";

open( CHR_FILE, ">$all_chr_file" );
my $all_chr_string;
if ( $names{CHR_file} ne "" ) {
	$names{CHR_file} =~ s/\r/\n/g;
	$names{CHR_file} =~ tr/A-Z/a-z/;
	if ( $names{CHR_file} =~ /\d/ ) {
		print CHR_FILE uc( $names{CHR_file} );
	}
	$all_chr_string = uc( $names{CHR_file} ), "\t";
}
if ( $names{CHR} ne "" ) {
	$names{CHR} =~ tr/A-Z/a-z/;
	$names{CHR} =~ s/^\n//;
	$names{CHR} =~ s/\r/\n/g;
	print CHR_FILE uc( $names{CHR} );
	$all_chr_string .= uc( $names{CHR} ), "\t";
}
close(CHR_FILE);

#check input validity and reformat chrfile
system("$reformat_chrfile $all_chr_file");

#Creat binned SNP files for SNP Classifier
system("$bin/map_choords_to_bin_indels.pl $bin_file $all_chr_file $pid");

#create index of old coords to new coords using pid_old_new_coords_map.txt
open( COORDS_MAP, "$tmp/$pid\_old_new_coords_map.txt" )
  || die("Cannot open old-new coords map file");
while (<COORDS_MAP>) {
	chomp;
	@elts       = split /\t/, $_;
	$old_coords = @elts[0];
	$new_coords = @elts[1];
	$index_old_new_coords{$new_coords} = $old_coords;
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
	$coding_info_file = "$coding_info_dir/$organism/chr$chromosome/@cols[3]";
	$chr_fasta_file   = "$coding_info_dir/$organism/chr$chromosome/@cols[4]";
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
	$coding_info_file = "$coding_info_dir/$organism/chr$chromosome/@cols[3]";
	$chr_fasta_file   = "$coding_info_dir/$organism/chr$chromosome/@cols[4]";
	$denormal_file    =
	  "$tmp/$pid\_snps_classified_chr$chromosome\_bin$bucket.denormal";
	open( DENORMAL, $denormal_file ) || die "(denormal file does not exist)";

	while (<DENORMAL>) {
		@elts        = split /\t/, $_;
		$coord       = @elts[1];
		@coord_elts  = split /\:/, $coord;
		$coord       = join( ',', @coord_elts );
		@coord_elts2 = split( '-', $coord, 2 );
		$coord       = join( ',', @coord_elts2 );
		$coord       = "$chromosome,$coord";
		#$coord = $index_old_new_coords{$binned_coord};
		#print "$binned_coord $coord<br>";
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
					$result.="\t$original_seq\t$modified_seq";
				}
				push @{ $coord_result_hash{$coord} }, $result;		#hash of arrays: multiple ENSTs per coord
                        }
                }
		elsif ( $_ !~ /subst/x ) {
                        $subst = "";
			$region = "INRON";
                        if ( $_ =~ /^(\d+).+?ENSG(\d+).+?ENST(\d+).+/) {
                                $hit_id = $1;
                                $en_transcript_id = "ENST$3";
                                $ensg = "ENSG$2";
                                $result = "$hit_id\t$ensg\t$en_transcript_id\t$subst\t$region";
                                #append seq (original and modified) to result
                                push @{ $coord_result_hash{$coord} }, $result;          #hash of arrays: multiple ENSTs per coord
                        }
                }

	}
}

#create original - modified protein files
foreach $coord (keys %coord_result_hash){
	my @elts;
	my @indel_result;
	@result_set = @{ $coord_result_hash{$coord} };
	foreach $result (@result_set){
		chomp $result;
		@elts = split /\t/,$result;
		my $hit_id = @elts[0];
		my $ensg = $elts[1];
		my $enst = $elts[2];
		my $subst = $elts[3];
		my $region = @elts[4];
		my $original_seq = $elts[5];
		my $modified_seq = $elts[6];
		my $outfile = "$tmp/$pid\_$enst\_$hit_id\_comparison.fa";
		if (-e $original_seq && -e $modified_seq){
			@indel_result = split /\t/, `$detect_indel $original_seq $modified_seq $outfile`;
			$aa_coords_original =  $indel_result[0];
			$aa_coords_modified =  $indel_result[1];
			$aa_change_original =  $indel_result[2];
			$aa_change_modified =  $indel_result[3];
			chomp $aa_coords_original,$aa_coords_modified,$aa_change_original,$aa_change_modified;	
			#print "@indel_result<br>";
			$result.="\t($aa_coords_original)-($aa_coords_modified)\t$aa_change_original-$aa_change_modified";
			#print "$result<br>";
		}
	}	
	$coord_result_hash{$coord} = [ @result_set ];
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

#print "$program_call $pid NCBI_blink $names{gi_number} $out 0 $subst_file 0 $seq_identity_filter $sequences_to_select $address > $comments 2>&1";
select(STDOUT);
$|       = 1;
$counter = -1;

print "<BR><BR>";

#build combined predictions file from individual predictions files.
$combined_siftresults_file = "$pid.combined_siftresults";
open (CLASSIFICATION_FILE,">$tmp/$pid.classification_file");
foreach $coord (keys %coord_result_hash) {
	@result_set = @{ $coord_result_hash{$coord} };
	foreach $result (@result_set){
		chomp $result;
		#print "RESULT: $result<br>";
		@elts             = split /\t/, $result;
		$comparison_file = "";
		$hit_id = @elts[0];
		$ensg = @elts[1];
		$enst = @elts[2];
		$subst = @elts[3];
		$region = @elts[4];
		$original_seq = @elts[5];
		$modified_seq = @elts[6];
		$coords_change = @elts[7];
		$aa_change = @elts[8];
		unless($original_seq eq ""){
			$comparison_file = "$tmp\/$pid\_$enst\_$hit_id\_comparison.fa";
		}
		print CLASSIFICATION_FILE "$coord\t$hit_id\t$ensg\t$enst\t$subst\t$coords_change\t$aa_change\t$region\t$comparison_file\n";
	}
}
close(CLASSIFICATION_FILE);

#read from combined file and print table
open( CLASSIFICATION_FILE, "$tmp/$pid.classification_file" )
  || die("Cannot open classification file");

$heading =
"Coordinates\tGene ID\tTranscript ID\tSubstitution Type\tRegion\tCoords change\tAmino acid change\tProtein Sequence Change\tRepeat Info\tRepeat Range";
open( OUTFILETABLE, ">$tmp/$combined_siftresults_file.predictions.table.html" );
open( OUTFILETSV,   ">$tmp/$combined_siftresults_file.predictions.tsv" );
print OUTFILETSV "$heading\n";
print OUTFILETABLE '<table border="1" cellspacing="0" cellpadding="4">';
print OUTFILETABLE '<thead>';
print OUTFILETABLE '<tr><th>';
$heading =~ s?\t?</th><th>?g;
print OUTFILETABLE "$heading</th></tr>\n";
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
	$aa_change = @elts[6];
	$region = @elts[7];
	$comparison_file = @elts[8];
	$table_row =
"$coord\t$ensg\t$enst\t$subst\t$region\t$coords_change\t$aa_change\t$comparison_file\n";

		print OUTFILETSV "$table_row";
		print OUTFILETABLE "<tr>\n";
		my @fields = split( '\t', $table_row );
		for $cell (@fields) {
			if ($cell =~ /comparison/i){
				print OUTFILETABLE "<td><a href=\"/sift-bin/catfile.csh?$cell+Alignment+PRE\" target=_blank> $enst change</td>";
			}
			elsif ($cell =~ /([a-z]|[A-Z])+\-([a-z]|[A-Z])+/){
				$cell =~ s/[a-z]+/<font color=red>$&<\/font>/g;
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
	

}
print OUTFILETABLE "</tr>\n</tbody>\n</table>\n<BR>";
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


