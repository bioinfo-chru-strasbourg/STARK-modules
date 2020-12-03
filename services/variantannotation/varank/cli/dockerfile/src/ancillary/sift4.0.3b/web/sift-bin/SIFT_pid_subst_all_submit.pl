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
  DBI->connect(
        "dbi:SQLite:dbname=/opt/www/sift/db/1.4/SIFT_DB_1.4.sqlite",
        "", "", { RaiseError => 1, AutoCommit => 1 } );
$sth_aln_pid_all =
  $db->prepare(
        "select address_alignment_all from pid_address where pid = ?");

$sth_aln_pid_best =
  $db->prepare(
        "select address_alignment_best from pid_address where pid = ?");

$sth_aln_gi_all =
  $db->prepare(
        "select address_alignment_all from pid_address where pid = (select pid from pid_gi where gi = ?)");

$sth_aln_gi_best =
  $db->prepare(
        "select address_alignment_best from pid_address where pid = (select pid from pid_gi where gi = ?)");

$sth_pid_subst =
  $db->prepare(
        "select substitution from dbsnp where pid = ?");
$sth_gi_subst =
  $db->prepare(
        "SELECT substitution FROM dbsnp where pid  = (select pid from pid_gi where gi = ?)");

$| = 1;
require 'SIFT_subroutines.pm';

#       Set file permissions to rw-rw----
system("umask 006");
my $bin            = "/opt/www/sift/htdocs/sift-bin";
my $tmp            = "/opt/www/sift/tmp";
my $pid            = $$;
my $datadir = "/opt/www/sift/data";
$program_call   = $bin . "/SIFTING2_for_subst_only.csh";
$return_address = "sift\@fhcrc.org";

# output the beginning text to be used on all pages
print "Content-type: text/html\n\n";
print "<body bgcolor=white>\n";

if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
	print "This script should be referenced with a METHOD of POST\n";
	exit;
}
#print "<A NAME=top><H1>Batch Protein Prediction</H1></A><BR>\n";

print "If your browser times out before results are shown, your results will be stored at  http://sift.jcvi.org/tmp/$pid\.combined_siftresults.predictions.table.html  for 1 hour.<BR>";

read( STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"} );
%names = &general_parse( $ENV{"CONTENT_TYPE"}, $QUERY_STRING );
my $address;
chomp($names{address});
$names{address} =~ s/\s+//;
if ( $names{address} ne "" ) {
        $address = $names{address};
}
#foreach (keys %names) {
#        print "The key $_ contains $names{$_}\n";
#}


### gi file####
my $all_gi_file = $tmp . "/$pid.gifile";
if ( $names{GI} eq "" && $names{GI_file} eq "" ) {
        print "<H1>Error</H1> Please enter some RefSeq accessions of GI numbers.<P>\n";
        exit;
}
open( GI_FILE, ">$all_gi_file" );
my $all_gi_string;
if ( $names{GI_file} ne "" ) {
        $names{GI_file} =~ s/\r/\n/g;
        $names{GI_file} =~ tr/A-Z/a-z/;
        print GI_FILE uc ($names{GI_file});
        $all_gi_string =uc( $names{GI_file}),"\t";
}
if ( $names{GI} ne "" ) {
        $names{GI} =~ tr/A-Z/a-z/;
        $names{GI} =~ s/\r/\n/g;
        print GI_FILE uc($names{GI});
        $all_gi_string .= uc($names{GI}),"\t";
}
close(GI_FILE);
my @fields = split( /[\t\n]/, $all_gi_string );
foreach $gi_line (@fields){
        chomp;
        if ($gi_line =~ /application/i || $gi_line !~ /\d/){
                next;
        }
        else{
                $gi_line =~ s/^\s+//g;
                $gi_line =~ s/\s+$//g;
                $gi_line =~ s/^\t+//g;
                $gi_line =~ s/\t+$//g;
		if ($gi_line !~ /[,]/){
			$gi_line=$gi_line."\,DBSNP";
		}
		push @gi_subst_list,$gi_line;
	}
}




### substitution files #######
foreach $gi_line (@gi_subst_list){
	@elts = split /\,/,$gi_line;
        $gi = @elts[0];
        $gi =~ s/^\s+//g;
        $gi =~ s/\s+$//g;
 	$gi_string.="\'$gi\',";
}
chop $gi_string;
$sth_pid_rsid =
  $db->prepare("select pid,substitution,rsid from dbsnp where pid in ($gi_string)");

$sth_gi_rsid =
  $db->prepare(
"select a.gi, b.substitution, b.rsid from dbsnp b, pid_gi a where b.pid in (select pid from pid_gi where gi in ($gi_string) ) and a.pid = b.pid"
  );
$sth_pid_rsid->execute();

while (@rows = $sth_pid_rsid->fetchrow_array()){
        $prot_id = @rows[0];
        $sub = @rows[1];
        $rsid = @rows[2];
        $pid_sub_rsid_hash{"$prot_id\t$sub"} = $rsid;
	if ($seen{"$prot_id\t$sub"} == 1){
		next;
	}
	else{	
		$pid_sub_string = $pid_sub_hash{$prot_id};
		$pid_sub_string.="$sub,";
		$pid_sub_hash{$prot_id} = $pid_sub_string;
		$seen{"$prot_id\t$sub"} = 1;
	}
}

$sth_gi_rsid->execute();
while (@rows = $sth_gi_rsid->fetchrow_array()){
        $prot_id = @rows[0];
        $sub = @rows[1];
        $rsid = @rows[2];
        $pid_sub_rsid_hash{"$prot_id\t$sub"} = $rsid;
	if ($seen{"$prot_id\t$sub"} == 1){
                next;
        }
        else{
		$pid_sub_string = $pid_sub_hash{$prot_id};
                $pid_sub_string.="$sub,";
                $pid_sub_hash{$prot_id} = $pid_sub_string;
                $seen{"$prot_id\t$sub"} = 1;
        }

}




my %gi_subst_file_hash;
foreach $gi_line (@gi_subst_list){
	@elts = split /\,/,$gi_line;
	$gi = @elts[0];
	$gi =~ s/^\s+//g;
	$gi =~ s/\s+$//g;
	push @gi_list,$gi;
	open (SUBST_FILE,">$tmp/$pid.$gi.substfile");
	for ($i = 1; $i < scalar @elts; $i++){
		@elts[$i] =~ s/^\s+//g;
 		@elts[$i] =~ s/\s+$//g;
		if (@elts[$i] =~ /dbsnp/i){
                       	$subst_str = $pid_sub_hash{$gi};
			chop $subst_str;
			@sub_arr = split /\,/, $subst_str;
			foreach $sub (@sub_arr){
				if ($subst_seen_hash{$sub} == 1){
					next;
				}
				else{
	                                print SUBST_FILE "$sub\n";
        	                        $subst_seen_hash{$sub} = 1;
                	        }

			}
			next;
		}
		
		if ($subst_seen_hash{@elts[$i]} == 1){
                	next;
                }
		else{
			print SUBST_FILE @elts[$i],"\n";
			$subst_seen_hash{@elts[$i]} = 1;
		}
	}
	close (SUBST_FILE);	
	$gi_subst_file_hash{$gi} = "$pid.$gi.substfile";
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
$| = 1;

foreach $gi (@gi_list){
	$subst_file = $gi_subst_file_hash{$gi};
	chomp;
	if ($gi !~ /XP/i && $gi !~ /NP/i){	#then its a gi number not a refseq id
		$sth_aln_best = $sth_aln_gi_best;	
		$sth_aln_all = $sth_aln_gi_all;
	}
	else{
		$sth_aln_best = $sth_aln_pid_best;
                $sth_aln_all = $sth_aln_pid_all;
	}
	$address_aln_best = "";
	$address_aln_all = "";
	if ($names{sequences_to_select} =~ /BEST/){
	
		$sth_aln_best->execute($gi);
        	$address_aln_best = "$datadir\/".$sth_aln_best->fetchrow_array();
		if ( $address_aln_best =~ /.+\/(.+)alignedfasta/){
#print "entered for gi $gi in here $address_aln_best\n<BR>";
			$giprefix = $1;
			$giprefix =~ s/\..+//g;
                        $alignedfasta = "$tmp\/$pid\.$gi.alignedfasta";
			system("cp $address_aln_best $alignedfasta");
                }
                else{
                        $alignedfasta = "";
			$unalignedfasta =  "$tmp\/$pid\.$gi.unaligned";
                }

	}
	else{
                $sth_aln_all->execute($gi);
                $address_aln_all = "$datadir\/".$sth_aln_all->fetchrow_array();
                $aln_all_hash{$gi} = $address_aln_all;
		if ( $address_aln_all =~ /.+\/(.+)alignedfasta/){
			$giprefix = $1;
			$giprefix =~ s/\..+//g;
			$alignedfasta = "$tmp\/$pid\.$gi.alignedfasta";
		}
		else{
			$unalignedfasta =  "$tmp/$pid\.$gi.unaligned";
			$alignedfasta = "";
		}
		system("cp $address_aln_all $alignedfasta");
	}
	$comments            = "$tmp/$pid.$gi.comments";
	$out                 = $tmp . "/$$.$gi.siftresults";

	if ($alignedfasta ne ""){ #aligned sequences obtained from databse so use option = alignedseqs for SIFTING2.csh
#print "entered alined fast and this exists $alignedfasta\n";
		system(
		   "$program_call $pid alignedseqs $alignedfasta $out 0 $subst_file 0 $seq_identity_filter $sequences_to_select $address > $comments 2>&1"
		);
		print "</PRE>";

		# SIFT successful -- OUTPUT RESULTS
		$commentscsh = $tmp . "/$pid.commentscsh";
		#$psiblastout = $tmp . "/$pid" . ".alignedfasta";
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
	
	else{	#aligned sequences could not be obtained from database so use option = NCBI_blink for SIFTING2.csh
#		print "$program_call $pid\.$gi NCBI_blink $gi $out 0 $tmp/$pid.$gi.substfile 0 $seq_identity_filter $sequences_to_select $address > $comments 2>&1 <BR>\n";
# PCN Pauline 01-29-2010
# this works, but I think it's going to saturate the server if we let
# people look up 1K sequences at NCBI. so just print out error message
	print "The sequence $gi is not found in our precomputed database. Please use <A HREF=\"http://sift.jcvi.org/www/SIFT_BLink_submit.html\">SIFT BLink</A HREF> for $gi.<BR>\n"; 
#		system(
#                   "$program_call $pid\.$gi NCBI_blink $gi $out 0 $tmp/$pid.$gi.substfile 0 $seq_identity_filter $sequences_to_select $address > $comments 2>&1" 
#                );

	}
	
}
$combined_siftresults_file = "$pid.combined_siftresults";
foreach $gi (@gi_list){
	system ("echo 'Protein Identifier: <b><u>$gi</u></b>' >> $tmp/$combined_siftresults_file.predictions");
	system ("cat $tmp/$pid.$gi.siftresults.predictions >> $tmp/$combined_siftresults_file.predictions");
}


open( FILE, "$tmp/$combined_siftresults_file.predictions" ) || die("Cannot open predictions file");

$heading =
"Protein ID\tSubstitution\tdbSNP ID\tPrediction\tScore\tMedian Info\tNumber of Seqs at position";
open( OUTFILETABLE, ">$tmp/$combined_siftresults_file.predictions.table.html" );
open( OUTFILETSV, ">$tmp/$pid\_predictions.tsv" );
print OUTFILETSV "$heading\n";
print OUTFILETABLE '<table border="1" cellspacing="0" cellpadding="4">';
print OUTFILETABLE '<thead>';
print OUTFILETABLE '<tr><th>';
$heading =~ s?\t?</th><th>?g;
print OUTFILETABLE "$heading</th></tr>\n";
my $warningflag = 0;
while (<FILE>) {
	chomp;
	if ( $_ =~ /Protein Identifier.+?\<u\>(.+?)\<\/u\>.+?/i ) {
		$proteinid = $1;
	}

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
		else {
			$pred = "DAMAGING";
		}
		$rsid = $pid_sub_rsid_hash{"$proteinid\t$subst"};
	}
	if ( $_ =~ /Median Sequence conservation: (.+)/i ) {		
		$median = $1;
		if ($median > 3.25 && $pred =~/DAMAGING/i){
			$warningflag = 1;
			$pred = $pred." *Warning! Low confidence.";
		}
	}
	if ( $_ =~ /Sequences represented at this position:(.+)/i ) {
		$numseqs   = $1;
		$table_row =
		  "$proteinid\t$subst\t$rsid\t$pred\t$score\t$median\t$numseqs";
		
		print OUTFILETSV "$table_row\n";
		print OUTFILETABLE "<tr>\n";
		my @fields = split('\t',$table_row);
		for $cell (@fields) {	
			if ($cell =~ /DAMAGING/ || $cell =~ /Warning/){
				print OUTFILETABLE "<td><font color=red>$cell</font></td>";		
			}	
			else{
				print OUTFILETABLE "<td>$cell</td>";		
			}	
				
		}
		print OUTFILETABLE "</tr>\n";
	}

}
print OUTFILETABLE "</tr>\n</tbody>\n</table>\n<BR>";
if ($warningflag ==1){

	print OUTFILETABLE "<font color=red>* Low confidence means that the protein alignment does not have enough sequence diversity. Because the position artifically appears to be conserved, an amino acid may incorrectly predicted to be damaging.</font><BR><BR>";
}
close(FILE);
close(OUTFILETABLE);
close(OUTFILETSV);
#print the table
open (OUTFILETABLE,"$tmp/$combined_siftresults_file.predictions.table.html") || die ("Cannot open predictions table");
while (<OUTFILETABLE>){
	print;	
}
print "Click <A HREF=\"\/tmp\/$pid\_predictions.tsv\">here</A> to download the following table in tab separated format. You can open it using excel with delimiter set as TAB";
print "<BR>";

print "<BR>";
print
   "<i>If you received a warning that the sequences were not diverse enough, you can have SIFT choose more diverse sequences <A HREF=\"/www/SIFT_seq_submit2.html\">here.</A></i><BR><BR>";

#email the results
if ($address ne ""){
	open (MESSAGE, ">$tmp/$pid.email_message.txt");
	print MESSAGE "Dear User\n\nThank you for using SIFT.\n\nPlease find the results of your recent query attached with this message.\nRemember this job id \"$pid\" for any future correspondance.\nDo not hesitate to contact us if you have any questions about SIFT.\n\nThanks\nSIFT Team\nJ Craig Venter Institute (West Coast Campus)\n10355 Science Center Drive\nSan Diego, CA 92121\nUSA";
	close (MESSAGE);

	system("\/usr\/bin\/mutt -F \/opt\/www\/sift\/htdocs\/.muttrc -a $tmp\/$pid\_predictions.tsv -s \"SIFT Results for Job ID $pid\" $address \< $tmp\/$pid.email_message.txt"); 
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

