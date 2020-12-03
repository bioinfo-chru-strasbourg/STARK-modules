#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

my $tmp = "/opt/www/sift/tmp";
my $pid = "2566";


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

my @predictions_files = `ls $tmp/$pid\.*.siftresults.predictions`;

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

