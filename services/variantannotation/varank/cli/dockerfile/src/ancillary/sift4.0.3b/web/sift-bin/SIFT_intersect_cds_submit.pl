#!/usr/local/bin/perl
#
#-------------------------------------------------------------------------

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use Tie::IxHash;
#use strict;

$| = 1;
require '/opt/www/sift/htdocs/sift-bin/SIFT_subroutines.pm';

#       Set file permissions to rw-rw----
system("umask 006");
my $bin             = "/opt/www/sift/htdocs/sift-bin";
my $tmp             = "/opt/www/sift/tmp";
my $pid             = $$;
my $coding_info_dir = "/opt/www/sift/coding_info";
my $coding_file;

print "Content-type: text/html\n\n";
print "<body bgcolor=white>\n";
my $outpage_url = "/sift-bin/catfile.csh?$tmp/$pid.outpage.html";
if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
        print "This script should be referenced with a METHOD of POST\n";
        exit;
}
read( STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"} );
%names = general_parse( $ENV{"CONTENT_TYPE"}, $QUERY_STRING );
my $address;
chomp( $names{address} );
$names{address} =~ s/\s+//;
if ( $names{address} ne "" ) {
        $address = $names{address};
}
## Check for validity of user inputs
my $all_chr_file =  $tmp . "/" . "$pid.allchrfile";
#print "all chr file called $all_chr_file\n";
#print "iiiii$organism<br>";
if ( $names{CHR} eq "" && $names{CHR_file} eq "" ) {
        print
"<H1>Error</H1> Please enter some chromosome coordinates with substitutions.<P>\n";
        exit;
}

my $organism = $names{organism};
$organism =~ s/^\s+//;
$organism =~ s/\s+$//;
if ( $organism =~ /Select Organism/i ) {
        print
"<H1>Error</H1> Please select organism after pressing back button on your browser.<P>\n";
        exit;
} elsif ($organism =~ /Homo_sapiens36/) { 
	$coding_file = $coding_info_dir . "/Homo_sapien/ens.hum.ncbi36.ver41.cds.merge.gff";
} elsif ($organism =~ /Homo_sapiens37/) {
	$coding_file = $coding_info_dir . "/Homo_sapien_37/ens.hum.ncbi37.ver55.cds.merge.gff";
}

#Read input list of chromosome coordinates and add to all_chr_string
my $all_chr_string;
if( $names{CHR_file} !~ /\d/){
        $names{CHR_file} = "";
}
if( $names{CHR} !~ /\d/){
        $names{CHR} = "";
}

if ($names{CHR_file} ne "" && $names{CHR} ne ""){
        print "<H1>Error</H1> Please choose only one of the following input methods after clicking the back button on your browser<BR>";
        print "1. Paste the input in the relevant textbox<BR>";
        print "2. Upload the text file containing input data";
        exit;
}
if ($names{CHR_file} eq "" && $names{CHR} eq ""){
        print "<H1>Error</H1> Please choose one of the following input methods after clicking the back button on your browser<BR>";
        print "1. Paste the input in the relevant textbox<BR>";
        print "2. Upload the text file containing input data";
        exit;
}

open( CHR_FILE, ">$all_chr_file" );
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
#print "finished reading $all_chr_file\n";
open (CHR_FILE,"$all_chr_file") || die ("Cannot open all chr file for validation");
#print "here 1234\n";
while (<CHR_FILE>){
        if ($_ =~ /\d/ && $_ =~ /\,/){
                $first_line = $_;
                last;
        }
}
close(CHR_FILE);
if ($first_line =~ /\d+,\d+,\d+,\-?1,/i){
        $COORD_SYSTEM = "SPACE";
}
elsif($first_line =~ /\d+,\d+,\-?1,/i){
        $COORD_SYSTEM = "RESIDUE";
}
else{
        print "<H1>Error</H1> Incorrect input format. Please see <A HREF=\"\/www\/chr_coords_example.html\">Sample format</A>";
        last;
}
print "<A NAME=top><H1>Variants In Coding</H1></A><BR>\n";

print
"Your input data has been recognized to use <A HREF=\/www\/chr_coords_example.html target=_blank>$COORD_SYSTEM based coordinate system</A>.<BR>\n";
print "If your browser times out before results are shown, your results will be stored at  http://sift.jcvi.org/tmp/$pid\.allchrfile.cds.txt for 1 hour.<BR>"; 
# if COORD SYSTEM is SPACE, then there are insertions with same coordinates

#if COORD SYSTEM is ABSOLUTE then convert chr file to space based.
my $new_chr_file = $all_chr_file . ".gff";
if ($COORD_SYSTEM eq "RESIDUE"){
	convert_residue_to_space_gff ($all_chr_file, $new_chr_file);
}
if ($COORD_SYSTEM eq "SPACE") {
	convert_space_with_insertions_gff ($all_chr_file, $new_chr_file);
}
	

 system ("sort -k1,1 -k4,4n -k5,5n $new_chr_file > $all_chr_file.sorted");
my @lines_to_get = `$bin/IntersectLocations.sh $all_chr_file.sorted gff $coding_file gff simple | grep -v Total | cut -f1 | uniq`;

my $all_chr_in_cds_file = $all_chr_file . ".cds.txt";

my ($cds_count, $non_cds_count) = get_lines ($new_chr_file, 
			$all_chr_in_cds_file, 
			@lines_to_get);

open (CHR_FILE,"$all_chr_in_cds_file") || die ("Cannot open all chr file");
open (OUTFILE,">$tmp/$new_pid.chrfile");
#print "<BR> Path for all_chr_in_cds_file $all_chr_in_cds_file <BR>";
print "<BR><BR>Coding: $cds_count\n<BR>Noncoding: $non_cds_count\n<BR><BR>\n";
print "1.  <b>Save</b> <A HREF=\"\/tmp\/$pid.allchrfile.cds.txt\">Coding Variants File</A><BR><BR>\n";
print "2.  Go to <A HREF=\"/www/SIFT_chr_coords_indels_submit.html\">SIFT Indels</A> or <A HREF=\"/www/SIFT_chr_coords_submit.html\">SIFT Genome</A> and upload coding variants file.<BR><BR>";

print "<b><u>List of Coding Variants</u></b><BR>\n";
while (<CHR_FILE>){
#        chomp;
	print "$_<BR>\n";
}
close(CHR_FILE);
close(OUTFILE);


print "<BR><BR> Problems? Please email us <A HREF=\"contact.pl\">us<A> the details of your submission.\n";
print "<BR>Your job id is $pid.\n";

#email the results
if ( $address ne "" ) {
	open( MESSAGE, ">$tmp/$pid.email_message.txt" );
	print MESSAGE
"Dear User\n\nThank you for using SIFT.\n\nPlease find the results of your recent query attached with this message.\nRemember this job id \"$pid\" for any future correspondance.\nDo not hesitate to contact us if you have any questions about SIFT.\n\nThanks\nSIFT Team\nJ Craig Venter Institute (West Coast Campus)\n10355 Science Center Drive\nSan Diego, CA 92121\nUSA";
	close(MESSAGE);
	system(
"mutt -F /opt/www/sift/htdocs/.muttrc -a $all_chr_in_cds_file -s \"SIFT Results for Job ID $pid\" $address <$tmp/$pid.email_message.txt"
	);
}

exit(0);

sub
get_lines
{
	my  ($all_chr_file, $all_chr_in_cds_file, @lines_to_get)  = @_;

	my %line_to_get_hash ;
	# more elegant way to do this by iterating through @lines_to_get, but oh well
	for (my $i=0; $i < @lines_to_get; $i++) {
		chomp ($lines_to_get[$i]);
		$lines_to_get_hash{$lines_to_get[$i]} = 1;
	}
	open (CHR_FILE, $all_chr_file) || die "can't open $all_chr_file";
	open (OUT_CHR_FILE, ">$all_chr_in_cds_file") || die "can't open $all_chr_in_cds_file";
	my $in_cds = 0;
	my $not_in_cds = 0;
        my $line;
	while ($line = <CHR_FILE>){
		chomp ($line);
		my @fields = split (/\t/, $line);
		my $line_num = $fields[1];
		if (exists ($lines_to_get_hash{$line_num})) {
			print OUT_CHR_FILE $fields[8] . "\n";
			$in_cds++;
		} else { 
			$not_in_cds++;
		}
	}
	close (CHR_FILE);
	close (OUT_CHR_FILE);
	return ($in_cds, $not_in_cds);
}

sub
convert_space_with_insertions_gff
{
       my ($all_chr_file, $new_chr_file) =  @_;

        open (CHR_FILE,"$all_chr_file") || die ("Cannot open all chr file");
        open (CHR_FILE_NEW,">$new_chr_file") || die ("Cannot open new all chr file");
        my $line_num = 0;
        my $line;
	while ($line = <CHR_FILE>){
                chomp ($line);
                if ($line !~ /\d/){next}
                @elts = split (/\,/, $line);
                $chr = @elts[0];
                $coord1 = @elts[1];
                $coord2 = @elts[2]; 
                $orn = @elts[3];
                $alleles = @elts[4];
                $comment = @elts[5];
		# for insertions
		if ($coord1 == $coord2) {
			$coord1 = $coord1 - 1;
		}
                print CHR_FILE_NEW "$chr\t$line_num\t$line_num\t$coord1\t$coord2\t.\t+\t.\t$line\n"; 
                $line_num++;
        }
        close(CHR_FILE);
        close (CHR_FILE_NEW);

} # end convert_space_with_insertions_gff

sub
convert_residue_to_space_gff
{
	my ($all_chr_file, $new_chr_file) =  @_;

   	open (CHR_FILE,"$all_chr_file") || die ("Cannot open all chr file");
        open (CHR_FILE_NEW,">$new_chr_file") || die ("Cannot open new all chr file");
	my $line_num = 0;
        my $line;
	while ($line = <CHR_FILE>){
                chomp ($line);
                if ($line !~ /\d/){next}
                @elts = split (/\,/, $line);
                $chr = @elts[0];
                $coord2 = @elts[1];
                $coord1 = $coord2-1;
                $orn = @elts[2];
                $alleles = @elts[3];
                $comment = @elts[4];
                print CHR_FILE_NEW "$chr\t$line_num\t$line_num\t$coord1\t$coord2\t.\t+\t.\t$line\n"; 
		$line_num++;
        }
        close(CHR_FILE);
        close (CHR_FILE_NEW);

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
check_ip_counts
{
	## Check that this IP address hasn't been used too much
	my $IP_address = $ENV{REMOTE_ADDR};

	my $remote_host = $ENV{REMOTE_HOST}; 
	my $ip_counts =
		`cat  /home/blocks/apache/logs/access_log  | grep POST | grep $IP_address | wc -l `;
	chomp($ip_counts);
	if ( $ip_counts == "" ) {
       	 $ip_counts =
		`cat /home/blocks/apache/logs/access_log  | grep POST | grep $remote_host | wc -l`;
        chomp($ip_counts);
	}
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

} # end check_counts 

sub print_outpage{
        open(OUTPAGE, "$tmp/$pid.outpage.html") || die("cannot open outpage");
        while (<OUTPAGE>){
                print;
        }
        close (OUTPAGE);
}




