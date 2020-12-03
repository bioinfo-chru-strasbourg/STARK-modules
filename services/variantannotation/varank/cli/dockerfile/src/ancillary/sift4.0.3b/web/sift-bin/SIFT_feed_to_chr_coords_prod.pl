#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use List::Util qw[min max];

$| = 1;
require 'SIFT_subroutines.pm';
system("umask 006");
my $bin             = "/opt/www/sift/htdocs/sift-bin";
my $tmp             = "/opt/www/sift/tmp";
my $pid             = $$;
my $num_coords_per_split = 1000;
#my $tmp             = "/opt/www/sift/tmp/$pid";
#mkdir $tmp;
#system("chmod 777 * -R $tmp");
# output the beginning text to be used on all pages
print "Content-type: text/html\n\n";
print "<body bgcolor=white>\n";
my $outpage_url = "/sift-bin/catfile.csh?$tmp/$pid.outpage.html";
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
## Check for validity of user inputs
my $all_chr_file = $tmp . "/$pid.allchrfile";
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
}
$seq_identity_filter = "90";

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
open (CHR_FILE,"$all_chr_file") || die ("Cannot open all chr file for validation");
while (<CHR_FILE>){
	if ($_ =~ /\d/ && $_ =~ /\,/){
		$first_line = $_;
		last;
	}
}
close(CHR_FILE);
if ($first_line =~ /\d+,\d+,\d+,-?1,[A,T,G,C]\/[A,T,G,C]/i){
        $COORD_SYSTEM = "SPACE";
}
elsif($first_line =~ /\d+,\d+,-?1,[A,T,G,C]\/[A,T,G,C]/i){
        $COORD_SYSTEM = "RESIDUE";
}
else{
        print "<H1>Error</H1> Incorrect input format. Please see <A HREF=\"\/www\/chr_coords_example.html\">Sample format</A>";
        last;
}

print "<A NAME=top><H1><center>S<font color=red>I</font>F<font
color=blue>T</font> Genome</center></H1></A><BR>\n";

print
"Your input data has been recognized to use <A HREF=\/www\/chr_coords_example.html target=_blank>$COORD_SYSTEM based coordinate system</A>. Your job id is $pid and is currently running.  Your job has been partitioned into datasets of $num_coords_per_split positions and the status of each job can be viewed in the <A HREF=$outpage_url target=\"_blank\">SIFT results status page<A>. Once the status of a job is <font color=green>'Complete'</font>, you may view or download the results. A partitioned job typically takes 6-7 min to complete.  <BR><BR>Proceed to <A HREF=$outpage_url target=\"_blank\">SIFT results status page.<A><BR><BR> Problems? Contact <A HREF=\"contact.pl\">us<A> with your job id.\n";

# Close the I/O handles
#close(STDin);
close(STDOUT);
#close(STDERR); 

#prepare output status page
open (OUTPAGE,">>$tmp/$pid.outpage.html") || die ("cannot open outpage");
print OUTPAGE "<A NAME=top><H1><center>S<font color=red>I</font>F<font
color=blue>T</font> Results Status</center></H1></A><BR>\n";
print OUTPAGE "<BR><BR>";
print OUTPAGE "<BR><a href=\"$outpage_url\">Refresh page</a>";
print OUTPAGE '<table border="1" cellspacing="2" cellpadding="14">';
print OUTPAGE '<thead>';
print OUTPAGE '<tr><th>';
$heading = "<thead><tr><th>Job</th><th>Job size</th><th>Job ID</th><th>Job status</th><th>View results</th><th>Download results</th></tr>";
print OUTPAGE "$heading</th></tr>\n";


#if COORD SYSTEM is ABSOLUTE then convert chr file to space based.
if ($COORD_SYSTEM eq "RESIDUE"){
	open (CHR_FILE,"$all_chr_file") || die ("Cannot open all chr file");
	open (CHR_FILE_NEW,">$all_chr_file.new") || die ("Cannot open new all chr file");
	while (<CHR_FILE>){
		chomp;
		if ($_ !~ /\d/){next}
		@elts = split /\,/, $_;
		$chr = @elts[0];
		$coord2 = @elts[1];
		$coord1 = $coord2-1;
		$orn = @elts[2];
		$alleles = @elts[3];
		$comment = @elts[4];
		print CHR_FILE_NEW "$chr,$coord1,$coord2,$orn,$alleles,$comment\n";
	}
	close(CHR_FILE);
	close (CHR_FILE_NEW);
	system("mv $all_chr_file.new $all_chr_file");
}


open (CHR_FILE,"$all_chr_file") || die ("Cannot open all chr file");
$count = 0;
$new_pid = $pid+1;
$input_set = 1;
open (OUTFILE,">$tmp/$new_pid.chrfile");

while (<CHR_FILE>){
	chomp;
	if ($_ !~ /\d/){next}
	$count++;
	if ($count % $num_coords_per_split== 0){
		push @pid_list,$new_pid;
		print OUTFILE "$_\n";
                $job_size = $num_coords_per_split;
		$start = $count-$job_size+1;
		$end = $count ;
                print OUTPAGE "<tr><td>Partitioned set $input_set</td>";
                print OUTPAGE "<td>Input rows $start to $end</td>";
                print OUTPAGE "<td>$new_pid</td>";
                print OUTPAGE "<td>Not started.</td><td>Not available</td></td><td>Not available</tr>";
                print OUTPAGE "</tr>\n";
		close(OUTFILE);
		$input_set++;
                $new_pid++;

		open (OUTFILE,">$tmp/$new_pid.chrfile");
	}
	else{
		print OUTFILE "$_\n";
		
	}
}
close(CHR_FILE);
close(OUTFILE);
$job_size = $count % $num_coords_per_split;
$start = $end+1;
$end = $start + $job_size-1;
if ($job_size == 0){
	system("rm -f $tmp/$new_pid.chrfile");
	print OUTPAGE "<tr BGCOLOR=\"\#DCDCDC\"><td>Complete set</td>";
        print OUTPAGE "<td>Input rows 1 to $end</td>";
        print OUTPAGE "<td>$pid</td>";
        print OUTPAGE "<td>Not started.</td></td><td>Not available</td></td><td>Not available</tr>";
        print OUTPAGE "</tr></table>\n";

}
else{
	push @pid_list, $new_pid;
	print OUTPAGE "<tr><td>Partitioned set $input_set</td>";
	print OUTPAGE "<td>Input rows $start to $end</td>";
	print OUTPAGE "<td>$new_pid</td>";
	print OUTPAGE "<td>Not started.</td></td><td>Not available</td></td><td>Not available</tr>";
	print OUTPAGE "</tr>\n";
	print OUTPAGE "<tr BGCOLOR=\"\#DCDCDC\"><td>Complete set</td>";
        print OUTPAGE "<td>Input rows 1 to $end</td>";
        print OUTPAGE "<td>$pid</td>";
        print OUTPAGE "<td>Not started.</td></td><td>Not available</td></td><td>Not available</tr>";
        print OUTPAGE "</tr></table><BR><BR>\n";

}
print OUTPAGE "<BR><BR>\n";
print OUTPAGE "<b>Batch Report</b><BR><BR>\n";
print OUTPAGE "Number of coding variants: <BR>\n";
print OUTPAGE "Coding variants predicted: <BR>\n";
print OUTPAGE "Tolerated: <BR>\n";
print OUTPAGE "Damaging: <BR>\n";
print OUTPAGE "Nonsynonymous: <BR>\n";
print OUTPAGE "Synonymous: <BR>\n";
print OUTPAGE "Novel: <BR>\n";

close (OUTPAGE);
#system("cp $tmp/$pid.outpage.html $tmp/$pid.outpage.swap.html");
#direct to the outpage url
#print "Location: $outpage_url\n\n";
#print_outpage();

open (BATCH_FILE,">$tmp/$pid.batchfile");
for ($i = 0; $i < scalar @pid_list; $i++){
	$new_pid = @pid_list[$i];
	chomp $new_pid;
	if ($i == scalar @pid_list -1){
		print BATCH_FILE "$pid\t$new_pid\t$tmp/$new_pid.chrfile\t$organism\t$seq_identity_filter\t$COORD_SYSTEM\tLAST\t$address\n";
	}
	else{
		print BATCH_FILE "$pid\t$new_pid\t$tmp/$new_pid.chrfile\t$organism\t$seq_identity_filter\t$COORD_SYSTEM\tNOT_LAST\t$address\n";
	}
}
system("$bin/SIFT_feed_to_chr_coords_batch.pl $tmp/$pid.batchfile");
close (BATCH_FILE);













sub print_outpage{
	open(OUTPAGE, "$tmp/$pid.outpage.html") || die("cannot open outpage");
	while (<OUTPAGE>){
        	print;
	}
	close (OUTPAGE);
}








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

