#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
use List::Util qw[min max];
use Template;

$| = 1;
require 'SIFT_subroutines.pm';
system("umask 006");
my $bin             = "/opt/www/sift/htdocs/sift-bin";
my $tmp             = "/opt/www/sift/tmp";
my $pid             = $$;
my $num_coords_per_split = 1000;
my $content = "";

print "Content-type: text/html\n\n";

my $outpage_url = "/sift-bin/catfile.csh?$tmp/$pid.outpage.html";

if($ENV{"REQUEST_METHOD"} ne "POST") {
  $content = "This script should be referenced with a METHOD of POST\n";
  &finish_script($content, -1);
} # end if($ENV{"REQUEST_METHOD"} ne "POST")

read( STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"} );
%names = general_parse( $ENV{"CONTENT_TYPE"}, $QUERY_STRING );

my $address;

chomp($names{address});

$names{address} =~ s/\s+//;

if($names{address} ne "") {
  $address = $names{address};
} # end if($names{address} ne "")

## Check for validity of user inputs
my $all_chr_file = $tmp . "/$pid.allchrfile";
#print "iiiii$organism<br>";

if($names{CHR} eq "" && $names{CHR_file} eq "") {
  $content .= "<p class=\"header1\">Error</p>\n";
  $content .= "<p>Please enter some chromosome coordinates with substitutions.</p>\n";
  &finish_script($content, -1);
} # end if($names{CHR} eq "" && $names{CHR_file} eq "")

my $organism = $names{organism};
$organism =~ s/^\s+//;
$organism =~ s/\s+$//;

if($organism =~ /Select Organism/i) {
  $content .= "<p class=\"header1\">Error</p>\n";
  $content .= "<p>Please select organism after pressing back button on your browser.</p>\n";
  &finish_script($content, -1);
} # end if($organism =~ /Select Organism/i)

$seq_identity_filter = "90";
#Read input list of chromosome coordinates and add to all_chr_string
my $all_chr_string;

if($names{CHR_file} !~ /\d/) {
  $names{CHR_file} = "";
} # end if($names{CHR_file} !~ /\d/)

if($names{CHR} !~ /\d/) {
  $names{CHR} = "";
} # end if($names{CHR} !~ /\d/)

my $oo1 = $names{oo1}==1 ? 1 : 0;
my $oo2 = $names{oo2}==1 ? 1 : 0;
my $oo3 = $names{oo3}==1 ? 1 : 0;
my $oo4 = $names{oo4}==1 ? 1 : 0;
my $oo5 = $names{oo5}==1 ? 1 : 0;
my $oo6 = $names{oo6}==1 ? 1 : 0;
my $oo7 = $names{oo7}==1 ? 1 : 0;
my $oo8 = $names{oo8}==1 ? 1 : 0;
my $oo9 = $names{oo9}==1 ? 1 : 0;
my $oo10 = $names{oo10}==1 ? 1 : 0;
my $oo11 = $names{oo11}==1 ? 1 : 0;
my $oo12 = $names{oo12}==1 ? 1 : 0;

my $output_options = "$oo1,$oo2,$oo3,$oo4,$oo5,$oo6,$oo7,$oo8,$oo9,$oo10,$oo11,$oo12";

if($names{CHR_file} ne "" && $names{CHR} ne "") {
  $content .= "<p class=\"header1\">Error</p>\n";
  $content .= "<p>Please choose only one of the following input methods after clicking the back button on your browser</p>\n";
  $content .= "<p>1. Paste the input in the relevant textbox</p>\n";
  $content .= "<p>2. Upload the text file containing input data</p>\n";
  &finish_script($content, -1);
} # end if($names{CHR_file} ne "" && $names{CHR} ne "")

if($names{CHR_file} eq "" && $names{CHR} eq "") {
  $content .= "<p class=\"header1\">Error</p>\n";
  $content .= "<p>Please choose one of the following input methods after clicking the back button on your browser</p>\n";
  $content .= "<p>1. Paste the input in the relevant textbox</p>\n";
  $content .= "<p>2. Upload the text file containing input data</p>\n";
  &finish_script($content, -1);
} # end if($names{CHR_file} eq "" && $names{CHR} eq "")

open( CHR_FILE, ">$all_chr_file" );

if($names{CHR_file} ne "") {
  $names{CHR_file} =~ s/\r/\n/g;
  $names{CHR_file} =~ tr/A-Z/a-z/;

  if($names{CHR_file} =~ /\d/) {
    print CHR_FILE uc( $names{CHR_file} );
  } # end if($names{CHR_file} =~ /\d/)

  $all_chr_string = uc($names{CHR_file}), "\t";
  $input_method = "FILE";
} # end if($names{CHR_file} ne "")
else {
  $names{CHR} =~ tr/A-Z/a-z/;
  $names{CHR} =~ s/^\n//;
  $names{CHR} =~ s/\r/\n/g;
  print CHR_FILE uc($names{CHR});
  $all_chr_string .= uc($names{CHR}), "\t";
  $input_method = "TEXT";
} # end else

close(CHR_FILE);
open (CHR_FILE,"$all_chr_file") || die ("Cannot open all chr file for validation");

while(<CHR_FILE>) {
  if($_ =~ /\d/ && $_ =~ /\,/) {
    $first_line = $_;
    last;
  } # end if($_ =~ /\d/ && $_ =~ /\,/)
} # end while(<CHR_FILE>)

close(CHR_FILE);

if($first_line =~ /[\d+,X,x,Y,y],\d+,\d+,\-?1,[A,T,G,C]\/[A,T,G,C]/i) {
  $COORD_SYSTEM = "SPACE";
} # end if($first_line =~ /[\d+,X,x,Y,y],\d+,\d+,\-?1,[A,T,G,C]\/[A,T,G,C]/i)
elsif($first_line =~ /[\d+,x,X,y,Y],\d+,\-?1,[A,T,G,C]\/[A,T,G,C]/i) {
  $COORD_SYSTEM = "RESIDUE";
} # end elsif($first_line =~ /[\d+,x,X,y,Y],\d+,\-?1,[A,T,G,C]\/[A,T,G,C]/i)
else {
  $content .= "<p class=\"header1\">Error</p>\n";
  $content .= "<p>Incorrect input format. Please see <a href=\"\/chr_coords_example.html\">Sample format</a></p>";
  last;	
}

$content .= "<p>Your input data has been recognized to use <a href=\"\/chr_coords_example.php\" rel=\"external\">$COORD_SYSTEM based coordinate system</a>. Your job id is $pid and is currently running.  Your job has been partitioned into datasets of $num_coords_per_split positions and the status of each job can be viewed in the <a href=\"$outpage_url\" rel=\"external\">SIFT results status page</a>. Once the status of a job is <span class=\"green\">'Complete'</span>, you may view or download the results. A partitioned job typically takes 6-7 min to complete.</p><p>Proceed to <a href=\"$outpage_url\" rel=\"external\">SIFT results status page.</a><p>Problems? Contact <a href=\"contact.pl\">us</a> with your job id.</p>\n";

&finish_script($content, 0); 

# Close the I/O handles
#close(STDin);
close(STDOUT);
#close(STDERR); 

#prepare output status page
my $outpage_content = '';
my $outpage_header = 'SIFT Results Status';
my $outpage_counter = 0;

open (OUTPAGE,">>$tmp/$pid.outpage.html") || die ("cannot open outpage");
$outpage_content .= "<p><a href=\"$outpage_url\">Refresh page</a></p>";
$outpage_content .= '<table>';
$heading = "<tr class=\"tableHeader\"><td><p>Job</p></td><td><p>Job size</p></td><td><p>Job ID</p></td><td><p>Job status</p></td><td><p>View results</p></td><td><p>Download results</p></td></tr>";
$outpage_content .= "$heading\n";


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
                $outpage_content .= "<tr class=\"tableRow".($outpage_counter++ % 2 == 0 ? 'Even' : 'Odd')."\"><td><p>Partitioned set $input_set</p></td>";
                $outpage_content .= "<td><p>Input rows $start to $end</p></td>";
                $outpage_content .= "<td><p>$new_pid</p></td>";
                $outpage_content .= "<td><p>Not started.</p></td><td><p>Not available</p></td><td><p>Not available</p></td>";
                $outpage_content .= "</tr>\n";
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
	$outpage_content .= "<tr class=\"tableRow".($outpage_counter++ % 2 == 0 ? 'Even' : 'Odd')."\"><td><p>Complete set</p></td>";
        $outpage_content .= "<td><p>Input rows 1 to $end</p></td>";
        $outpage_content .= "<td><p>$pid</p></td>";
        $outpage_content .= "<td><p>Not started.</p></td><td><p>Not available</p></td><td><p>Not available</p></td>";
        $outpage_content .= "</tr></table>\n";

}
else{
	push @pid_list, $new_pid;
	$outpage_content .= "<tr class=\"tableRow".($outpage_counter++ % 2 == 0 ? 'Even' : 'Odd')."\"><td><p>Partitioned set $input_set</p></td>";
	$outpage_content .= "<td>Input rows $start to $end</p></td>";
	$outpage_content .= "<td><p>$new_pid</p></td>";
	$outpage_content .= "<td><p>Not started.</p></td><td><p>Not available</p></td><td><p>Not available</p></td>";
	$outpage_content .= "</tr>\n";
	$outpage_content .= "<tr class=\"tableRow".($outpage_counter++ % 2 == 0 ? 'Even' : 'Odd')."\"><td><p>Complete set</p></td>";
        $outpage_content .= "<td><p>Input rows 1 to $end</p></td>";
        $outpage_content .= "<td><p>$pid</p></td>";
        $outpage_content .= "<td><p>Not started.</p></td><td><p>Not available</p></td><td><p>Not available</p></td>";
        $outpage_content .= "</tr></table>\n";

}

$outpage_content .= "<p class=\"header1\">Batch Report</p>\n";
$outpage_content .= "<p>Number of input (non-intronic) variants: <br />\n";
$outpage_content .= "Coding variants: <br />\n";
$outpage_content .= "Coding variants predicted: <br />\n";
$outpage_content .= "Tolerated: <br />\n";
$outpage_content .= "Damaging: <br />\n";
$outpage_content .= "Nonsynonymous: <br />\n";
$outpage_content .= "Synonymous: <br />\n";
$outpage_content .= "Novel: <br />\n";

  my $outpage_title = "SIFT: Feed to CHR Coords";
  my @stylesheets = qw(/stylesheets/main.css);
 
  my $outpage_tt = Template->new({
    ABSOLUTE => 1,
  });

  my $template_file = '/usr/local/common/web'.lc($ENV{"WEBTIER"}).'/templates/3_column_fixed_width.tpl';
  my $vars = {
    main_content => $outpage_content,
    page_header => $outpage_header,
    title => $outpage_title,
    stylesheets => \@stylesheets,
  };

  my $template_output = '';
  
  $outpage_tt->process($template_file, $vars, \$template_output) || die $oupage_tt->error();
  print OUTPAGE $template_output;

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
		print BATCH_FILE "$pid\t$new_pid\t$tmp/$new_pid.chrfile\t$organism\t$seq_identity_filter\t$COORD_SYSTEM\tLAST\t$output_options\t$address\n";
	}
	else{
		print BATCH_FILE "$pid\t$new_pid\t$tmp/$new_pid.chrfile\t$organism\t$seq_identity_filter\t$COORD_SYSTEM\tNOT_LAST\t$output_options\t$address\n";
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
sub hextodec {
        unpack( "N", pack( "H8", substr( "0" x 8 . shift, -8 ) ) );
}

  ###########################################################################
  # This function calls the template object and injects the HTML variables
  # into the template. 
  #
  # @param : the HTML content
  # @param : the exit code
  # @return void
  ###########################################################################
  sub finish_script {
    my($content, $code) = @_;

    my $title = "SIFT: Genome Center";
    my $header = "SIFT Genome Center";
    my @stylesheets = qw(/stylesheets/main.css);

    my $tt = Template->new({
      ABSOLUTE => 1,
    });

    my $template_file = '/usr/local/common/web'.lc($ENV{"WEBTIER"}).'/templates/3_column_fixed_width.tpl';
    my $vars = {
      main_content => $content,
      title => $title,
      page_header => $header,
      stylesheets => \@stylesheets,
    };

    $tt->process($template_file, $vars) || die $tt->error();

    if($code < 0) {  
      exit($code);
    } # end if($code < 0)
  } # end finish_script
  ###########################################################################