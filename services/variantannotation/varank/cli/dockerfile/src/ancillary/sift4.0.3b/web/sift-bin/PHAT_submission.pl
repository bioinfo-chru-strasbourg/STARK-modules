#!/usr/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

# PHAT_submission.pl
# Execute phat_submission and display the output
# (1) take sequence and turn into profile using pro_to_prf.pl
# (2) execute swat against given database sequences
# for PHAT paper (Ng, Henikoff, Henikoff 2000)

#	Set file permissions to rw-rw----
system("umask 006");
$bin = ".";
$data = "../data";
$tmp = "../tmp";

# output the beginning text to be used on all pages
print "Content-type: text/html\n\n";
print "<TITLE>SWAT Results</TITLE>\n";
print "<body bgcolor=white>\n";
print "<basefont size=4>\n";


if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
   print "This script should be referenced with a METHOD of POST\n";
   exit;
}

read (STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"});
%names = &general_parse($ENV{"CONTENT_TYPE"}, $QUERY_STRING);

# check that there's a query sequence
if ($names{queryfile} eq "" && $names{queryseq} eq "") {
	print "<H1>Error</H1>Please enter a sequence.<P>\n";
	exit;
}

#create a file of the query sequence
#
 $query = $tmp . "/$$.fasta";

open (QUERY, ">$query");
if ($names{queryfile} ne "")
{ # SWAT wants profile in capital letters, so put query sequence in caps
  $names{queryfile} =~ tr/a-z/A-Z/;
  print QUERY $names{queryfile}; 
}
if ($names{queryseq} ne "")
{ 
  $names{queryseq} =~ tr/a-z/A-Z/;
  print QUERY $names{queryseq}; 
}
print QUERY "\n";
close (QUERY);

# check that transmembrane regions are specificed
if ($names{querytxmemfile} eq "" && $names{querytxmemseq} eq "") {
        print "<H1>Error</H1>Please enter a sequence that designates transmembrane regions.<P>\n";
        exit;
}

#create a file of the transmembrane region sequence

$txmem = $tmp . "/$$.txmem";
open (TXMEM, ">$txmem");
if ($names{querytxmemfile} ne "")
{ print TXMEM $names{querytxmemfile}; }
if ($names{querytxmemseq} ne "")
{ print TXMEM $names{querytxmemseq}; }
print TXMEM "\n";
close (TXMEM);       

# check that there's a database sequence
if ($names{dbfile} eq "" && $names{database_seq} eq "") {
        print "<H1>Error</H1>Please enter database sequences.<P>\n";
        exit;
}

#create a file of the query sequence


$dbase = $tmp . "/$$.dbase";

open (DB, ">$dbase");
if ($names{dbfile} ne "")
{ print DB $names{dbfile}; }
if ($names{database_seq} ne "")
{ print DB $names{database_seq}; }
print DB "\n";
close (DB);       

#=========================================================================
#	Run c-shell now
#	Following ought to work, but perl is unreliable
#open(OUT, "$bin/sample.csh $names{program} $seq $blk 2>&1 |");
#wait;

$nontxmem_matrix = $data . "/" . $names{nontxmem_matrix};
$txmem_matrix = $data . "/" . $names{txmem_matrix};

$out = $query . ".out";
 
$program_call = $bin . "/PHAT_submission.csh";
$queryclean = $query . ".clean";
$dbaseclean = $dbase . ".clean";

system("$bin/fastaseqs $query $queryclean > /dev/null");
system ("$bin/fastaseqs $dbase $dbaseclean > /dev/null");


#system("$bin/PHAT_submission.csh $queryclean $txmem $dbaseclean $nontxmem_matrix $txmem_matrix > $out 2>&1");

system("$bin/PHAT_submission.csh $queryclean $txmem $dbaseclean $nontxmem_matrix $txmem_matrix > $out 2>/dev/null");

print "</PRE>";

print "<A NAME=top><H1><center>SWAT Results</center></H1></A>";
print "<PRE>";
print "<b>Search and alignment results start at *********************** </b><BR><BR>\n";

print "Matrix used on nontransmembrane regions:\n";
print "<b>";
if ($names{nontxmem_matrix} eq "BLOSUM62") {
	print "BLOSUM 62";
} elsif ($names{nontxmem_matrix} eq "jones_250.bla") {
	print "Jones, Taylor, Thornton PAM 250";
} elsif ($names{nontxmem_matrix} eq "BLOSUM55") {
	print "BLOSUM 55";
}

print "</b><BR>\n";

print "Matrix used on transmembrane regions:\n";
print "<b>";
if ($names{txmem_matrix} eq "BLOSUM62") {
        print "BLOSUM 62";
} elsif ($names{txmem_matrix} eq "jones_170.bla") {
	print "Jones, Taylor, Thornton PAM 170";
} elsif ($names{txmem_matrix} eq "BLOSUM55") {
	print "BLOSUM 55";
} elsif ($names{txmem_matrix} eq "xmem23_170.bla") {
	print "Jones, Taylor Thornton transmembrane 170";
} elsif ($names{txmem_matrix} eq "PHAT_T75_B73.bla") {
	print "PHAT 75/73";
}
print "</b><BR><BR>\n";


open(OUT, "<$out");
while ($_ = <OUT>) { print; }
close(OUT);

system ("rm $out");

print "</PRE>";
print "</PRE><A HREF=\"#top\">[return to top]</A><P>";
#-------------------------------------------------------------------------
exit (0);

#-------------------------------------------------------------------------
#
# parameter: a string that is the html QUERY_STRING environment variable
# returns: an associative array of name/value pairs.  The name is the key.
sub parse_query {
  local($query_string) = @_;
  local(%ans, @q, $pair);
#print $query_string;
  # break up into individual name/value lines
  @q = split(/&/, $query_string);

  foreach $pair (@q) {
    # break the name/value pairs up
    # use split rather than regular expressions because the value may have
    #  newlines in it
    split(/=/, $pair, 2);                        

    # change '+' to ' '
    $_[1] =~ s/\+/ /g;

    # change the escaped characters (has to be after the split on '&' and '=')
    $_[1] =~ s/%(..)/pack("c",&hextodec("$1"))/eg;

    $ans{$_[0]} = $_[1];
  }

  return %ans;
}                                     
#-------------------------------------------------------------------------
# parameter: a hex representation of a number (doesn't need to be a string)
# returns: the decimal representation of the number
sub hextodec {
  unpack("N", pack("H8", substr("0" x 8 . shift, -8)));
}

#-------------------------------------------------------------------------
# $names = &general_parse($ENV{CONTENT_TYPE}, $QUERY_STRING);
# parameters:   CONTENT_TYPE
#               QUERY_STRING
# returns: an associative array of name/value pairs.  The name is the key.

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
  local($content_type, $query_string) = @_;
  local(%ans, @q, $pair, $loc, $boundary, $temp, $temp1);

  if ($content_type eq "application/x-www-form-urlencoded")
  {
     # break up into individual name/value lines
     @q = split(/&/, $query_string);

     foreach $pair (@q) {
       # break the name/value pairs up
       # use split rather than regular expressions because the value may have
       #  newlines in it
       split(/=/, $pair, 2);

       # change '+' to ' '
       $_[1] =~ s/\+/ /g;

       # change the escaped characters (must be after the split on '&' and '=')
       $_[1] =~ s/%(..)/pack("c",&hextodec("$1"))/eg;

       $ans{$_[0]} = $_[1];
     } #end of foreach $pair 

  } #end of if ($content_type) 
  else
  {
     $loc = index($content_type, "boundary=");
     if ($loc > 0)
     {
        $temp = substr($content_type, $loc+9);
#               Why is this necessary? (boundary= doesn't match actual)
        $boundary = "--".$temp;
        # break up into individual name/value lines
        @q = split(/$boundary/, $query_string);


        foreach $pair (@q) {
	 # break the name/value pairs up
          $loc = index($pair, "name=");
          $temp = substr($pair, $loc+5);
#         $loc = index($temp, "\n\n");
          $loc = index($temp, "\n");
          $temp1 = substr($temp, $loc+2);
#   Get rid of stuff after the name; including semicolon if any
          $loc_semi = index($temp, ";");
          $loc_eol = index($temp, "\n");
          $loc = $loc_eol;
          if ($loc_semi > 0 && $loc_semi < $loc) {$loc = $loc_semi; }
          if ($loc > 0) { $temp = substr($temp, 0, $loc); }
#               Get rid of quotes around the name
          $temp =~ s/\"//g;
#               Still has a trailing whitespace character ...
          $temp =~ s/\s//g;
#               Substitute spaces with nothing
#		Need to strip leading/ending whitespace off of $temp1,
#               but be careful not to strip off internal CRs
#               MAC file lines end in just \r, no \n, so makelis won't find all
#               of the sequences; DOS file lines end in \r\n, UNIX in \n.
#               Change \r\n to \n and then \r to \n
#######PROGRAM -SPECIFIC!!!!!!!######################
#In my case, I want to keep the newlines in fields which have "file" or 'seq"
# and remove newlines everywhere else.
	   if  ($temp =~ /file/ || $temp =~ /seq/) 
 	   { $temp1 =~ s/\r\n/\n/g; $temp1 =~ s/\r/\n/g; }
	   # for other variables that are NOT files or file-like, remove extra
	   #whitespace
	    else { $temp1 =~ s/\s//g;} 
	if ($temp ne "") { $ans{$temp} = $temp1; }
     }# end of foreach
     } #end of if loc > 0
     else
     {  print "Cannot parse\n";
        print "content_type=$content_type\n";
        print "query_string=$query_string\n";
     }
  }
  return %ans;
#print "</PRE>";
}   # end of general_parse

