#!/usr/bin/perl
#  SIFT_seq_submit2.pl
# SIFT_seq_submit2.pl
# version 2
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

$| = 1;
require 'SIFT_subroutines.pm';
#	Set file permissions to rw-rw----
system("umask 006");
$bin = ".";
$tmp = "../tmp";
$pid = $$;
$program_call = $bin . "/SIFTING2.csh";
#$program_call = $bin . "/TESTING2.csh";
$return_address = "sift\@fhcrc.org";

# output the beginning text to be used on all pages
print "Content-type: text/html\n\n";
#print "<TITLE>SIFTING Results</TITLE>\n";
print "<body bgcolor=white>\n";

if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
   print "This script should be referenced with a METHOD of POST\n";
   exit;
}

read (STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"});
%names = &general_parse($ENV{"CONTENT_TYPE"}, $QUERY_STRING);

if ($names{address} ne "") {
	$address = $names{address};
}

# check that there's a query sequence
if ($names{query_file} eq "" && $names{queryseq} eq "") {
	print "<H1>Error</H1>Please enter your query sequence.<P>\n";
	exit;
}

$queryseqfile = $tmp . "/$pid.seq";
$out = $tmp . "/$pid.siftresults"; # where sifting results are printed
$comments = $tmp . "/$pid.comments";
open (QUERYSEQFP, ">$queryseqfile");
if ($names{query_file} ne "")
{ 
  $names{query_file} =~ tr/a-z/A-Z/;
  print QUERYSEQFP $names{query_file}; 
}
if ($names{queryseq} ne "")
{ 
 $names{queryseq} =~ tr/a-z/A-Z/;
  print QUERYSEQFP $names{queryseq};
}
print QUERYSEQFP "\n";
close (QUERYSEQFP);

### substitution file #######
$subst_file = $tmp . "/$pid.substfile";
if ($names{substitution_file} eq "" && $names{substitutions} eq "") {
        $subst_file = "-";
}
if ($subst_file ne "-") { # have some substitutions
        open (SUBSTITUTION, ">$subst_file");
        if ($names{substitution_file} ne "") {
                $names{substitution_file} =~ tr/a-z/A-Z/;
                print SUBSTITUTION $names{substitution_file};
        }
        if ($names{substitutions} ne "") {
                $names{substitutions} =~ tr/a-z/A-Z/;
                print SUBSTITUTION $names{substitutions};
        }
        close (SUBSTITUTION);
}

## SIFT prediction operations

#if ($names{exp_option} eq "TRUE") {
        $exp_option = 1;
#} else { $exp_option = 0; }
$info = $names{info};


$seq_identity_filter = $names{seq_identity_filter};
$seq_identity_filter = trim ($seq_identity_filter);

#if ($names{gap_option} ne "TRUE") {
#        $gap_option = 0;
#} else {$gap_option = 1; }

######  Calling the program #########

print "<A NAME=top><H1><center>S<font color=red>I</font>F<font color=blue>T</font> results</center></H1></A><BR>\n";

print "Your job id is $pid<BR>\n";
print "<A HREF=\"http://blocks.fhcrc.org/sift-bin/format.pl?$pid\">Check here</A> for results, they will be available for approximately 24 hours.<P>\n";
#       Make sure STDOUT gets flushed before the exec()
select(STDOUT); $| = 1;

if ($address ne "") {
	print "Results will be sent to $address<BR><BR>\n";
	print "<b>Use web browser to view these files.</b><BR>";
	print "<i>Please note that tables will take some time to load.</i>";
	# Put it in the queue
        #NOTE: Results are actually mailed by $program_call
        system("./add_queue_entry.pl SIFT_queue $program_call $pid query_seq $queryseqfile $out $info $subst_file 0 $seq_identity_filter $names{database} $address > $comments 2>/dev/null");
}
else
{
	system("$program_call $pid query_seq $queryseqfile $out $info $subst_file 0  $seq_identity_filter $names{database} $address > $comments 2>&1 &");
}

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
	   if  ($temp =~ /file/ || $temp =~ /seq/ || $temp =~ /subst/) 
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

