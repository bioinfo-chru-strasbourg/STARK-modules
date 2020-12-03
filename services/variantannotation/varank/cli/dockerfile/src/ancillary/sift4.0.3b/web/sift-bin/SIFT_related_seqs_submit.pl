#!/usr/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
require 'SIFT_subroutines.pm';

# SIFT_related_seqs_submit.pl
# Execute SIFTING.csh with first field related_seqs

#	Set file permissions to rw-rw----
system("umask 006");
$bin = "/opt/www/sift/htdocs/sift-bin";
$tmp = "/opt/www/sift/tmp";

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

# check that there's a query sequence
if ($names{seqs_file} eq "" && $names{seqs} eq "") {
	print "<H1>Error</H1>Please enter aligned sequences.<P>\n";
	exit;
}

#create a file for the aligned sequences
#
$unalignedseq = $tmp . "/$$.unalignedfasta";
$out = $tmp . "/$$.siftresults"; # where sifting results are printed
$comments = $tmp . "/$$.comments";
open (UNALIGNEDSEQ, ">$unalignedseq");
if ($names{seqs_file} ne "")
{ 
  $names{seqs_file} =~ tr/a-z/A-Z/;
  print UNALIGNEDSEQ $names{seqs_file}; 
}
if ($names{seqs} ne "")
{ 
  $names{seqs} =~ tr/a-z/A-Z/;
  print UNALIGNEDSEQ $names{seqs};
}
print UNALIGNEDSEQ "\n";
close (UNALIGNEDSEQ);


$subst_file = $tmp . "/$$.substfile";
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

## SIFT prediction options

        $exp_option = 1;

if ($names{gap_option} ne "TRUE") {
        $gap_option = 0;
} else {$gap_option = 1; }


$seq_identity_filter = $names{seq_identity_filter};
$seq_identity_filter = trim ($seq_identity_filter);

######  Calling the program and printing results #########

print "<A NAME=top><H1><center>S<font color=red>I</font>F<font color=blue>T</font> results</center></H1></A><BR>\n";

$program_call = $bin . "/SIFTING.csh";
$pid = $$;
print "<PRE>";
system("$program_call $pid related_seqs $unalignedseq $out 27 $subst_file $gap_option $seq_identity_filter > $comments 2>/dev/null");

@error_files = `ls $tmp/$$.*.error`;
if (@error_files > 0) {
        # there's an error
        #print "<b>ERROR!</b>  Please email <A HREF=\"/blocks-bin/contact.pl\">us</A> with your sequence if you do not understand the error.<BR>"; 
        foreach $errorfile (@error_files) {
                open (OUT, "<$errorfile");
                while ($_ = <OUT>) {print;}
                close (OUT);
        }

#        system ("rm -f $tmp/$$.*");
#        exit (-1);
}

print "</PRE>";

$commentscsh = $tmp . "/$pid.commentscsh"; 
$psiblastout = $tmp . "/$pid" . ".alignedfasta";

#comments
#open (OUT, "<$comments");
#while ($_ = <OUT>) {print;}
#close (OUT);

print "<A HREF=\"/sift-bin/catfile.csh?$tmp/$pid.msf+Alignment+PRE\">PSIBLAST alignment of submitted sequences</A><BR>\n";

print "<A HREF=\"/sift-bin/catfile_notitle.csh?$tmp/$pid.alignedfasta+PRE\" TARGET=_blank>Alignment in FASTA format </A> (for modification)<BR>";

print_psiblast_alignment_desc();
print "<BR>";

$no_table_files = `ls $tmp/$$.aatable*|grep -v html | wc -l`;

tolerated_nottolerated_aminoacidtable($tmp,$no_table_files, $$);

$commentscsh = $tmp . "/$$.commentscsh";
$loc = rindex ($out, "/");
$file = substr ($out, $loc+1, 100);
print "<BR>\n";

print_score_desc ($tmp,$file);
print "<BR>\n";


print_predictions ($tmp,$file);
print "<BR>\n";


#errors that are wrong but don't disrupt program
open (OUT, "<$tmp/$pid.err");
while ($_ = <OUT>) {
# don't want to print Jorja's messages, only my own.  so my error
# messages all have ERROR in the line
      if (/ERROR/ || /WARNING/ || /\*\*\*/) {
            print;
      }
}
close(OUT);

system ("rm -f $tmp/$$.comments*");
system ("rm -f $tmp/$$.tmp");

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

