#!/usr/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.
require 'SIFT_subroutines.pm';

# seq_to_subfamily_block.pl
# Execute seq_to_subfamily.csh and display the output

#	Set file permissions to rw-rw----
system("umask 002");
$bin = "/opt/www/sift/htdocs/sift-bin";
$tmp = "/opt/www/sift/tmp";

$pid = $$;


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
if ($names{alignment_file} eq "" && $names{aligned_seqs} eq "") {
	print "<H1>Error</H1>Please enter aligned sequences.<P>\n";
	exit;
}

#create a file for the aligned sequences
#
system ("umask 000");
$alignedseq = $tmp . "/$pid.alignedfasta";
$out = $tmp . "/$pid.siftresults"; # where sifting results are printed
$comments = $tmp . "/$pid.comments";
open (ALIGNEDSEQ, ">$alignedseq");
if ($names{alignment_file} ne "")
{ 
  $names{alignment_file} =~ tr/a-z/A-Z/;
  print ALIGNEDSEQ $names{alignment_file}; 
}
if ($names{aligned_seqs} ne "")
{ 
  $names{aligned_seqs} =~ tr/a-z/A-Z/;
  print ALIGNEDSEQ $names{aligned_seqs};
}
print ALIGNEDSEQ "\n";
close (ALIGNEDSEQ);

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

## SIFT prediction options

	$exp_option = 1;

if ($names{gap_option} ne "TRUE") {
	$gap_option = 0;
} else {$gap_option = 1; }


$seq_identity_filter = $names{seq_identity_filter};
$seq_identity_filter = trim ($seq_identity_filter);

######  Calling the program and printing results #########

print "<A NAME=top><H1><center>S<font color=red>I</font>F<font color=blue>T</font> results</center></H1></A><BR>\n";

#  ERROR HANDLING
#  $$.err will contain errors that do not stop program from running
#  like max. # of sequences exceeded or mistake in substitution
# .error files contain error messages that explain why program was stopped
#

$program_call = $bin . "/SIFTING2.csh";
system("$program_call $pid alignedseqs $alignedseq $out 27 $subst_file $gap_option $seq_identity_filter > $comments 2>/dev/null");

#print "the seeq idenityt filter is $seq_identity_filter kdfawe<BR>"; 
#print"OUT file <BR>";
#print "<PRE>";
#open (OUT, "$alignedseq") ;
#while (<OUT>) {print $_;}
#close (OUT);
#print "SUBST file\n";
#$file = $tmp . "/" . $pid . ".commentscsh";
#print "trying to get comments from $file";
#open (OUT, "$file") || die "can't open $tmp/$pid.commentscsh";
#while (<OUT>) {print $_;}
#close (OUT);
#print "</PRE>";
@error_files = `ls $tmp/$pid.*.error`;
if (@error_files > 0) {
        # there's an error
        print "<b>ERROR!</b>  Please email <A HREF=\"/sift-bin/contact.pl\">us</A> with";
	print " your sequence if you do not understand the error.<BR>";
print "ERROR files @error_files<BR>";
print "$pid\n";
        foreach $errorfile (@error_files) {
print "$errorfile name<BR>";
                open (OUT, "<$errorfile");
                while ($_ = <OUT>) {print;}
                close (OUT);
        }
#       system ("rm -f $tmp/$$.*");
        exit (-1);
}

$no_table_files = `ls $tmp/$pid.aatable* |grep -v html| wc -l`;

print "<A HREF=\"/sift-bin/catfile_notitle.csh?$tmp/$pid.alignedfasta.gapsremoved+PRE\" TARGET=_blank>Alignment</A> (The first sequence is taken as the reference sequence.  Any position with a gap in the first sequence is removed.  Please check that the alignment is correct.)<BR>\n";

@table_files = `ls $tmp/$pid.aatable*`;
$subject = "\"SIFT results\"";  # literal double quotes for mailx
if ($no_table_files > 0) {
	tolerated_nottolerated_aminoacidtable($tmp,$no_table_files, $pid);
}

$commentscsh = $tmp . "/$pid.commentscsh"; 
$loc = rindex ($out, "/");
$file = substr ($out, $loc+1, 100);

print "<BR>\n";

print_score_desc ($tmp,$file);

print "<BR>\n";
print_predictions ($tmp,$file);
print "<BR>\n";

print "<BR>\n";

open (OUT, "<$tmp/$pid.err");
while ($_ = <OUT>) {
# don't want to print Jorja's messages, only my own.  so my error
# messages all have ERROR in the line
      if (/ERROR/ || /WARNING/ || /\*\*\*/) {
            print;
      }
}
close(OUT);



#system ("rm -f $tmp/$$.comments*");
#system ("rm -f $tmp/$$.tmp");

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
# and substitutions file and remove newlines everywhere else.
	   if  ($temp =~ /file/ || $temp =~ /seq/ || $temp =~ /subst/) 
 	   { 
		$temp1 =~ s/\r\n/\n/g; $temp1 =~ s/\r/\n/g; 
           }
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

