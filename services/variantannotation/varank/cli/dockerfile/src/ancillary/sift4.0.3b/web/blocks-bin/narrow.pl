#!/usr/bin/perl

#	pssm.pl
# Execute pssm.csh & display the output

#	Set file permissions to rw-rw----
system("umask 006");
#	Don't like to hard-code file names this way ...
$tmpdir = "../tmp/$$";
system("mkdir $tmpdir");

$bin = ".";

# output the beginning text to be used on all pages
print "Content-type: text/html\n\n";
print "<TITLE>Narrow PSSM Results</TITLE>\n";


if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
   print "This script should be referenced with a METHOD of POST\n";
   exit;
}

read (STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"});
#%names = &parse_query($QUERY_STRING);
%names = &general_parse($ENV{"CONTENT_TYPE"}, $QUERY_STRING);

$err = "$tmpdir/err";
$fin = "$tmpdir/fin";
#	$mast & $blimps are created by jcooper
$mast = "$fin.mast";
$blimps = "$fin.blimps";
$bout = "$fin.bout";
$mout = "$fin.mout";

#	Get the program name
#if ($names{program} eq "") {
#   print "<H1>Error</H1> Please select a searching program.<P>\n";
#   exit;
#}

#	Problem here; $names{qfile} always has some \n's, never empty ...
$qftemp = $names{qfile};
$qftemp =~ s/\s//g;
$qutemp = $names{query};
$qutemp =~ s/\s//g;
#	Get the query & write it to a file
if ($qutemp eq "" && $qftemp eq "")
{
         print "<H1>Error</H1> Please enter a query.<P>\n";
         exit;
}

#  create a file of the input
open(FIN, ">$fin");
#	Problem here; $names{qfile} always has some \n's, never empty ...
if ($qftemp ne "")
{ print FIN $names{qfile}; }
elsif ($qutemp ne "")
{ print FIN $names{query}; }
print FIN "\n";
close(FIN);

#=========================================================================
#	Run shell now

#print "<B>Your search is running, please wait ...</B><BR>";
#	Try to flush stdout so they see the above message
select(STDOUT); $| = 1;

system("$bin/narrow.csh $names{type} $fin > $err 2>&1");

print "<A NAME=top><H1>Narrow PSSM Results</H1></A>";

print "<PRE>";
open(ERR, "<$err");
while ($_ = <ERR>)
{
  print;
}
close(ERR);
print "</PRE>";


print "[<A HREF=\"/blocks-bin/catfile.sh?$fin\">Input</A>] ";
print "<P>";

print "[<A HREF=\"/blocks-bin/catfile.sh?$blimps\">BLIMPS PSSM</A>] ";
print "[<A HREF=\"/blocks-bin/matrix.csh?$blimps\">BLIMPS Search</A>]  ";
print "<BR>BLIMPS search vs Swiss-Prot will take several minutes to complete";
print "<P>";

print "[<A HREF=\"/blocks-bin/catfile.sh?$mast\">MAST PSSM</A>] ";
$group = "NarrowBlock";
print "[<A HREF=\"/blocks-bin/mastpssm.sh?$mast+$group\">MAST Search at SDSC</A>]  ";
print "[<A HREF=\"http://meme.sdsc.edu/meme/website/mast-intro.html\">About MAST</A>] ";
print "<BR>MAST search results vs several popular databases are returned via email from SDSC";
print "<BR>Use a large E-value cutoff for narrow PSSMs";
#print "<P>";


#print "<A HREF=\"#top\">[return to top]</A><P>";
#-------------------------------------------------------------------------
exit (0);

#=========================================================================

#-------------------------------------------------------------------------
#  NOTE: Not really general, does special stuff for names{query}
#	 and for names{sequences}
# $names = &general_parse($ENV{CONTENT_TYPE}, $QUERY_STRING);
# parameters:	CONTENT_TYPE
#		QUERY_STRING
# returns: an associative array of name/value pairs.  The name is the key.

# CONTENT_TYPE: application/x-www-form-urlencoded
# QUERY_STRING: key1=val1&key2=val2

# CONTENT_TYPE: multipart/form-data; boundary=<boundary>
# QUERY_STRING: <boundary>
#		Content-Disposition: form-data; name="qfile"; filename="file.name"
#		<blank line>
#		<content of file.pros>
#		<boundary>
#		Content-Disposition: form-data; name="sequences"
#		<blank line>
#		<data from input form>
#		<boundary>

sub general_parse {
  local($content_type, $query_string) = @_;
  local(%ans, @q, $pair, $loc, $boundary, $temp, $temp1);

#print "<PRE>$content_type\n\n";
#print "query_string:\n$query_string\n\n";

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
     }

  }
  else
  {
     $loc = index($content_type, "boundary=");
     if ($loc > 0)
     { 
        $temp = substr($content_type, $loc+9);
#	Why is this necessary? (boundary= sometimes doesn't match actual?)
 	$boundary = "--".$temp;
#       $boundary = $temp;
#print "boundary=$boundary\n\n";
        # break up into individual name/value lines
        @q = split(/$boundary/, $query_string);
        foreach $pair (@q) {
          # break the name/value pairs up
#print "pair=$pair\n\n";
          $loc = index($pair, "name=");
	  $temp = substr($pair, $loc+5);
#	  $loc = index($temp, "\n\n");
 	  $loc = index($temp, "\n");
	  $temp1 = substr($temp, $loc+2);
#print "1 temp=$temp\n";
#		Get rid of stuff after the name
#	Need to look for a ; before the first \n ... <<< find another way!
	  $loc1 = index($temp, ";");
          $loc = index($temp, "\n");
          if ($loc1 > 0 && $loc1 < $loc) { $loc = $loc1; }
 	  if ($loc > 0) { $temp = substr($temp, 0, $loc); }
#print "2 loc=$loc temp=$temp\n";
#		Get rid of quotes around the name
          $temp =~ s/\"//g;
#print "3 temp=$temp\n";
#		Still has a trailing whitespace character ...
          $temp =~ s/\s//g;

#		Need to strip leading/ending whitespace off of $temp1,
#		but be careful not to strip off internal CRs in "query"
#		and "sequences"
#print "4 temp=$temp\ntemp1=$temp1\n";
	  if ($temp ne "qfile" && $temp ne "query")
 	  { $temp1 =~ s/\s//g; }

#		MAC file lines end in just \r, no \n;
#		DOS file lines end in \r\n; UNIX in \n.
#		mablock uses fgets() which requires \n.
#		Change \r\n to \n, then change \r to \n
 	  if ($temp eq "qfile" || $temp eq "query")
  	  { $temp1 =~ s/\r\n/\n/g;  $temp1 =~ s/\r/\n/g; }
#print "temp=$temp\ntemp1=$temp1\n";
	  if ($temp ne "") { $ans{$temp} = $temp1; }
        }
     }
     else
     {  print "Cannot parse $content_type\n"; }
  }
#print "</PRE>";
  return %ans;
}

# parameter: a hex representation of a number (doesn't need to be a string)
# returns: the decimal representation of the number
sub hextodec {
  unpack("N", pack("H8", substr("0" x 8 . shift, -8)));
}

