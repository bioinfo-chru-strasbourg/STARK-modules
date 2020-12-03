#!/usr/bin/perl

#	impala.pl
# Execute impala & display the output
# 4/13/07 Execute getblock.pl instead of getblock.sh

#	Set file permissions to rw-rw----
system("umask 006");
$bin = "./";

# output the beginning text to be used on all pages
print "Content-type: text/html\n\n";
print "<TITLE>IMPALA Results</TITLE>\n";
if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
   print "This script should be referenced with a METHOD of POST\n";
   exit;
}

read (STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"});
%names = &general_parse($ENV{"CONTENT_TYPE"}, $QUERY_STRING);

#print $names{database};

if ($names{database} =~ /minus/) { $db = "blminus"; }
elsif ($names{database} =~ /plus/) { $db = "blplus"; }
elsif ($names{database} =~ /prints/) { $db = "prints"; }
else { $db = "blplus";  }

#print "$names{database} $db\n";
if ($db eq "") { $db = "blplus"; }
#$db = $names{database};
#$db =~ /\s+//g;

if ($names{ex} eq "") { $names{ex} = "5"; }
$ex = $names{ex};

#	Get the sequence & write it to a file
if ($names{sequence} eq "") {
   print "<H1>Error</H1> Please enter a query sequence.<P>\n";
   exit;
}
$seq = "../tmp/$$.seq";
open(SEQ, ">$seq");
print SEQ $names{sequence};
print SEQ "\n";
close(SEQ);

#	Filter flag, default is not to filter
$filt = "F";
if ($names{filter} eq "T") { $filt = "T"; }

#=========================================================================
#	Run shell now

$out = "../tmp/$$.out";
$arguments = $seq.' '.$db.' '.$ex.' '.$filt;
#print $arguments;
#	Only the first 2 arguments are getting passed this way
system("$bin/impala.csh $seq $db $ex $filt > $out 2>&1");

#	why doesn't umask take care of this?
system("chmod 660 ../tmp/$$.*");

print "<HTML><BASE HREF=\"http://blocks.fhcrc.org\">\n";
print "<TITLE>Block Searcher Impala Results</TITLE>\n";
print "<H1>Block Searcher Impala Results</H1>\n";

print"<PRE>\n";
open(OUT, "<$out");
while ($_ = <OUT>)
{
  # to replace the html special characters with escape sequences
# s/&/&amp/g ;
# s/>/&gt/g ;
# s/</&lt/g ;
# s/"/&quot/g ;

  if (/^IPB\d\d\d\d\d\d\S*/) {
    s|(^IPB\d\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1\2">\1\2</A>|g
  }
  elsif (/^BL\d\d\d\d\d\S*/) {
    s|(^BL\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1\2">\1\2</A>|g
  }
  elsif (/^DM\d\d\d\d\d\S*/) {
    s|(^DM\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1\2">\1\2</A>|g
  }
  elsif (/^PD\d\d\d\d\d\S*/) {
    s|(^PD\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1\2">\1\2</A>|g
  }
  elsif (/^BP\d\d\d\d\d\S*/) {
    s|(^BP\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1\2">\1\2</A>|g
  }
  elsif (/^PR\d\d\d\d\d\S*/) {
    s|(^PR\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1\2">\1\2</A>|g
  }
  elsif (/^PF\d\d\d\d\d\S*/) {
    s|(^PF\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1\2">\1\2</A>|g
  }
  print;
}

close(OUT);

print "</PRE>";
print "</PRE><A HREF=\"#top\">[return to top]</A><P>";
#-------------------------------------------------------------------------
exit (0);


#-------------------------------------------------------------------------
# $names = &general_parse($ENV{CONTENT_TYPE}, $QUERY_STRING);
# parameters:	CONTENT_TYPE
#		QUERY_STRING
# returns: an associative array of name/value pairs.  The name is the key.

# WARNING:  Some of this routine is program-dependent!!!

# CONTENT_TYPE: application/x-www-form-urlencoded
# QUERY_STRING: key1=val1&key2=val2

# CONTENT_TYPE: multipart/form-data; boundary=<boundary>
# QUERY_STRING: <boundary>
#		Content-Disposition: form-data; name="key1"
#		<blank line>
#		val1
#		<boundary>
#		Content-Disposition: form-data; name="key2"
#		<blank line>
#		val2
#		<boundary>

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
     }

  }
  else
  {
     $loc = index($content_type, "boundary=");
     if ($loc > 0)
     { 
        $temp = substr($content_type, $loc+9);
#		Why is this necessary? (boundary= doesn't match actual)
	$boundary = "--".$temp;
#print "<PRE>";
#print "$query_string\n";
#print "boundary=$boundary\n";
        # break up into individual name/value lines
        @q = split(/$boundary/, $query_string);
        foreach $pair (@q) {
          # break the name/value pairs up
#print "pair=$pair\n";
          $loc = index($pair, "name=");
	  $temp = substr($pair, $loc+5);
#	  $loc = index($temp, "\n\n");
 	  $loc = index($temp, "\n");
	  $temp1 = substr($temp, $loc+2);
#print "1temp=$temp\ntemp1=$temp1\n";
#		Get rid of stuff after the name; including semicolon if any
	  $loc_semi = index($temp, ";");
	  $loc_eol = index($temp, "\n");
	  $loc = $loc_eol;
          if ($loc_semi > 0 && $loc_semi < $loc) {$loc = $loc_semi; }
	  if ($loc > 0) { $temp = substr($temp, 0, $loc); }
#		Get rid of quotes around the name
          $temp =~ s/\"//g;
#		Still has a trailing whitespace character ...
          $temp =~ s/\s//g;
#		Need to strip leading/ending whitespace off of $temp1,
#		but be careful not to strip off internal CRs
#>>>>>       Following line is program-dependent !!!
	  if ($temp ne "dbfile" && $temp ne "sequence" && $temp ne "database")
 	  { $temp1 =~ s/\s//g; }
#		MAC file lines end in just \r, no \n, so makelis won't find all
#		of the sequences; DOS file lines end in \r\n, UNIX in \n.
	  if ($temp eq "seqfile")
#	  { $temp1 =~ s/\r/\n/g; $temp1 =~ s/\n\n/\n/g; }
#		Change \r\n to \n and then \r to \n
 	  { $temp1 =~ s/\r\n/\n/g; $temp1 =~ s/\r/\n/g; }
#print "2temp=$temp\ntemp1=$temp1\n";
	  if ($temp ne "") { $ans{$temp} = $temp1; }
        }
     }
     else
     {  print "Cannot parse\n";  
        print "content_type=$content_type\n";
        print "query_string=$query_string\n";
     }
  }
  return %ans;
#print "</PRE>";
}   # end of general_parse
#-------------------------------------------------------------------------
# parameter: a hex representation of a number (doesn't need to be a string)
# returns: the decimal representation of the number
sub hextodec {
  unpack("N", pack("H8", substr("0" x 8 . shift, -8)));
}

