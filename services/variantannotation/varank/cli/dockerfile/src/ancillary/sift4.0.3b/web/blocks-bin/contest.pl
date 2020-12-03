#!/usr/bin/perl

# contact.pl

#	Set file permissions to rw-rw----
system("umask 006");

$mailprog = "/usr/bin/mailx";
#$webmaster = "jorja\@fhcrc.org";
#$siftmaster = "sift\@fhcrc.org";
$webmaster = "jorja";
$siftmaster = "sift";
$lamamaster = "lama";
$tmp = "../tmp";
$maxsecs = 300;		# 5 minutes
@colors = ("red","blue","green","maroon","purple","lime");
@spam_list = ("ku.name");

select(STDOUT); $| = 1;
print "Content-type: text/html\n\n";
print "<HTML><TITLE>Blocks Contact</TITLE>\n";
print "<H1 ALIGN=CENTER><IMG ALIGN=MIDDLE SRC=\"/blocks/icons/hutch_logo.gif\">
Blocks WWW Server Contact Form
<IMG ALIGN=MIDDLE SRC=\"/blocks/icons/blocks.xbm\"></H1><HR>\n";

#----------------------------------------------------------------------
#Process the form
$x = $ENV{"REQUEST_METHOD"};
if ( $ENV{"REQUEST_METHOD"} ne "POST" )      #GET
{
   $QUERY_STRING = $ENV{"QUERY_STRING"};
}
else
{
   read (STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"});
}
%names = &parse_query($QUERY_STRING);
$form = $names{form};
#print "<PRE>$x form=$names{form} randnum=$names{randnum}\n</PRE>";
if ($form eq "") { $form = "start"; }

#Try to detect probably spammers
if ( $form ne "start" && $ENV{"REQUEST_METHOD"} ne "POST" )
{   exit(-1);  }

#----------------------------------------------------------------
#Send an input form with a random number embedded in it
if ($form eq "start")
{
   $x = int(rand(5));
   $color = $colors[$x];
   if ($x < 0 || $x > 5) { $x = 0; }
   $randnum = int(rand(100000)) + 1;
   $gmtsecs = time();
   #  Send input form
   print "<FORM METHOD=\"POST\" ACTION=\"/blocks-bin/contest.pl\">";
   print "<INPUT TYPE=hidden name=form VALUE=\"form1\">";
   print "<INPUT TYPE=hidden name=randnum VALUE=\"$randnum\">";
   print "<INPUT TYPE=hidden name=gmtsecs VALUE=\"$gmtsecs\">\n";
   print "<P>We will try to respond to your questions within one working day.<BR>\n";
   print "<font size=5 color=$color>$randnum</font>\n";
   print "<UL><P><LI>To discourage spammers, please type in the <font color=$color>$color</font> number printed above: <INPUT SIZE=5 NAME=\"randchk\"><P>\n";
   print "<LI>Please enter your email address: <INPUT SIZE=50 NAME=\"email\"><P>\n";
   print "<LI>Please enter your name: <INPUT SIZE=50 NAME=\"name\">\n";
   print "<LI><P>Select a topic: <SELECT NAME=\"topic\" SIZE=1>\n";
   print "<OPTION VALUE=\"db\" SELECTED>Blocks Database\n";
   print "<OPTION VALUE=\"bs\">Block Searcher\n";
   print "<OPTION VALUE=\"bm\">Block Maker\n";
   print "<OPTION VALUE=\"ch\">CODEHOP\n";
   print "<OPTION VALUE=\"si\">SIFT\n";
   print "<OPTION VALUE=\"la\">LAMA\n";
   print "<OPTION VALUE=\"bl\">BLIMPS Software\n";
   print "<OPTION VALUE=\"ot\">Other\n";
   print "</SELECT>\n";
   print "<LI><P>Enter your question or comment:\n";
   print "<TEXTAREA NAME=\"comment\" ROWS=10 COLS=60></TEXTAREA>\n";
   print "<LI><INPUT TYPE=submit VALUE=\"Send comment\">\n";
   print "<INPUT TYPE=reset VALUE=\"Reset & Clear\">\n";
   print "</UL></FORM><HR>\n";
   print "<A HREF=\"/blocks/\">[Blocks Home]</A>\n";
   print "<A HREF=\"/blocks-bin/getblock.sh\">[Get Blocks]</A>\n";
   print "<A HREF=\"/blocks/blocks_search.html\">[Block Searcher]</A>\n";
   print "<A HREF=\"/blocks/make_blocks.html\">[Block Maker]</A>\n";
   print "<A HREF=\"/blocks/codehop.html\">[Codehop]</A>\n";
   print "<HR>Page last modified <MODIFICATION_DATE>Feb 2008</MODIFICATION_DATE></HTML>\n";
   exit(0);
} # end of start form

#-------------------------------------------------------------------------
#	Edit
if ($names{randchk} ne $names{randnum}) 
{   
#  print "<P><B>ERROR:</B> Incorrect spam number.<P>\n";
   print "<A HREF=\"/blocks/\">[Blocks Home]</A>\n";
   exit(-1);    
}
#If they waited too long, exit
if ((time() - $names{gmtsecs}) > $maxsecs)
{   
   print "<P><B>Timeout</B><P>\n";
   exit(-1);    
}

$nerr = 0;
if ($names{email} eq "") {
   print "<H1>Error</H1> Please enter your email address.<P>\n";
   $nerr++;
}
if ($names{comment} eq "") {
   print "<H1>Error</H1> Please enter your comment.<P>\n";
   $nerr++;
}
if ($nerr > 0) { exit(-1); }
#Check for known spammers
foreach $spam (@spam_list)
{
   if ($names{email} =~ m/$spam/) { exit(-1); }
}

$com = "$tmp/$$.com";
open(BLK, ">$com");
print BLK "REMOTE_HOST=$ENV{\"REMOTE_HOST\"}\n";
print BLK "REMOTE_USER=$ENV{\"REMOTE_USER\"}\n";
print BLK "REMOTE_ADDR=$ENV{\"REMOTE_ADDR\"}\n";
print BLK "REMOTE_IDENT=$ENV{\"REMOTE_IDENT\"}\n";
print BLK "email=$names{email}\n";
print BLK "name=$names{name}\n";
print BLK "topic=$names{topic}\n";
print BLK "$names{comment}\n";
close(BLK);

if ($names{topic} eq "si")
{
   system("$mailprog -s \"comment $names{topic}\" -r $names{email} $siftmaster < $com");
}
elsif ($names{topic} eq "la")
{
   system("$mailprog -s \"comment $names{topic}\" -r $names{email} $lamamaster < $com");
}
else
{
   system("$mailprog -s \"comment $names{topic}\" -r $names{email} $webmaster < $com");
}
print "<P>Thank you for your comment.<P>\n";
print "<A HREF=\"/blocks/\">[Blocks Home]</A>\n";
print "<A HREF=\"/blocks-bin/getblock.sh\">[Get Blocks]</A>\n";
print "<A HREF=\"/blocks/blocks_search.html\">[Block Searcher]</A>\n";
print "<A HREF=\"/blocks/make_blocks.html\">[Block Maker]</A>\n";
print "<A HREF=\"/blocks/codehop.html\">[Codehop]</A>\n";


#-------------------------------------------------------------------------
exit (0);

#-------------------------------------------------------------------------
#
# parameter: a string that is the html QUERY_STRING environment variable
# returns: an associative array of name/value pairs.  The name is the key.
sub parse_query {
  local($query_string) = @_;
  local(%ans, @q, $pair);
#print "<PRE>qs=$query_string\n</PRE>";
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

# parameter: a hex representation of a number (doesn't need to be a string)
# returns: the decimal representation of the number
sub hextodec {
  unpack("N", pack("H8", substr("0" x 8 . shift, -8)));
}



