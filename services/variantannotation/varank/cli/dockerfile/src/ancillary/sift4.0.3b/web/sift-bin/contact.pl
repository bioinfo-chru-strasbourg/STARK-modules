#!/usr/local/bin/perl

# This program is licensed to you under the Fred
# Hutchinos Cancer Research Center (FHCRC)
# NONCOMMERICAL LICENSE.  A copy of the license may be found at
# http://blocks.fhcrc.org/sift/license.html and should be attached
# to this software.

# contact.pl

#	Set file permissions to rw-rw----
system("umask 006");
$siftmaster1="smurphy\@jcvi.org";
$siftmaster2="png\@jcvi.org";
$tmp = "/opt/www/sift/tmp";
$maxsecs = 300;		# 5 minutes
@colors = ("red","blue","green","maroon","purple","lime");
@spam_list = ("ku.name");

select(STDOUT); $| = 1;
print "Content-type: text/html\n\n";
print "<HTML><TITLE>SIFT Contact</TITLE>\n";
print "<H1 ALIGN=CENTER>
SIFT Contact Form
</H1><HR>\n";

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
   print "<FORM METHOD=\"POST\" ACTION=\"/sift-bin/contact.pl\">";
   print "<INPUT TYPE=hidden name=form VALUE=\"form1\">";
   print "<INPUT TYPE=hidden name=randnum VALUE=\"$randnum\">";
   print "<INPUT TYPE=hidden name=gmtsecs VALUE=\"$gmtsecs\">\n";
   print "<P>We will try to respond to your questions within one working day.<BR>\n";
   print "<font size=5 color=$color>$randnum</font>\n";
   print "<UL><P><LI>To discourage spammers, please type in the <font color=$color>$color</font> number printed above <font color=red>(required)</font>: <INPUT SIZE=5 NAME=\"randchk\"><P>\n";
   print "<LI>Please enter your email address <font color=red>(required)</font>: <INPUT SIZE=50 NAME=\"email\"><P>\n";
   print "<LI>Please enter your name: <INPUT SIZE=50 NAME=\"name\">\n";
   print "<LI><P>Enter your question or comment <font color=red>(required)</font>:\n";
   print "<TEXTAREA NAME=\"comment\" ROWS=10 COLS=60></TEXTAREA>\n";
   print "<LI><INPUT TYPE=submit VALUE=\"Send comment\">\n";
   print "<INPUT TYPE=reset VALUE=\"Reset & Clear\">\n";
   print "</UL></FORM><HR>\n";
   print "<A HREF=\"/\">[SIFT Home]</A>\n";
   print "<HR>Page last modified <MODIFICATION_DATE>Nov 2008</MODIFICATION_DATE></HTML>\n";
   exit(0);
} # end of start form

#-------------------------------------------------------------------------
#	Edit
if ($names{randchk} ne $names{randnum}) 
{   
#  print "<P><B>ERROR:</B> Incorrect spam number.<P>\n";
   print "<A HREF=\"/\">[SIFT Home]</A>\n";
   exit(-1);    
}
#If they waited too long, exit
if ((time() - $names{gmtsecs}) > $maxsecs)
{   
   print "<P><B>Timeout</B><P>\n";
   exit(-1);    
}

$nerr = 0;
if ($names{email} eq "" || $names{comment} eq ""){
	print "<H1>Errors</H1> Please press the back button on your browser and resolve the following errors.<P>\n";
}
if ($names{email} eq "") {
   print "<H3>Error</H3> Email address field cannot be left empty.<P>\n";
   $nerr++;
}
if ($names{comment} eq "") {
   print "<H3>Error</H3> Comment field cannot be neft empty.<P>\n";
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
print BLK "randchk=$names{randchk}\n";
print BLK "REMOTE_HOST=$ENV{\"REMOTE_HOST\"}\n";
print BLK "REMOTE_USER=$ENV{\"REMOTE_USER\"}\n";
print BLK "REMOTE_ADDR=$ENV{\"REMOTE_ADDR\"}\n";
print BLK "REMOTE_IDENT=$ENV{\"REMOTE_IDENT\"}\n";
print BLK "email=$names{email}\n";
print BLK "name=$names{name}\n";
print BLK "topic=$names{topic}\n";
print BLK "$names{comment}\n";
close(BLK);
system("mutt -F /opt/www/sift/htdocs/.muttrc -s \"SIFT COMMENT for $$\" $siftmaster1 $siftmaster2 < $com");
print "<P>Thank you for your comment.<P>\n";
print "<A HREF=\"/\">[SIFT Home]</A>\n";

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



