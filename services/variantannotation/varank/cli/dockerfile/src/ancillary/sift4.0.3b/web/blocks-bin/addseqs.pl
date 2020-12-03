#!/usr/bin/perl

# Execute addseqs.csh & process the output
#	Try to add new sequences to existing blocks

#	Set file permissions to rw-rw----
system("umask 006");

$program = "./addseqs.csh";
$tmp = "../tmp";

select(STDOUT); $| = 1;

print "Content-type: text/html\n\n";
print "<TITLE>Do-It-Yourself Block Search Results</TITLE>\n";
print "<A NAME=top><H1>Do-It-Yourself Block Search Results</H1></A>";
print "<PRE>";


if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
   print "This script should be referenced with a METHOD of POST\n";
   exit;
}

read (STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"});
%names = &parse_query($QUERY_STRING);


#	Get the blocks & write it to a file
if ($names{blocks} eq "") {
   print "<H1>Error</H1> Please enter some blocks.<P>\n";
   exit;
}
$blk = "$tmp/$$.blk";
open(BLK, ">$blk");
print BLK $names{blocks};
close(BLK);

#	Get the sequence & write it to a file
if ($names{sequence} eq "") {
   print "<H1>Error</H1> Please enter some protein sequences.<P>\n";
   exit;
}
$seq = "$tmp/$$.seq";
open(SEQ, ">$seq");
print SEQ $names{sequence};
close(SEQ);

#=========================================================================
#	Run shell now

$addseqs = "$tmp/$$.addseqs";
open(ERR, "$program $blk $seq $addseqs 2>&1 |");
wait;
while ($_ = <ERR>) { print; }
close(ERR);

select(STDOUT); $| = 1;

print "</PRE><A HREF=\"/blocks/tmp/$$.addseqs\">[Addseqs Output]</A><P>\n";
print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/diy.sh?$blk.add\">DIY Search</A>] ";
print "[<A TARGET=\"blocks\" HREF=\"/blocks/help/about_diy.html\">About DIY</A>]<P>\n";
print "[<A TARGET=\"blocks\" HREF=\"/blocks-bin/addseqs.sh?$blk.add\">Add Sequences</A>] ";
print "[<A TARGET=\"blocks\" HREF=\"/blocks/help/about_addseqs.html\">About Add Sequences</A>]<P>\n";
print "<HR><PRE>\n";


open(ADD, "<$blk.add");
while ($_ = <ADD>) { print; }
close(ADD);


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

# parameter: a hex representation of a number (doesn't need to be a string)
# returns: the decimal representation of the number
sub hextodec {
  unpack("N", pack("H8", substr("0" x 8 . shift, -8)));
}



