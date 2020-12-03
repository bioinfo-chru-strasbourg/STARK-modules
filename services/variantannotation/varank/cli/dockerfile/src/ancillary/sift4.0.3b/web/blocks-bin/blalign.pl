#!/usr/bin/perl

# Execute blalign.csh & process the output

#	Set file permissions to rw-rw----
system("umask 006");

$program = "./blalign.csh";
$tmp = "../tmp";
$bin = "./blimps-bin";

select(STDOUT); $| = 1;

print "Content-type: text/html\n\n";
print "<TITLE>Re-format Blocks as Alignment</TITLE>\n";
print "<A NAME=top><H1>Re-format Blocks as Alignment</H1></A>";
print "<PRE>";


if ( $ENV{"REQUEST_METHOD"} ne "POST" ) {
   print "This script should be referenced with a METHOD of POST\n";
   exit;
}

read (STDIN, $QUERY_STRING, $ENV{"CONTENT_LENGTH"});
%names = &parse_query($QUERY_STRING);


#	Get the blocks & write them to a file
if ($names{blocks} eq "") {
   print "<H1>Error</H1> Please enter some blocks.<P>\n";
   exit;
}
$blk = "$tmp/$$.blk";
open(BLK, ">$blk");
print BLK $names{blocks};
close(BLK);

#	Get the output format type
if ($names{style} eq "") {
   print "<H1>Error</H1> Please select an output format.<P>\n";
   exit;
}
#print "style=$names{style}\n";

#=========================================================================
#	Run blalign now

system("$bin/blalign $blk $names{style} > $tmp/$$.out");

select(STDOUT); $| = 1;

open(OUT, "< $tmp/$$.out");
while ($_ = <OUT>) { print; }
close(OUT);


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



