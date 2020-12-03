#!/usr/bin/perl

# Execute cobbler.csh & process the output

#	Set file permissions to rw-rw----
system("umask 006");

$program = "./cobbler.csh";
$tmp = "../tmp";

select(STDOUT); $| = 1;

print "Content-type: text/html\n\n";
print "<TITLE>Cobbler Results</TITLE>\n";
print "<A NAME=top><H1>Cobbler Results</H1></A>";
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
   print "<H1>Error</H1> Please enter a protein sequence.<P>\n";
   exit;
}
$seq = "$tmp/$$.seq";
open(SEQ, ">$seq");
print SEQ $names{sequence};
close(SEQ);

#=========================================================================
#	Run shell now


$cob = "$tmp/$$.cob";
open(ERR, "$program $blk $seq $cob 2>&1 |");
wait;
while ($_ = <ERR>) { print; }
close(ERR);

select(STDOUT); $| = 1;

open(COB, "<$cob");
@input = <COB>;

print "</PRE>To do a Blast search, copy the cobbler sequence below ";
print "then click on a Blast link\n<P><PRE>";

#	Print sequence for blast search, skip title line
#print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/cgi-bin/BLAST/nph-blast?PROGRAM=blastp&DESCRIPTIONS=500&DATALIB=nr&SEQUENCE=";
#print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/blast.cgi?Jform=0&PROGRAM=blastp&DESCRIPTIONS=500&DATALIB=nr&SEQUENCE=";
#$i = 1;
#while ($i < @input) {
#  print @input[$i];
#  ++$i;
#}
#print "\" TARGET=\"getblock\">Blast Search</A>]</B>\n";

#	Print sequence again for gap-blast search
#print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/cgi-bin/BLAST/nph-newblast?PROGRAM=blastp&GAPPED_ALIGNMENT=is_set&DESCRIPTIONS=500&ALIGNMENTS=100&DATALIB=nr&SEQUENCE=";
#print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/blast.cgi?Jform=1&PROGRAM=blastp&GAPPED_ALIGNMENT=is_set&DESCRIPTIONS=500&ALIGNMENTS=100&DATALIB=nr&SEQUENCE=";
#$i = 1;
#while ($i < @input) {
#  print @input[$i];
#  ++$i;
#}
#print "\" TARGET=\"getblock\">Gap-Blast Search</A>]</B>\n";

#	Print sequence again for psi-blast search
#print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/cgi-bin/BLAST/nph-psi_blast?PROGRAM=blastp&GAPPED_ALIGNMENT=is_set&DESCRIPTIONS=500&ALIGNMENTS=100&DATALIB=nr&SEQUENCE=";
#$i = 1;
#while ($i < @input) {
#  print @input[$i];
#  ++$i;
#}
#print "\" TARGET=\"getblock\">Psi-Blast Search</A>]</B>\n";

#	Print it to the screen, including title line
print "\n\n";
$i = 0;
while ($i < @input) {
  print @input[$i];
  ++$i;
}
close(COB);
print "</PRE>";

print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/blast.cgi?Jform=0";
print "\" TARGET=\"getblock\">Blast Search</A>]</B>\n";
print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/blast.cgi?Jform=1";
print "\" TARGET=\"getblock\">Gap-Blast Search</A>]</B>\n";
print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/psiblast.cgi";
print "\" TARGET=\"getblock\">Psi-Blast Search</A>]</B>\n";

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



