#!/usr/bin/perl
#	htmlize-swiss.pl <blimps or mast output file>

$infile = NULL;
if (@ARGV == 0) {
  $infile = STDIN;
}
elsif (@ARGV > 1) {
  print "Usage:\nhtmlize-swiss.pl [filename]\n";
}
else {
  unless (open($infile, @ARGV[0])) {
    print STDERR "Can't open @ARGV[0]: $!\n";
    return;
  }
}

print "Content-type: text/html\n\n";
print "<HTML><BASE HREF=\"http://blocks.fhcrc.org\">\n";
print "<TITLE>Raw Blimps Results</TITLE>\n";
print "<H1>Raw Blimps Results: PSSM vs SWISS/TREMBL</H1>\n";
print"<PRE>\n";

#	Print until the search results start
$_ = <$infile>;
until (eof($infile) || /^AC#/ || /^----------/)
{
   print;
   $_ = <$infile>;
}
print;

#	This puts links on some odd stuff in the mast output ...
while ($_ = <$infile>) 
{
  @bar = split(/\s+/, $_);
  @ids = split(/\|/, $bar[0]);
  if ($_ =~ m/^\w+\|\w+/) 
  {
    s!(^\w+\|)(\w+)!\1<A HREF="http://expasy.cbr.nrc.ca/cgi-bin/niceprot.pl?\2">\2</A>!
  }
  elsif ($_ =~ m/\w+\s+/)
  {
    s!(^\w+)!<A HREF="http://expasy.cbr.nrc.ca/cgi-bin/niceprot.pl?\1">\1</A>!
  }
  print;
}
#print "</PRE><A HREF=\"#top\">[return to top]</A>   ";
print "</PRE><P><A HREF=\"/blocks/\">[Blocks Home]</A>\n";
print "</HTML>";
exit(0);

