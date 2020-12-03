#!/usr/bin/perl

#	htmlize-blkprob.pl
#  Note: only executed by the email server, WWW server executes
#	htmlize-blksort.pl; this version uses absolute addresses
#  5/ 3/03 Changes to work with perl5
#  5/ 3/03 html changes to allow independent execution
#  4/13/07 Execute getblock.pl instead of getblock.sh

print "Content-type: text/html\n\n";
print "<HTML>";
print "<TITLE>Block Searcher Results</TITLE>\n";
print "<H1>Block Searcher Results</H1>\n";

$infile = STDIN;
if (@ARGV > 0) 
{
  unless (open($infile, @ARGV[0]))
  {
    print "Cannot open @ARGV[0]: $!\n<BR>";
    print "Output file has probably expired\n";
    exit(-1);
  }
}

#
print "<A HREF=\"#bkhits\">Go to hits</A>\n";
print "<BASE HREF=\"http://blocks.fhcrc.org\">\n";

# put the marker for the blksort introduction
print "<A NAME=\"bkintro\"><H3>Introduction</H3></A><P><PRE>\n";

# output the introduction until the Query= line
while ($_ = <$infile>) {
  last if (/^Query=/);
  # to replace the html special characters with escape sequences
  s/&/&amp;/g ;
  s/>/&gt;/g ;
  s/</&lt;/g ;
  s/"/&quot;/g ;
  print;
}

# put the marker for the blksort results
print "</PRE><A HREF=\"#top\">[Return to top]</A><P><PRE>\n";
print "<A NAME=\"bkhits\"><H3>Hits</H3></A><P>\n" ;

# to replace the html special characters with escape sequences
s/&/&amp;/g ;
s/>/&gt;/g ;
s/</&lt;/g ;
s/"/&quot;/g ;
print;  # must get the line missed

while ($_ = <$infile>) {
  # to replace the html special characters with escape sequences
  s/&/&amp;/g ;
  s/>/&gt;/g ;
  s/</&lt;/g ;
  s/"/&quot;/g ;
  if (/IPB\d\d\d\d\d\d\S*/) {
    s|(IPB\d\d\d\d\d\d)(\S*)|<A HREF="http://blocks.fhcrc.org/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  elsif (/BL\d\d\d\d\d\S*/) {
    s|(BL\d\d\d\d\d)(\S*)|<A HREF="http://blocks.fhcrc.org/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  elsif (/PR\d\d\d\d\d\S*/) {
    s|(PR\d\d\d\d\d)(\S*)|<A HREF="http://blocks.fhcrc.org/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  elsif (/PD\d\d\d\d\d\S*/) {
    s|(PD\d\d\d\d\d)(\S*)|<A HREF="http://blocks.fhcrc.org/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  elsif (/BP\d\d\d\d\d\S*/) {
    s|(BP\d\d\d\d\d)(\S*)|<A HREF="http://blocks.fhcrc.org/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  elsif (/DM\d\d\d\d\d\S*/) {
    s|(DM\d\d\d\d\d)(\S*)|<A HREF="http://blocks.fhcrc.org/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  elsif (/PF\d\d\d\d\d\S*/) {
    s|(PF\d\d\d\d\d)(\S*)|<A HREF="http://blocks.fhcrc.org/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  print;
}  # end of $infile

print "</PRE><A HREF=\"#top\">[Return to top]</A>   ";
print "<A HREF=\"http://blocks.fhcrc.org/blocks\">[Blocks Home]</A>\n";
print "</HTML>";
exit(0);
