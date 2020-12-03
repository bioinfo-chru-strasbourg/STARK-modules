#!/usr/bin/perl
#
#  5/ 3/03 changes for perl5
#  4/13/07 execute getblock.pl instead of getblock.sh

$infile = STDIN;
if (@ARGV > 0) {
  unless (open($infile, @ARGV[0])) {
    print STDERR "Can't open @ARGV[0]: $!\n";
    exit(-1);
  }
}

print "Content-type: text/html\n\n";
print "<HTML>\n";
print "<TITLE>Block Searcher Results: Old Format</TITLE>\n";
print "<H1>Block Searcher Results: Old Format</H1>\n";

printf("<P><UL>\n");
printf("<LI> <A HREF=\"#bkintro\">Introduction</A>\n");
printf("<LI> <A HREF=\"#bkhits\">Go To Hits</A>\n");
printf("</UL><P>\n");
print "<BASE HREF=\"http://blocks.fhcrc.org\">\n";

# put the marker for the blksort introduction
print "<A NAME=\"bkintro\"><H3>Introduction</H3></A><P>\n" ;

# output the introduction until the Query= line
print "<PRE>\n";
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
    s|(IPB\d\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  elsif (/BL\d\d\d\d\d\S*/) {
    s|(BL\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  elsif (/PR\d\d\d\d\d\S*/) {
    s|(PR\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  elsif (/PD\d\d\d\d\d\S*/) {
    s|(PD\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  elsif (/BP\d\d\d\d\d\S*/) {
    s|(BP\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  elsif (/DM\d\d\d\d\d\S*/) {
    s|(DM\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  elsif (/PF\d\d\d\d\d\S*/) {
    s|(PF\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  print;
}  # end of $infile

print "<A HREF=\"/blocks\">[Blocks Home]</A>\n";
print "</HTML>";
exit(0);
