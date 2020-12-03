#!/usr/bin/perl

#4/13/07 Execute getblock.pl instead of getblock.sh

print "Content-type: text/html\n\n";
print "<HTML><BASE HREF=\"http://blocks.fhcrc.org\">\n";
print "<TITLE>Block Searcher Raw Blimps Results</TITLE>\n";
print "<H1>Block Searcher Raw Blimps Results</H1>\n";

$infile = STDIN;
if (@ARGV > 0) {
  unless (open($infile, @ARGV[0])) {
    print "Cannot open @ARGV[0]: $!\n<BR>";
    print "Output file has probably expired.\n";
    exit(-1);
  }
}

print"<PRE>\n";
while ($_ = <$infile>) {
  # to replace the html special characters with escape sequences
  s/&/&amp/g ;
  s/>/&gt/g ;
  s/</&lt/g ;
  s/"/&quot/g ;
  if (/IPB\d\d\d\d\d\d\S*/) {
    s|(IPB\d\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  if (/BL\d\d\d\d\d\S*/) {
    s|(BL\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  if (/DM\d\d\d\d\d\S*/) {
    s|(DM\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  if (/PD\d\d\d\d\d\S*/) {
    s|(PD\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  if (/BP\d\d\d\d\d\S*/) {
    s|(BP\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  if (/PR\d\d\d\d\d\S*/) {
    s|(PR\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  if (/PF\d\d\d\d\d\S*/) {
    s|(PF\d\d\d\d\d)(\S*)|<A HREF="/blocks-bin/getblock.pl?\1#\1\2">\1\2</A>|g
  }
  print;
}
print "</PRE><P><A HREF=\"/blocks\">[Blocks Home]</A>\n";
print "</HTML>";
exit(0);
