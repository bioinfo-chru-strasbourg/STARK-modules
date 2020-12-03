#!/usr/bin/perl
#

print "Content-type: text/html\n\n";

#	Print sequence for blocks search
#	Seems to only work with pure sequence; no title, no spaces, etc.

print "<P>These work.<BR>\n";
print "[<A HREF=\"http://www.proweb.org/proweb-bin/blockmap.cgi?type=FAM&dbtype=IPB&name=IPB000072&condensed=YES\">blockmap FAM</A>]\n";

print "[<A HREF=\"http://www.proweb.org/proweb-bin/blockmap.cgi?type=SEQ&dbtype=IPB&name=VEGH_ORFN7&condensed=YES\">blockmap SEQ</A>]\n";

print "<P>Shouldn't this work?<BR>\n";
print "[<A HREF=\"http://www.proweb.org/proweb-bin/blockmap.cgi?type=FAM&dbtype=USER&condensed=YES&map=
>IPB001525 6 222 1761 C5_DNA_meth
C-5 cytosine-specific DNA methylase
user-query 519 6
A 11 24
B 76 91
C 100 127
D 154 180
E 437 454
F 465 477
\">blockmap USER FAM</A>]\n";

print "<P>This doesn't work yet<BR>\n";
print "[<A HREF=\"http://www.proweb.org/proweb-bin/blockmap.cgi?type=SEQ&dbtype=USER&name=user-query&condensed=YES&map=>IPB000953 1 284 5322 Chromo
Chromo domain
user-query 519 1
A 12 49
>IPB001525 6 222 1761 C5_DNA_meth
C-5 cytosine-specific DNA methylase
user-query 519 6
A 11 24
B 76 91
C 100 127
D 154 180
E 437 454
F 465 477
\">blockmap USER SEQ</A>]\n";

exit;

