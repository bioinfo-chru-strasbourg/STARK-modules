#!/usr/bin/perl
#
#	bm_format.pl <id>
#	Formats Block Maker output in <id>
#-----------------------------------------------------------------------
#  5/16/03 Separated from makeblocks.pl
#  8/25/03 Added blalign.sh links
#  8/27/03 Added diy.sh links
#  3/26/04 Run LAMA & SIFT in a blank window
#  2/18/07 blockmap here instead of proweb
#  4/ 9/07 treeviewer here instead of proweb
#-----------------------------------------------------------------------
$tmp = "../tmp/bm";

if (@ARGV >= 0) {  $id = $ARGV[0]; }
else { exit(-1); }

$prefix = "$tmp/$id/$id";

#=========================================================================
# htmlize the output, which is mostly written to separate files

print "Content-type: text/html\n\n";
print "<TITLE>Block Maker results</TITLE>\n";

#print "prefix=$prefix\n";

if (-s "$prefix.err")
{
   print "<PRE>";
   open(ERR, "<$prefix.err");
   while ($err = <ERR>) { print "$err"; $nerr++; }
   close(ERR);
   exit(-1);
}

if (-s "$prefix.warn")
{   print "Check for Warnings before re-submitting"; }

print "<A NAME=top><H1>Block Maker Results</H1></A><P>\n<UL>
<LI> <A TARGET=\"blockmaker\" HREF=\"/blocks-bin/catfile.sh?blockmaker.stp\">Introduction</A>";
print "<LI> <A TARGET=\"blockmaker\" HREF=\"/blocks/help/bm/hints.html\">Hints</A> on saving these results for future use";

if (-s "$prefix.warn")
{
print "<LI><A TARGET=\"blockmaker\" HREF=\"/blocks-bin/catfile.sh?$prefix.warn\">Warnings</A>";
}

print "<P>
<LI> <A HREF=\"#motif\"> BLOCKS from MOTIF</A><BR>
<LI> <A HREF=\"#gibbs\"> BLOCKS from GIBBS</A>";

print "<P>
<LI> <A TARGET=\"blockmaker\" HREF=\"/blocks-bin/bm_map.csh?$prefix.mblks.mapfile+$prefix.gblks.mapfile\">BLOCK Maps</A>
     [<A TARGET=\"blockmaker\" ALIGN=RIGHT HREF=\"/blocks/help/about_maps.html\">  About Block Maps</A>]";

print "</UL><HR>\n\n";

#=========================================================================
print "<A NAME=\"motif\"><H3>BLOCKS from MOTIF</H3></A><P><PRE>\n";

#	Print the blalign output 
open(ALN, "<$prefix.maln");
while ($alnrec = <ALN>)
{ print $alnrec; }
close(ALN);

if (-s "$prefix.mblks") 
{
# print "</PRE><TABLE border cellspacing=0>";
  print "</PRE><TABLE>";

  print "<TR><TH><A TARGET=\"blockmaker\" HREF=\"/blocks/help/blocks_format.html\">Formatted BLOCKS</A></TH>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/catfile.sh?$prefix.mblks\">BLOCKS format</A>]</TD>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/catfile.sh?$prefix.mmast\" TARGET=\"blockmaker\">MAST Searchable format</A>]</TD>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/catfile.sh?$prefix.mcob\" TARGET=\"blockmaker\">COBBLER Sequence</A>]</TD>";
  print "<TD>[<A TARGET=\"mblks\" ALIGN=RIGHT HREF=\"/blocks/help/about_cobbler.html\">About COBBLER</A>]</TD>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/blalign.sh?$prefix.mblks\">Re-format</A>]</TD>";
  print "</TR>\n";

  print "<TR><TH><A TARGET=\"blockmaker\" HREF=\"/blocks/help/about_logos.html\">Logos</A></TH>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/logo.csh?$prefix.mblks+ps\">Postscript</A>]</TD>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/logo.csh?$prefix.mblks+pdf\">PDF</A>]</TD>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/logo.csh?$prefix.mblks+gif\">GIF</A>]</TD>";
  print "</TR>\n";
#
#<FORM METHOD=POST ACTION=\"http://www.proweb.org/proweb-bin/blockmap.cgi\" TARGET=\"mblks\">
  # The blockmap name must match what's in the mapfile
  open(GREP, "grep \"^>\" $prefix.mblks.mapfile |");
  $map_rec = <GREP>;
  ($mapname) = $map_rec =~ m/^>(\S+)/;
  print "<TR><TH><A HREF=\"/blocks/help/about_maps.html\"TARGET=\"blockmaker\">Maps</A></TH>";
  print "<TD VALIGN=bottom>
<FORM METHOD=POST ACTION=\"/blocks-bin/blockmap.pl\" TARGET=\"mblks\">
<INPUT TYPE=hidden NAME=\"name\" VALUE=\"$mapname\">
<INPUT TYPE=hidden NAME=\"type\" VALUE=\"FAM\">
<INPUT TYPE=hidden NAME=\"condensed\" VALUE=\"YES\">
<INPUT TYPE=hidden NAME=\"dbtype\" VALUE=\"USER\">
<INPUT TYPE=hidden NAME=\"map\" VALUE=\"";
  open(MAP, "<$prefix.mblks.mapfile");
  while ($maprec = <MAP>)
  { print $maprec; }
  close(MAP);
  print "\"><INPUT TYPE=submit VALUE=\"Graphical Map\"></FORM>\n";
  print "</TD>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/map.csh?$prefix.mblks.mapfile\">Text Map</A>]</TD>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks/tmp/bm/$id/$id.mblks.mapfile\">Map Positions</A>]</TD>";
  print "</TR>\n";

#
  print "<TR><TH><A TARGET=\"blockmaker\" HREF=\"/blocks/help/about_trees.html\">Tree</A></TH>";
  print "<TD VALIGN=bottom>";
# print "<FORM METHOD=POST ACTION=\"http://www.proweb.org/proweb-bin/trees.cgi\"
  print "<FORM METHOD=POST ACTION=\"/blocks-bin/treeviewer.pl\"
 TARGET=\"mblks\">
<INPUT TYPE=hidden NAME=\"format\" VALUE=\"png\">
<INPUT TYPE=hidden NAME=\"length\" VALUE=\"yes\">
<INPUT TYPE=hidden NAME=\"tree\" VALUE=\"";
  open(TRE, "<$prefix.mtree");
  while ($trerec = <TRE>)
  {
#    $trerec =~ s/(\W)/sprintf("%%%02X", ord($1))/eg;
     print $trerec;
  }
  close(TRE);
  print "\"><INPUT TYPE=hidden NAME=\"blocks\" VALUE=\"";
  open(BLK, "<$prefix.mblks");
  while ($blkrec = <BLK>)
  { print $blkrec; }
  close(BLK);
  print "\">";
  print "<INPUT TYPE=hidden NAME=\"seqs\" VALUE=\"";
  open(SEQ, "<$prefix.pros");
  while ($seqrec = <SEQ>)
  { print $seqrec; }
  close(SEQ);
  print "\"><INPUT TYPE=submit VALUE=\"ProWeb TreeViewer\"></FORM>\n";
  print "</TD>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/drawgram.csh?$prefix.mtree+dat\">Data</A>]</TD>";
# print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/drawgram.csh?$prefix.mtree+xbm\">XBitmap</A>]</TD>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/drawgram.csh?$prefix.mtree+ps\">Postscript</A>]</TD>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/drawgram.csh?$prefix.mtree+pdf\">PDF</A>]</TD>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/drawgram.csh?$prefix.mtree+gif\">GIF</A>]</TD>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/drawgram.csh?$prefix.mtree+new\">Newick</A>]</TD>";
  print "</TR>\n";

#
  print "<TR><TH><A TARGET=\"blockmaker\" HREF=\"http://www.proweb.org/3dblocks_intro.html\">Structures</A></TH>";
  print "<TD VALIGN=bottom>
<FORM METHOD=POST ACTION=\"http://www.proweb.org/proweb-bin/3dblocks.cgi\"
 TARGET=\"mblks\">
<INPUT TYPE=hidden NAME=\"CUTOFF\" VALUE=\"0.0001\">
<INPUT TYPE=hidden NAME=\"SAVE\" VALUE=\"NO\">
<INPUT TYPE=hidden NAME=\"blocks_data\" VALUE=\"\n";
   open(BLK, "<$prefix.mblks");
   while ($blkrec = <BLK>)
   { print $blkrec; }
   close(BLK);
   print "\"><INPUT TYPE=hidden NAME=\"SEQ\" VALUE=\"";
   open(SEQ, "<$prefix.mcob");
   while ($seqrec = <SEQ>)
   { print $seqrec; }
   close(SEQ);
   print "\"><INPUT TYPE=submit VALUE=\"3D Blocks\"></FORM></TD>";
   print "<TD>(takes several minutes)</TD>";
   print "</TR>\n";

#
  print "<TR><TH>Search using BLOCKS as query</TH>";
  print "<TD>[<A TARGET=\"_blank\" HREF=\"/blocks-bin/LAMA_search.sh?$prefix.mblks\">LAMA</A>]</TD>";
  print "<TD>[<A TARGET=\"blockmaker\" HREF=\"/blocks/help/LAMA_help.html\">About LAMA</A>]</TD>";
  $group = "MotifBlocks";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/bm_mast.sh?$prefix.mmast+$group\">MAST</A>]</TD>";
  print "<TD>[<A TARGET=\"blockmaker\" HREF=\"http://www.sdsc.edu/MEME/meme/website/mast-intro.html\">About MAST</A>]</TD>";
  print "</TR>\n";

#
  print "<TR><TH>Search using COBBLER sequence as query</TH>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&AUTO_FORMAT=Semiauto&PROGRAM=blastp&DATABASE=nr&QUERY=";
  open(SEQ, "<$prefix.mcob");
  while ($seqrec = <SEQ>)
  { 
     #   don't print the title line
     if (!($seqrec =~ m/^>/))  { print $seqrec; }
  }
  close(SEQ);
  print "\">BLAST</A>]</TD>";
#
  print "<TD>[<A TARGET=\"mblks\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?AUTO_FORMAT=Semiauto-&PROGRAM=blastp&RUN_PSIBLAST=on&CDD_SEARCH=on&COMPOSITION_BASED_STATISTICS=on&DESCRIPTIONS=250&ALIGNMENTS=100&ALIGNMENT_VIEW=Pairwise&DATABASE=nr&END_OF_HTTPGET=YES&QUERY=";
  open(SEQ, "<$prefix.mcob");
  while ($seqrec = <SEQ>)
  {
     #   don't print the title line
     if (!($seqrec =~ m/^>/))  { print $seqrec; }
  }
  close(SEQ);
  print "\">PSI-BLAST</A>]</TD>";
  print "</TR>\n";

  print "<TR><TH>Search these Blocks</TH>\n";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/diy.sh?$prefix.mblks\">DIY Search</A>]</TD></TR>\n";
#
  print "<TR><TH><A TARGET=\"blockmaker\" HREF=\"/blocks/help/CODEHOP/CODEHOP_help.html\">Primers</A></TH>";
  print "<TD>[<A TARGET=\"mblks\" HREF=\"/blocks-bin/codehop.sh?$prefix.mblks\">CODEHOP</A>]</TD>";
#
  print "<TR><TH><A TARGET=\"blockmaker\" HREF=\"/sift/SIFT_on_blocks.html\">Substitutions in Blocks</A></TH>";
  print "<TD>[<A TARGET=\"_blank\" HREF=\"/blocks-bin/sift.sh?$prefix.mblks\">SIFT</A>]</TD>";
  print "</TR>\n";
#

  print "</TABLE>";
}  # end of mblks

print "<P>\n";

#===========================================================================
print "<HR>";
print "<A NAME=\"gibbs\"><H3>BLOCKS from GIBBS</H3></A><P><PRE>\n" ;

#	Print the blalign output 
open(ALN, "<$prefix.galn");
while ($alnrec = <ALN>)
{ print $alnrec; }
close(ALN);

if (-s "$prefix.gblks") 
{
# print "</PRE><TABLE border cellspacing=0>";
  print "</PRE><TABLE>";

  print "<TR><TH><A TARGET=\"blockmaker\" HREF=\"/blocks/help/blocks_format.html\">Formatted BLOCKS</A></TH>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/catfile.sh?$prefix.gblks\">BLOCKS format</A>]</TD>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/catfile.sh?$prefix.gmast\" TARGET=\"blockmaker\">MAST Searchable format</A>]</TD>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/catfile.sh?$prefix.gcob\" TARGET=\"blockmaker\">COBBLER Sequence</A>]</TD>";
  print "<TD>[<A TARGET=\"gblks\" ALIGN=RIGHT HREF=\"/blocks/help/about_cobbler.html\">About COBBLER</A>]</TD>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/blalign.sh?$prefix.gblks\">Re-format</A>]</TD>";
  print "</TR>\n";

  print "<TR><TH><A TARGET=\"blockmaker\" HREF=\"/blocks/help/about_logos.html\">Logos</A></TH>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/logo.csh?$prefix.gblks+ps\">Postscript</A>]</TD>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/logo.csh?$prefix.gblks+pdf\">PDF</A>]</TD>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/logo.csh?$prefix.gblks+gif\">GIF</A>]</TD>";
  print "</TR>\n";
#
  # The blockmap name must match what's in the mapfile
  open(GREP, "grep \"^>\" $prefix.gblks.mapfile |");
  $map_rec = <GREP>;
  ($mapname) = $map_rec =~ m/^>(\S+)/;
  print "<TR><TH><A HREF=\"/blocks/help/about_maps.html\"TARGET=\"blockmaker\">Maps</A></TH>";
  print "<TD VALIGN=bottom>
<FORM METHOD=POST ACTION=\"/blocks-bin/blockmap.pl\" TARGET=\"gblks\">
<INPUT TYPE=hidden NAME=\"name\" VALUE=\"$mapname\">
<INPUT TYPE=hidden NAME=\"type\" VALUE=\"FAM\">
<INPUT TYPE=hidden NAME=\"condensed\" VALUE=\"YES\">
<INPUT TYPE=hidden NAME=\"dbtype\" VALUE=\"USER\">
<INPUT TYPE=hidden NAME=\"map\" VALUE=\"";
  open(MAP, "<$prefix.gblks.mapfile");
  while ($maprec = <MAP>)
  { print $maprec; }
  close(MAP);
  print "\"><INPUT TYPE=submit VALUE=\"Graphical Map\"></FORM>\n";
  print "</TD>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/map.csh?$prefix.gblks.mapfile\">Text Map</A>]</TD>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks/tmp/bm/$id/$id.gblks.mapfile\">Map Positions</A>]</TD>";
  print "</TR>\n";

#
  print "<TR><TH><A TARGET=\"blockmaker\" HREF=\"/blocks/help/about_trees.html\">Tree</A></TH>";
  print "<TD VALIGN=bottom>";
  print "<FORM METHOD=POST ACTION=\"/blocks-bin/treeviewer.pl\"
 TARGET=\"gblks\">
<INPUT TYPE=hidden NAME=\"format\" VALUE=\"png\">
<INPUT TYPE=hidden NAME=\"length\" VALUE=\"yes\">
<INPUT TYPE=hidden NAME=\"tree\" VALUE=\"";
  open(TRE, "<$prefix.gtree");
  while ($trerec = <TRE>)
  {
#    $trerec =~ s/(\W)/sprintf("%%%02X", ord($1))/eg;
     print $trerec;
  }
  close(TRE);
  print "\"><INPUT TYPE=hidden NAME=\"blocks\" VALUE=\"";
  open(BLK, "<$prefix.gblks");
  while ($blkrec = <BLK>)
  { print $blkrec; }
  close(BLK);
  print "\">";
  print "<INPUT TYPE=hidden NAME=\"seqs\" VALUE=\"";
  open(SEQ, "<$prefix.pros");
  while ($seqrec = <SEQ>)
  { print $seqrec; }
  close(SEQ);
  print "\"><INPUT TYPE=submit VALUE=\"ProWeb TreeViewer\"></FORM>\n";
  print "</TD>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/drawgram.csh?$prefix.gtree+dat\">Data</A>]</TD>";
# print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/drawgram.csh?$prefix.gtree+xbm\">XBitmap</A>]</TD>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/drawgram.csh?$prefix.gtree+ps\">Postscript</A>]</TD>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/drawgram.csh?$prefix.gtree+pdf\">PDF</A>]</TD>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/drawgram.csh?$prefix.gtree+gif\">GIF</A>]</TD>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/drawgram.csh?$prefix.gtree+new\">Newick</A>]</TD>";
  print "</TR>\n";

#
  print "<TR><TH><A TARGET=\"blockmaker\" HREF=\"http://www.proweb.org/3dblocks_intro.html\">Structures</A></TH>";
  print "<TD VALIGN=bottom>
<FORM METHOD=POST ACTION=\"http://www.proweb.org/proweb-bin/3dblocks.cgi\"
 TARGET=\"gblks\">
<INPUT TYPE=hidden NAME=\"CUTOFF\" VALUE=\"0.0001\">
<INPUT TYPE=hidden NAME=\"SAVE\" VALUE=\"NO\">
<INPUT TYPE=hidden NAME=\"blocks_data\" VALUE=\"\n";
   open(BLK, "<$prefix.gblks");
   while ($blkrec = <BLK>)
   { print $blkrec; }
   close(BLK);
   print "\"><INPUT TYPE=hidden NAME=\"SEQ\" VALUE=\"";
   open(SEQ, "<$prefix.gcob");
   while ($seqrec = <SEQ>)
   { print $seqrec; }
   close(SEQ);
   print "\"><INPUT TYPE=submit VALUE=\"3D Blocks\"></FORM></TD>";
   print "<TD>(takes several minutes)</TD>";
   print "</TR>\n";

#
  print "<TR><TH>Search using BLOCKS as query</TH>";
  print "<TD>[<A TARGET=\"_blank\" HREF=\"/blocks-bin/LAMA_search.sh?$prefix.gblks\">LAMA</A>]</TD>";
  print "<TD>[<A TARGET=\"blockmaker\" HREF=\"/blocks/help/LAMA_help.html\">About LAMA</A>]</TD>";
  $group = "MotifBlocks";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/mast.sh?$prefix.gmast+$group\">MAST</A>]</TD>";
  print "<TD>[<A TARGET=\"blockmaker\" HREF=\"http://www.sdsc.edu/MEME/meme/website/mast-intro.html\">About MAST</A>]</TD>";
  print "</TR>\n";

#
  print "<TR><TH>Search using COBBLER sequence as query</TH>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Web&AUTO_FORMAT=Semiauto&PROGRAM=blastp&DATABASE=nr&QUERY=";
  open(SEQ, "<$prefix.gcob");
  while ($seqrec = <SEQ>)
  { 
     #   don't print the title line
     if (!($seqrec =~ m/^>/))  { print $seqrec; }
  }
  close(SEQ);
  print "\">BLAST</A>]</TD>";
#
  print "<TD>[<A TARGET=\"gblks\" HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?AUTO_FORMAT=Semiauto-&PROGRAM=blastp&RUN_PSIBLAST=on&CDD_SEARCH=on&COMPOSITION_BASED_STATISTICS=on&DESCRIPTIONS=250&ALIGNMENTS=100&ALIGNMENT_VIEW=Pairwise&DATABASE=nr&END_OF_HTTPGET=YES&QUERY=";
  open(SEQ, "<$prefix.gcob");
  while ($seqrec = <SEQ>)
  {
     #   don't print the title line
     if (!($seqrec =~ m/^>/))  { print $seqrec; }
  }
  close(SEQ);
  print "\">PSI-BLAST</A>]</TD>";
  print "</TR>\n";

  print "<TR><TH>Search these Blocks</TH>\n";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/diy.sh?$prefix.gblks\">DIY Search</A>]</TD></TR>\n";
#
#
  print "<TR><TH><A TARGET=\"blockmaker\" HREF=\"/blocks/help/CODEHOP/CODEHOP_help.html\">Primers</A></TH>";
  print "<TD>[<A TARGET=\"gblks\" HREF=\"/blocks-bin/codehop.sh?$prefix.gblks\">CODEHOP</A>]</TD>";
#
  print "<TR><TH><A TARGET=\"blockmaker\" HREF=\"/sift/SIFT_on_blocks.html\">Substitutions in Blocks</A></TH>";
  print "<TD>[<A TARGET=\"_blank\" HREF=\"/blocks-bin/sift.sh?$prefix.gblks\">SIFT</A>]</TD>";
  print "</TR>\n";
#

  print "</TABLE>";

}  # end of gblks

print "<P>\n";

#-------------------------------------------------------------------------
print "<HR><A HREF=\"/blocks/\">[BLOCKS home]</A>\n";
exit(0);

