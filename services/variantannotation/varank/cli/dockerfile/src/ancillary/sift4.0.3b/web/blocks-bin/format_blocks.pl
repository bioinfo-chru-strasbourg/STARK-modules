#!/usr/bin/perl
#		format_blocks.pl ID
#	Format a file of raw blocks and display them with links
#	Looks for blocks in $tmp/ID.blks and seqs in $tmp/ID.seqs
#--------------------------------------------------------------------------
#  2/25/07 From process_blocks.pl
#--------------------------------------------------------------------------
#
if (@ARGV >= 0) {  $ID = $ARGV[0]; }
else { exit(-1); }

$bin = ".";
$tmp = "../tmp";

#		These files must exist
$blks = "$tmp/$ID.blks";
$seqs = "$tmp/$ID.seqs";

$cblks = "$tmp/$ID.cblks";
$err = "$tmp/$ID.err";
$out = "$tmp/$ID.out";
#		These files are left
$wblks = "$tmp/$ID.wblks";
$tree = "$tmp/$ID.treefile";
#	htmlize-map.pl expects this format!
$map = "$wblks.map";
$pssm = "$tmp/$ID.pssm";
$mast = "$tmp/$ID.mast";

#	Be sure files and directories are created with rw-rw---- permissions
#system("umask 006");

# output the beginning text to be used on all pages
print "Content-type: text/html\n\n";
print "<TITLE>Block Information</TITLE>\n";

#	Check that $blks has non-zero size
if (-s "$blks") 
{
   #  calibrate the blocks for diy searching
   system("$bin/calibrate.csh $blks $cblks > /dev/null");
   #  create a file of PB-weighted blocks
   # NOTE: Should give them the option of using existing weights (if any)
   #       or other weighting schemes
   system("$bin/blimps-bin/blweight $cblks $wblks P M > /dev/null");

   #  create the PSSM files == $pssm & $mast
   system("$bin/blimps-bin/blk2pssm $wblks $pssm B > /dev/null");
   system("$bin/blimps-bin/blk2pssm $wblks $mast M > /dev/null");

   #  create the treefile == $tree
   system("$bin/maketree.csh $wblks > /dev/null");

   #  create the map file == $map
   system("$bin/blimps-bin/makeblockmap $wblks $map > /dev/null 2>&1");

   #  remove unnecessary files
#  system("rm $cblks");
}

#=========================================================================
#	Assuming here that all the other files got made if $wblks has
#	non-zero size
if (-s "$wblks") 
{
  print "<H1>Block Information</H1>\n";
# print "<STRONG>";
  print "<PRE>";
  open(OUT, "cat $out |");
  while ($_ = <OUT>) { print; }
  close(OUT);	
  print "</PRE>";

  #---------------------------------blocks ------------------------------
  print "<P><I>Sequence-Weighted Blocks:</I><BR>";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/catfile.sh?$wblks\">Blocks Format</A>] ";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/catfile.sh?$pssm\">Blimps PSSM</A>] ";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/catfile.sh?$mast\">MAST PSSM</A>] ";
  print "[<A  TARGET=\"help\"HREF=\"/blocks/help/PSSM_def.html\">About PSSMs</A>] ";
  print "<BR>";

  #---------------------------------sequences----------------------------
  if (-s "$seqs") 
  {
     print "<P><I>Sequences:</I><BR>";
     print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/catfile.sh?$seqs\">Sequences</A>] ";
     print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/make_blocks.sh?$seqs\">Block Maker</A>] ";
     print "[<A TARGET=\"help\" HREF=\"/blocks/help/about_blocks.html\">About Block Maker</A>]\n";
     print "<BR>";
  }

  #------------------------------------map-------------------------------
  # The blockmap name must match what's in the mapfile
  if (-s "$map") 
  {
  open(GREP, "grep \"^>\" $map |");
  $map_rec = <GREP>;
  ($mapname) = $map_rec =~ m/^>(\S+)/;
  print "<P><I>Map:     </I><BR>";
  print "
<FORM METHOD=POST ACTION=\"/blocks-bin/blockmap.pl\">
<INPUT TYPE=hidden NAME=\"name\" VALUE=\"$mapname\">
<INPUT TYPE=hidden NAME=\"type\" VALUE=\"FAM\">
<INPUT TYPE=hidden NAME=\"condensed\" VALUE=\"YES\">
<INPUT TYPE=hidden NAME=\"dbtype\" VALUE=\"USER\">
<INPUT TYPE=hidden NAME=\"map\" VALUE=\"";
  open(MAP, "<$map");
  while ($maprec = <MAP>)
  { print $maprec; }
  close(MAP);
  print "\"><INPUT TYPE=submit VALUE=\"Graphical Map\"></FORM>  \n";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/map.csh?$map\">Text Map</A>]  \n";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks/tmp/$wblks.map\">Map Positions</A>]  \n";
  print "[<A HREF=\"/blocks/help/about_maps.html\"TARGET=\"help\">About Maps</A>]";
  }

  #-----------------------------------logos------------------------------
  print "<P><I>Logos:     </I><BR>";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/logo.csh?$wblks+ps\">Postscript</A>] ";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/logo.csh?$wblks+pdf\">PDF</A>] ";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/logo.csh?$wblks+gif\">GIF</A>] ";
  print "[<A TARGET=\"help\" HREF=\"/blocks/help/about_logos.html\">About logos</A>]\n ";

  #-----------------------------------tree------------------------------
  print "<P><I>Tree:     </I><BR>";
  if (-s "$tree") 
  {
  print "
<FORM METHOD=POST ACTION=\"http://www.proweb.org/proweb-bin/trees.cgi\"
 TARGET=\"_blank\">
<INPUT TYPE=hidden NAME=\"format\" VALUE=\"png\">
<INPUT TYPE=hidden NAME=\"length\" VALUE=\"yes\">
<INPUT TYPE=hidden NAME=\"tree\" VALUE=\"";
  open(TRE, "<$tree");
  while ($trerec = <TRE>)
  { print $trerec; }
  close(TRE);
  print "\"><INPUT TYPE=hidden NAME=\"blocks\" VALUE=\"";
  open(BLK, "<$wblks");
  while ($blkrec = <BLK>)
  { print $blkrec; }
  close(BLK);
  print "\"><INPUT TYPE=hidden NAME=\"seqs\" VALUE=\"";
  open(SEQ, "<$seqs");
  while ($seqrec = <SEQ>)
  { print $seqrec; }
  close(SEQ);
  print "\"><INPUT TYPE=submit VALUE=\"ProWeb TreeViewer\"></FORM>\n";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/maketree.csh?$wblks+xbm\">XBitmap</A>] ";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/maketree.csh?$wblks+ps\">Postscript</A>] ";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/maketree.csh?$wblks+pdf\">PDF</A>] ";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/maketree.csh?$wblks+gif\">GIF</A>] ";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/maketree.csh?$wblks+new\">NEW</A>] ";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/maketree.csh?$wblks+dat\">Data</A>] ";
  print "[<A TARGET=\"help\" HREF=\"/blocks/help/about_trees.html\">About trees</A>]\n";
  }   # end of if tree file 

  #-----------------------------------structures------------------------
#	Don't have a cobbler sequence to send to 3dblocks.cgi
# print "&SEQ=";
  print "<P><I>Structures (takes several minutes):</I><BR>";
  print "
<FORM METHOD=POST ACTION=\"http://www.proweb.org/proweb-bin/3dblocks.cgi\"
 TARGET=\"3dblocks\">
<INPUT TYPE=hidden NAME=\"METHOD\" VALUE=\"Mast\">
<INPUT TYPE=hidden NAME=\"CUTOFF\" VALUE=\"0.0001\">
<INPUT TYPE=hidden NAME=\"SAVE\" VALUE=\"NO\">
<INPUT TYPE=hidden NAME=\"blocks_data\" VALUE=\"";
  open(BLK, "<$wblks");
  while ($blkrec = <BLK>)
  { print $blkrec; }
  close(BLK);
  print "\"><INPUT TYPE=submit VALUE=\"3D Blocks\"></FORM>\n";
  print "  [<A TARGET=\"help\" HREF=\"http://www.proweb.org/3dblocks_intro.html\">About 3D Blocks</A>]</FORM>\n";

  print "<P><I>Multiple alignment search:</I><BR>";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/LAMA_search.sh?$wblks\">LAMA</A>] ";
  print "[<A TARGET=\"help\" HREF=\"/blocks/help/LAMA_help.html\">About LAMA</A>]   ";
  print "<SPACER type=horizontal size=10>";
  $group = "UserBlocks";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/mast.sh?$wblks+$group\">MAST</A>] ";
  print "[<A TARGET=\"help\" HREF=\"http://meme.sdsc.edu/meme/website/mast-intro.html\">About MAST</A>]\n";

  print "<P><I>COBBLER sequence and single sequence search:</I><BR>";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/cobbler.sh?$wblks\">COBBLER</A>] ";
  print "[<A TARGET=\"help\" HREF=\"/blocks/help/about_cobbler.html\">About COBBLER</A>]\n";

  #-----------------------------------searches------------------------
  print "<P><I>Search these blocks:</I><BR>";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/diy.sh?$wblks\">DIY Search</A>] ";
  print "[<A TARGET=\"help\" HREF=\"/blocks/help/about_diy.html\">About DIY</A>]\n";

  print "<P><I>Primers:</I><BR>";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/codehop.sh?$wblks\">CODEHOP</A>] ";
  print "[<A TARGET=\"help\" HREF=\"/blocks/help/CODEHOP/CODEHOP_help.html\">About CODEHOP</A>]\n";

  print "<P><I>Substitutions:</I><BR>";
  print "[<A TARGET=\"_blank\" HREF=\"/blocks-bin/sift.sh?$wblks\">SIFT</A>] ";
  print "[<A TARGET=\"help\" HREF=\"/sift/SIFT_on_blocks.html\">About SIFT on Blocks</A>]\n";

}
else 
{ 
  print "<PRE>Problem processing input:<BR>";
  open(ERR, "cat $err |");
  while ($_ = <ERR>) { print; }
  close(ERR);	
  print "\n<P><\PRE>Please refer to "; 
  print "[<A TARGET=\"_blank\" HREF=\"/blocks/blocks_format.html\">acceptable formats</A>].";
}

print "<P>\n";

exit(0);
