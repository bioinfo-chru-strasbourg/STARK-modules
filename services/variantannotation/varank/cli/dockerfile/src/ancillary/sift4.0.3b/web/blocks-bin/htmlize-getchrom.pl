#!/usr/bin/perl
#	htmlize-getchrom.pl
#  getblock output is assumed to consist of:
#	* Introduction from blksort.stp or prints.stp
#	* The blocks, followed by "%d blocks processed"
#	* The cobbler.pros sequence, followed by "end cobbler"
#	* The block map lines, followed by "end map"
#	* The tree lines, followed by "end tree"
# 	* The proweb.dat entries, followed by "end proweb"
# 	* The blinks.dat entries, followed by "end blinks"
# 	* The blk2pdb.dat entries, followed by "end pdb"
# 	* The cyrca.dat entries, followed by "end cyrca"
#	* prosite entry (prosite.dat, prosite.doc)
#----------------------------------------------------------------
#  2/25/96  Modified for getblock with COBBLER sequence  JGH
#           Got rid of while (1) statements
#  5/29/96  Modified for getblock with PROWEB data.  JGH
#  7/11/96  Added limit on cputime, fixed error in while loop at
#           line 286 which caused looping.  JGH
#  7/12/96  Put eof checks in all while statements JGH
#  8/ 7/96  Revised internal prosite links.  JGH
#  9/ 3/96  Added tree link.  JGH
#  9/11/96  Added frames. JGH
# 11/ 6/96  Added MAST link.  JGH
#  3/ 6/97  Moved sequence links from SwissProt to NCBI to generalize
#	    for Prints sequence IDs. Check for PR or BL block_ACs.
#		logo.sh & LAMA_search still need work ...
#  3/ 9/97  logo.sh just once
#  7/31/97  Added maps (map.sh)
#  8/11/97  Moved logo, lama & mast links to the top
#  8/13/97  Removed frames, use target window for documentation instead.
#	    Use target window for Blast search.
#  9/ 4/97  Added about_trees, about_maps
#  9/ 8/97  Added "Mime" type for trees, newick.sh
#  9/18/97  Added link to Psi-Blast
# 10/ 1/97  Changed Blastp & Psi-Blastp for DESCRIPTIONS=500
# 10/ 9/97  Added link to Emotif
# 10/12/97  Renamed logo.sh logo.ps.sh; moved trees to top
#  5/ 8/98  Added pdb link. 
#  7/ 7/98  Added ProDom.
# 10/ 8/98  Added Domo.
# 12/ 1/98  Pre-made gif logos
# 12/ 6/98  Revised tree display shells (drawgram.csh)
# 12/ 8/98  logo.csh replaces logo.*.csh|sh
# 12/28/98  Don't edit ">", "<", etc. Causes problems with some browsers
#		when <PRE>
#  1/ 7/99  Expect blk2pdb.dat for other than BL; format changes
#  3/22/99  Use sprot links for all but prints
#  3/27/99  Add internal link for group number at top
#  6/11/99  New sequence names: id|ac
#  7/23/99  Changed pdbfile to call pdbfile.pl
#  8/21/99  Prints links to swiss now (23.1); may not be a tree
# 11/ 6/99  Blast links are broken ...
#  4/29/00  Blocks 12 (IPB...)
#  2/ 2/01  Link to proweb tree stuff
#  2/21/01  cyrca stuff
# 11/29/01  Fixed psi-blast link
#  1/ 6/02  Re-arranged tree links
#  1/15/02  Link to blocks_format.html; rearrange links
# 10/ 4/02  ChromDB stuff (dropped pfam, domo, prodom)
# 10/13/02  Changed proweb-bin/trees.cgi to form
# 11/18/02  Fixed sequence links
#  1/ 9/03  TreeViewer changes : trees are now on bateson
#  3/24/03  ProWeb block map
#  3/29/03  New chromdb family names
#  4/10/03  New chromdb links
#  4/11/03  Separate links for chromdb|trembl ids
#  4/30/03  Fix link to Chromdb (requires ../www/chromdb/groups.list)
#----------------------------------------------------------------
#
#link URLs and expressions
$BLKS_LINK = '<A HREF="http://blocks.fhcrc.org/blocks-bin/getblock.pl?' ;
$EMTF_LINK = '<A HREF="http://dna.stanford.edu/cgi-bin/3motif/nph-blocks?blocks=' ;
$PDB_LINK = '<A HREF="http://molbio.info.nih.gov/cgi-bin/moldraw?' ;
#$SEQ_LINK = '<A HREF="http://www3.ncbi.nlm.nih.gov:80/htbin-post/Entrez/query?form=6&db=p&Dopt=g&uid=' ;


#	Prevent this script from running forever - 30 seconds
system("limit cputime 30");

if (@ARGV == 0) {
  $infile = STDIN;
}
elsif (@ARGV > 1) {
  print "Usage:\nhtmlize-block [filename]\n";
}
else {
  unless (open($infile, @ARGV[0])) {
    print STDERR "Can't open @ARGV[0]: $!\n";
    return;
  }
}

@input = <$infile>;

#     If a <FILEHANDLE> is used in a context that is  looking  for
#     an  array,  an  array  consisting  of all the input lines is
#     returned, one line per array element.  It's easy to  make  a
#     LARGE data space this way, so use with care.


#----------------------------------------------------------------
#	Read entire input file, saving some lines in arrays

$past_blocks = 0;
$past_intro = 0;
$past_cobb = 0;
$past_map = 0;
$past_tree = 0;
$past_web = 0;
$nweb = 0;
$past_lnk = 0;
$nlnk = 0;
$past_pdb = 0;
$past_cyrca = 0;
$npdb = 0;
$ncyrca = 0;
$ntree = 0;

#	Write the blocks for blk2pssm to use later for MAST
open(FBLK, ">../tmp/$$.blocks");
open(FMAP, ">../tmp/$$.blocks.mapfile");
open(FTREE, ">../tmp/$$.treefile");
open(FPDB, ">../tmp/$$.pdbfile");
open(FCYR, ">../tmp/$$.cyrca");
system("chmod -f 660 ../tmp/$$.*");
#
for ($j = 0; $j < @input; $j++)
{
  # beyond the blocks yet?
  if (@input[$j] =~ /^\d+ blocks processed/ ) 
  {   $past_blocks = 1; close(FBLK);   }

  # beyond the cobbler sequence?
  if (@input[$j] =~ /end cobbler/ ) 
  {   $past_cobb = 1;   }

  # beyond the block map?
  if (@input[$j] =~ /end map/ ) 
  {   $past_map = 1; close(FMAP);  }

  # beyond the tree?
  if (@input[$j] =~ /end tree/ ) 
  {   $past_tree = 1; close(FTREE); }

  # beyond proweb yet?
  $past_web = 1 if @input[$j] =~ /end proweb/ ;

  # beyond blinks yet?
  $past_lnk = 1 if @input[$j] =~ /end blinks/ ;

  # beyond the pdb?
  if (@input[$j] =~ /end pdb/ ) 
  {   $past_pdb = 1; close(FPDB); }

  # beyond cyrca?
  if (@input[$j] =~ /end cyrca/ ) 
  {   $past_cyrca = 1; close(FCYR); }

  # read the line numbers that the blocks start at
  if (@input[$j] =~ /^ID   / && !$past_blocks)
  {  @block_lines = (@block_lines, $j); 
     @block_ID = (@block_ID, (@input[$j] =~ /^ID\s+(.*);/)) ;
     $past_intro = 1; 
  }

  # read the block numbers for use in the TOC
  @block_AC = (@block_AC, (@input[$j] =~ /^AC\s+(.*);/))  
    if @input[$j] =~ /^AC   / && !$past_blocks ;

  # 	Description
  @block_DE = (@block_DE, (@input[$j] =~ /^DE\s+(.*)$/)) 
    if @input[$j] =~ /^DE   / && !$past_blocks ;

  #	Prepare the blocks
  if ($past_intro && !$past_blocks)
  {  print FBLK @input[$j];  }

  #	Prepare the block map file, don't write "end cobbler" line
  if ($past_cobb && !$past_map && !(@input[$j] =~ /end cobbler/))
  {  print FMAP @input[$j];  }

  #	Prepare the tree file, don't write "end map" or ">" title line
  if ($past_map && !$past_tree && !(@input[$j] =~ /end map/) &&
      !(@input[$j] =~ /^>/) )
  {  print FTREE @input[$j]; $ntree +=1;  }

  #	Keep track of proweb lines
  if ($past_tree && ! $past_web && !(@input[$j] =~ /end tree/))
  { $nweb +=1; }
  #	Keep track of prodom lines
  if ($past_web && ! $past_lnk && !(@input[$j] =~ /end proweb/))
  { $nlnk +=1; }
  #	Keep track of pdb lines
  if ($past_lnk && ! $past_pdb && !(@input[$j] =~ /end blinks/))
  { $npdb +=1; }

# Double subscript blink? Save line numbers instead of array?
# if ($past_proweb && !past_lnk && !(@input[$j] =~ /end blinks/) )
#    @blink_AC = (@blink_AC, (@input[$j] =~ /^AC\s+(.*);/))  
#	blink[0] = AC, blink[1] = link, blink[2] = description
#         @blink = split(/!/, @input[$j]);

  #	Prepare the pdb file, don't write "end blinks" or ">" title line
  if ($past_lnk && !$past_pdb && !(@input[$j] =~ /end blinks/) )
  {
     #	Just pass the info to pdbfile.pl now
     print FPDB @input[$j];
  }

  #	Prepare the cyrca file, don't write "end pdb"
  if ($past_pdb && !$past_cyrca && !(@input[$j] =~ /end pdb/) )
  {
     print FCYR @input[$j];
#    print FCYR "http://blocks.fhcrc.org/~cyrca/block2set.cgi?@input[$j]"
     $ncyrca += 1;
  }

#>>>>>>> prosite stuff is currently missing; see getblock.sh
  # read the prosite ID for use in the title and TOC
  @prosite_ID = (@prosite_ID, (@input[$j] =~ /^ID\s+(.*);/)) 
    if @input[$j] =~ /^ID   / && $past_blocks && $past_web ;

  # read the prosite AC for use in the TOC
  @prosite_AC = (@prosite_AC, (@input[$j] =~ /^AC\s+(.*);/)) 
    if @input[$j] =~ /^AC   / && $past_blocks && $past_web ;

  # read the prosite DE for use in the TOC and title
  @prosite_DE = (@prosite_DE, (@input[$j] =~ /^DE\s+(.*)$/)) 
    if @input[$j] =~ /^DE   / && $past_blocks && $past_web ;

  # read the prosite DO for use in the TOC and title
  @prosite_DO = (@prosite_DO, (@input[$j] =~ /^DO\s+(.*);/)) 
    if @input[$j] =~ /^DO   / && $past_blocks && $past_web ;

  # read in the prosite sequence numbers and names
  # this reads a string like this (minus the '#'): 
#DR   P02820, OSTC_BOVIN, T; P02822, OSTC_CHICK, T; P15504, OSTC_DRONO, T;
#  $names{OSTC_BOVIN} = P02820 and $prosite{P02820} = OSTC_BOVIN

  if (@input[$j] =~ /^DR   /) {
    @bar = split(/; */, @input[$j]);

    foreach $bar (@bar) {
      ($bar =~ /(\S+)\s*,\s+(\S+)\s*,/) 
       && ($names{$2} = $1)
       && ($prosite{$1} = $2);
    }
  }
}     #  end of input[j]
#=========================================================================
#   Now entire input file has been read & stored
#   Check to see if any blocks were found
if (@block_lines == 0) {
  print "</PRE>

<TITLE>Getblock Error</TITLE>
<H1>Getblock Error</H1>
Invalid block number or nonexistent block.<P>
 
The correct format for the block number should be IPC??????.<P>
 
For example:  IPC001525, where '0' is the number zero.<P>
  
";
  exit;
}

#--------------------------------------------------------------------------

#	Get the number part of the block name
$interpro_flag = 1;
$prints_flag = 0;
$chrom_flag = 0;

(@block_AC[0] =~ /^IPC(\d+)/) && ($block_group_num = $1);
if (! $block_group_num)
{
   #		Try for chrom; doesn't have numbers in it
   (@block_AC[0] =~ /^CH(\d+)/) && ($block_group_num = $1);
   if ($block_group_num) 
   { $chrom_flag = 1; $interpro_flag = 0; 
     $group = substr(@block_AC[0], 0, 8); 
   }
   #		Try for Prints
   else
   {
      (@block_AC[0] =~ /^PC(\d+)/) && ($block_group_num = $1);
      if ($block_group_num) 
      { $prints_flag = 1; $interpro_flag = 0;
        $group = substr(@block_AC[0], 0, 7); 
        $pgroup = $group;
        $pgroup  =~ s/^PC/PR/;
      }
   }
}     #  end if not interpro
else
{   
      $group = substr(@block_AC[0], 0, 9); 
      $ipr =  "IPR".$block_group_num;
      $ipb =  "IPB".$block_group_num;
}
#print "AC=@block_AC[0] group=$group prints_flag=$prints_flag pgroup=$pgroup\n";
#exit;

if ($chrom_flag)
{   $SEQ_LINK = '<A HREF="http://chromdb.biosci.arizona.edu/cgi-bin/v3/query2.cgi?list=gene&id=';  }
else
{   $SEQ_LINK = '<A HREF="http://www.expasy.ch/cgi-bin/get-sprot-entry?';  }

#============================================================================
   
#  Output title 
print "</PRE>\n";
print "<A NAME=\"",$group,"\"></A>\n";
print "<TITLE>Blocks for $group</TITLE>\n";
print "<A NAME=top><H1>$group: @block_ID[0]</H1></A>\n";
print "@block_DE[0]<P>\n";

#print "<A NAME=top><H2>TABLE OF CONTENTS</H2></A>\n";
print "<MENU>\n";
print "<UL>\n";
#print "<LI> <A HREF=\"#bintro\">Introduction</A>\n";
print "<LI> <A HREF=\"/blocks/help/blocks_format.html\">Introduction</A>\n";
foreach $blk (@block_AC) {
  print "<LI> <A HREF=\"#", $blk, "\">Block number ", $blk, "</A>\n";
}

print "<P>\n";
if ( $prints_flag )
{
   print "<LI>";
   print "PRINTS entry <A HREF=\"http://www.bioinf.man.ac.uk/cgi-bin/dbbrowser/PRINTS/DoPRINTS.pl?cmd_a=Display&qua_a=/Brief&fun_a=Text&qst_a=$pgroup\"TARGET=\"getblock\"> $pgroup</A>";
   print " (source of blocks)\n";
}
elsif ( $chrom_flag )
{
   print "<LI>";
#>>>>>>>>>>>>this is wrong, need a cross-ref file<<<<<<<<<<<<<<<<<<
   open(GRP, "grep '^$group' ../www/chromdb/groups.list |");
   $line = <GRP>;
   ($x1, $chrom) = split(/\s+/, $line);
   close(GRP);
   print "ChromDB entry <A HREF=\"http://chromdb.biosci.arizona.edu/cgi-bin/v3/homoutput.cgi?class=$chrom&org=all&homotype=one\"TARGET=\"getblock\"> $chrom</A>" ;
   print " (source of sequences used to make blocks)\n";
}
if ( $interpro_flag )
{
   print "<LI>";
   print "InterPro entry <A HREF=\"http://www.ebi.ac.uk/interpro/IEntry?ac=$ipr\"TARGET=\"getblock\"> $ipr</A>" ;
   print " (source of sequences used to make blocks)\n";
}

print "<P>\n";
print "<LI>Block Maps.";
print "[<A HREF=\"http://www.proweb.org/test-bin/blockmap.cgi?name=$group&type=FAM&condensed=YES&dbtype=CH\">Graphical Map</A>]";
print "     [<A HREF=\"/blocks-bin/map.csh?../tmp/$$.blocks.mapfile\">Text Map</A>]";
print "     [<A HREF=\"../tmp/$$.blocks.mapfile\">Map Positions</A>]";
print "     [<A HREF=\"/blocks/help/about_maps.html\">About Maps</A>]";

#---------------------------------------------------------------------------
print "<P>\n";
print "<LI>Logos.";
print "[<A HREF=\"/blocks/help/about_logos.html\">About Logos</A>]";
print "<BR>\nSelect display format: ";
if (-e "../www/logos/$group.gif")
{  
   print "<A HREF=\"/blocks/logos/$group.gif\">[GIF]</A> \n"; 
   print "<A HREF=\"/blocks-bin/convert.csh?../www/logos/$group.gif+pdf\">[PDF]</A> ";
}
else
{  
   print "<A HREF=\"/blocks-bin/logo.csh?../tmp/$$.blocks+gif\">[GIF]</A> ";
   print "<A HREF=\"/blocks-bin/logo.csh?../tmp/$$.blocks+pdf\">[PDF]</A> ";
}
print "<A HREF=\"/blocks-bin/logo.csh?../tmp/$$.blocks+ps\">[Postscript]</A> ";

#---------------------------------------------------------------------------
#	Some groups with < 4 seqs don't have trees
if ( -s "../tmp/$$.treefile" )	# $ntree > 0 ???
{
print "<P>\n";
print "<LI>Tree from blocks alignment. ";
print "[<A HREF=\"/blocks/help/about_trees.html\" TARGET=\"tree\">About Trees</A>]  ";
print "[<A HREF=\"http://www.proweb.org/treeinfo.html\" TARGET=\"tree\">About ProWeb TreeViewer</A>]<BR>\n";

#print "<FORM METHOD=POST ACTION=\"http://www.proweb.org/proweb-bin/trees.cgi\"
# TARGET=\"tree\">
#<INPUT TYPE=hidden NAME=\"tree\" VALUE=\"\n";
#open(TRE, "<../tmp/$$.treefile");
#while ($trerec = <TRE>)
#{ print $trerec; }
#close(TRE);
#print "\"><INPUT TYPE=submit VALUE=\"ProWebTreeViewer\">\n";
#print "</FORM>\n";

#---------------------------------------------------------------------------
print "[<A HREF=\"/blocks-bin/drawgram.csh?../tmp/$$.treefile+dat\" TARGET=\"tree\">Data</A>]  ";
print " [<A HREF=\"http://www.proweb.org/proweb-bin/trees.cgi?family=$group\" TARGET=\"proweb\">ProWeb TreeViewer</A>] ";
print "[<A HREF=\"/blocks-bin/drawgram.csh?../tmp/$$.treefile+xbm\" TARGET=\"tree\">XBitmap</A>]  ";
print "[<A HREF=\"/blocks-bin/drawgram.csh?../tmp/$$.treefile+gif\" TARGET=\"tree\">GIF</A>] ";
print "[<A HREF=\"/blocks-bin/drawgram.csh?../tmp/$$.treefile+pdf\" TARGET=\"tree\">PDF</A>] ";
print "[<A HREF=\"/blocks-bin/drawgram.csh?../tmp/$$.treefile+ps\" TARGET=\"tree\">Postscript</A>]  ";
print "[<A HREF=\"/blocks-bin/drawgram.csh?../tmp/$$.treefile+new\" TARGET=\"tree\">Newick</A>]";
}
print "<P>\n";

#---------------------------------------------------------------------------
#         Send 3dblocks the blocks and cobbler sequence
print "<LI>
<FORM METHOD=POST ACTION=\"http://www.proweb.org/proweb-bin/3dblocks.cgi\"
 TARGET=\"3dblock\">
<INPUT TYPE=hidden NAME=\"SAVE\" VALUE=\"NO\">
<INPUT TYPE=hidden NAME=\"CUTOFF\" VALUE=\"0.0001\">
<INPUT TYPE=hidden NAME=\"blocks_data\" VALUE=\"\n";
open(BLK, "<../data-chromdb/blks/$group.blks");
while ($blkrec = <BLK>)
{ print $blkrec; }
close(BLK);
print "\"><INPUT TYPE=hidden NAME=\"SEQ\" VALUE=\"";
open(SEQ, "<../data-chromdb/cobs/$group.pro");
while ($seqrec = <SEQ>)
{ print $seqrec; }
close(SEQ);
print "\"><INPUT TYPE=submit VALUE=\"3D Blocks\">";
print " to find structures using COBBLER sequence (takes several minutes) ";
print "  [<A HREF=\"http://www.proweb.org/3dblocks_intro.html\" TARGET=\"3dblock\">About 3D Blocks</A>]";
print "</FORM>\n";

#if ($npdb > 0)
#
#Unfortunately, pdbfile.pl requires the full path
#  $fullblk = "/howard/servers/blocks/tmp/$$.blocks";
#  print "<P>\n<LI><A HREF=\"/blocks-bin/pdbfile.pl?../tmp/$$.pdbfile+$fullblk\">PDB entries</A>\n";
#

if ($ncyrca > 0)
{
    print "<P>\n<LI><A HREF=\"/blocks-bin/catfile.sh?../tmp/$$.cyrca\">CYRCA entries</A>";
#  print "  [<A HREF=\"http://blocks.fhcrc.org/~cyrca/index.html\">About CYRCA</A>]";
}

print "<P>\n";
print "<LI>Search Blocks vs Other Databases:<BR><UL>\n";
print "<LI><A HREF=\"#cobbler\">COBBLER sequence</A>\n";
print " and BLAST searches";
print "     [<A HREF=\"/blocks/help/about_cobbler.html\">About COBBLER</A>]";

print "<LI>";
print "<A HREF=\"/blocks-bin/mast.sh?../tmp/$$.blocks+$group\" TARGET=\"getblock\">MAST Search</A>";
print " of all blocks vs a sequence database ";
print "     [<A HREF=\"http://meme.sdsc.edu/meme/website/mast-intro.html\">About MAST</A>]";

print "<LI>";
print "<A HREF=\"/blocks-bin/LAMA_search.sh?../tmp/$$.blocks\" TARGET=\"_blank\">LAMA search</A>";
print " of all blocks vs a blocks database  ";
print "     [<A HREF=\"/blocks/help/LAMA_help.html\">About LAMA</A>]";
print "</UL>";

print "<P>\n";
print "<LI>";
print "<A HREF=\"/blocks-bin/codehop.sh?../tmp/$$.blocks\">CODEHOP</A>";
print " to design PCR primers from blocks";
print "     [<A HREF=\"/blocks/help/CODEHOP/tips.html\">About CODEHOP</A>]";

print "<P>\n";
print "<LI>";
print "<A HREF=\"/blocks-bin/sift.sh?../tmp/$$.blocks\" TARGET=\"_blank\">SIFT</A>";
print " to predict amino acid substitutions in blocks";
print "     [<A HREF=\"/sift/SIFT_on_blocks.html\">About SIFT</A>]";
print "<P>\n";

#
if ($nweb > 0 || $nlnk >0)
{print "<LI> <A HREF=\"#proweb\">Additional Links</A>\n"};
print "</UL>\n";
print "</MENU>\n";
print "<P><HR>\n";

print "<PRE>\n";

#------------------------------------------------------------------------
# Output Intro (header & body) (to first block)
#print "<A NAME=\"bintro\"><H3>Introduction</H3></A><P>\n";

for ($i=0; $i<@block_lines[0]; ++$i) {
  print @input[$i];
}

print "</PRE><A HREF=\"#top\">[Return to top]</A><P><PRE>\n";


#-----------------------------------------------------------------------
# Output Blocks headers and bodies with links inserted
# Trembl accession numbers are 6 chars:
#	1=[O,P,Q] 2=[0-9] 3,4,5=[A-Z,0-9] 6=[0-9]
# Swiss ids have and "_" in the middle
# Chromdb accession numbers are 9 chars:
#	1,2,3=[A-Z] 4,5,6,6,8,9=[0-9]
#>>> attach blinks to appropriate block somehow

#	$id is the full seq id in the block, $id1 is left, $id2 is right of |
for ($j=0; $j<@block_lines ; ++$j) {
  print "<H3><A NAME=\"", @block_AC[$j],"\">Block ", @block_AC[$j], "</A>";
  print "</H3>";
  print "<P>\n";
  $i = @block_lines[$j];
  #	could just use $id here instead of @names{$id} ...
  while ($i < @input) {
    if (@input[$i] =~ /\s*\S+\ +\(/) 
    {
      (@input[$i] =~ /^\s*(\S+)/) && ($id = $1) if @input[$i] =~ /\s*\S+\ +\(/;
      #	Assume we have a sequence name == $id; could be "id|ac"
      ($id1, $id2) = split(/\|/, $id);
      #  Figure out the id type of each part
      $id1type = $id2type = 0;
      if    ($id1 =~ m/\_/) { $id1type = 1; }    #swiss
      elsif ($id1 =~ m/[O-Q][0-9][A-Z,0-9]{3,3}[0-9]/) { $id1type = 2; } #trembl
      elsif ($id1 =~ m/[A-Z]{3,3}[0-9]{6,6}/) { $id1type = 3; }    #chromdb
      if ($id2 ne "")
      {
      if    ($id2 =~ m/\_/) { $id2type = 1; }    #swiss (shouldn't happen)
      elsif ($id2 =~ m/[O-Q][0-9][A-Z,0-9]{3,3}[0-9]/) { $id2type = 2; } #trembl
      elsif ($id2 =~ m/[A-Z]{3,3}[0-9]{6,6}/) { $id2type = 3; }    #chromdb
      }
      #  If swiss and/or trembl, just link to trembl
      if (($id1type == 2) || ($id1type == 1 && $id2type != 2)) 
      {
         @input[$i] =~ s|(\s*)(\S+)(.*)|\1<A HREF="http://www.expasy.ch/cgi-bin/get-sprot-entry?$id1">\2</A>\3|;
      }
      elsif ($id1type < 3 && $id2type == 2)
      {
         @input[$i] =~ s|(\s*)(\S+)(.*)|\1<A HREF="http://www.expasy.ch/cgi-bin/get-sprot-entry?$id2">\2</A>\3|;
      }
      #   Swiss only (shouldn't happen any more)
      elsif ($id1type == 1 && $id2type == 0)
      {
         @input[$i] =~ s|(\s*)(\S+)(.*)|\1<A HREF="http://www.expasy.ch/cgi-bin/get-sprot-entry?$id1">\2</A>\3|;
      }
      #	  ChromDB only
      elsif ($id1type == 3 && $id2type == 0)
      {
         @input[$i] =~ s|(\s*)(\S+)(.*)|\1<A HREF="http://chromdb.biosci.arizona.edu/cgi-bin/v3/query.cgi?list=gene&id=$id1">\2</A>\3|;
      }
      #	  ChromDB and Trembl (2 links)
      elsif ($id1type == 3 && $id2type == 2)
      {
         @input[$i] =~ s#(\s*)(\S+)\|(\S+)(.*)#\1<A HREF="http://chromdb.biosci.arizona.edu/cgi-bin/v3/query.cgi?list=gene&id=$id1">\2</A>\|<A HREF="http://www.expasy.ch/cgi-bin/get-sprot-entry?$id2">\3</A>\4#;
      }
      #  Other combos shouldn't be encountered


      print @input[$i];
    }   # end of sequence line in block

    else { print @input[$i]; }

    last if @input[$i] =~ /^\/\//;
    ++$i;
  }
  print "</PRE><A HREF=\"#top\">[Return to top]</A><P><PRE>\n";
}

# Output lines following blocks ?

++$i;
while ($i < @input) {
  last if @input[$i] =~ /^\d* blocks processed/;
  print @input[$i];
  ++$i;
}
#		Don't print the "blocks processed" line
++$i;
 
#-------------------------------------------------------------
# Output the COBBLER sequence
print "<A NAME=\"cobbler\"><H3>COBBLER sequence (region containing Blocks only)</H3></A>";

print "</PRE>To do a BLAST search, copy the cobbler sequence";
print " below then click on a BLAST link:<PRE><BR>\n";

#	Print sequence for blast search
print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/blast.cgi?Jform=0&PROGRAM=blastp&DESCRIPTIONS=500&ALIGNMENTS=100&DATALIB=nr&SEQUENCE=";
#   Save the position, put the sequence in the blast link
$j = $i;
$i = $j + 3;				#  First three lines are headings
while ($i < @input) {
  last if @input[$i] =~ /end cobbler/;
  print @input[$i];
  ++$i;
}
print "\" TARGET=\"getblock\">Blast Search</A>]</B>";
#print " (copy the COBBLER sequence first)\n";

#	Print sequence again for gap-blast search
print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/blast.cgi?Jform=1&PROGRAM=blastp&GAPPED_ALIGNMENT=is_set&DESCRIPTIONS=500&ALIGNMENTS=100&DATALIB=nr&SEQUENCE=";
$i = $j + 3;
while ($i < @input) {
  last if @input[$i] =~ /end cobbler/;
  print @input[$i];
  ++$i;
}
print "\" TARGET=\"getblock\">Gap-Blast Search</A>]</B>";
#print " (copy the COBBLER sequence first)\n";

#	Print sequence again for psi-blast search
#print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/psiblast.cgi?PROGRAM=blastp&GAPPED_ALIGNMENT=is_set&DESCRIPTIONS=500&ALIGNMENTS=100&DATALIB=nr&SEQUENCE=";
print "<B>[<A HREF=\"http://www.ncbi.nlm.nih.gov/blast/Blast.cgi?AUTO_FORMAT=Semiauto&PROGRAM=blastp&RUN_PSIBLAST=on&CDD_SEARCH=on&COMPOSITION_BASED_STATISTICS=on&DESCRIPTIONS=250&ALIGNMENTS=100&ALIGNMENT_VIEW=Pairwise&DATABASE=nr&END_OF_HTTPGET=YES&QUERY=";
$i = $j + 3;
while ($i < @input) {
  last if @input[$i] =~ /end cobbler/;
  print @input[$i];
  ++$i;
}
print "\" TARGET=\"getblock\">PSI-Blast Search</A>]</B>";
#print " using default parameters\n";
#print " (copy the COBBLER sequence first)\n";
print "\n";

#	Print sequence to the screen
$i = $j;
while ($i < @input) {
  last if @input[$i] =~ /end cobbler/;
  print @input[$i];
  ++$i;
}
#	Don't print the "end cobbler" line
++$i;
print "</PRE><A HREF=\"#top\">[Return to top]</A><P><PRE>\n";

#-------------------------------------------------------------------------
#	Skip to the proweb entries

while ($i < @input) {
  last if @input[$i] =~ /end tree/;
  ++$i;
}
#		Skip the "end tree" line
++$i;

#-------------------------------------------------------------
if ($nweb > 0 || $nlnk > 0)
{
print "<A NAME=\"proweb\"><H3>Additional Links (separate browser window)</H3></A><P>\n";
}

# Output the PROWEB data
while ($i < @input) {
  last if @input[$i] =~ /end proweb/;
  if (!(@input[$i] =~ /^ID   /) && !(@input[$i] =~ /^AC   /) &&
      !(@input[$i] =~ /^\/\//) ) {
          @bar = split(/!/, @input[$i]);
          print "<A HREF=\"", @bar[0], "\"TARGET=\"getblock\">",  @bar[1], "</A>";
          print "<BR>\n";
  }
  ++$i;
}
#	Don't print the "end proweb" line
++$i;

#	Metafam link
if ($prints_flag)
{
   print "<A HREF=\"http://metafam.ahc.umn.edu/superSet.jsp?SearchType=ProteinFamily&SearchString=$pgroup&ShowApplet=TRUE&ShowKeys=TRUE\"TARGET=\"getblock\">MetaFam $pgroup</A>";
}
elsif ($interpro_flag)
{
   print "<A HREF=\"http://metafam.ahc.umn.edu/superSet.jsp?SearchType=ProteinFamily&SearchString=$ipb&ShowApplet=TRUE&ShowKeys=TRUE\"TARGET=\"getblock\">MetaFam $ipb</A>";
}

#	Output a blank line before starting more links
print "<P>\n";

# Output the ProDom links in blinks.dat - note @bar[] subscripts!
while ($i < @input) {
  last if @input[$i] =~ /end blinks/;
  if (!(@input[$i] =~ /^ID   /) && !(@input[$i] =~ /^AC   /) &&
      !(@input[$i] =~ /^\/\//) ) {
          @bar = split(/!/, @input[$i]);
          print "@bar[0]: <A HREF=\"", @bar[1], "\"TARGET=\"getblock\">",  @bar[2], "</A>";
          print "<BR>\n";
  }
  ++$i;
}
#	Don't print the "end blinks" line
++$i;

#--------------------------------------------------------------

exit;
